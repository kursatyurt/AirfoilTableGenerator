#!/usr/bin/env python
"""Digitise NACA 23015 lift data from figure 28 of NACA Report 832.

Report 832 (Graham, Nitzberg & Olson, 1945), Ames 1x3.5 ft high-speed tunnel.
Figure 28 (top of PDF page 48) plots section lift coefficient against angle of
attack for M = 0.300 to 0.825 as a staggered carpet, each curve offset 4 deg
from the previous one.

    python reference/digitize_tr832_23015.py       # fetches the PDF if needed

Method, and why each step is there -- identical in structure to the NACA 0015
extractor in reference/digitize_tr832.py, recalibrated for this page:
  * page 48 is rendered at 400 dpi and deskewed (this page is tilted -0.46 deg,
    which is enough to defeat row/column gridline detection)
  * gridlines are removed as long straight runs; this fragments the curves where
    they cross, so a morphological closing reconnects them
  * for each Mach the straight part of the curve is fitted iteratively inside a
    band, which is robust to the Mach labels drawn along the curves

Calibration was re-derived from this page and not copied: the 36 vertical and 23
horizontal gridlines were located directly, giving 79.20 px per grid square,
0.1 cl and 2 deg per square, cl = 0 on the 15th row from the top, and the 14
curve origins on the 14 axis labels reading "0" (the first at x = 487.5 px).

Self-check, replacing the symmetry test used for the symmetric 0015: the 23015
is cambered, so thin-airfoil theory on the 230 mean line (p = 0.15, r = 0.2025,
k1 = 15.957) is used instead -- it gives alpha_L0 = -1.09 deg, and the published
experimental value is about -1.2 deg. ZERO_LIFT_REF/TOL below assert that every
one of the fourteen fitted curves lands there. Because Prandtl-Glauert scales cl
but not alpha_L0, the fourteen values must also agree with EACH OTHER; the
spread is printed and asserted. An axis error of one gridline is 2 deg and fails
both tests loudly.

Scope: lift slopes only. No CLmax is emitted at all, unlike the 0015. The 23015
carries much more lift, so its curves peak higher and overlap heavily, and the
same peak tracer that worked on figure 27 was found to jump between neighbours
here: the M = 0.400 trace merges onto the M = 0.300 curve just past their
crossing near cl = 1.18, and reports its maximum on the transition rather than
on the real peak (1.19 instead of the ~1.31 that a direct pixel read of the scan
gives). Only M = 0.300 came out right, and one point is not a CLmax-versus-Mach
table, so the peak stage was removed rather than shipped with caveats.
"""
import subprocess, urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PDF = HERE / "naca_tr832.pdf"
URL = "https://ntrs.nasa.gov/api/citations/19930091909/downloads/19930091909.pdf"

MACH = [0.300, 0.400, 0.500, 0.550, 0.600, 0.625, 0.650, 0.675,
        0.700, 0.725, 0.750, 0.775, 0.800, 0.825]
SKEW_DEG = -0.46
# calibration of the deskewed 400 dpi crop of page 48, in pixels
Y_CL0, PX_PER_CL, PX_PER_DEG, STAGGER, X_M030 = 1442.5, 792.0, 39.60, 158.40, 487.5
# thin-airfoil zero-lift angle of the 230 mean line, from mean_line_alpha_L0()
ZERO_LIFT_REF, ZERO_LIFT_TOL = -1.09, 0.9


def mean_line_alpha_L0(r=0.2025, k1=15.957):
    """Thin-airfoil alpha_L0 (deg) of the NACA 5-digit 230 mean line."""
    th = np.linspace(0.0, np.pi, 200001)
    x = (1 - np.cos(th)) / 2
    dyc = np.where(x < r, k1 / 6 * (3 * x**2 - 6 * r * x + r**2 * (3 - r)),
                   -k1 * r**3 / 6)
    return np.degrees(-np.trapezoid(dyc * (np.cos(th) - 1), th) / np.pi)


def load_figure():
    from PIL import Image, ImageOps
    from scipy import ndimage
    if not PDF.exists():
        print(f"fetching {URL}")
        urllib.request.urlretrieve(URL, PDF)
    png = HERE / "_tr832_p48.png"
    if not png.exists():
        subprocess.run(["pdftoppm", "-f", "48", "-l", "48", "-r", "400", "-png",
                        "-singlefile", str(PDF), str(png.with_suffix(""))], check=True)
    im = Image.open(png).convert("L")
    W, H = im.size
    im = ImageOps.autocontrast(im.crop((int(W * .06), int(H * .02), int(W * .97), int(H * .50))))
    im = im.rotate(SKEW_DEG, resample=Image.BICUBIC, fillcolor=255)
    mask = np.array(im) < 190
    for axis in (1, 0):                                   # strip gridlines
        mask = _strip(mask, axis, 60)
    return ndimage.binary_closing(mask, structure=np.ones((7, 7)))


def _strip(mask, axis, maxlen):
    """Blank straight runs longer than maxlen -- gridlines, not curve segments."""
    m = mask.copy()
    for k in range(m.shape[1] if axis == 0 else m.shape[0]):
        idx = np.nonzero(m[:, k] if axis == 0 else m[k, :])[0]
        if not len(idx):
            continue
        start = prev = idx[0]
        for i in list(idx[1:]) + [10 ** 9]:
            if i != prev + 1:
                if prev - start + 1 > maxlen:
                    if axis == 0:
                        m[start:prev + 1, k] = False
                    else:
                        m[k, start:prev + 1] = False
                start = i
            prev = i
    return m


def fit_linear(ink_x, ink_y, x0):
    """Iteratively fit the straight part of one curve. Returns (slope/deg, alpha0, line)."""
    slope = -2.0                                  # 0.1 cl/deg in this pixel scale
    inter = Y_CL0 + 0.11 * PX_PER_CL * -1 + 2.0 * x0   # camber puts cl(0) near +0.11
    sel = None
    for _ in range(6):
        # cambered section: the straight run sits high, roughly -0.3 <= cl <= +0.7
        band = (ink_y > Y_CL0 - 0.70 * PX_PER_CL) & (ink_y < Y_CL0 + 0.30 * PX_PER_CL)
        sel = band & (np.abs(ink_y - (inter + slope * ink_x)) < 28)
        if sel.sum() < 200:
            return None
        slope, inter = np.polyfit(ink_x[sel], ink_y[sel], 1)
    return (-slope * PX_PER_DEG / PX_PER_CL,
            ((Y_CL0 - inter) / slope - x0) / PX_PER_DEG,
            (slope, inter), int(sel.sum()))


def main():
    a_l0 = mean_line_alpha_L0()
    print(f"thin-airfoil alpha_L0 of the 230 mean line: {a_l0:+.2f} deg")
    assert abs(a_l0 - ZERO_LIFT_REF) < 0.05, "mean-line integration changed"

    mask = load_figure()
    iy, ix = np.nonzero(mask)
    ix, iy = ix.astype(float), iy.astype(float)

    slopes = []
    print(f"{'M':<7}{'slope /deg':>11}{'zero-lift deg':>15}{'px':>7}")
    for k, M in enumerate(MACH):
        x0 = X_M030 + STAGGER * k
        fit = fit_linear(ix, iy, x0)
        if fit is None:
            print(f"{M:<7}  linear fit failed")
            continue
        slope, a0, _line, npix = fit
        assert abs(a0 - a_l0) < ZERO_LIFT_TOL, \
            f"M={M}: zero-lift {a0:.2f} deg vs theory {a_l0:.2f} -- calibration is wrong"
        slopes.append((M, slope, a0))
        print(f"{M:<7}{slope:>11.4f}{a0:>15.2f}{npix:>7}")

    a0s = [a for _, _, a in slopes]
    spread = max(a0s) - min(a0s)
    print(f"\nzero-lift angle spread across {len(a0s)} curves: {spread:.2f} deg "
          f"(mean {np.mean(a0s):+.2f}, theory {a_l0:+.2f})")
    assert spread < 0.5, "zero-lift angle drifts with Mach -- the stagger is wrong"

    hdr = ("# NACA 23015, digitised from figure 28 of NACA Report 832 (Graham, Nitzberg\n"
           "# & Olson, 1945), Ames 1x3.5 ft high-speed tunnel, 1e6 <= Re <= 2e6.\n"
           "# Produced by reference/digitize_tr832_23015.py -- see that file for method\n"
           "# and limits.\n"
           "# zero_lift_deg is a self-check, not data: the section is cambered, so instead\n"
           "# of symmetry the check is thin-airfoil theory on the 230 mean line, which gives\n"
           f"# alpha_L0 = {a_l0:+.2f} deg (published experiment: about -1.2 deg). Prandtl-Glauert\n"
           "# leaves alpha_L0 Mach-independent, so the fourteen values must also agree with\n"
           f"# each other; they span {spread:.2f} deg. That validates the axis calibration,\n"
           "# the 4 deg stagger, and each curve's identity.\n")
    (HERE / "naca23015_lift_slope_vs_mach.csv").write_text(
        hdr + "mach,lift_slope_per_deg,zero_lift_deg\n"
        + "".join(f"{M},{s:.4f},{a:+.2f}\n" for M, s, a in slopes))
    print(f"\nwrote naca23015_lift_slope_vs_mach.csv ({len(slopes)} rows); "
          "no CLmax file -- see the scope note in this file's docstring")


if __name__ == "__main__":
    main()
