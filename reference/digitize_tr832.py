#!/usr/bin/env python
"""Digitise NACA 0015 lift data from figure 27 of NACA Report 832.

Report 832 (Graham, Nitzberg & Olson, 1945) is the original source behind the
NACA 0015 sheet in the Army Helicopter Design Datcom, and NASA's scan of it is
far cleaner than the Datcom photocopy: crisp line art on white rather than a
halftone screen. Figure 27 plots section lift coefficient against angle of
attack for M = 0.300 to 0.825 as a staggered carpet, each curve offset 4 deg
from the previous one.

    python reference/digitize_tr832.py            # fetches the PDF if needed

Method, and why each step is there:
  * page 47 is rendered at 400 dpi and deskewed (the scan is tilted 0.51 deg,
    which is enough to defeat row/column gridline detection)
  * gridlines are removed as long straight runs; this fragments the curves where
    they cross, so a morphological closing reconnects them
  * for each Mach the straight part of the curve is fitted iteratively inside a
    band, which is robust to the Mach labels drawn along the curves
  * the peak is then traced upward from the top of that straight segment

Self-check: the section is symmetric, so every curve must pass through CL = 0 at
its own alpha = 0. The fitted zero-lift angle is reported for all fourteen and is
within +/-0.21 deg -- an independent confirmation of the axis calibration, the
4 deg stagger, and each curve's identity. An earlier calibration attempt that
was one gridline out produced a 2 deg zero-lift angle and was rejected by this
same test.

Scope: the peak trace is only trustworthy for M = 0.300, 0.400 and 0.500. Above
that the curves crowd together and the tracer jumps to a neighbour, so no CLmax
is emitted for M >= 0.55. Lift slopes are emitted for all fourteen.
"""
import subprocess, sys, urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PDF = HERE / "naca_tr832.pdf"
URL = "https://ntrs.nasa.gov/api/citations/19930091909/downloads/19930091909.pdf"

MACH = [0.300, 0.400, 0.500, 0.550, 0.600, 0.625, 0.650, 0.675,
        0.700, 0.725, 0.750, 0.775, 0.800, 0.825]
PEAK_OK = {0.300, 0.400, 0.500}          # see the scope note above
SKEW_DEG = 0.512
# calibration of the deskewed 400 dpi crop, in pixels
Y_CL0, PX_PER_CL, PX_PER_DEG, STAGGER, X_M030 = 1140.0, 790.0, 39.55, 156.06, 546.0


def load_figure():
    from PIL import Image, ImageOps
    from scipy import ndimage
    if not PDF.exists():
        print(f"fetching {URL}")
        urllib.request.urlretrieve(URL, PDF)
    png = HERE / "_tr832_p47.png"
    if not png.exists():
        subprocess.run(["pdftoppm", "-f", "47", "-l", "47", "-r", "400", "-png",
                        "-singlefile", str(PDF), str(png.with_suffix(""))], check=True)
    im = Image.open(png).convert("L")
    W, H = im.size
    im = ImageOps.autocontrast(im.crop((int(W * .06), int(H * .51), int(W * .97), int(H * .95))))
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
    slope, inter = -2.0, Y_CL0 + 2.0 * x0
    sel = None
    for _ in range(6):
        band = (ink_y > Y_CL0 - 0.33 * PX_PER_CL) & (ink_y < Y_CL0 + 0.55 * PX_PER_CL)
        sel = band & (np.abs(ink_y - (inter + slope * ink_x)) < 28)
        if sel.sum() < 200:
            return None
        slope, inter = np.polyfit(ink_x[sel], ink_y[sel], 1)
    return (-slope * PX_PER_DEG / PX_PER_CL,
            ((Y_CL0 - inter) / slope - x0) / PX_PER_DEG,
            (slope, inter), int(sel.sum()))


def centres(mask, x, ylo, yhi):
    x, ylo, yhi = int(round(x)), max(0, int(ylo)), min(mask.shape[0], int(yhi))
    col, out, s = mask[ylo:yhi, x], [], None
    for i, v in enumerate(list(col) + [False]):
        if v and s is None:
            s = i
        elif not v and s is not None:
            if i - s <= 24:
                out.append(ylo + (s + i - 1) / 2)
            s = None
    return out


def trace_peak(mask, xs, ys):
    y, dy, pts, miss, x = float(ys), -2.0, {}, 0, int(round(xs))
    while 40 < x < mask.shape[1] - 40 and 25 < y < 1755 and miss < 60:
        x += 1
        pred = y + dy
        cand = [c for c in centres(mask, x, pred - 13, pred + 14) if abs(c - pred) < 11]
        if cand:
            c = min(cand, key=lambda v: abs(v - pred))
            dy = 0.6 * dy + 0.4 * (c - y)
            y, miss, pts[x] = c, 0, c
        else:
            y, miss = pred, miss + 1
    return pts


def main():
    mask = load_figure()
    iy, ix = np.nonzero(mask)
    ix, iy = ix.astype(float), iy.astype(float)

    slopes, peaks = [], []
    print(f"{'M':<7}{'slope /deg':>11}{'zero-lift deg':>15}{'px':>7}{'CLmax':>9}{'a_stall':>9}")
    for k, M in enumerate(MACH):
        x0 = X_M030 + STAGGER * k
        fit = fit_linear(ix, iy, x0)
        if fit is None:
            print(f"{M:<7}  linear fit failed")
            continue
        slope, a0, (spx, inter), npix = fit
        assert abs(a0) < 0.5, f"M={M}: zero-lift {a0:.2f} deg -- calibration is wrong"
        slopes.append((M, slope, a0))

        ytop = Y_CL0 - 0.36 * PX_PER_CL
        xs = (ytop - inter) / spx
        c = centres(mask, xs, ytop - 30, ytop + 31)
        pts = trace_peak(mask, xs, min(c, key=lambda v: abs(v - ytop)) if c else ytop)
        cl_max = a_max = None
        if M in PEAK_OK and pts:
            dd = sorted(((x - x0) / PX_PER_DEG, (Y_CL0 - y) / PX_PER_CL) for x, y in pts.items())
            a_max, cl_max = max(dd, key=lambda p: p[1])
            peaks.append((M, cl_max, a_max))
        print(f"{M:<7}{slope:>11.4f}{a0:>15.2f}{npix:>7}"
              f"{('%9.3f' % cl_max) if cl_max else '       --'}"
              f"{('%9.1f' % a_max) if a_max else '       --'}")

    hdr = ("# NACA 0015, digitised from figure 27 of NACA Report 832 (Graham, Nitzberg\n"
           "# & Olson, 1945), Ames 1x3.5 ft high-speed tunnel, 1e6 <= Re <= 2e6.\n"
           "# Produced by reference/digitize_tr832.py -- see that file for method and limits.\n"
           "# zero_lift_deg is a self-check, not data: symmetry demands 0, and every curve\n"
           "# lands within 0.21 deg, which validates the axis calibration and curve identity.\n")
    (HERE / "naca0015_lift_slope_vs_mach.csv").write_text(
        hdr + "mach,lift_slope_per_deg,zero_lift_deg\n"
        + "".join(f"{M},{s:.4f},{a:+.2f}\n" for M, s, a in slopes))
    (HERE / "naca0015_clmax_vs_mach.csv").write_text(
        hdr + "# Only M <= 0.500: above that the curves crowd and the peak tracer jumps\n"
        "# to a neighbouring curve, so no value is claimed.\n"
        "mach,cl_max,alpha_stall_deg\n"
        + "".join(f"{M},{c:.3f},{a:.1f}\n" for M, c, a in peaks))
    print(f"\nwrote naca0015_lift_slope_vs_mach.csv ({len(slopes)} rows) "
          f"and naca0015_clmax_vs_mach.csv ({len(peaks)} rows)")


if __name__ == "__main__":
    main()
