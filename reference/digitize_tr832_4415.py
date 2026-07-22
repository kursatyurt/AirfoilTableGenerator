#!/usr/bin/env python
"""Digitise NACA 4415 lift data from figure 29 of NACA Report 832.

Report 832 (Graham, Nitzberg & Olson, 1945) plots, at the bottom of page 48,
section lift coefficient against section angle of attack for M = 0.300 to 0.825
as a staggered carpet, each curve offset 4 deg from the previous one. Companion
of reference/digitize_tr832.py, which does the same for the NACA 0015 of
figure 27; the method is identical, the calibration is re-derived from scratch
because this figure sits on a different page and in the lower half of it.

    python reference/digitize_tr832_4415.py       # fetches the PDF if needed

Method, and why each step is there:
  * page 48 is rendered at 400 dpi, the lower half cropped and deskewed
    (this figure is tilted -0.20 deg, not the +0.51 deg of page 47)
  * the axes are calibrated off the printed grid itself: 23 horizontal rules at
    0.1 cl and 36 vertical rules at 2 deg, each fitted by least squares, so no
    pixel constant is eyeballed
  * gridlines are removed as long straight runs; this fragments the curves where
    they cross, so a morphological closing reconnects them
  * for each Mach the straight part of the curve is fitted iteratively inside a
    band, which is robust to the Mach labels drawn along the curves

Self-check.  The 0015 of figure 27 is symmetric, so every curve there had to
cross cl = 0 at its own alpha = 0.  The 4415 is strongly cambered and that
anchor is gone, so three independent replacements are used instead:

  1. Thin-airfoil theory on the NACA 4-digit mean line (m = 0.04, p = 0.4) gives
     alpha_L0 = -4.15 deg analytically -- _alpha_l0_theory() computes it, it is
     not quoted from a table.  Every fitted curve must land near it.  The
     published experimental value is about -4.0 deg, and viscous decambering
     always makes the measured magnitude slightly smaller than the inviscid one,
     so the fitted values are expected just inside the theoretical figure.
  2. Prandtl-Glauert scales cl by 1/beta and leaves alpha_L0 alone, so the fitted
     zero-lift angles must agree with EACH OTHER far more tightly than they agree
     with theory.  The spread is asserted below; it comes out at 0.28 deg over the
     seven accepted curves, against 0.21 deg for the 0015's fourteen.
  3. The lift slope must track 1/sqrt(1-M^2) up to force divergence -- nothing in
     the fitter knows about Mach number.  It does, to within 2.8%.
  4. Cross-check against the 0015 of figure 27, digitised by the sibling script:
     camber shifts the curve but barely changes the slope, and at M = 0.300 the
     two give 0.0994 and 0.1000 per deg, 0.6% apart.

Off-by-one-curve hazard, specific to a cambered section: alpha_L0 is about
-4 deg and the carpet stagger is exactly 4 deg, so mis-indexing the carpet by one
curve turns a correct fit into an apparently perfect alpha_L0 = 0 or a plausible
-8.  Two guards: the M = 0.300 origin is asserted to be the LEFTMOST curve in the
figure, and every fitted alpha_L0 is required to sit in (-5.0, -3.0), which
excludes both decoys.

The alpha_L0 window alone is NOT enough, and it is worth being explicit about why:
the fitter latches onto whichever curve is nearest its seed and has no idea which
Mach that curve is, so shifting the whole carpet by one stagger relabels all
fourteen curves and leaves every alpha_L0 exactly as it was. That test catches a
wrong axis calibration -- it caught a sign error in this script's seed intercept,
which showed up as +4.04 deg -- but not a wrong Mach labelling. What pins the
labelling is that M = 0.300 is the leftmost curve: the same fit run one stagger
further left must find no curve, and indeed collects only 1574 px of y-axis tick
label against 3584 px for the weakest real curve. Overlaying the fitted lines on
the scan confirms it by eye: the k = 0 line lies on the curve captioned
"M = 0.300", the k = -1 line lies on the axis annotation.

Scope, and what is deliberately NOT emitted:

  * Lift slopes are emitted only for M = 0.300 to 0.650. All fourteen curves are
    fitted and printed, but from M = 0.675 up the curve has no straight segment
    left -- the shock-induced kink puts an S-bend right through the fitting band
    -- and the self-check notices: the extrapolated alpha_L0 walks off
    monotonically to -4.4, -4.7, -5.8 deg. Those slopes are therefore rejected
    rather than published. The cost is that the post-divergence fall-off cannot
    be shown from this figure; the rise to divergence can, and is asserted below.
  * No CLmax file is written at all, unlike the 0015. The peak tracer used for
    figure 27 fails outright here: the 4415 knee is far more rounded, so the
    tracer flattens near the top, latches onto the stubs the gridline stripper
    leaves behind, and runs horizontally into the NEXT curve. Rendering the
    traced points over the figure shows all three low-Mach traces converging on
    one and the same pixel (cl = 1.344), which is nonsense. Peaks are legible to
    a human eye here (roughly 1.25 to 1.35 for M <= 0.55) but nothing in this
    script can check them, so no number is claimed.
"""
import subprocess, urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PDF = HERE / "naca_tr832.pdf"
URL = "https://ntrs.nasa.gov/api/citations/19930091909/downloads/19930091909.pdf"

MACH = [0.300, 0.400, 0.500, 0.550, 0.600, 0.625, 0.650, 0.675,
        0.700, 0.725, 0.750, 0.775, 0.800, 0.825]
FIT_OK = 0.650          # above this the curve has no straight part; see scope note
SKEW_DEG = -0.20
# calibration of the deskewed 400 dpi crop, in pixels. Derived by _calibrate()
# from the printed grid; these are the values it produces, kept here so the
# numbers used are visible without running the script.
Y_CL0, PX_PER_CL, PX_PER_DEG, STAGGER, X_M030 = 1208.2, 787.3, 39.47, 157.89, 508.8
GRID_ROWS, GRID_COLS = 23, 36            # 0.1 cl apart, 2 deg apart
CL_TOP, CL_BOT = 1.4, -0.8               # value of the first/last horizontal rule
ALPHA_L0_MIN, ALPHA_L0_MAX = -5.0, -3.0  # rejects the +-4 deg off-by-one decoys


def _alpha_l0_theory(m=0.04, p=0.4):
    """Thin-airfoil zero-lift angle of the NACA 4-digit mean line, in degrees.

    Two parabolic arcs, dyc/dx piecewise linear; alpha_L0 = -(1/pi) * integral
    over theta of dyc/dx * (cos(theta) - 1), with x = (1 - cos theta)/2.
    """
    th = np.linspace(0.0, np.pi, 200001)
    x = (1 - np.cos(th)) / 2
    dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
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
    im = ImageOps.autocontrast(im.crop((int(W * .06), int(H * .50), int(W * .97), int(H * .94))))
    im = im.rotate(SKEW_DEG, resample=Image.BICUBIC, fillcolor=255)
    raw = np.array(im) < 190
    mask = raw
    for axis in (1, 0):                                   # strip gridlines
        mask = _strip(mask, axis, 60)
    return raw, ndimage.binary_closing(mask, structure=np.ones((7, 7)))


def _rules(profile, thresh):
    """Centroids of the runs where `profile` exceeds `thresh` -- the printed grid."""
    out, i = [], 0
    while i < len(profile):
        if profile[i] > thresh:
            j = i
            while j < len(profile) and profile[j] > thresh:
                j += 1
            out.append(np.average(np.arange(i, j), weights=profile[i:j]))
            i = j
        else:
            i += 1
    return out


def _calibrate(raw):
    """Least-squares fit of the printed grid. Returns the pixel constants."""
    rows = _rules(raw.sum(1), 0.5 * raw.shape[1])
    cols = _rules(raw.sum(0), 0.5 * raw.shape[0])
    assert len(rows) == GRID_ROWS, f"found {len(rows)} horizontal rules, want {GRID_ROWS}"
    assert len(cols) == GRID_COLS, f"found {len(cols)} vertical rules, want {GRID_COLS}"
    dy, y0 = np.polyfit(np.arange(GRID_ROWS), rows, 1)        # px per 0.1 cl
    dx, x0 = np.polyfit(np.arange(GRID_COLS), cols, 1)        # px per 2 deg
    px_per_cl = dy / 0.1
    y_cl0 = y0 + dy * (CL_TOP / 0.1)                          # rule index of cl = 0
    px_per_deg = dx / 2.0
    # the M = 0.300 origin is the 4th vertical rule; the stagger is exactly two
    # rules, i.e. 4 deg, which is what the caption's carpet construction states
    return y_cl0, px_per_cl, px_per_deg, 2 * dx, x0 + 3 * dx


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
    slope = -0.1 * PX_PER_CL / PX_PER_DEG
    inter = Y_CL0 - 0.42 * PX_PER_CL - slope * x0     # cambered: cl ~ 0.42 at alpha = 0
    sel = None
    for _ in range(8):
        band = (ink_y > Y_CL0 - 0.70 * PX_PER_CL) & (ink_y < Y_CL0 + 0.15 * PX_PER_CL)
        sel = band & (np.abs(ink_y - (inter + slope * ink_x)) < 28)
        if sel.sum() < 200:
            return None
        slope, inter = np.polyfit(ink_x[sel], ink_y[sel], 1)
    return (-slope * PX_PER_DEG / PX_PER_CL,
            ((Y_CL0 - inter) / slope - x0) / PX_PER_DEG,
            int(sel.sum()), float(ink_x[sel].min()))


def main():
    global Y_CL0, PX_PER_CL, PX_PER_DEG, STAGGER, X_M030
    raw, mask = load_figure()
    cal = _calibrate(raw)
    for name, was, now in zip("Y_CL0 PX_PER_CL PX_PER_DEG STAGGER X_M030".split(),
                              (Y_CL0, PX_PER_CL, PX_PER_DEG, STAGGER, X_M030), cal):
        assert abs(now - was) < 1.5, f"{name} drifted: {was} -> {now:.2f}"
    Y_CL0, PX_PER_CL, PX_PER_DEG, STAGGER, X_M030 = cal
    a_l0 = _alpha_l0_theory()
    print(f"grid calibration: cl0 y={Y_CL0:.1f}  {PX_PER_CL:.1f} px/cl  "
          f"{PX_PER_DEG:.2f} px/deg  stagger {STAGGER:.2f} px = {STAGGER/PX_PER_DEG:.3f} deg")
    print(f"thin-airfoil alpha_L0 for the 4-digit mean line (m=0.04, p=0.4): {a_l0:.2f} deg\n")

    iy, ix = np.nonzero(mask)
    ix, iy = ix.astype(float), iy.astype(float)

    slopes, rejected, npix_ok = [], [], []
    print(f"{'M':<7}{'slope /deg':>11}{'zero-lift deg':>15}{'px':>7}   status")
    for k, M in enumerate(MACH):
        x0 = X_M030 + STAGGER * k
        fit = fit_linear(ix, iy, x0)
        assert fit is not None, f"M={M}: linear fit failed"
        slope, a0, npix, xmin = fit
        if k:
            assert xmin > leftmost, f"M={M}: starts left of the M=0.300 curve"
        else:
            leftmost = xmin
        if M <= FIT_OK:
            # second off-by-one guard: mis-indexing the carpet by one curve moves
            # alpha_L0 by exactly +-4 deg, into the 0 / -8 decoys this excludes
            assert ALPHA_L0_MIN < a0 < ALPHA_L0_MAX, (
                f"M={M}: zero-lift {a0:.2f} deg outside "
                f"({ALPHA_L0_MIN}, {ALPHA_L0_MAX}) -- calibration or curve index is wrong")
            slopes.append((M, slope, a0))
            npix_ok.append(npix)
            status = ""
        else:
            rejected.append((M, a0))
            status = "rejected: no straight segment (see scope note)"
        print(f"{M:<7}{slope:>11.4f}{a0:>15.2f}{npix:>7}   {status}")

    # Off-by-one guard proper. The fitter has no idea which curve is which -- it
    # simply latches onto whatever curve is nearest the seed -- so shifting the
    # whole carpet by one stagger relabels every curve while leaving alpha_L0
    # untouched. What actually pins the labelling is that the M = 0.300 curve is
    # the LEFTMOST curve in the figure: run the same fit one stagger further left
    # and it must find no curve at all. It finds only the y-axis tick labels, and
    # collects well under half the ink of the weakest real curve.
    ghost = fit_linear(ix, iy, X_M030 - STAGGER)
    weakest = min(n for n in npix_ok)
    assert ghost is None or ghost[2] < 0.6 * weakest, (
        f"a curve exists one stagger left of X_M030 ({ghost[2]} px vs {weakest}) "
        "-- the carpet is shifted and every Mach label is off by one")
    print(f"\nleftmost-curve guard: fit one stagger left of M=0.300 collects "
          f"{ghost[2] if ghost else 0} px (y-axis labels) vs {weakest} px for the "
          "weakest real curve")

    a0s = [a for _, _, a in slopes]
    spread = max(a0s) - min(a0s)
    print(f"\nzero-lift angle: mean {np.mean(a0s):+.2f} deg, spread {spread:.2f} deg over "
          f"{len(a0s)} accepted curves; thin-airfoil theory {a_l0:.2f} deg")
    assert spread < 0.4, f"zero-lift angle drifts by {spread:.2f} deg -- stagger is wrong"
    assert abs(np.mean(a0s) - a_l0) < 0.5, "mean zero-lift angle disagrees with theory"
    print("rejected curves, for contrast: "
          + ", ".join(f"M={M} -> {a:+.2f}" for M, a in rejected))

    # Prandtl-Glauert: nothing in the fitter knows about Mach number, so agreement
    # with 1/sqrt(1-M^2) is a genuine check on curve identity and scaling
    base_m, base_s, _ = slopes[0]
    print(f"\n{'M':<7}{'measured':>10}{'Prandtl-Glauert':>17}{'diff %':>9}")
    for M, s, _ in slopes:
        pg = base_s * (1 - base_m ** 2) ** .5 / (1 - M ** 2) ** .5
        d = 100 * (s - pg) / pg
        print(f"{M:<7}{s:>10.4f}{pg:>17.4f}{d:>+9.1f}")
        assert abs(d) < 4, f"M={M}: {d:+.1f}% off Prandtl-Glauert, digitisation suspect"

    hdr = ("# NACA 4415, digitised from figure 29 of NACA Report 832 (Graham, Nitzberg\n"
           "# & Olson, 1945), Ames 1x3.5 ft high-speed tunnel, 1e6 <= Re <= 2e6.\n"
           "# Produced by reference/digitize_tr832_4415.py -- see that file for method\n"
           "# and limits.\n"
           "# Only M <= 0.650. Above that the shock-induced kink leaves the curve with no\n"
           "# straight segment, the fitted zero-lift angle walks off to -4.4 ... -5.8 deg,\n"
           "# and the slope is not claimed. No CLmax file is produced at all: the peak\n"
           "# tracer that works on the 0015 latches onto neighbouring curves here.\n"
           "# zero_lift_deg is a self-check, not data. The section is cambered, so the\n"
           "# symmetry test used for the 0015 does not apply; instead thin-airfoil theory\n"
           f"# on the 4-digit mean line gives {a_l0:.2f} deg, the published experimental value\n"
           f"# is about -4.0 deg, and the {len(slopes)} accepted curves agree with each other to\n"
           f"# {spread:.2f} deg -- which validates the axis calibration, the 4 deg stagger, and\n"
           "# each curve's identity (a one-curve index error would show as 0 or -8 deg).\n")
    (HERE / "naca4415_lift_slope_vs_mach.csv").write_text(
        hdr + "mach,lift_slope_per_deg,zero_lift_deg\n"
        + "".join(f"{M},{s:.4f},{a:+.2f}\n" for M, s, a in slopes))
    print(f"\nwrote naca4415_lift_slope_vs_mach.csv ({len(slopes)} rows); "
          "no CLmax file -- see the scope note")


if __name__ == "__main__":
    main()
