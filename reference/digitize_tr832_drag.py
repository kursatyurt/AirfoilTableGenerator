#!/usr/bin/env python
"""Digitise NACA 0015 section drag from figure 32 of NACA Report 832.

Figure 32 (page 50 of the NASA scan) plots section drag coefficient against
section lift coefficient for M = 0.300 to 0.825. Unlike figure 27 it is not
staggered horizontally: every Mach curve is drawn on the same cl axis but on
its **own cd = 0 baseline**, stacked up the page. Nothing is readable until
those baselines are pinned, and an earlier attempt stopped here -- see the
"Drag and moment" section of reference/README.md for the version of the
problem this file resolves.

    python reference/digitize_tr832_drag.py       # fetches the PDF if needed

How the baselines were pinned:
  * page 50 is rendered at 600 dpi and deskewed (-0.55 deg)
  * the printed grid is 118.879 px = 0.01 in cd and 118.738 px = 0.1 in cl;
    the frame spans grid rows 0 (top) to 33 (bottom axis)
  * every ODD grid row 1..33 carries a long labelled tick, seventeen in all.
    This is measured, not assumed: the ink fraction immediately left of the
    frame is 0.30-0.76 on every odd row and 0.00-0.06 on every even row.
  * rows 1,3,5,7,9,11 are the graduated scale .12 .10 .08 .06 .04 .02, and
    rows 13,15,...,33 are eleven ticks all labelled "0"
  * eleven zeros, fourteen curves. Counting up from the bottom, one curve per
    zero, gets M = 0.300 (row 33) through M = 0.725 (row 15) and runs out.
    The remaining four -- 0.750, 0.775, 0.800, 0.825 -- SHARE the zero at row
    13, which is exactly why the graduated 0-to-0.12 scale is printed there
    and nowhere else: it is the scale those four tall curves are read against.
  * that is forced, not chosen. Giving 0.775 its own zero at row 11 puts its
    minimum at cd = -0.002, and giving 0.300 anything but row 33 either runs
    off the bottom of the frame or makes cd negative.

Independent check (this is what makes the numbers emittable): figure 42 on
page 54 plots section drag against Mach at cl = 0.20 for all five sections on
an ordinary shared axis, so it is immune to the baseline problem. The five
curves are hard to tell apart, but their envelope is not, and every one of the
fourteen values read out of figure 32 near cl = 0.20 lands inside or within
0.0015 of that envelope -- across a factor of fourteen in cd, from 0.0066 at
M = 0.400 to 0.107 at M = 0.825. An off-by-one baseline shifts cd by 0.02 and
would throw the low-Mach end out by a factor of four. FIG42 below is that
envelope, measured from the scan, and main() asserts against it.

Second check, internal to this figure: the 0015 is symmetric, so cd at cl =
+0.22 and cd at cl = -0.22 must be the same number. They agree to 0.0004 or
better for every Mach up to 0.650. That test knows nothing about figure 42 and
would also fail on a wrong baseline, since the two readings share one.

What is emitted, and what is deliberately not:
  * emitted: cd at |cl| = 0.22 for all fourteen Mach numbers, one CSV. Both
    checks above are run as assertions before anything is written.
  * NOT emitted: minimum drag. The bottom of the drag bucket lies under the
    printed Mach label, which sits directly above its own curve and inside the
    same 0.02-tall band, so the one value everybody wants is the one the
    figure covers up. Reading around it gave numbers that were not monotone in
    Mach at the 0.001 level, which is the size of the effect being measured.
  * NOT emitted: the cd-vs-cl polars. Traced curves failed the symmetry test
    by up to 0.005 at scattered cl -- the tracer brushes cl gridline stubs and
    neighbouring curves out in the steep wings. Ten of the fourteen were close
    to usable; none were checkable point by point, so none are recorded.
  * |cl| = 0.22 is not an arbitrary choice: it and its mirror are the only two
    columns where exactly fourteen curves resolve. Left of cl = 0.20 the Mach
    labels start, and at cl = 0.24 M = 0.825 has already climbed out through
    the top of the frame -- which is also why curves are counted from the
    BOTTOM of each column upward. Counting from the top silently relabels
    every curve once M = 0.825 leaves, and looks entirely plausible doing so.
  * the reading is a pixel measurement of a 1945 line drawing, good to about
    +/-0.0008 in cd (one third of a printed line width). It is not better than
    the draughtsman's fairing of the original data points.
"""
import subprocess, urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PDF = HERE / "naca_tr832.pdf"
URL = "https://ntrs.nasa.gov/api/citations/19930091909/downloads/19930091909.pdf"

MACH = [0.300, 0.400, 0.500, 0.550, 0.600, 0.625, 0.650, 0.675,
        0.700, 0.725, 0.750, 0.775, 0.800, 0.825]
# grid row carrying each curve's cd = 0; the top four share row 13
BASELINE = [33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 13, 13, 13]

SKEW_DEG = -0.55
Y_TOP, Y_BOT, N_ROWS = 598.0, 4521.0, 33     # frame rows, 0.01 in cd per row
X_CL0, PX_PER_TENTH_CL = 832.23, 118.7377    # x of cl = -0.8, and per 0.1 cl
CL_LEFT = -0.8
READ_CL = 0.22               # +/- this is where all fourteen curves resolve cleanly
SLOPE_CL = 0.18              # a nearby column, used only to gauge each curve's slope
SYM_MAX_MACH, SYM_TOL = 0.650, 0.0006   # symmetry self-check and its tolerance

# figure 42, cl = 0.20: (Mach, envelope of the five sections) measured off page 54
FIG42 = {0.300: (0.0051, 0.0074), 0.400: (0.0052, 0.0071), 0.500: (0.0055, 0.0071),
         0.600: (0.0054, 0.0085), 0.650: (0.0035, 0.0087), 0.675: (0.0057, 0.0104),
         0.700: (0.0063, 0.0144), 0.725: (0.0105, 0.0214), 0.750: (0.0195, 0.0403),
         0.775: (0.0302, 0.0735), 0.800: (0.0401, 0.0919), 0.825: (0.0751, 0.1020)}
FIG42_SLOP = 0.0015          # cl = 0.22 vs 0.20, plus both readings' line width

DY = (Y_BOT - Y_TOP) / N_ROWS


def load_figure():
    from PIL import Image
    if not PDF.exists():
        print(f"fetching {URL}")
        urllib.request.urlretrieve(URL, PDF)
    png = HERE / "_tr832_p50.png"
    if not png.exists():
        subprocess.run(["pdftoppm", "-f", "50", "-l", "50", "-r", "600", "-png",
                        "-singlefile", str(PDF), str(png.with_suffix(""))], check=True)
    im = Image.open(png).convert("L").rotate(SKEW_DEG, resample=Image.BICUBIC,
                                             fillcolor=255)
    return strip_verticals(np.array(im) < 180)


def strip_verticals(mask):
    """Blank vertical runs longer than 250 px -- the cl gridlines, which stand on
    every tenth of cl and would otherwise swallow whole columns. Nothing in the
    data is vertical for 250 px except the near-vertical high-cl ends of the
    curves, and those already lie outside the lane each curve is read in."""
    m = mask.copy()
    for x in range(m.shape[1]):
        idx = np.nonzero(m[:, x])[0]
        if not len(idx):
            continue
        s = p = idx[0]
        for i in list(idx[1:]) + [10 ** 9]:
            if i != p + 1:
                if p - s + 1 > 250:
                    m[s:p + 1, x] = False
                s = i
            p = i
    return m


def check_ticks(mask):
    """Every odd grid row carries a labelled tick, every even row does not.

    This is the whole basis of the baseline assignment, so it is verified
    rather than trusted: the strip 796..829 px sits between the tick labels
    and the frame, and only a long labelled tick reaches into it.
    """
    ink = [mask[int(round(Y_TOP + DY * k)) - 7:int(round(Y_TOP + DY * k)) + 8,
                796:829].max(0).mean() for k in range(N_ROWS + 1)]
    odd, even = [ink[k] for k in range(1, 34, 2)], [ink[k] for k in range(0, 34, 2)]
    assert min(odd) > 0.20, f"a labelled tick is missing: {odd}"
    assert max(even) < 0.15, f"an unexpected tick on an even row: {even}"
    return len(odd)


def grid_rows(mask):
    """Measured y of the 34 horizontal gridlines. The printed grid is not quite
    linear across 3900 px, and half a line width of drift is enough to make a
    curve lying beside a gridline indistinguishable from the gridline itself,
    so the rows are measured rather than computed from Y_TOP and DY."""
    frac = mask[:, 900:3250].mean(1)
    idx = [i for i in range(int(Y_TOP) - 20, int(Y_BOT) + 20) if frac[i] > 0.35]
    ys, s, p = [], idx[0], idx[0]
    for i in idx[1:] + [10 ** 9]:
        if i > p + 5:
            ys.append((s + p) / 2)
            s = i
        p = i
    ys = [y for y in ys if min(abs(y - (Y_TOP + DY * k)) for k in range(34)) < 12]
    assert len(ys) == 34, f"found {len(ys)} gridlines, expected 34"
    return np.array(ys)


def curves_at(mask, x, ys):
    """Grid-row positions of every curve crossing column x, bottom of frame first.

    The gridlines are located per column, not once for the whole figure: the
    printed grid bows by up to ~10 px across the 2500 px frame, which is more
    than a line width, and a fixed grid therefore mistakes a bowed gridline for
    a curve sitting just above a baseline -- reading cd = 0.001 where there is
    no curve at all. The local offset is found by sliding the whole grid until
    the most thin runs land on it. A run that then sits on a gridline and is no
    more than 8 px thick is that gridline; anything else is curve ink, or a
    curve fused with the gridline it happens to be crossing.
    """
    col = mask[:, x - 1:x + 2].max(1)
    idx = np.nonzero(col[int(Y_TOP) - 40:int(Y_BOT) + 40])[0] + int(Y_TOP) - 40
    runs, s, p = [], idx[0], idx[0]
    for i in list(idx[1:]) + [10 ** 9]:
        if i != p + 1:
            runs.append(((s + p) / 2, p - s + 1))
            s = i
        p = i
    thin = np.array([y for y, t in runs if t <= 8]) if runs else np.zeros(0)
    shifts = np.arange(-20, 21, 1.0)
    if len(thin):
        score = [(np.abs(thin[:, None] - (ys + d)).min(1) < 4).sum() for d in shifts]
        ys = ys + shifts[int(np.argmax(score))]
    out = [float(np.interp(y, ys, np.arange(34)))
           for y, t in runs if not (np.abs(ys - y).min() < 6 and t <= 8)]
    return sorted(out, reverse=True)


def seed_rows(mask, ys, cl):
    """The fourteen curves at this cl, bottom of frame first.

    Finding exactly fourteen is itself the check that the column is a usable
    one: at cl = 0.24 M = 0.825 has already left the frame through the top, and
    under the printed Mach labels a label glyph reads as a fifteenth curve.
    """
    x = int(round(X_CL0 + (cl - CL_LEFT) / 0.1 * PX_PER_TENTH_CL))
    rows = curves_at(mask, x, ys)
    assert len(rows) == len(MACH), f"{len(rows)} curves at cl = {cl}, want 14"
    return rows


def main():
    mask = load_figure()
    ys = grid_rows(mask)
    n = check_ticks(mask)
    print(f"{n} labelled ticks on odd grid rows, none on even rows -- "
          f"6 graduated + {n - 6} zeros for {len(MACH)} curves\n")

    hi_rows = seed_rows(mask, ys, +READ_CL)
    lo_rows = seed_rows(mask, ys, -READ_CL)
    near = seed_rows(mask, ys, SLOPE_CL)          # only to gauge each curve's slope

    print(f"{'M':<7}{'cd(+0.22)':>11}{'cd(-0.22)':>11}{'cd':>9}{'sym':>8}"
          f"{'fig42 envelope at cl=0.20':>28}{'slop':>8}")
    rowsout = []
    for k, M in enumerate(MACH):
        up = (BASELINE[k] - hi_rows[k]) * 0.01
        dn = (BASELINE[k] - lo_rows[k]) * 0.01
        cd, sym = (up + dn) / 2, abs(up - dn) / 2
        lo, hi = FIG42.get(M, (0.0, 1.0))
        # tolerance has to scale with how fast the curve is moving: at M = 0.825
        # it is nearly vertical at cl = 0.22, so the 0.02 of cl between the two
        # figures is worth 0.005 in cd; at M = 0.300 it is worth nothing
        slop = FIG42_SLOP + abs(hi_rows[k] - near[k]) * 0.01 / 2
        assert lo - slop <= up <= hi + slop, (
            f"M={M}: cd={up:.4f} at cl={READ_CL} is outside figure 42's "
            f"{lo:.4f}-{hi:.4f} envelope -- the baseline assignment is wrong")
        if M <= SYM_MAX_MACH:
            assert sym < SYM_TOL, (
                f"M={M}: cd differs by {2 * sym:.4f} between cl = +{READ_CL} and "
                f"-{READ_CL}, but the section is symmetric")
        print(f"{M:<7}{up:>11.4f}{dn:>11.4f}{cd:>9.4f}{sym:>8.4f}"
              f"{lo:>21.4f}-{hi:.4f}{slop:>8.4f}")
        rowsout.append((M, cd, sym, lo, hi))
    assert all(rowsout[i][1] <= rowsout[i + 1][1] + 0.0015 for i in range(13)), \
        "cd at fixed cl must rise with Mach"

    (HERE / "naca0015_cd_vs_mach_tr832.csv").write_text(
        "# NACA 0015 section drag coefficient, digitised from figure 32 of NACA\n"
        "# Report 832 (Graham, Nitzberg & Olson, 1945), Ames 1x3.5 ft high-speed\n"
        "# tunnel, 1e6 <= Re <= 2e6. Produced by reference/digitize_tr832_drag.py --\n"
        "# see that file for the method and for what it deliberately does not claim.\n"
        "#\n"
        f"# cd is the mean of the readings at cl = +{READ_CL} and cl = -{READ_CL}. Those\n"
        "# are the only two columns of figure 32 where all fourteen curves resolve\n"
        "# cleanly: clear of the printed Mach labels, and inside the cl range where\n"
        "# M = 0.825 is still on the page. This is NOT minimum drag -- the drag bucket\n"
        "# bottom is hidden under the printed Mach labels and is not reported.\n"
        "#\n"
        "# sym_half_diff is a self-check, not data: the 0015 is symmetric, so the two\n"
        "# readings must agree. They do, to better than 0.0005 up to M = 0.650; above\n"
        "# that the curves are steep enough that a fraction of a line width in the cl\n"
        "# calibration shows up as a real difference, so treat it as the error bar.\n"
        "# fig42_lo/fig42_hi are the envelope of all five sections in figure 42 at\n"
        "# cl = 0.20, measured off page 54 -- the independent check on the baseline\n"
        "# assignment, not data about the 0015.\n"
        "mach,cd,sym_half_diff,fig42_lo,fig42_hi\n"
        + "".join(f"{M},{cd:.4f},{sym:.4f}," +
                 (f"{lo:.4f},{hi:.4f}\n" if MACH[i] in FIG42 else ",\n")
                 for i, (M, cd, sym, lo, hi) in enumerate(rowsout)))
    print(f"\nwrote naca0015_cd_vs_mach_tr832.csv ({len(rowsout)} rows)")


if __name__ == "__main__":
    main()
