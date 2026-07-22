#!/usr/bin/env python
"""Digitise figures 40, 41 and 42 of NACA Report 832 (Graham, Nitzberg & Olson,
1945) -- lift-curve slope, angle of zero lift and section drag coefficient
against Mach number.

    .venv/bin/python reference/digitize_tr832_fig40_42.py

WHAT IS ACTUALLY ON THESE FIGURES
---------------------------------
Not the 0006/0009/0012/0015/0018 thickness family.  Report 832's "five
representative sections" are five different cambers all at 15 per cent
thickness: NACA 65_2-215 (a=0.5), 66,2-215 (a=0.6), 0015, 23015 and 4415.  Only
the 0015 is symmetric.  All three figures are on printed page 49 = PDF page 54,
and each plots all five sections on one set of ordinary axes as five dash
patterns:

    solid                65_2-215 (a=0.5)
    short dashes         66,2-215 (a=0.6)
    long dash + short    0015
    long dash + 2 short  23015
    long dashes          4415

Figures 40 and 42 are taken at c_l = 0.20; figure 41 is the angle of zero lift.

METHOD
------
  * PDF page 54 is rendered at 400 dpi with pdftoppm.  Unlike figure 27 this
    page needs no deskew -- check_grid() confirms the ruling is uniform to about
    1.4 px rms over 1400 px, so any skew is under 0.1 deg.
  * Each figure is cropped to a fixed box and calibrated from its ruling: a
    uniform square grid of CELL = 87.3 px.  FIGS gives, per figure, the pixel
    row and column of a known value and the data units per grid cell.  The
    pitch and the straightness are asserted against the ink at run time.
  * The ruling is removed by position rather than by run length.  It is printed
    much paler than the curves and wanders a pixel or two, so no row holds one
    long unbroken run and length-based stripping does not find it.  Instead ink
    within GRID_HALF of a rule row is deleted where it is thin (GRID_THICK), and
    rule columns are declared unreadable outright.  The thinness test is what
    preserves a curve drawn along a rule -- the 4415 of figure 41 rides the -4
    deg line for most of its length.
  * Curves are then followed column by column by a slope-predicting tracer with
    a gate that widens while it coasts a dash gap and stops rather than
    extrapolate once it has coasted too long.

WHAT IS EMITTED, AND WHAT IS NOT
--------------------------------
Telling the five dash patterns apart is the whole difficulty, and this script
refuses to guess.  A named curve is emitted only where its identity rests on
evidence independent of the dash pattern:

  fig 41  NACA 0015, M = 0.30..0.83.  The only symmetric section of the five and
          therefore the only curve that can sit at zero.  It is the sole ink
          within 0.6 deg of the zero rule at every column of a whole grid cell
          around the seed, which is asserted.
  fig 41  NACA 4415, M = 0.31..0.75.  Alone in the entire lower half of the
          plot, near its published low-speed zero-lift angle of -4.0 deg.  Its
          trace stops in the middle of the post-divergence plunge, where it
          crosses the crowd of rules below -4 deg.
  fig 41  65_2-215, 66,2-215 and 23015 all lie inside a 0.5 deg band around
          -1.3 deg.  Their mutual order cannot be recovered.  NOT emitted.
  fig 40  the five form an upper pair and a lower trio; the members of each
          group merge to within a line width below M ~ 0.66 and cross one
          another above M ~ 0.75.  No individual curve could be pinned by
          anything independent, so no per-airfoil slope is emitted.  What is
          emitted is the bundle envelope -- the highest and lowest of the five
          at each Mach -- over M = 0.305..0.795, which is well defined and, as
          the checks below show, directly testable.
  fig 42  all five lie inside one band about 0.002 wide up to M ~ 0.6 and then
          fan out and cross.  Envelope again, and only over M = 0.305..0.570:
          past that the two edges merge into a single drawn line and the tracker
          cannot separate them again without guessing, so the post-divergence
          drag rise -- which is plainly there on the figure -- is NOT emitted.

SELF-CHECKS (all asserts; the numbers are printed on every run)
  * grid pitch and straightness, all three figures        -> calibration
  * fig 41 NACA 0015 zero-lift angle is 0 within 0.3 deg.  Measured +0.193 deg.
    Symmetry demands zero, so this tests the y calibration and the identity at
    once -- the same test that validated digitize_tr832.py.
  * fig 41 NACA 4415 zero-lift angle is -4.0 within 0.3 deg.  Measured -4.046.
  * fig 40 the slope rises with Mach and tracks Prandtl-Glauert: the M = 0.30
    value predicts the M = 0.60 value to within 15 per cent.  Measured +2.6 %.
  * fig 40 the peak of the envelope (force divergence) falls between M = 0.64
    and 0.76.  Measured M = 0.680.  The lower edge then collapses to 69 per cent
    of its peak by M = 0.80.
  * fig 40 the envelope brackets reference/naca0015_lift_slope_vs_mach.csv,
    which was extracted independently from figure 27.  All 7 of its points below
    M = 0.68 fall inside, with zero slack needed.  Above divergence the two are
    not comparable -- figure 27's slope is fitted over a wide alpha range while
    figure 40's is the local slope at c_l = 0.20 -- so the check stops there.
  * fig 42 both edges of the band lie in 0.004..0.012 below M = 0.60 (measured
    0.0048..0.0086) and the band is flat to better than a factor 1.6.
"""
import subprocess, urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PDF = HERE / "naca_tr832.pdf"
URL = "https://ntrs.nasa.gov/api/citations/19930091909/downloads/19930091909.pdf"
PAGE, DPI = 54, 400

# --- calibration ----------------------------------------------------------
# crop: (left, top, right, bottom) in 400 dpi page pixels
# x0px/x0val/y0px/y0val: pixel of a known gridline and its data value
# per_cell: data units per grid cell, x then y
FIGS = {
    40: dict(crop=(350, 350, 1800, 2000), x0px=176.5, x0val=0.30,
             y0px=1512.5, y0val=0.0, per_cell=(0.05, 0.02),
             band=(820, 1400), seed=(1035, 1090)),
    41: dict(crop=(2075, 350, 3525, 2000), x0px=109.5, x0val=0.30,
             y0px=897.0, y0val=0.0, per_cell=(0.05, 1.0),
             band=(700, 1400)),
    42: dict(crop=(1200, 2450, 2650, 4075), x0px=160.5, x0val=0.30,
             y0px=1503.5, y0val=0.0, per_cell=(0.05, 0.01),
             band=(600, 1500), seed=(1428, 1468)),
}
CELL = 87.3          # nominal grid pitch at 400 dpi, both axes, all three figs
GRID_TOL = 3.0       # px, allowed rms departure of detected grid from uniform
INK = 170            # grey level counted as curve ink (the curves are black)
GRID_INK = 215       # grey level that also catches the pale ruled grid
GRID_HALF = 5        # px, half-width of the band deleted around a lattice line
GRID_THICK = 3       # px, ink this thin across a lattice line is ruling, not data
TRACE_TOL = 7.0      # px, prediction window half-width
MAXMISS = 40         # columns a tracer may coast across a dash gap
ENV_WIN = 12         # px, how far each envelope edge may move in one column
ENV_GROW = 0.50      # px per missed column the gate opens outwards
ENV_GROW_IN = 0.10   # ... and inwards, much slower, or a neighbouring curve
                     #     captures the edge and silently narrows the envelope
ENV_REOPEN = 12      # px an edge may jump back out to when the bundle,
                     #     having merged into one line, fans out again
ENV_MISS = 45        # columns an envelope edge may lose its curve before the walk stops
BIN = 0.005          # Mach, width of the bin the envelope is averaged into
SEED_M = 0.45        # Mach at which the fig 41 tracers are seeded
CROSS_TOL = 0.004    # per deg, slack in the fig 40 vs figure-27 cross-check
CROSS_MAX = 0.68     # Mach, above which that cross-check is not meaningful


# --- image ----------------------------------------------------------------
def page():
    from PIL import Image
    if not PDF.exists():
        print(f"fetching {URL}")
        urllib.request.urlretrieve(URL, PDF)
    png = HERE / f"_tr832_p{PAGE}.png"
    if not png.exists():
        subprocess.run(["pdftoppm", "-f", str(PAGE), "-l", str(PAGE), "-r", str(DPI),
                        "-png", "-singlefile", str(PDF), str(png.with_suffix(""))],
                       check=True)
    return Image.open(png).convert("L")


def _runlen(mask, axis):
    """Length of the run of True containing each pixel, measured along `axis`."""
    from scipy import ndimage
    st = np.zeros((3, 3), bool)
    st[1, 1] = True
    st[0, 1] = st[2, 1] = axis == 0
    st[1, 0] = st[1, 2] = axis == 1
    lab, n = ndimage.label(mask, structure=st)
    size = np.bincount(lab.ravel())
    size[0] = 0
    return size[lab]


def _band(n, anchor, half):
    """Boolean mask of the pixels within `half` of anchor + k*CELL."""
    b = np.zeros(n, bool)
    for k in range(int(np.floor(-anchor / CELL)) - 1, int((n - anchor) / CELL) + 2):
        c = int(round(anchor + k * CELL))
        if -half <= c < n + half:                 # a negative stop would slice
            b[max(0, c - half):c + half + 1] = True  # from the far end instead
    return b


def figure(im, spec):
    """Return (loose mask for grid checking, curve mask, columns to ignore).

    Run-length stripping does not work on this scan: the ruling is pale (only
    visible at the loose GRID_INK threshold) and wanders a pixel or two, so no
    single row holds one long unbroken run.  The ruling is instead removed by
    position, using the very calibration lattice declared in FIGS -- which
    check_grid() has already confirmed against the ink.

      * horizontal rules: ink within GRID_HALF of a rule row is deleted only
        where it is thin (a run of at most GRID_THICK px across the rule).  The
        thinness test is what saves a curve that rides along a rule -- NACA 4415
        in figure 41 is drawn on top of the -4 deg line for most of its length,
        and unconditional blanking would erase it.
      * vertical rules: a curve crossing one is nearly parallel to nothing, so
        there is no thinness test that separates them.  Those columns are simply
        declared unreadable and skipped; they are 11 px in every 87, and both
        the tracer (which coasts) and the envelope (which bins) step over them.
    """
    g = np.array(im.crop(spec["crop"]))
    loose = g < GRID_INK
    strict = g < INK
    rows = _band(strict.shape[0], spec["y0px"], GRID_HALF)
    cols = _band(strict.shape[1], spec["x0px"], GRID_HALF)
    curve = strict & ~(rows[:, None] & (_runlen(strict, 0) <= GRID_THICK))
    curve[:, cols] = False
    return loose, curve, cols


def gridlines(raw, axis):
    """Candidate straight-line positions along `axis` (0 = rows).  These are not
    all grid: a dead-flat data curve, such as NACA 4415 in figure 41, shows up
    here too.  lattice() is what separates the two."""
    s = raw.sum(1 - axis).astype(float)
    hit = np.nonzero(s > 0.40 * s.max())[0]
    out, g = [], [hit[0]]
    for i in hit[1:]:
        if i - g[-1] <= 4:
            g.append(i)
        else:
            out.append(sum(g) / len(g))
            g = [i]
    out.append(sum(g) / len(g))
    return np.array(out)


def lattice(raw, axis):
    """Fit the uniform ruling to the candidates and return (positions, pitch,
    rms).  The anchor is chosen by consensus, so candidates that are really data
    curves (fig 41's 4415 rides -4 deg, fig 42's drag bundle is flat) fall off
    the lattice and are excluded instead of being blanked as grid."""
    p = gridlines(raw, axis)
    best = max(p, key=lambda a: np.sum(np.abs((p - a) - np.round((p - a) / CELL) * CELL) < 4))
    k = np.round((p - best) / CELL)
    keep = np.abs((p - best) - k * CELL) < 4
    pk, kk = p[keep], k[keep]
    pitch, off = np.polyfit(kk, pk, 1)
    rms = float(np.sqrt(np.mean((pk - (pitch * kk + off)) ** 2)))
    n = raw.shape[axis]
    ks = np.arange(np.floor(-off / pitch) - 1, np.ceil((n - off) / pitch) + 1)
    return pitch * ks + off, float(pitch), rms, int(keep.sum())


def check_grid(raw, name):
    """Assert the ruling really is a uniform CELL-pitch grid: wrong pitch means
    the wrong render dpi, a large rms means the page is skewed and the whole
    calibration below would be wrong."""
    worst = 0.0
    for axis in (0, 1):
        _, pitch, rms, n = lattice(raw, axis)
        print(f"  grid {name} axis{axis}: pitch {pitch:.2f} px, {n} lines, rms {rms:.2f} px")
        assert abs(pitch - CELL) < 1.5, f"{name}: grid pitch {pitch:.2f} != {CELL}"
        worst = max(worst, rms)
    assert worst < GRID_TOL, f"{name}: grid rms {worst:.2f} px -- page is skewed"
    return worst


# --- tracing --------------------------------------------------------------
def runs(mask, x, lo, hi):
    col = mask[lo:hi, x]
    idx = np.nonzero(col)[0]
    if not len(idx):
        return []
    out, s, p = [], idx[0], idx[0]
    for i in list(idx[1:]) + [10 ** 9]:
        if i != p + 1:
            out.append(lo + (s + p) / 2)
            s = i
        p = i
    return out


def trace(mask, x0, y0, band, step, xlim):
    """Follow one curve from (x0, y0).  Returns {x: y} for columns where exactly
    one ink run sits in the prediction window, and the set of ambiguous columns."""
    y, dy, pts, amb, miss = float(y0), 0.0, {}, set(), 0
    x = int(x0)
    while xlim[0] <= x <= xlim[1] and miss < MAXMISS:
        x += step
        pred = y + dy * step
        if not (band[0] < pred < band[1]):
            break
        c = [v for v in runs(mask, x, band[0], band[1]) if abs(v - pred) < TRACE_TOL]
        if len(c) == 1:
            dy = 0.7 * dy + 0.3 * (c[0] - y) / step
            y, miss, pts[x] = c[0], 0, c[0]
        else:
            if len(c) > 1:
                amb.add(x)
            y, miss = pred, miss + 1
    return pts, amb


def envelope(mask, spec, xs):
    """Top and bottom of the five-curve bundle at each column.

    Robust exactly where the individual curves are not: it needs only that some
    curve is inked, never which one.  The accepted rows are windowed on the
    previous column's answer (ENV_WIN) so that ruling left over beyond the plot
    frame cannot be mistaken for data; the walk is seeded from `seed`, the band
    the five curves occupy at the left-hand end of the plot."""
    # the two edges are tracked independently -- tying them together lets the
    # envelope collapse onto whichever curve happens to be inked in a column
    # where all the others are in a dash gap, and it can never reopen.  Each edge
    # carries a slope so that the window still finds it on the near-vertical
    # drag rise past divergence.
    y = list(spec["seed"])
    dy, miss, out = [0.0, 0.0], [0, 0], {}
    for x in xs:
        r = runs(mask, x, *spec["band"])
        for i, pick in enumerate((min, max)):
            # coast with a widening gate: while the edge is in a dash gap its
            # position is held and the search window grows, so it re-acquires the
            # same curve instead of drifting off it
            n = miss[i] + 1
            pred = y[i] + dy[i] * n
            # ... and the gate only widens outwards.  Widening it inwards lets a
            # neighbouring curve capture the edge whenever the outermost one is
            # in a long dash gap, which silently narrows the envelope.
            out_, in_ = ENV_GROW * miss[i], ENV_GROW_IN * miss[i]
            up, dn = (out_, in_) if i == 0 else (in_, out_)
            near = [v for v in r if -(ENV_WIN + up) <= v - pred <= ENV_WIN + dn]
            if near:
                v = pick(near)
                dy[i], y[i], miss[i] = 0.6 * dy[i] + 0.4 * (v - y[i]) / n, v, 0
            else:
                miss[i] += 1
        # reopen: where the five curves merge (fig 42 is one line around M=0.6)
        # both edges land on the same ink, and without this they could never
        # separate again as the bundle fans out past divergence
        ext = [v for v in r if y[0] - ENV_REOPEN <= v <= y[1] + ENV_REOPEN]
        if ext and not max(miss):
            y[0], y[1] = min(ext), max(ext)
        # once an edge has lost its curve for longer than a dash gap the walk is
        # extrapolating, not reading: stop rather than invent numbers
        if max(miss) > ENV_MISS:
            break
        if r:
            out[x] = (y[0], y[1])
    return out


# --- units ----------------------------------------------------------------
def to_x(spec, px):
    return spec["x0val"] + (px - spec["x0px"]) / CELL * spec["per_cell"][0]


def to_y(spec, py):
    return spec["y0val"] - (py - spec["y0px"]) / CELL * spec["per_cell"][1]


def from_x(spec, v):
    return spec["x0px"] + (v - spec["x0val"]) / spec["per_cell"][0] * CELL


HDR = ("# NACA Report 832 (Graham, Nitzberg & Olson, 1945), Ames 1x3.5 ft\n"
       "# high-speed tunnel, 1e6 <= Re <= 2e6.  Figure {n} on printed page 49.\n"
       "# The five sections of Report 832 are NACA 65_2-215 (a=0.5), 66,2-215\n"
       "# (a=0.6), 0015, 23015 and 4415 -- five cambers at 15% thickness, NOT a\n"
       "# 0006..0018 thickness family.\n"
       "# Produced by reference/digitize_tr832_fig40_42.py -- see that file for\n"
       "# the method, the self-checks and the limits.\n")


def main():
    im = page()
    out = {}

    # ---------------- figure 41: the one figure with identifiable curves ----
    spec = FIGS[41]
    raw41, m41, _ = figure(im, spec)
    print("figure 41 -- angle of zero lift")
    check_grid(raw41, "fig41")
    xseed = int(from_x(spec, SEED_M))

    curves41, level = {}, {}
    for name, target in (("NACA 0015", 0.0), ("NACA 4415", -4.0)):
        want = spec["y0px"] - target / spec["per_cell"][1] * CELL
        # count what sits near the target over a whole grid cell of columns: the
        # identity claim is that exactly one curve does, and the per-column count
        # is 0 (a dash gap or a blanked rule column) or 1, never more
        rule = _band(m41.shape[0], spec["y0px"], 2)     # rows the ruling occupies

        def cand(x):
            """Ink within 0.6 deg of the target.  A heavy rule that survived the
            grid removal (the zero line of fig 41 is drawn bolder than the rest)
            is dropped -- but only if something else is left, because the 4415 is
            itself drawn along the -4 deg rule."""
            c = [v for v in runs(m41, x, *spec["band"]) if abs(v - want) < 0.6 * CELL]
            off = [v for v in c if not rule[int(round(v))]]
            return off or c

        counts = [len(cand(x)) for x in range(xseed, xseed + int(CELL))]
        assert max(counts) == 1, (f"fig41 {name}: up to {max(counts)} curves within 0.6 deg "
                                  f"of {target} deg near M={SEED_M} -- identity not unique")
        x0 = xseed + counts.index(1)
        y0 = cand(x0)[0]
        xlim = (int(from_x(spec, 0.30)), int(from_x(spec, 0.90)))
        fwd, a1 = trace(m41, x0, y0, spec["band"], +1, xlim)
        bwd, a2 = trace(m41, x0, y0, spec["band"], -1, xlim)
        pts = {**bwd, **fwd, x0: y0}
        curves41[name] = pts
        flat = [to_y(spec, y) for x, y in pts.items() if to_x(spec, x) < 0.60]
        level[name] = float(np.mean(flat))
        print(f"  {name}: {len(pts)} columns, M {to_x(spec, min(pts)):.3f}"
              f"..{to_x(spec, max(pts)):.3f}, level below M=0.6 {level[name]:+.3f} deg,"
              f" {len(a1 | a2)} ambiguous columns")

    # the two physics checks that fix both the calibration and the identities
    assert abs(level["NACA 0015"]) < 0.3, (
        f"fig41: NACA 0015 zero-lift {level['NACA 0015']:+.2f} deg; the section is "
        "symmetric so this must be 0 -- calibration or curve identity is wrong")
    assert abs(level["NACA 4415"] + 4.0) < 0.3, (
        f"fig41: NACA 4415 zero-lift {level['NACA 4415']:+.2f} deg, published -4.0 "
        "-- curve misidentified")

    binned41 = {}
    for n, pts in curves41.items():
        for x, y in pts.items():
            M = to_x(spec, x)
            if 0.295 <= M <= 0.90:
                binned41.setdefault((round(M / BIN) * BIN, n), []).append(to_y(spec, y))
    rows41 = sorted((M, n, float(np.mean(v)), len(v)) for (M, n), v in binned41.items())
    out[41] = (HDR.format(n=41) +
               "# Only the NACA 0015 and the NACA 4415 are given.  The 0015 is the only\n"
               "# symmetric section of the five and so the only curve that can sit at\n"
               "# zero; the 4415 is alone in the whole lower half of the plot near -4\n"
               "# deg.  Both identities are fixed by physics, not by reading the dash\n"
               "# pattern off the scan.  The other three (65_2-215, 66,2-215, 23015) lie\n"
               "# inside a 0.5 deg band around -1.3 deg; their order is not recoverable\n"
               "# from this scan and they are deliberately omitted.\n"
               f"# Self-check: 0015 level below M=0.6 is {level['NACA 0015']:+.3f} deg\n"
               "#             (symmetry demands 0.000 -- this validates the y calibration)\n"
               f"#             4415 level below M=0.6 is {level['NACA 4415']:+.3f} deg\n"
               "#             (published low-speed value -4.0)\n"
               "# columns = readable image columns averaged into the Mach bin.\n"
               "mach,airfoil,angle_of_zero_lift_deg,columns\n" +
               "".join(f"{M:.3f},{n},{v:+.3f},{c}\n" for M, n, v, c in rows41))

    # ---------------- figures 40 and 42: bundle envelope only ---------------
    env = {}
    for n, title in ((40, "lift-curve slope at c_l=0.20"),
                     (42, "section drag coefficient at c_l=0.20")):
        spec = FIGS[n]
        raw, m, _ = figure(im, spec)
        print(f"figure {n} -- {title}")
        check_grid(raw, f"fig{n}")
        e = envelope(m, spec, range(int(from_x(spec, 0.30)), int(from_x(spec, 0.86))))
        binned = {}
        for x, (top, bot) in e.items():
            binned.setdefault(round(to_x(spec, x) / BIN) * BIN, []).append((top, bot))
        env[n] = sorted((M, to_y(spec, np.mean([t for t, _ in v])),
                         to_y(spec, np.mean([b for _, b in v])), len(v))
                        for M, v in binned.items() if len(v) >= 3 and 0.295 <= M <= 0.855)
        print(f"  envelope over {len(env[n])} Mach stations, "
              f"M {env[n][0][0]:.3f}..{env[n][-1][0]:.3f}")

    # --- figure 40 checks --------------------------------------------------
    M40 = np.array([r[0] for r in env[40]])
    hi40 = np.array([r[1] for r in env[40]])
    lo40 = np.array([r[2] for r in env[40]])
    peak = float(M40[int(np.argmax(hi40))])
    print(f"  fig40 upper edge {hi40[0]:.4f} at M={M40[0]:.2f}, peak {hi40.max():.4f} "
          f"at M={peak:.3f}, {hi40[-1]:.4f} at M={M40[-1]:.2f}")
    assert 0.64 < peak < 0.76, f"fig40: force divergence at M={peak:.2f}, expected ~0.70"
    # the lower edge is the one that shows the collapse inside the traced range
    fall = float(lo40[-1] / lo40.max())
    print(f"  fig40 lower edge falls to {100 * fall:.0f} % of its peak by M={M40[-1]:.2f}")
    assert fall < 0.80, "fig40: no post-divergence collapse"
    i3, i6 = int(np.argmin(abs(M40 - 0.30))), int(np.argmin(abs(M40 - 0.60)))
    pg = hi40[i3] * np.sqrt(1 - 0.30 ** 2) / np.sqrt(1 - 0.60 ** 2)
    print(f"  fig40 Prandtl-Glauert: {hi40[i3]:.4f} at M=0.30 predicts {pg:.4f} at "
          f"M=0.60, measured {hi40[i6]:.4f} ({100 * (hi40[i6] / pg - 1):+.1f} %)")
    assert hi40[i6] > hi40[i3], "fig40: slope does not rise with Mach"
    assert abs(hi40[i6] / pg - 1) < 0.15, "fig40: departs from Prandtl-Glauert"

    ref = np.array([[float(v) for v in ln.split(",")[:2]] for ln in
                    (HERE / "naca0015_lift_slope_vs_mach.csv").read_text().splitlines()
                    if ln and not ln.startswith(("#", "mach"))])
    # only below force divergence: figure 27's slope is fitted over a wide alpha
    # range, figure 40's is the local slope at c_l = 0.20, and past divergence the
    # curve is too bent for the two to mean the same thing
    ref = ref[(ref[:, 0] <= CROSS_MAX) & (ref[:, 0] >= M40[0])]
    worst, inside = 0.0, 0
    for M, s in ref:
        b, u = np.interp(M, M40, lo40), np.interp(M, M40, hi40)
        inside += b - CROSS_TOL <= s <= u + CROSS_TOL
        worst = max(worst, max(0.0, b - s, s - u))
    print(f"  fig40 NACA 0015 cross-check (M<={CROSS_MAX}): {inside}/{len(ref)} points of "
          f"naca0015_lift_slope_vs_mach.csv inside the envelope, worst miss {worst:.4f}")
    assert inside == len(ref), "fig40: envelope disagrees with the fig-27 0015 slopes"

    # --- figure 42 checks --------------------------------------------------
    M42 = np.array([r[0] for r in env[42]])
    hi42 = np.array([r[1] for r in env[42]])
    lo42 = np.array([r[2] for r in env[42]])
    sub = M42 < 0.60
    print(f"  fig42 band below M=0.60: {lo42[sub].min():.4f}..{hi42[sub].max():.4f}, "
          f"peak {hi42.max():.4f} at M={M42[int(np.argmax(hi42))]:.3f}")
    assert 0.004 <= hi42[sub].max() <= 0.012, \
        f"fig42: top of the band {hi42[sub].max():.4f} outside 0.004..0.012"
    assert 0.004 <= lo42[sub].min() <= 0.012, \
        f"fig42: bottom of the band {lo42[sub].min():.4f} outside 0.004..0.012"
    # below divergence the drag of all five is flat; that flatness is the check
    assert hi42[sub].max() / hi42[sub].min() < 1.6, "fig42: pre-divergence band not flat"

    for n, ylab in ((40, "lift_curve_slope_per_deg"), (42, "section_drag_coefficient")):
        out[n] = (HDR.format(n=n) +
                  "# NO per-airfoil curve is given for this figure.  All five dash\n"
                  "# patterns lie within a line width of one another below M ~ 0.66 and\n"
                  "# cross each other above M ~ 0.75, and nothing independent of the dash\n"
                  "# pattern distinguishes them, so naming any single trace would be a\n"
                  "# guess.  What is given instead is the envelope of the bundle -- the\n"
                  "# highest and the lowest of the five curves at each Mach number.  All\n"
                  "# five sections lie between the two value columns.\n"
                  "# columns = readable image columns averaged into the Mach bin.\n"
                  f"mach,{ylab}_max,{ylab}_min,columns\n" +
                  "".join(f"{M:.3f},{a:.5f},{b:.5f},{c}\n" for M, a, b, c in env[n]))

    for n, stem in ((40, "lift_slope"), (41, "zero_lift_angle"), (42, "drag")):
        f = HERE / f"tr832_fig{n}_{stem}_vs_mach.csv"
        f.write_text(out[n])
        print(f"wrote {f.name} "
              f"({sum(1 for l in out[n].splitlines() if not l.startswith('#')) - 1} rows)")


if __name__ == "__main__":
    main()
