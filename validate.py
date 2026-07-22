#!/usr/bin/env python
"""Compare a computed polar against a digitised experimental reference.

    python validate.py runs/naca0015_M030/polar.csv reference/naca0015_cl_alpha_M030.csv

Reports the lift-slope and zero-lift angle from both, the point-by-point CL
difference over the linear range, and CLmax with the angle it occurs at.
"""
import csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_csv(path, xcol, ycol):
    with open(path) as f:
        rows = list(csv.DictReader(l for l in f if not l.startswith("#")))
    return [(float(r[xcol]), float(r[ycol])) for r in rows]


def lift_slope(pts, lo=1.0, hi=8.0):
    """Least-squares fit of CL = a*alpha + b over the linear range."""
    s = [(x, y) for x, y in pts if lo <= x <= hi]
    if len(s) < 2:
        return None, None
    n = len(s)
    mx = sum(x for x, _ in s) / n
    my = sum(y for _, y in s) / n
    sxx = sum((x - mx) ** 2 for x, _ in s)
    sxy = sum((x - mx) * (y - my) for x, y in s)
    a = sxy / sxx
    return a, -(my - a * mx) / a          # slope per deg, zero-lift alpha


def interp(pts, x):
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0) if x1 != x0 else y0
    return None


def main(computed, reference):
    cfd = read_csv(computed, "aoa", "cl")
    exp = read_csv(reference, "alpha_deg", "cl")
    print(f"computed:  {computed}  ({len(cfd)} points)")
    print(f"reference: {reference} ({len(exp)} points)\n")

    for label, pts in (("computed", cfd), ("reference", exp)):
        a, a0 = lift_slope(pts)
        i = max(range(len(pts)), key=lambda k: pts[k][1])
        print(f"{label:>10}: slope {a:.4f}/deg ({a * 57.3:.2f}/rad), "
              f"zero-lift {a0:+.2f} deg, CLmax {pts[i][1]:.3f} at {pts[i][0]:.1f} deg")

    print(f"\n{'alpha':>7} {'CL cfd':>9} {'CL exp':>9} {'diff':>8} {'%':>7}")
    diffs = []
    for x, y in cfd:
        ye = interp(exp, x)
        if ye is None:
            continue
        d = y - ye
        diffs.append(abs(d))
        pct = 100 * d / ye if abs(ye) > 0.05 else float("nan")
        print(f"{x:7.1f} {y:9.4f} {ye:9.4f} {d:+8.4f} {pct:+7.1f}")
    if diffs:
        print(f"\nmean |dCL| = {sum(diffs) / len(diffs):.4f}, max |dCL| = {max(diffs):.4f}")
        print("reference is digitised from a plot: roughly +/-0.02 CL of its own uncertainty")


def _selftest():
    pts = [(0.0, 0.0), (2.0, 0.2), (4.0, 0.4), (8.0, 0.8)]
    a, a0 = lift_slope(pts, lo=0, hi=8)
    assert abs(a - 0.1) < 1e-9 and abs(a0) < 1e-9, (a, a0)
    assert abs(interp(pts, 3.0) - 0.3) < 1e-9
    assert interp(pts, 99) is None
    shifted = [(x, y + 0.1) for x, y in pts]        # +0.1 offset -> -1 deg zero-lift
    a, a0 = lift_slope(shifted, lo=0, hi=8)
    assert abs(a0 + 1.0) < 1e-9, a0
    print("validate selftest ok")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        _selftest()
    else:
        main(sys.argv[1], sys.argv[2])
