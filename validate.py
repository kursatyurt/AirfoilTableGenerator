#!/usr/bin/env python
"""Compare a computed polar against a digitised experimental reference.

    python validate.py runs/naca0015_M030/polar.csv reference/naca0015_cl_alpha_M030.csv 0.30

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


def check_cd(computed, mach, cl_ref=0.22):
    """CD against TR-832 figure 32, which reports cd only at |cl| = 0.22.

    The experiment gives one number per Mach, not a polar, so the computed
    drag is interpolated to the same cl rather than compared angle by angle.
    """
    rows = read_csv(computed, "cl", "cd")
    cd = interp(sorted(rows), cl_ref)
    # fig42_lo/hi are blank for the Mach numbers figure 42 does not cover
    num = lambda v: float(v) if v.strip() else None
    ref = {float(r["mach"]): (float(r["cd"]), num(r["fig42_lo"]), num(r["fig42_hi"]))
           for r in _rows(ROOT / "reference" / "naca0015_cd_vs_mach_tr832.csv")}
    exp = ref.get(round(mach, 3))
    print(f"\nCD at cl = {cl_ref}")
    if cd is None:
        print(f"  computed polar does not reach cl = {cl_ref}")
        return
    if exp is None:
        print(f"  computed {cd:.5f}; no TR-832 row for M = {mach}")
        return
    print(f"  computed {cd:.5f}   TR-832 fig 32 {exp[0]:.5f}   "
          f"{100 * (cd - exp[0]) / exp[0]:+.1f}%")
    if exp[1] is not None:
        print(f"  fig 42 envelope for the five sections at cl=0.20: {exp[1]:.4f}-{exp[2]:.4f}")
    print("  fully turbulent SU2 vs a partly laminar 1945 tunnel: computed should run high")


def check_cm(computed):
    """CM has no experimental reference -- see reference/README.md.

    Figure 37's baselines were resolved, but the digitised curves failed the
    cm(cl=0) = 0 symmetry test by about -0.04 per curve, which no common
    offset removes, so no moment data was recorded. Only symmetry is checkable.
    """
    pts = read_csv(computed, "aoa", "cm")
    cm0 = interp(pts, 0.0) if pts[0][0] != 0.0 else pts[0][1]
    print(f"\nCM: no experimental reference exists (TR-832 fig 37 rejected, see "
          f"reference/README.md)")
    print(f"  symmetry check, cm at alpha=0: {cm0:+.5f}  (must be 0)")
    assert abs(cm0) < 0.005, f"cm({0}) = {cm0}, symmetry broken -- mesh or markers"


def _rows(path):
    with open(path) as f:
        return list(csv.DictReader(l for l in f if not l.startswith("#")))


def main(computed, reference, mach=None):
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
    if mach is not None:
        check_cd(computed, mach)
    check_cm(computed)


def _selftest():
    pts = [(0.0, 0.0), (2.0, 0.2), (4.0, 0.4), (8.0, 0.8)]
    a, a0 = lift_slope(pts, lo=0, hi=8)
    assert abs(a - 0.1) < 1e-9 and abs(a0) < 1e-9, (a, a0)
    assert abs(interp(pts, 3.0) - 0.3) < 1e-9
    assert interp(pts, 99) is None
    shifted = [(x, y + 0.1) for x, y in pts]        # +0.1 offset -> -1 deg zero-lift
    a, a0 = lift_slope(shifted, lo=0, hi=8)
    assert abs(a0 + 1.0) < 1e-9, a0
    # the TR-832 drag row must survive being read back, since check_cd keys on Mach
    # every row must parse, not just M=0.30: fig42_lo/hi are blank for 0.550 and 0.625
    r = {float(x["mach"]): x for x in _rows(ROOT / "reference" / "naca0015_cd_vs_mach_tr832.csv")}
    assert abs(float(r[0.3]["cd"]) - 0.0078) < 1e-9, r[0.3]
    assert not r[0.55]["fig42_lo"].strip(), "expected a blank fig42 column to exercise"
    import tempfile
    p = Path(tempfile.mkdtemp()) / "polar.csv"
    p.write_text("aoa,cl,cd,cm,converged\n0,0.0,0.010,0.0,1\n4,0.44,0.012,-0.004,1\n")
    check_cm(p)                                    # symmetric case must pass
    p.write_text("aoa,cl,cd,cm,converged\n0,0.0,0.010,-0.05,1\n4,0.44,0.012,-0.06,1\n")
    try:
        check_cm(p)
    except AssertionError:
        pass
    else:
        raise AssertionError("check_cm accepted a non-zero cm at alpha=0")
    print("validate selftest ok")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        _selftest()
    else:
        main(sys.argv[1], sys.argv[2],
             float(sys.argv[3]) if len(sys.argv) > 3 else None)
