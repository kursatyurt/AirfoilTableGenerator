#!/usr/bin/env python
"""Verify the transcribed Datcom coordinate tables against the analytic NACA 00xx
thickness distribution. A mistyped or misread digit shows up immediately as an
outlier; genuine measured-model deviations stay small.

    python reference/check_reference.py
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent


def naca_thickness(x, t):
    """Standard NACA 4-digit half-thickness, open trailing edge."""
    return (t / 0.2) * (0.2969 * x**0.5 - 0.1260 * x - 0.3516 * x**2
                        + 0.2843 * x**3 - 0.1015 * x**4)


def load(name):
    with open(HERE / name) as f:
        rows = [r for r in csv.reader(l for l in f if not l.startswith("#"))]
    assert rows[0] == ["x_c", "y_c"], rows[0]
    return [(float(a), float(b)) for a, b in rows[1:]]


def check(name, t, tol):
    pts = load(name)
    worst = max(((abs(y - naca_thickness(x, t)), x, y) for x, y in pts))
    assert pts[0] == (0.0, 0.0), "leading edge must sit at the origin"
    assert pts[-1][0] == 1.0, "table must run to the trailing edge"
    assert all(b[0] > a[0] for a, b in zip(pts, pts[1:])), "x/c must increase"
    assert abs(max(y for _, y in pts) - t / 2) < 0.001, "max half-thickness != t/2"
    print(f"{name}: {len(pts)} points, worst deviation {worst[0]:.5f} at x/c={worst[1]}")
    assert worst[0] < tol, f"{name}: {worst} exceeds {tol} -- likely a transcription error"


def check_compressibility():
    """Lift slope must follow Prandtl-Glauert until shocks break it.

    Nothing in the digitiser knows about compressibility -- each Mach curve is
    fitted independently -- so agreement with 1/sqrt(1-M^2) is a genuine check
    that the curves were identified and scaled correctly.
    """
    with open(HERE / "naca0015_lift_slope_vs_mach.csv") as f:
        rows = [r for r in csv.DictReader(l for l in f if not l.startswith("#"))]
    data = [(float(r["mach"]), float(r["lift_slope_per_deg"]), float(r["zero_lift_deg"]))
            for r in rows]
    base_m, base_s, _ = data[0]
    print(f"\n{'M':<7}{'measured':>10}{'Prandtl-Glauert':>17}{'diff %':>9}")
    for m, s, a0 in data:
        assert abs(a0) < 0.5, f"M={m}: zero-lift {a0} deg, symmetry demands 0"
        pg = base_s * (1 - base_m**2) ** 0.5 / (1 - m**2) ** 0.5
        d = 100 * (s - pg) / pg
        print(f"{m:<7}{s:>10.4f}{pg:>17.4f}{d:>+9.1f}" + ("" if m <= 0.675 else "   (past divergence)"))
        if m <= 0.675:
            assert abs(d) < 5, f"M={m}: {d:+.1f}% off Prandtl-Glauert, digitisation suspect"
    # beyond force divergence the real airfoil must fall below the inviscid trend
    tail = [(m, s) for m, s, _ in data if m >= 0.70]
    for m, s in tail:
        pg = base_s * (1 - base_m**2) ** 0.5 / (1 - m**2) ** 0.5
        assert s < pg, f"M={m}: slope {s} exceeds Prandtl-Glauert past divergence"
    print("lift slope follows Prandtl-Glauert to M=0.675, then falls away as it must")


def check_cross_source():
    """The M=0.30 curve was digitised twice, from two different scans."""
    with open(HERE / "naca0015_cl_alpha_M030.csv") as f:
        datcom = [(float(r["alpha_deg"]), float(r["cl"]))
                  for r in csv.DictReader(l for l in f if not l.startswith("#"))]
    lin = [(a, c) for a, c in datcom if 2 <= a <= 8]
    n = len(lin)
    mx = sum(a for a, _ in lin) / n
    my = sum(c for _, c in lin) / n
    slope = (sum((a - mx) * (c - my) for a, c in lin)
             / sum((a - mx) ** 2 for a, _ in lin))
    with open(HERE / "naca0015_lift_slope_vs_mach.csv") as f:
        tr832 = float(next(r for r in csv.DictReader(l for l in f if not l.startswith("#"))
                           if float(r["mach"]) == 0.3)["lift_slope_per_deg"])
    d = 100 * (slope - tr832) / tr832
    print(f"\nM=0.30 lift slope: Datcom scan {slope:.4f}/deg, TR-832 scan {tr832:.4f}/deg "
          f"({d:+.1f}%)")
    assert abs(d) < 3, "independent digitisations of the same experiment disagree"


if __name__ == "__main__":
    # tolerance is loose enough for real measured ordinates, tight enough that a
    # wrong digit (>= 0.001 in the 4th place) fails
    check("naca0012_coordinates.csv", 0.12, 0.0008)
    check("naca0015_coordinates.csv", 0.15, 0.0008)
    check_compressibility()
    check_cross_source()
    print("\nreference data ok")
