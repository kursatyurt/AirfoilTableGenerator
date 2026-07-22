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


if __name__ == "__main__":
    # tolerance is loose enough for real measured ordinates, tight enough that a
    # wrong digit (>= 0.001 in the 4th place) fails
    check("naca0012_coordinates.csv", 0.12, 0.0008)
    check("naca0015_coordinates.csv", 0.15, 0.0008)
    print("reference tables ok")
