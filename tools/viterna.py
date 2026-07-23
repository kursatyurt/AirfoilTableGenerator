#!/usr/bin/env python3
"""Viterna 360-degree airfoil polar extrapolation.

Extends a CFD/experimental polar (available over a limited alpha range) out to
the full +-180 deg range needed by rotor codes, using the classic Viterna
flat-plate model anchored at the measured stall points. CL/CD follow Viterna;
CM is a bounded placeholder (see cm_extrap).

Reusable: `extrapolate_column(a, cl, cd, cm, out_alpha, cdmax)` returns CL/CD/CM
on `out_alpha`, interpolating CFD inside the data range and Viterna outside.

Run `python tools/viterna.py` for a self-check.
"""
import numpy as np

sind = lambda d: np.sin(np.radians(d))
cosd = lambda d: np.cos(np.radians(d))


def viterna_coeffs(a_anchor, cl_a, cd_a, cdmax):
    """Viterna A1,A2,B1,B2 anchored at stall angle a_anchor (deg, >0)."""
    ar = np.radians(a_anchor)
    B1 = cdmax
    A1 = B1 / 2.0
    A2 = (cl_a - cdmax * np.sin(ar) * np.cos(ar)) * np.sin(ar) / np.cos(ar) ** 2
    B2 = (cd_a - cdmax * np.sin(ar) ** 2) / np.cos(ar)
    return A1, A2, B1, B2


def extrap_pos(a, coeffs, a_high, cl_adj=0.7):
    """CL,CD for a in (a_high, 180] deg. Guards the sin->0 blow-up near 180."""
    A1, A2, B1, B2 = coeffs
    if a <= 90:
        cl = A1 * sind(2 * a) + A2 * cosd(a) ** 2 / sind(a)
        cd = B1 * sind(a) ** 2 + B2 * cosd(a)
    elif a <= 180 - a_high:
        r = 180 - a
        cl = -cl_adj * (A1 * sind(2 * r) + A2 * cosd(r) ** 2 / sind(r))
        cd = B1 * sind(r) ** 2 + B2 * cosd(r)
    else:  # (180-a_high, 180]: linear ramp CL->0, hold CD flat at the seam
        r = a_high
        cl_seam = -cl_adj * (A1 * sind(2 * r) + A2 * cosd(r) ** 2 / sind(r))
        cd = B1 * sind(r) ** 2 + B2 * cosd(r)
        t = (a - (180 - a_high)) / a_high
        cl = cl_seam * (1 - t)
    return cl, cd


def cm_extrap(a, a_anchor, cm_anchor):
    """Hold near-stall CM through deep stall, ramp to 0 at +-180.
    ponytail: flat-plate CM placeholder -- CM barely drives rotor power/thrust;
    upgrade to a cp-travel model if pitch-link loads ever matter."""
    if a >= 180:
        return 0.0
    t = (a - a_anchor) / (180 - a_anchor)
    return float(cm_anchor * (1 - t))


def extrapolate_column(a, cl, cd, cm, out_alpha, cdmax=2.05):
    """Full-range CL/CD/CM on out_alpha. a ascending (deg); interp CFD inside,
    Viterna outside; the negative side mirrors the positive extrapolation."""
    a = np.asarray(a, float)
    a_high, a_low = a[-1], a[0]
    cp = viterna_coeffs(a_high, cl[-1], cd[-1], cdmax)
    cn = viterna_coeffs(-a_low, cl[0], cd[0], cdmax)   # mirror anchor
    CL, CD, CM = [], [], []
    for aa in out_alpha:
        if a_low <= aa <= a_high:
            CL.append(float(np.interp(aa, a, cl)))
            CD.append(float(np.interp(aa, a, cd)))
            CM.append(float(np.interp(aa, a, cm)))
        elif aa > a_high:
            l, d = extrap_pos(aa, cp, a_high)
            CL.append(l); CD.append(d); CM.append(cm_extrap(aa, a_high, cm[-1]))
        else:
            l, d = extrap_pos(-aa, cn, -a_low)
            CL.append(-l); CD.append(d); CM.append(cm_extrap(-aa, -a_low, cm[0]))
    return CL, CD, CM


def _selfcheck():
    # synthetic attached polar over [-10,16]; extrapolate to +-180
    a = np.arange(-10, 17, 2.0)
    cl = 0.11 * (a + 1.0)                 # slope 0.11/deg, zero-lift ~-1 deg
    cd = 0.008 + 0.0004 * a ** 2
    cm = np.full_like(a, -0.02)
    grid = [-180, -90, -20, -10, 0, 8, 16, 45, 90, 135, 180]
    CL, CD, CM = extrapolate_column(a, cl, cd, cm, grid, cdmax=2.05)
    i = {g: k for k, g in enumerate(grid)}
    assert abs(CL[i[0]] - 0.11) < 0.02, CL[i[0]]           # slope preserved at 0
    assert abs(CL[i[-180]]) < 0.05 and abs(CL[i[180]]) < 0.05, "CL~0 at +-180"
    assert abs(CD[i[90]] - 2.05) < 0.1, CD[i[90]]          # flat plate at 90
    assert CD[i[45]] > CD[i[0]], "CD rises into deep stall"
    print("viterna self-check OK:  CL@0=%.3f  CD@90=%.2f  CL@180=%.3f"
          % (CL[i[0]], CD[i[90]], CL[i[180]]))


if __name__ == "__main__":
    _selfcheck()
