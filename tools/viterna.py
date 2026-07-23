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


def extrapolate_column(a, cl, cd, cm, out_alpha, cdmax=2.05, warn=None):
    """Full-range CL/CD/CM on out_alpha. a ascending (deg); interp CFD up to the
    stall anchors, Viterna beyond.

    Anchors at the CLmax / CLmin points, not the sweep endpoints: steady-RANS
    post-stall points are unreliable (per the tool README), so Viterna past the
    lift peak is more trustworthy than keeping them, and anchoring there is C0-
    continuous by construction. `warn(msg)`, if given, is called for suspect
    inputs (stall not captured -> anchor is pre-stall; large interior gaps that
    np.interp silently bridges)."""
    a = np.asarray(a, float)
    cl = np.asarray(cl, float); cd = np.asarray(cd, float); cm = np.asarray(cm, float)
    ip, im = int(np.argmax(cl)), int(np.argmin(cl))   # + / - stall anchors
    a_high, a_low = a[ip], a[im]
    cp = viterna_coeffs(a_high, cl[ip], cd[ip], cdmax)
    cn = viterna_coeffs(-a_low, -cl[im], cd[im], cdmax)   # reflect angle AND lift
    if warn is not None:
        if ip == len(a) - 1:
            warn(f"+ stall not captured: CLmax at sweep top {a_high:g} deg; Viterna anchored pre-stall")
        if im == 0:
            warn(f"- stall not captured: CLmin at sweep bottom {a_low:g} deg; Viterna anchored pre-stall")
        gaps = np.diff(a)
        for k in np.where(gaps > 3.0)[0]:
            warn(f"interp bridges {gaps[k]:g} deg gap ({a[k]:g} -> {a[k+1]:g} deg)")
    CL, CD, CM = [], [], []
    for aa in out_alpha:
        if a_low <= aa <= a_high:
            CL.append(float(np.interp(aa, a, cl)))
            CD.append(float(np.interp(aa, a, cd)))
            CM.append(float(np.interp(aa, a, cm)))
        elif aa > a_high:
            l, d = extrap_pos(aa, cp, a_high)
            CL.append(l); CD.append(d); CM.append(cm_extrap(aa, a_high, cm[ip]))
        else:  # reflect back: negative deep stall mirrors the positive curve
            l, d = extrap_pos(-aa, cn, -a_low)
            CL.append(-l); CD.append(d); CM.append(cm_extrap(-aa, -a_low, cm[im]))
    return CL, CD, CM


def _selfcheck():
    # synthetic cambered polar with a real stall peak at +14 / -12, so the
    # anchors are interior and post-stall points get replaced by Viterna.
    a = np.arange(-16, 21, 2.0)
    cl = np.clip(0.11 * (a + 1.0), -1.35, 1.5)   # slope 0.11/deg, zero-lift ~-1
    cl[a > 14] = 1.5 - 0.03 * (a[a > 14] - 14)   # droop past +stall
    cl[a < -12] = -1.35 + 0.03 * (-12 - a[a < -12])   # recover toward 0 past -stall
    cd = 0.008 + 0.0004 * a ** 2
    cm = np.full_like(a, -0.02)
    grid = [-180, -90, -60, -45, -20, -14, -10, 0, 8, 14, 20, 45, 90, 135, 180]
    msgs = []
    CL, CD, CM = extrapolate_column(a, cl, cd, cm, grid, cdmax=2.05, warn=msgs.append)
    i = {g: k for k, g in enumerate(grid)}
    assert abs(CL[i[0]] - 0.11) < 0.03, CL[i[0]]           # slope preserved at 0
    assert abs(CL[i[-180]]) < 0.05 and abs(CL[i[180]]) < 0.05, "CL~0 at +-180"
    assert abs(CD[i[90]] - 2.05) < 0.1, CD[i[90]]          # flat plate at 90
    assert CD[i[45]] > CD[i[0]], "CD rises into deep stall"
    # negative deep stall: correct sign (no flip), then monotonic recovery to 0
    assert CL[i[-45]] < 0 and CL[i[-60]] < 0, "neg-stall CL stays negative"
    assert CL[i[-90]] > CL[i[-60]] > CL[i[-45]], "neg CL recovers to 0 past the 45deg bump"
    assert abs(CL[i[-90]]) < 0.05, "CL~0 at -90"
    assert not msgs, f"unexpected warnings: {msgs}"        # stall IS captured here
    # a truncated (pre-stall) sweep must warn on both ends
    warn2 = []
    extrapolate_column(a[:10], cl[:10], cd[:10], cm[:10], grid, warn=warn2.append)
    assert any("stall not captured" in m for m in warn2), "should warn on truncated sweep"
    print("viterna self-check OK:  CL@0=%.3f  CD@90=%.2f  CL@-45=%.3f  CL@180=%.3f"
          % (CL[i[0]], CD[i[90]], CL[i[-45]], CL[i[180]]))


if __name__ == "__main__":
    _selfcheck()
