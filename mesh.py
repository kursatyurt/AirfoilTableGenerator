"""Structured C-mesh around an airfoil, written as SU2 with markers Airfoil/Inlet/Outlet.

Forked from FALCON's meshing.py (Prisha22/FALCON). Changes from upstream:
  * the airfoil is split into upper / leading-edge / lower blocks by INDEX around
    the leading-edge point, not by a stateful scan with duplicate-point pops.
    Upstream's scan gave a symmetric airfoil an asymmetric surface discretisation
    (484 upper vs 495 lower nodes on NACA 0012) and hence CL != 0 at zero incidence.
  * the lower far-field line's progression sign now mirrors the upper one.
  * the output path is an argument instead of 'airfoil.su2' in the cwd, and the
    block counts are keyword arguments.
"""
import math

import numpy as np

RHO, MU, A_SOUND = 1.225, 1.813e-5, 341.348  # sea-level air


def first_cell_thickness(Re, M, y_plus=1.0):
    """Wall spacing for a target y+, from a flat-plate skin-friction correlation."""
    cf = (2 * np.log10(float(Re)) - 0.65) ** -2.3
    tau_w = cf * 0.5 * RHO * (float(M) * A_SOUND) ** 2
    u_star = np.sqrt(tau_w / RHO)
    return y_plus * MU / (RHO * u_star)


def _split(x, y, x_split):
    """TE->upper->LE->lower->TE coordinates -> (upper_aft, leading_edge, lower_aft).

    Split by index about the LE, so a symmetric airfoil yields mirror-image blocks.
    Blocks share their junction points. An open trailing edge (Selig files usually
    end at y = +/-t/2) is closed onto the MIDPOINT of the two TE points -- closing it
    onto the upper one, as upstream did, skews an otherwise symmetric airfoil.
    """
    pts = [[float(a), float(b), 0.0] for a, b in zip(x, y)]
    te = [(pts[0][0] + pts[-1][0]) / 2, (pts[0][1] + pts[-1][1]) / 2, 0.0]
    pts[0] = pts[-1] = te
    le = int(np.argmin(x))                                    # leading-edge point
    iu = max(i for i in range(le + 1) if x[i] > x_split)       # last upper point aft of the split
    il = min(i for i in range(le, len(x)) if x[i] > x_split)   # first lower point aft of the split
    return pts[: iu + 1], pts[iu : il + 1], pts[il :]


def generate_mesh(xcoords, ycoords, Re, M, y_plus=1.0, path="airfoil.su2",
                  n_airfoil=401, n_wake=301, n_le=180, inlet_radius=15.0,
                  downstream=25.0, growth=1.2, x_split=0.1, show_graphics=False,
                  hide_output=True):
    import gmsh

    te_thickness = 1e-3
    te_growth = 1.015

    dy = first_cell_thickness(Re, M, y_plus)
    log_arg = (inlet_radius * (growth - 1) / dy) + 1
    n_volume = int(math.log(log_arg) / math.log(growth)) if log_arg > 1 else 120
    print(f"y+={y_plus} -> first cell {dy:.4e} m, {n_volume} normal layers at growth {growth}")

    upper_pts, le_pts, lower_pts = _split(np.asarray(xcoords), np.asarray(ycoords), x_split)

    gmsh.initialize()
    gmsh.model.add("airfoil")
    if hide_output:
        gmsh.option.setNumber("General.Terminal", 0)
    geo = gmsh.model.geo

    upper = [geo.addPoint(*p) for p in upper_pts]
    le = [geo.addPoint(*p) for p in le_pts[1:-1]]
    lower = [geo.addPoint(*p) for p in lower_pts[1:]]
    af_te, af_top, af_bottom = upper[0], upper[-1], lower[0]
    le = [af_top] + le + [af_bottom]
    lower = lower + [af_te]

    af_upper = geo.addBSpline(upper)          # TE -> top junction
    af_le = geo.addBSpline(le)                # top junction -> LE -> bottom junction
    af_lower = geo.addBSpline(lower)          # bottom junction -> TE

    center = geo.addPoint(x_split, 0, 0)
    inlet_top = geo.addPoint(x_split, inlet_radius, 0)
    inlet_bottom = geo.addPoint(x_split, -inlet_radius, 0)
    front = geo.addCircleArc(inlet_top, center, inlet_bottom)
    afTop_inletTop = geo.addLine(af_top, inlet_top)
    inletBottom_afBottom = geo.addLine(inlet_bottom, af_bottom)
    inlet_sec = geo.addPlaneSurface([geo.addCurveLoop(
        [front, inletBottom_afBottom, -af_le, afTop_inletTop])])

    top_te = geo.addPoint(1, inlet_radius, 0)
    top_line = geo.addLine(top_te, inlet_top)
    topTe_afTe = geo.addLine(top_te, af_te)
    top_sec = geo.addPlaneSurface([geo.addCurveLoop(
        [-af_upper, -afTop_inletTop, top_line, -topTe_afTe])])

    bottom_te = geo.addPoint(1, -inlet_radius, 0)
    bottom_line = geo.addLine(inlet_bottom, bottom_te)
    afTe_bottomTe = geo.addLine(af_te, bottom_te)
    bottom_sec = geo.addPlaneSurface([geo.addCurveLoop(
        [-inletBottom_afBottom, -af_lower, -afTe_bottomTe, bottom_line])])

    top_wake_pt = geo.addPoint(downstream, inlet_radius, 0)
    center_wake_pt = geo.addPoint(downstream, 0, 0)
    bottom_wake_pt = geo.addPoint(downstream, -inlet_radius, 0)
    top_wake_line = geo.addLine(top_wake_pt, top_te)
    center_wake_line = geo.addLine(af_te, center_wake_pt)
    outlet_top = geo.addLine(center_wake_pt, top_wake_pt)
    top_wake_sec = geo.addPlaneSurface([geo.addCurveLoop(
        [topTe_afTe, center_wake_line, outlet_top, top_wake_line])])
    outlet_bottom = geo.addLine(bottom_wake_pt, center_wake_pt)
    bottom_wake_line = geo.addLine(bottom_te, bottom_wake_pt)
    bottom_wake_sec = geo.addPlaneSurface([geo.addCurveLoop(
        [-center_wake_line, outlet_bottom, bottom_wake_line, afTe_bottomTe])])

    tc = geo.mesh.setTransfiniteCurve
    tc(front, n_le, "Bump", coef=-0.1)
    tc(af_le, n_le)
    # normal direction: sign follows whether the curve points away from the wall
    tc(afTop_inletTop, n_volume, "Progression", growth)
    tc(inletBottom_afBottom, n_volume, "Progression", -growth)
    tc(topTe_afTe, n_volume, "Progression", -growth)
    tc(afTe_bottomTe, n_volume, "Progression", growth)
    tc(outlet_top, n_volume, "Progression", growth)
    tc(outlet_bottom, n_volume, "Progression", -growth)
    # chordwise: each far-field line matches the airfoil curve it faces
    tc(af_upper, n_airfoil, "Progression", -te_growth)
    tc(top_line, n_airfoil, "Progression", -te_growth)
    tc(af_lower, n_airfoil, "Progression", te_growth)
    tc(bottom_line, n_airfoil, "Progression", te_growth)
    # wake: start at the TE thickness, stretch to a uniform far-wake spacing
    wake_growth = min(((downstream - 1.0) / (n_wake * 0.5) / (te_thickness * 2.0))
                      ** (1.0 / (n_wake - 1)), 1.05)
    tc(top_wake_line, n_wake, "Progression", -wake_growth)
    tc(center_wake_line, n_wake, "Progression", wake_growth)
    tc(bottom_wake_line, n_wake, "Progression", wake_growth)

    sections = [inlet_sec, top_sec, bottom_sec, top_wake_sec, bottom_wake_sec]
    for s in sections:
        geo.mesh.setTransfiniteSurface(s)
        geo.mesh.setRecombine(2, s)

    geo.addPhysicalGroup(1, [af_upper, af_le, af_lower], name="Airfoil")
    geo.addPhysicalGroup(1, [outlet_top, outlet_bottom, top_line, top_wake_line,
                             bottom_line, bottom_wake_line], name="Outlet")
    geo.addPhysicalGroup(1, [front], name="Inlet")
    geo.addPhysicalGroup(2, sections, name="FlowDomain")
    geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write(str(path))
    if show_graphics:
        gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
        gmsh.fltk.run()
    gmsh.finalize()

    # ponytail: gmsh counts FlowDomain in NMARK; SU2 wants only the 3 boundaries.
    text = open(path).read().splitlines(keepends=True)
    with open(path, "w") as f:
        f.writelines("NMARK= 3\n" if l.startswith("NMARK") else l for l in text)
    return path


def _selftest():
    """Split must be mirror-symmetric for a symmetric airfoil."""
    th = np.linspace(0, math.pi, 61)
    xs = (1 - np.cos(th)) / 2                      # cosine-clustered chordwise
    t = 0.6 * (0.2969 * np.sqrt(xs) - 0.1260 * xs - 0.3516 * xs**2
               + 0.2843 * xs**3 - 0.1015 * xs**4)  # NACA 0012 thickness
    x = np.concatenate([xs[::-1], xs[1:]])
    y = np.concatenate([t[::-1], -t[1:]])
    up, le, lo = _split(x, y, 0.1)
    assert len(up) == len(lo), (len(up), len(lo))
    assert np.allclose([p[1] for p in up], [-p[1] for p in lo[::-1]])
    assert np.allclose([p[1] for p in le], [-p[1] for p in le[::-1]])
    assert up[-1] == le[0] and le[-1] == lo[0]     # blocks share junctions
    assert 1e-6 < first_cell_thickness(1e6, 0.15) < 1e-4
    print("mesh selftest ok")


if __name__ == "__main__":
    _selftest()
