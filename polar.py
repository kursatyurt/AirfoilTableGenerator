#!/usr/bin/env python
"""Airfoil polar: FALCON's gmsh C-mesh + SU2 RANS-SA sweep over AoA.

    source env.sh
    python polar.py --airfoil naca0012 --re 1e6 --mach 0.15 --aoa -4:16:2 --np 8
"""
import argparse, math, os, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FALCON = ROOT / "opt" / "FALCON"
sys.path.insert(0, str(FALCON))

RHO, A_SOUND = 1.225, 341.348  # sea level, matches FALCON's meshing.py


def auto_farfield(mach):
    """Return a conservative C-mesh radius in chords for the Mach number.

    A fixed-outlet convergence sweep accepted 25c at M=0.1 and M=0.3, then
    35c at M=0.5, using two lift counts and one drag count as the limits.
    Grow from that high-Mach anchor rather than repeating a domain study for
    each angle in a production table.
    """
    return max(25.0, 25.0 + 50.0 * (mach - 0.3))


def auto_downstream(farfield, mach):
    """Return the outlet location in chords, ramping the wake only when needed."""
    wake_target = min(75.0, 35.0 + 200.0 * max(0.0, mach - 0.1))
    return max(wake_target, farfield + 10.0)


COMMON = """\
MATH_PROBLEM= DIRECT
KIND_TURB_MODEL= SA
MARKER_FAR= ( Inlet, Outlet )
MARKER_PLOTTING= ( Airfoil )
MARKER_MONITORING= ( Airfoil )
REF_ORIGIN_MOMENT_X= 0.25
REF_ORIGIN_MOMENT_Y= 0.00
REF_ORIGIN_MOMENT_Z= 0.00
REF_LENGTH= 1.0
REF_AREA= 1.0
NUM_METHOD_GRAD= GREEN_GAUSS
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
VENKAT_LIMITER_COEFF= 0.03
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
SLOPE_LIMITER_TURB= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT
TIME_DISCRE_TURB= EULER_IMPLICIT
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-6
LINEAR_SOLVER_ITER= 15
% Coefficients must be Cauchy-stable and the residual must reach its floor.
% The residual tail damps the remaining CD oscillation; the tighter Cauchy test
% prevents a numerically quiet residual from accepting a drifting polar value.
CONV_CAUCHY_ELEMS= 100
CONV_CAUCHY_EPS= 1E-5
CONV_RESIDUAL_MINVAL= -6
CONV_STARTITER= 100
MESH_FILENAME= airfoil.su2
MESH_FORMAT= SU2
SOLUTION_FILENAME= solution_flow.dat
RESTART_FILENAME= restart_flow.dat
TABULAR_FORMAT= CSV
OUTPUT_FILES= ( RESTART, PARAVIEW, SURFACE_CSV )
"""


TRANSITION = {
    # LM rides on top of KIND_TURB_MODEL= SA; it does not replace it
    "none": "",                      # fully turbulent from the leading edge
    "lm": "KIND_TRANS_MODEL= LM\n",  # Langtry-Menter gamma-Re_theta, two extra equations
}


def make_cfg(regime, aoa, re, mach, iters, restart, transition="none", tu=0.001,
             output_suffix="", unsteady=False, urans_steps_per_chord=200,
             urans_convective_times=10, urans_inner=30):
    v = mach * A_SOUND
    output_suffix = output_suffix or ("_urans" if unsteady else "")
    tag = f"{aoa:+.2f}{output_suffix}"
    res = "RMS_PRESSURE" if regime == "inc" else "RMS_DENSITY"
    # LM's two extra transport equations make the transition front hunt at high
    # CFL and stall convergence. Capping the adaptive CFL ceiling trades wall time
    # for coefficients that actually settle. Fully-turbulent SA keeps the fast ceiling.
    cfl_max = 15.0 if transition != "none" else 75.0
    if unsteady:
        inc_cfl = "CFL_NUMBER= 2.0\nCFL_ADAPT= NO"
    elif transition == "none":
        # Start gently, then ramp to CFL 100 as in SU2's incompressible NACA case.
        inc_cfl = "CFL_NUMBER= 10.0\nCFL_ADAPT= YES\nCFL_ADAPT_PARAM= ( 0.8, 1.1, 1.0, 100.0 )"
    else:
        inc_cfl = "CFL_NUMBER= 10.0\nCFL_ADAPT= NO"
    if regime == "inc":
        # mu from Re with chord = 1 m
        mu = RHO * v * 1.0 / re
        head = f"""\
SOLVER= INC_RANS
INC_DENSITY_INIT= {RHO}
INC_VELOCITY_INIT= ( {v * math.cos(math.radians(aoa)):.8f}, {v * math.sin(math.radians(aoa)):.8f}, 0.0 )
INC_NONDIM= INITIAL_VALUES
INC_DENSITY_REF= {RHO}
VISCOSITY_MODEL= CONSTANT_VISCOSITY
MU_CONSTANT= {mu:.10e}
FREESTREAM_NU_FACTOR= 4.0
MARKER_HEATFLUX= ( Airfoil, 0.0 )
AOA= {aoa}
{inc_cfl}
CONV_NUM_METHOD_FLOW= {"LD2" if unsteady else "FDS"}
"""
    else:
        head = f"""\
SOLVER= RANS
MACH_NUMBER= {mach}
REYNOLDS_NUMBER= {re:.6g}
REYNOLDS_LENGTH= 1.0
FREESTREAM_TEMPERATURE= 288.15
FREESTREAM_PRESSURE= 101325.0
MARKER_HEATFLUX= ( Airfoil, 0.0 )
AOA= {aoa}
CFL_NUMBER= {20.0 if unsteady else 5.0}
CFL_ADAPT= {"NO" if unsteady else "YES"}
CFL_ADAPT_PARAM= ( 0.1, 2.0, 5.0, {cfl_max} )
CONV_NUM_METHOD_FLOW= {"JST" if unsteady else "ROE"}
ENTROPY_FIX_COEFF= 0.05
LOW_MACH_PREC= YES
"""
    common = COMMON
    if unsteady:
        # Match SU2's unsteady NACA RANS practice: centered low-dissipation
        # flow flux, no MUSCL reconstruction, and inexpensive implicit inner solves.
        common = common.replace("MUSCL_FLOW= YES\nSLOPE_LIMITER_FLOW= VENKATAKRISHNAN\nVENKAT_LIMITER_COEFF= 0.03",
                                "MUSCL_FLOW= NO\nVENKAT_LIMITER_COEFF= 0.03")
        urans_steps = math.ceil(urans_steps_per_chord * urans_convective_times)
        # Discard the first 40% before averaging the physical force signal.
        unsteady_cfg = f"""\
TIME_DOMAIN= YES
TIME_MARCHING= DUAL_TIME_STEPPING-2ND_ORDER
TIME_STEP= {1.0 / (v * urans_steps_per_chord):.10e}
TIME_ITER= {urans_steps}
INNER_ITER= {urans_inner}
MAX_TIME= {urans_convective_times / v * 1.01:.10e}
WINDOW_START_ITER= {int(urans_steps * 0.4)}
"""
        history = "HISTORY_OUTPUT= ( ITER, RMS_RES, AERO_COEFF, TAVG_AERO_COEFF, AOA )"
        screen = f"SCREEN_OUTPUT= ( TIME_ITER, INNER_ITER, {res}, LIFT, DRAG, MOMENT_Z, TAVG_LIFT, TAVG_DRAG, TAVG_MOMENT_Z )"
        limit = ""
    else:
        unsteady_cfg = ""
        history = "HISTORY_OUTPUT= ( ITER, RMS_RES, AERO_COEFF, AOA )"
        screen = f"SCREEN_OUTPUT= ( INNER_ITER, {res}, RMS_MOMENTUM-X, LIFT, DRAG, MOMENT_Z )"
        limit = f"ITER= {iters}"
    return head + common + TRANSITION[transition] + unsteady_cfg + f"""\
CONV_FIELD= ( LIFT, DRAG, MOMENT_Z, {res} )
{screen}
{history}
FREESTREAM_TURBULENCEINTENSITY= {tu}
RESTART_SOL= {"YES" if restart else "NO"}
{limit}
CONV_FILENAME= history_{tag}
VOLUME_FILENAME= flow_{tag}
SURFACE_FILENAME= surface_{tag}
"""


def parse_range(spec):
    """'-4:16:2' -> inclusive range; '0,2,4' -> explicit list. Used for AoA and Mach."""
    if ":" in spec:
        lo, hi, step = (float(x) for x in spec.split(":"))
        n = int(round((hi - lo) / step))
        return [lo + i * step for i in range(n + 1)]
    return [float(x) for x in spec.split(",")]


parse_aoa = parse_range  # back-compat alias


def past_stall(extreme_cl, cl, direction, drop):
    """True after lift reverses by ``drop`` on an AoA march branch."""
    return drop > 0 and (cl < extreme_cl - drop if direction > 0 else cl > extreme_cl + drop)


def read_history(path):
    """Last row of an SU2 history CSV -> {CL, CD, CMz, ...}. Headers are quoted+padded."""
    import pandas as pd
    df = pd.read_csv(path)
    df.columns = [c.strip().strip('"').strip() for c in df.columns]
    return df.iloc[-1].to_dict()


def find_dat(name):
    p = Path(name)
    if p.is_file():
        return p
    p = FALCON / "Airfoil_DAT_Selig" / (name if name.endswith(".dat") else name + ".dat")
    if p.is_file():
        return p
    db = FALCON / "Airfoil_DAT_Selig"
    import difflib
    near = difflib.get_close_matches(Path(name).stem, [f.stem for f in db.glob("*.dat")], n=8, cutoff=0.5)
    sys.exit(f"airfoil not found: {name} (not a file, not in {db})"
             + (f"\ndid you mean: {', '.join(near)}" if near else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--airfoil", required=True, help="name in FALCON's Selig database, or path to a .dat")
    ap.add_argument("--re", type=float, default=1e6)
    ap.add_argument("--mach", default="0.15", help="single value, or lo:hi:step / comma list "
                    "to sweep Mach; each Mach gets its own mesh and subdir")
    ap.add_argument("--aoa", default="-4:16:2", help="lo:hi:step or comma list")
    ap.add_argument("--regime", choices=["inc", "comp"], default="comp",
                    help="comp (default) = compressible RANS with Roe + low-Mach "
                         "preconditioning, consistent across the whole Mach sweep; "
                         "inc = INC_RANS, only for strictly incompressible cases")
    ap.add_argument("--np", type=int, default=None,
                    help="MPI ranks; defaults to machine.conf from tune_np.py, else half the cores")
    ap.add_argument("--iters", type=int, default=10000,
                    help="max solver iterations per angle (default 10000); the run stops "
                         "earlier when the LIFT/DRAG/MOMENT_Z Cauchy criterion is met")
    ap.add_argument("--yplus", type=float, default=1.0)
    ap.add_argument("--farfield", type=float,
                    help="optional farfield radius in chords; by default it is sized from Mach "
                         "(25c minimum, 35c at M=0.5)")
    ap.add_argument("--transition", choices=list(TRANSITION), default="none",
                    help="laminar-turbulent transition on top of SA: 'lm' = Langtry-Menter "
                         "gamma-Re_theta (two extra transport equations)")
    ap.add_argument("--tu", type=float, default=0.001,
                    help="freestream turbulence intensity for the transition models "
                         "(0.001 = 0.1%%, a low-turbulence wind tunnel)")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--stall-drop", type=float, default=0.02,
                    help="switch a branch to URANS after CL drops this far below its peak; "
                         "set 0 to keep the branch steady (default 0.02)")
    ap.add_argument("--urans-steps-per-chord", type=int, default=200,
                    help="physical URANS time resolution after stall (default 200)")
    ap.add_argument("--urans-convective-times", type=float, default=10.0,
                    help="physical URANS duration after stall, in chord transit times (default 10)")
    ap.add_argument("--urans-inner-iters", type=int, default=30,
                    help="dual-time iterations per physical URANS step (default 30)")
    # ponytail: argparse reads a leading '-' as an option, so "--aoa -4:16:2" fails.
    argv, rest = [], list(sys.argv[1:])
    while rest:
        v = rest.pop(0)
        argv.append(f"--aoa={rest.pop(0)}" if v == "--aoa" and rest else v)
    a = ap.parse_args(argv)
    if a.stall_drop < 0:
        ap.error("--stall-drop must be non-negative")
    if a.farfield is not None and a.farfield <= 0:
        ap.error("--farfield must be positive")
    if a.urans_steps_per_chord < 1 or a.urans_convective_times <= 0 or a.urans_inner_iters < 1:
        ap.error("URANS resolution, duration, and inner iterations must be positive")

    if a.np is None:
        from tune_np import stored_np
        tuned = stored_np()
        a.np = tuned or max(1, os.cpu_count() // 2)
        if tuned is None:
            print(f"no machine.conf; guessing {a.np} ranks. Run 'python tune_np.py' once "
                  f"to measure the right number for this machine.")

    dat = find_dat(a.airfoil)
    base = Path(a.outdir or ROOT / "runs" / dat.stem)

    from read_airfoil import read_airfoil_coordinates
    x, y = read_airfoil_coordinates(str(dat.parent), dat.name)

    machs = parse_range(a.mach)
    for mach in machs:
        # one mesh + sweep per Mach; the wall spacing depends on Mach, so each
        # needs its own dir. Single Mach keeps the flat runs/<stem>/ layout.
        case = base if len(machs) == 1 else base / f"M{mach:g}"
        run_sweep(a, dat, x, y, mach, case)


def run_sweep(a, dat, x, y, mach, case):
    case.mkdir(parents=True, exist_ok=True)
    farfield = a.farfield if a.farfield is not None else auto_farfield(mach)

    mesh = case / "airfoil.su2"
    if not mesh.exists():
        from mesh import generate_mesh
        generate_mesh(x, y, a.re, mach, y_plus=a.yplus, path=mesh,
                      inlet_radius=farfield, downstream=auto_downstream(farfield, mach))
    print(f"mesh: {mesh} (farfield {farfield:g}c)")

    sol, restart_dat = case / "solution_flow.dat", case / "restart_flow.dat"

    def su2(cfg_name, log_name):
        with open(case / log_name, "w") as log:
            return subprocess.call(["mpirun", "-n", str(a.np), "SU2_CFD", cfg_name],
                                   cwd=case, stdout=log, stderr=subprocess.STDOUT)

    def solve(aoa, restart, transition, suffix="", tu=None, unsteady=False):
        """Run one angle; return (CL, CD, CMz, converged) or None. Leaves the
        converged field in solution_flow.dat for the next (warm-started) angle."""
        tag = f"{aoa:+.2f}{'_urans' if unsteady else ''}"
        cfg = case / f"aoa_{tag}{suffix}.cfg"
        cfg.write_text(make_cfg(a.regime, aoa, a.re, mach, a.iters, restart, transition,
                                a.tu if tu is None else tu, "_urans" if unsteady else "", unsteady,
                                a.urans_steps_per_chord, a.urans_convective_times, a.urans_inner_iters))
        label = " URANS" if unsteady else (" turbulent seed" if suffix == "_seed" else "")
        print(f"--- M {mach:g} AoA {aoa:g}{label} ({a.np} ranks) ...", end=" ", flush=True)
        rc = su2(cfg.name, f"aoa_{tag}{suffix}.log")
        hist = case / f"history_{tag}.csv"
        if rc != 0 or not hist.exists():
            print(f"FAILED (rc={rc}), see aoa_{tag}{suffix}.log")
            return None
        h = read_history(hist)
        n_iter = sum(1 for _ in open(hist)) - 1  # history rows = iterations run
        converged = unsteady or n_iter < a.iters
        if converged:                            # only a converged field warm-starts the
            shutil.copy(restart_dat, sol)        # next angle; a diverged one (e.g. post-
                                                 # stall) would poison the whole march
        if unsteady:
            return h["tavg[CL]"], h["tavg[CD]"], h.get("tavg[CMz]", float("nan")), converged
        return h["CL"], h["CD"], h.get("CMz", float("nan")), converged

    # Fan outward from the angle nearest zero lift. The near-zero angles are the hardest
    # for LM (symmetric loading -> transition front hunts); doing them first from a
    # fully-turbulent seed, then marching into progressively more asymmetric (and stabler)
    # incidences with 1-deg warm-started steps, avoids the cold -4 deg start and the long
    # cold-transient. Each branch (up from the pivot, then down) is monotonic in AoA.
    angles = parse_aoa(a.aoa)
    pivot = min(angles, key=abs)
    up = sorted(x for x in angles if x >= pivot)
    down = sorted((x for x in angles if x < pivot), reverse=True)

    rows = {}

    def record(aoa, r):
        if r:
            rows[aoa] = r
            print(f"CL={r[0]:.4f} CD={r[1]:.5f}" + ("" if r[3] else "  [not converged]"))

    # pivot: seed with the SAME transition solver but at high freestream turbulence
    # (SEED_TU), which drives gamma->1 so it converges fast and stable like a fully-
    # turbulent solve. Restarting the real low-Tu LM run from it is variable-compatible
    # (same solver) -- unlike seeding from plain SA, which lacks LM's transition
    # variables and makes SU2 diverge (NaN) on the mismatched restart.
    SEED_TU = 0.10
    seeded = False
    if a.transition != "none":
        if solve(pivot, False, a.transition, "_seed", tu=SEED_TU):
            seeded = True
    record(pivot, solve(pivot, seeded, a.transition))
    pivot_sol = case / "solution_flow.pivot.dat"
    if sol.exists():
        shutil.copy(sol, pivot_sol)  # save the pivot field to reseed the down branch

    peak_cl = rows[pivot][0] if pivot in rows else -float("inf")
    post_stall = False
    for aoa in up[1:]:               # ascending above the pivot, warm-started
        r = solve(aoa, True, a.transition, unsteady=post_stall)
        record(aoa, r)
        if r:
            if not post_stall and past_stall(peak_cl, r[0], 1, a.stall_drop):
                post_stall = True
                print(f"  stall detected at {aoa:g} deg (CL {r[0]:.4f}, peak {peak_cl:.4f}); "
                      "switching the remaining up branch to URANS")
            peak_cl = max(peak_cl, r[0])
    if down and pivot_sol.exists():
        shutil.copy(pivot_sol, sol)  # rewind to the pivot field before fanning down
    trough_cl = rows[pivot][0] if pivot in rows else float("inf")
    post_stall = False
    for aoa in down:                 # descending below the pivot, warm-started
        r = solve(aoa, True, a.transition, unsteady=post_stall)
        record(aoa, r)
        if r:
            if not post_stall and past_stall(trough_cl, r[0], -1, a.stall_drop):
                post_stall = True
                print(f"  stall detected at {aoa:g} deg (CL {r[0]:.4f}, trough {trough_cl:.4f}); "
                      "switching the remaining down branch to URANS")
            trough_cl = min(trough_cl, r[0])

    rows = [(aoa, *rows[aoa]) for aoa in sorted(rows)]
    if not rows:
        print(f"M {mach:g}: no converged runs")
        return

    csv = case / "polar.csv"
    csv.write_text("aoa,cl,cd,cm,converged\n" +
                   "".join(f"{r[0]:g},{r[1]:.6f},{r[2]:.6f},{r[3]:.6f},{int(r[4])}\n" for r in rows))
    print(f"\n{csv}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    al, cl, cd = [r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for i, (xs, ys, xl, yl) in enumerate([(al, cl, "AoA [deg]", "CL"),
                                          (al, cd, "AoA [deg]", "CD"),
                                          (cd, cl, "CD", "CL")]):
        ax[i].plot(xs, ys, "o-")
        ax[i].set_xlabel(xl); ax[i].set_ylabel(yl); ax[i].grid(True)
    fig.suptitle(f"{dat.stem}  Re={a.re:.3g}  M={mach:g}  ({a.regime})")
    fig.tight_layout()
    fig.savefig(case / "polar.png", dpi=130)
    print(case / "polar.png")


def selftest():
    assert parse_aoa("-4:16:2") == [-4 + 2 * i for i in range(11)], parse_aoa("-4:16:2")
    assert parse_aoa("2,4,6") == [2.0, 4.0, 6.0]
    assert parse_aoa("0:0:1") == [0.0]
    ms = parse_range("0.15:0.6:0.15")  # Mach sweep uses the same parser
    assert len(ms) == 4 and abs(ms[0] - 0.15) < 1e-9 and abs(ms[-1] - 0.6) < 1e-9, ms
    assert parse_range("0.3") == [0.3]  # single Mach still parses to a one-element list
    assert auto_farfield(0.1) == 25.0 and auto_farfield(0.5) == 35.0
    assert auto_farfield(0.8) == 50.0
    assert auto_downstream(25.0, 0.1) == 35.0 and auto_downstream(25.0, 0.3) == 75.0
    assert auto_downstream(75.0, 0.5) == 85.0
    assert past_stall(1.20, 1.17, 1, 0.02)
    assert not past_stall(1.20, 1.19, 1, 0.02)
    assert past_stall(-1.20, -1.17, -1, 0.02)
    import tempfile
    p = Path(tempfile.mkdtemp()) / "h.csv"
    p.write_text('"Inner_Iter",       "rms[Rho]",           "CL",           "CD",          "CMz"\n'
                 '0, -1.0, 0.1, 0.02, -0.01\n'
                 '1, -8.5, 0.4412, 0.00931, -0.00123\n')
    h = read_history(p)
    assert abs(h["CL"] - 0.4412) < 1e-9 and abs(h["CD"] - 0.00931) < 1e-9, h
    # Coefficients and the regime-specific residual must all converge.
    assert "CONV_CAUCHY_ELEMS= 100" in COMMON and "CONV_CAUCHY_EPS= 1E-5" in COMMON
    assert "CONV_RESIDUAL_MINVAL= -6" in COMMON
    cfg =make_cfg("inc", 2.0, 1e6, 0.15, 500, True)
    assert "MU_CONSTANT= 6.2722695000e-05" in cfg and "RESTART_SOL= YES" in cfg, cfg
    assert "CONV_FIELD= ( LIFT, DRAG, MOMENT_Z, RMS_PRESSURE )" in cfg
    comp0 = make_cfg("comp", 0.0, 1e6, 0.8, 500, False)
    assert "CONV_FIELD= ( LIFT, DRAG, MOMENT_Z, RMS_DENSITY )" in comp0
    assert "SOLVER= RANS" in comp0
    # the convective scheme is solver-specific: FDS is incompressible-only, ROE
    # compressible-only, and SU2 rejects the wrong one at startup
    assert "CONV_NUM_METHOD_FLOW= FDS" in cfg and "ROE" not in cfg
    comp = make_cfg("comp", 0.0, 1e6, 0.8, 500, False)
    assert "CONV_NUM_METHOD_FLOW= ROE" in comp and "FDS" not in comp
    assert "LOW_MACH_PREC= YES" in comp  # comp is the default; must converge at low Mach too
    # transition rides on top of SA, it never replaces the turbulence model
    assert list(TRANSITION) == ["none", "lm"]  # BCM removed; LM is the only transition model
    for t in TRANSITION:
        c = make_cfg("inc", 0.0, 1e6, 0.15, 500, False, transition=t, tu=0.002)
        assert "KIND_TURB_MODEL= SA" in c and "FREESTREAM_TURBULENCEINTENSITY= 0.002" in c
    assert "KIND_TRANS_MODEL= LM" in make_cfg("inc", 0.0, 1e6, 0.15, 500, False, "lm")
    # transition must run at a gentler CFL than fully-turbulent SA, or the front
    # hunts and the coefficients stall short of convergence (as on the first LM
    # 0015 sweep). Turbulent keeps the fast ceiling.
    assert "( 0.1, 2.0, 5.0, 75.0 )" in make_cfg("comp", 0.0, 1e6, 0.3, 500, False, "none")
    assert "( 0.1, 2.0, 5.0, 15.0 )" in make_cfg("comp", 0.0, 1e6, 0.3, 500, False, "lm")
    inc_steady = make_cfg("inc", 0.0, 1e6, 0.15, 500, False, "none")
    assert "CFL_NUMBER= 10.0" in inc_steady and "CFL_ADAPT= YES" in inc_steady
    assert "CONV_NUM_METHOD_FLOW= FDS" in inc_steady and "MUSCL_FLOW= YES" in inc_steady
    assert "MUSCL_TURB= NO" in inc_steady and "LINEAR_SOLVER_ITER= 15" in inc_steady
    assert "CFL_NUMBER= 10.0" in make_cfg("inc", 0.0, 1e6, 0.15, 500, False, "lm")
    urans = make_cfg("inc", 14.0, 1e6, 0.15, 500, True, unsteady=True)
    assert "TIME_DOMAIN= YES" in urans and "TIME_MARCHING= DUAL_TIME_STEPPING-2ND_ORDER" in urans
    assert "TIME_ITER= 2000" in urans and "INNER_ITER= 30" in urans
    assert "WINDOW_START_ITER= 800" in urans and "TAVG_AERO_COEFF" in urans
    assert "CONV_FILENAME= history_+14.00_urans" in urans and "ITER= 500" not in urans
    assert "CONV_NUM_METHOD_FLOW= LD2" in urans and "MUSCL_FLOW= NO" in urans
    assert "MUSCL_TURB= NO" in urans and "LINEAR_SOLVER_ITER= 15" in urans
    comp_urans = make_cfg("comp", 14.0, 1e6, 0.15, 500, True, unsteady=True)
    assert "CONV_NUM_METHOD_FLOW= JST" in comp_urans and "CFL_NUMBER= 20.0" in comp_urans
    print("selftest ok")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        main()
