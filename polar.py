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
VENKAT_LIMITER_COEFF= 0.05
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= YES
TIME_DISCRE_FLOW= EULER_IMPLICIT
TIME_DISCRE_TURB= EULER_IMPLICIT
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-6
LINEAR_SOLVER_ITER= 20
% CONV_FIELD/SCREEN_OUTPUT are regime-specific (RMS_DENSITY comp, RMS_PRESSURE
% inc) and appended per solver below. LIFT/DRAG/MOMENT_Z are SU2 COEFFICIENT
% fields (CL/CD/CMz), non-dimensional -- not dimensional forces. Cauchy EPS
% 1E-4 = one drag count in CD; RESIDUAL_MINVAL floors the residual so an angle
% cannot report converged while CD is still oscillating on a live residual.
CONV_CAUCHY_ELEMS= 100
CONV_CAUCHY_EPS= 1E-4
CONV_RESIDUAL_MINVAL= -6
CONV_STARTITER= 10
MESH_FILENAME= airfoil.su2
MESH_FORMAT= SU2
SOLUTION_FILENAME= solution_flow.dat
RESTART_FILENAME= restart_flow.dat
TABULAR_FORMAT= CSV
OUTPUT_FILES= ( RESTART, PARAVIEW, SURFACE_CSV )
HISTORY_OUTPUT= ( ITER, RMS_RES, AERO_COEFF, AOA )
"""


TRANSITION = {
    # LM rides on top of KIND_TURB_MODEL= SA; it does not replace it
    "none": "",                      # fully turbulent from the leading edge
    "lm": "KIND_TRANS_MODEL= LM\n",  # Langtry-Menter gamma-Re_theta, two extra equations
}


def make_cfg(regime, aoa, re, mach, iters, restart, transition="none", tu=0.001):
    v = mach * A_SOUND
    tag = f"{aoa:+.2f}"
    # Residual field is solver-specific: the compressible solver tracks density,
    # the incompressible pressure. It joins the Cauchy coefficient fields so an
    # angle converges only once CL/CD/CMz are Cauchy-stable AND the residual has
    # floored (CONV_RESIDUAL_MINVAL), killing the CD oscillation.
    res = "RMS_PRESSURE" if regime == "inc" else "RMS_DENSITY"
    # LM's two extra transport equations make the transition front hunt at high
    # CFL and stall convergence. Capping the adaptive CFL ceiling trades wall time
    # for coefficients that actually settle. Fully-turbulent SA keeps the fast ceiling.
    cfl_max = 15.0 if transition != "none" else 50.0
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
CFL_NUMBER= {25.0 if transition == "none" else 10.0}
CONV_NUM_METHOD_FLOW= FDS
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
CFL_NUMBER= 5.0
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.1, 2.0, 5.0, {cfl_max} )
CONV_NUM_METHOD_FLOW= ROE
ENTROPY_FIX_COEFF= 0.05
LOW_MACH_PREC= YES
"""
    return head + COMMON + TRANSITION[transition] + f"""\
CONV_FIELD= ( LIFT, DRAG, MOMENT_Z, {res} )
SCREEN_OUTPUT= ( INNER_ITER, {res}, RMS_MOMENTUM-X, LIFT, DRAG )
FREESTREAM_TURBULENCEINTENSITY= {tu}
RESTART_SOL= {"YES" if restart else "NO"}
ITER= {iters}
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
    ap.add_argument("--farfield", type=float, default=15.0,
                    help="farfield radius in chords; go well past 15 for transonic runs")
    ap.add_argument("--transition", choices=list(TRANSITION), default="none",
                    help="laminar-turbulent transition on top of SA: 'lm' = Langtry-Menter "
                         "gamma-Re_theta (two extra transport equations)")
    ap.add_argument("--tu", type=float, default=0.001,
                    help="freestream turbulence intensity for the transition models "
                         "(0.001 = 0.1%%, a low-turbulence wind tunnel)")
    ap.add_argument("--outdir", default=None)
    # ponytail: argparse reads a leading '-' as an option, so "--aoa -4:16:2" fails.
    argv, rest = [], list(sys.argv[1:])
    while rest:
        v = rest.pop(0)
        argv.append(f"--aoa={rest.pop(0)}" if v == "--aoa" and rest else v)
    a = ap.parse_args(argv)

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

    mesh = case / "airfoil.su2"
    if not mesh.exists():
        from mesh import generate_mesh
        generate_mesh(x, y, a.re, mach, y_plus=a.yplus, path=mesh,
                      inlet_radius=a.farfield, downstream=max(25.0, a.farfield + 10))
    print(f"mesh: {mesh}")

    sol, restart_dat = case / "solution_flow.dat", case / "restart_flow.dat"

    def su2(cfg_name, log_name):
        with open(case / log_name, "w") as log:
            return subprocess.call(["mpirun", "-n", str(a.np), "SU2_CFD", cfg_name],
                                   cwd=case, stdout=log, stderr=subprocess.STDOUT)

    def solve(aoa, restart, transition, suffix="", tu=None):
        """Run one angle; return (CL, CD, CMz, converged) or None. Leaves the
        converged field in solution_flow.dat for the next (warm-started) angle."""
        tag = f"{aoa:+.2f}"
        cfg = case / f"aoa_{tag}{suffix}.cfg"
        cfg.write_text(make_cfg(a.regime, aoa, a.re, mach, a.iters, restart, transition,
                                a.tu if tu is None else tu))
        label = " turbulent seed" if suffix == "_seed" else ""
        print(f"--- M {mach:g} AoA {aoa:g}{label} ({a.np} ranks) ...", end=" ", flush=True)
        rc = su2(cfg.name, f"aoa_{tag}{suffix}.log")
        hist = case / f"history_{tag}.csv"
        if rc != 0 or not hist.exists():
            print(f"FAILED (rc={rc}), see aoa_{tag}{suffix}.log")
            return None
        h = read_history(hist)
        n_iter = sum(1 for _ in open(hist)) - 1  # history rows = iterations run
        converged = n_iter < a.iters
        if converged:                            # only a converged field warm-starts the
            shutil.copy(restart_dat, sol)        # next angle; a diverged one (e.g. post-
                                                 # stall) would poison the whole march
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

    for aoa in up[1:]:               # ascending above the pivot, warm-started
        record(aoa, solve(aoa, True, a.transition))
    if down and pivot_sol.exists():
        shutil.copy(pivot_sol, sol)  # rewind to the pivot field before fanning down
    for aoa in down:                 # descending below the pivot, warm-started
        record(aoa, solve(aoa, True, a.transition))

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
    import tempfile
    p = Path(tempfile.mkdtemp()) / "h.csv"
    p.write_text('"Inner_Iter",       "rms[Rho]",           "CL",           "CD",          "CMz"\n'
                 '0, -1.0, 0.1, 0.02, -0.01\n'
                 '1, -8.5, 0.4412, 0.00931, -0.00123\n')
    h = read_history(p)
    assert abs(h["CL"] - 0.4412) < 1e-9 and abs(h["CD"] - 0.00931) < 1e-9, h
    # convergence needs Cauchy-stable coefficients AND a floored residual, so the
    # residual field joins CONV_FIELD -- but its name is solver-specific.
    assert "CONV_CAUCHY_EPS= 1E-4" in COMMON and "CONV_RESIDUAL_MINVAL= -6" in COMMON
    cfg =make_cfg("inc", 2.0, 1e6, 0.15, 500, True)
    assert "MU_CONSTANT= 6.2722695000e-05" in cfg and "RESTART_SOL= YES" in cfg, cfg
    # inc solver tracks pressure, never density -- RMS_DENSITY would abort SU2
    assert "CONV_FIELD= ( LIFT, DRAG, MOMENT_Z, RMS_PRESSURE )" in cfg
    assert "CONV_FIELD= ( LIFT, DRAG, MOMENT_Z, RMS_DENSITY )" not in cfg
    comp0 = make_cfg("comp", 0.0, 1e6, 0.8, 500, False)
    assert "CONV_FIELD= ( LIFT, DRAG, MOMENT_Z, RMS_DENSITY )" in comp0
    assert "CONV_FIELD= ( LIFT, DRAG, MOMENT_Z, RMS_PRESSURE )" not in comp0
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
    assert "( 0.1, 2.0, 5.0, 50.0 )" in make_cfg("comp", 0.0, 1e6, 0.3, 500, False, "none")
    assert "( 0.1, 2.0, 5.0, 15.0 )" in make_cfg("comp", 0.0, 1e6, 0.3, 500, False, "lm")
    assert "CFL_NUMBER= 25.0" in make_cfg("inc", 0.0, 1e6, 0.15, 500, False, "none")
    assert "CFL_NUMBER= 10.0" in make_cfg("inc", 0.0, 1e6, 0.15, 500, False, "lm")
    print("selftest ok")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        main()
