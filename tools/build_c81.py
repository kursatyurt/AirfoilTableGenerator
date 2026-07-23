#!/usr/bin/env python3
"""Assemble a whole ARCS C81 table from SU2 polar.csv columns, any airfoil.

Reads every runs/<airfoil>_m0NN/polar.csv (NN = Mach*100), keeps the converged
CFD points, and Viterna-extrapolates each column to the full +-180 deg range so
the ARCS C81 interpolator (which THROWS outside the tabulated range) is always
in bounds. Emits CL, CD and CM on a common alpha grid x the Mach grid.

The Mach grid gains a 0.0 column: the incompressible-limit copy of the lowest
Mach run (the polar is ~Mach-independent as M->0, and inboard low-speed sections
carry negligible power) -- required because ARCS queries M->0 inboard and would
otherwise throw.

Usage: python tools/build_c81.py --airfoil vr12 [--func VR_12_su2] [--cdmax 2.05]
       writes runs/<func>.cpp.txt (paste-ready) and prints a summary.
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from viterna import extrapolate_column

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(HERE, "runs")

# Common output alpha grid: fine 1-deg core over the CFD range, coarser to +-180.
CORE = list(range(-20, 21))
WING = [22, 24, 27, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
        150, 160, 170, 180]
OUT_ALPHA = sorted(set([-a for a in WING] + CORE + WING))   # -180..180


def load(path):
    a, cl, cd, cm = [], [], [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            if int(float(r["converged"])) != 1:
                continue
            a.append(float(r["aoa"])); cl.append(float(r["cl"]))
            cd.append(float(r["cd"])); cm.append(float(r["cm"]))
    o = np.argsort(a)
    return (np.array(a)[o], np.array(cl)[o], np.array(cd)[o], np.array(cm)[o])


def fmt(cols, machs):
    rows = []
    for i in range(len(OUT_ALPHA)):
        rows.append("      {" + ", ".join(f"{cols[m][i]:.6f}" for m in machs) + "}")
    return ",\n".join(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--airfoil", default="vr12", help="run-dir prefix: runs/<airfoil>_m0NN")
    ap.add_argument("--func", default=None, help="C81 function name (default <AIRFOIL>_su2)")
    ap.add_argument("--cdmax", type=float, default=2.05,
                    help="flat-plate CD at 90 deg for Viterna (1.11+0.018*AR; ~2 for 2D)")
    args = ap.parse_args()
    func = args.func or f"{args.airfoil.upper()}_su2"

    data = {}
    for d in sorted(glob.glob(os.path.join(RUNS, f"{args.airfoil}_m0*"))):
        m = re.search(rf"{re.escape(args.airfoil)}_m0(\d+)", d)
        pol = os.path.join(d, "polar.csv")
        if not m or not os.path.exists(pol):
            continue
        a, cl, cd, cmv = load(pol)
        if a.size >= 4:
            data[int(m.group(1)) / 100.0] = (a, cl, cd, cmv)
    if not data:
        sys.exit(f"no converged polar.csv under runs/{args.airfoil}_m0*")

    warns = []
    warn = lambda msg: warns.append(msg)

    machs_cfd = sorted(data)
    CL, CD, CM = {}, {}, {}
    for m in machs_cfd:
        CL[m], CD[m], CM[m] = extrapolate_column(*data[m], OUT_ALPHA, args.cdmax,
                                                 warn=lambda s, m=m: warn(f"M{m:g}: {s}"))
    m0 = machs_cfd[0]                       # M=0 incompressible-limit copy
    CL[0.0], CD[0.0], CM[0.0] = CL[m0], CD[m0], CM[m0]
    machs = [0.0] + [m for m in machs_cfd if m > 0.0]

    # Plausibility checks. WARN (don't abort): the tool is airfoil-generic, and a
    # legitimately thick/thin section can sit outside these VR12-tuned bands.
    i0, i90 = OUT_ALPHA.index(0), OUT_ALPHA.index(90)
    op = [m for m in machs if 0.1 <= m <= 0.6] or machs[1:]
    min_cd0 = min(CD[m][i0] for m in op)
    if not 0.006 <= min_cd0 <= 0.016:
        warn(f"min CD@0 = {min_cd0:.4f} outside typical [0.006,0.016]")
    slope = (CL[machs[1]][OUT_ALPHA.index(4)] - CL[machs[1]][OUT_ALPHA.index(-4)]) / 8.0
    if not 0.08 <= slope <= 0.13:
        warn(f"CL slope@0 = {slope:.3f}/deg outside typical [0.08,0.13]")
    for m in machs:                        # Viterna invariants -- these should always hold
        if not (abs(CL[m][0]) < 0.05 and abs(CL[m][-1]) < 0.05):
            warn(f"M{m:g}: CL not ~0 at +-180 (extrapolation broken)")
        if CD[m][i90] <= 1.0:
            warn(f"M{m:g}: CD@90 = {CD[m][i90]:.2f} not flat-plate-like")

    m_lit = ", ".join(f"{m:g}" for m in machs)
    a_lit = ",".join(f"{a:g}" for a in OUT_ALPHA)
    body = lambda tag, C: f"""  static constexpr const std::size_t num_mach_{tag}  = {len(machs)};
  static constexpr const std::size_t num_alpha_{tag} = {len(OUT_ALPHA)};
  static const std::vector<float> vals_mach_{tag}{{{m_lit}}};
  static const std::vector<float> vals_alpha_{tag}{{{a_lit}}};
  static const std::vector<std::vector<float>> coeffs_{tag}{{
{fmt(C, machs)}}};"""
    out = os.path.join(RUNS, f"{func}.cpp.txt")
    with open(out, "w") as f:
        f.write(f"""const C81& {func}() {{
  static std::string const name = "{func}";
{body('CL', CL)}

{body('CD', CD)}

{body('CM', CM)}

  static auto const c81 = C81(name, vals_alpha_CL, vals_mach_CL, coeffs_CL,
                              vals_alpha_CD, vals_mach_CD, coeffs_CD,
                              vals_alpha_CM, vals_mach_CM, coeffs_CM);
  return c81;
}}
""")
    print(f"{func}: columns {machs}   ({len(OUT_ALPHA)} alpha x {len(machs)} mach)")
    print(f"min CD@0 = {min_cd0:.4f}   CL slope@0 = {slope:.3f}/deg   cdmax = {args.cdmax}")
    print(f"wrote {out}")
    if warns:
        print(f"\n{len(warns)} warning(s) -- review before trusting the table:", file=sys.stderr)
        for w in warns:
            print(f"  ! {w}", file=sys.stderr)


if __name__ == "__main__":
    main()
