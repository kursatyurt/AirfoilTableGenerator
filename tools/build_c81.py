#!/usr/bin/env python3
"""Assemble a whole C81 airfoil table from SU2 polar.csv columns, any airfoil.

Reads every <case>/<airfoil>_mNNN/polar.csv (NNN = round(Mach*100), the layout
run_rotor_table.sh writes under runs/<airfoil>_c<chord> or runs/<airfoil>_Re<re>),
keeps the converged CFD points, and Viterna-extrapolates each column to the full
+-180 deg range so the C81 interpolator is always in bounds. Emits CL, CD and CM
on a common alpha grid x the Mach grid.

The Mach grid gains a 0.0 column: the incompressible-limit copy of the lowest
Mach run (the polar is ~Mach-independent as M->0, and inboard low-speed sections
carry negligible power) -- required because queries at M->0 inboard would
otherwise be out of bounds.

Output format follows the C81 specification:
  Line 1: 30-char padded name, then six space-separated counts
  Lines 2+: CL block, CD block, CM block (each with Mach row + alpha rows)

Usage: python tools/build_c81.py --airfoil VR12 --case runs/VR12_c0.08 \
           --name VR_12 --cdmax 2.05
       writes <case>/<name>.c81 and prints a summary.  Every option is required;
       nothing is guessed.
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

# Common output alpha grid: fine 1-deg core over the CFD range, coarser to +-180.
CORE = list(range(-20, 21))
WING = [22, 24, 27, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
        150, 160, 170, 180]
OUT_ALPHA = sorted(set([-a for a in WING] + CORE + WING))   # -180..180


def resolve_case(airfoil, case):
    """Return (case_dir, airfoil_as_spelled_on_disk).

    run_rotor_table.sh derives OUTROOT as runs/<AIRFOIL>_c<CHORD> (or _Re<RE>)
    and writes <OUTROOT>/<AIRFOIL>_mNNN.  The airfoil is spelled inconsistently
    across the runs (n0012, NACA0015, VR12), so match the Mach dirs
    case-insensitively and then use the on-disk spelling.
    """
    case_dir = case if os.path.isabs(case) else os.path.join(HERE, case)
    if not os.path.isdir(case_dir):
        sys.exit(f"no such case directory: {case_dir}")
    for d in sorted(glob.glob(os.path.join(case_dir, "*_m*"))):
        stem = re.match(r"(.+)_m\d+$", os.path.basename(d))
        if os.path.isdir(d) and stem and stem.group(1).lower() == airfoil.lower():
            return case_dir, stem.group(1)
    sys.exit(f"no {airfoil}_mNNN column directories under {case_dir}")


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


def format_c81_block(cols, machs, alphas):
    """Format a single coefficient block (CL, CD, or CM) in C81 format.

    Returns a list of lines:
      - First line: Mach values
      - Subsequent lines: alpha value followed by coefficients at each Mach
    """
    lines = []
    # Mach row
    lines.append("  ".join(f"{m:.6f}" for m in machs))
    # Alpha rows: each row is alpha followed by values at each Mach
    for i, alpha in enumerate(alphas):
        row_vals = [f"{alpha:.6f}"] + [f"{cols[m][i]:.6f}" for m in machs]
        lines.append("  ".join(row_vals))
    return lines


def main():
    ap = argparse.ArgumentParser(
        description="Assemble a C81 airfoil table from the polar.csv files of a "
                    "finished run_rotor_table.sh sweep.",
        epilog="""\
expected input layout (what run_rotor_table.sh produces):

  runs/VR12_c0.08/            <- pass this path as --case
    VR12_m030/polar.csv       <- one directory per Mach, NNN = round(Mach*100)
    VR12_m050/polar.csv          the "VR12" part is what you pass as --airfoil
    VR12_m065/polar.csv

output:

  runs/VR12_c0.08/VR_12.c81   <- <case>/<name>.c81, i.e. written *into* --case

example:

  python tools/build_c81.py \\
      --case    runs/VR12_c0.08 \\
      --airfoil VR12 \\
      --name    VR_12 \\
      --cdmax   2.05

Every option is required; nothing is guessed.""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, metavar="DIR",
                    help="INPUT DIRECTORY: the sweep root holding the "
                         "<airfoil>_mNNN/polar.csv subdirectories "
                         "(e.g. runs/VR12_c0.08). Relative paths are resolved "
                         "against the repo root. The .c81 file is written here too.")
    ap.add_argument("--airfoil", required=True, metavar="NAME",
                    help="airfoil name as spelled in the <airfoil>_mNNN subdirectory "
                         "names, e.g. VR12, n0012, NACA0015 (case-insensitive)")
    ap.add_argument("--name", required=True, metavar="NAME",
                    help="name of the table itself: goes in the C81 header and "
                         "becomes the output filename <case>/<name>.c81 "
                         "(max 30 chars, e.g. VR_12)")
    ap.add_argument("--cdmax", type=float, required=True, metavar="CD",
                    help="flat-plate CD at 90 deg used by the Viterna "
                         "extrapolation to +-180 deg; 1.11+0.018*AR for a finite "
                         "wing, ~2.0 for a 2D section (e.g. 2.05)")
    args = ap.parse_args()
    case_dir, airfoil = resolve_case(args.airfoil, args.case)
    table_name = args.name

    data = {}
    for d in sorted(glob.glob(os.path.join(case_dir, f"{airfoil}_m*"))):
        m = re.match(rf"{re.escape(airfoil)}_m(\d+)$", os.path.basename(d))
        pol = os.path.join(d, "polar.csv")
        if not m or not os.path.exists(pol):
            continue
        a, cl, cd, cmv = load(pol)
        if a.size >= 4:
            data[int(m.group(1)) / 100.0] = (a, cl, cd, cmv)
    if not data:
        sys.exit(f"no converged polar.csv under {case_dir}/{airfoil}_mNNN")

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

    # Build C81 file content
    n_mach = len(machs)
    n_alpha = len(OUT_ALPHA)

    # Line 1: 30-char padded name + six counts (all blocks use same grid)
    name_field = table_name[:30].ljust(30)
    counts = f"{n_mach} {n_alpha} {n_mach} {n_alpha} {n_mach} {n_alpha}"
    header_line = name_field + counts

    # Format the three coefficient blocks
    cl_lines = format_c81_block(CL, machs, OUT_ALPHA)
    cd_lines = format_c81_block(CD, machs, OUT_ALPHA)
    cm_lines = format_c81_block(CM, machs, OUT_ALPHA)

    # Write C81 file
    out = os.path.join(case_dir, f"{table_name}.c81")
    with open(out, "w") as f:
        f.write(header_line + "\n")
        f.write("# CL block\n")
        f.write("\n".join(cl_lines) + "\n")
        f.write("# CD block\n")
        f.write("\n".join(cd_lines) + "\n")
        f.write("# CM block\n")
        f.write("\n".join(cm_lines) + "\n")

    print(f"{table_name}: {len(OUT_ALPHA)} alpha x {len(machs)} mach")
    print(f"Mach values: {machs}")
    print(f"min CD@0 = {min_cd0:.4f}   CL slope@0 = {slope:.3f}/deg   cdmax = {args.cdmax}")
    print(f"wrote {out}")
    if warns:
        print(f"\n{len(warns)} warning(s) -- review before trusting the table:", file=sys.stderr)
        for w in warns:
            print(f"  ! {w}", file=sys.stderr)


if __name__ == "__main__":
    main()
