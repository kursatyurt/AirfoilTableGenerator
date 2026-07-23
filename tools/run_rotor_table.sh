#!/usr/bin/env bash
# Whole-table CFD sweep for a rotor section: one SU2 lm polar per Mach column,
# each at its rotor-coupled Reynolds number, then build_c81.py + Viterna assemble
# the +-180 C81 table. Generic over airfoil and blade chord.
#
#   Re(M) = M * SOUND * CHORD / NU     (fixed-RPM rotor: Re and Mach are coupled)
#
# Configure via env vars (defaults = VR12 on the UT Austin rotor):
#   AIRFOIL=vr12  CHORD=0.08  MACHS="0.1 0.2 0.3 0.4 0.5 0.6"
#   SOUND=340.3  NU=1.46e-5  AOA=-14:20:1  ITERS=10000  TU=0.001
#
# Mach columns are independent, so on a big box run them concurrently:
#   CONCURRENT=1  -> background every column and wait (NP defaults to 12, so
#                    6 columns fill ~72 cores; warm-starts stay intact within a
#                    column). A single 2D solve stops scaling past ~1-2 dozen
#                    ranks, so throughput comes from concurrent columns, not
#                    bigger NP -- see the machine.conf study in the README.
#   default (CONCURRENT=0) -> one column at a time (NP defaults to 4); safe for
#                    a laptop or a background-job wall-clock budget.
# Override NP to change ranks-per-column.  Usage:  bash tools/run_rotor_table.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source env.sh

AIRFOIL="${AIRFOIL:-vr12}"
CHORD="${CHORD:-0.08}"
MACHS="${MACHS:-0.1 0.2 0.3 0.4 0.5 0.6}"
SOUND="${SOUND:-340.3}"
NU="${NU:-1.46e-5}"
AOA="${AOA:--14:20:1}"
ITERS="${ITERS:-10000}"
TU="${TU:-0.001}"
CONCURRENT="${CONCURRENT:-0}"
NP="${NP:-$([ "$CONCURRENT" = 1 ] && echo 12 || echo 4)}"

run_column() {  # $1 = Mach
  local M="$1" RE REGIME FF NN OUT LOG
  RE=$(python -c "print(f'{$M*$SOUND*$CHORD/$NU:.4g}')")
  # inc for M<0.3, compressible otherwise; transonic wants a bigger farfield
  if python -c "import sys; sys.exit(0 if $M<0.3 else 1)"; then
    REGIME=inc;  FF=15
  else
    REGIME=comp; FF=25
  fi
  NN=$(python -c "print(f'{int(round($M*100)):03d}')")
  OUT="runs/${AIRFOIL}_m${NN}"
  LOG="runs/${AIRFOIL}_m${NN}.log"
  {
    echo "=== $OUT  M$M Re$RE $REGIME np$NP  $(date) ==="
    python polar.py --airfoil "$AIRFOIL" --mach "$M" --re "$RE" --aoa "$AOA" \
      --transition lm --regime "$REGIME" --farfield "$FF" --np "$NP" \
      --iters "$ITERS" --tu "$TU" --outdir "$OUT"
    echo "=== $OUT DONE $(date) exit=$? ==="
  } &>> "$LOG"
}

for M in $MACHS; do
  if [ "$CONCURRENT" = 1 ]; then
    run_column "$M" &
  else
    run_column "$M"
  fi
done
[ "$CONCURRENT" = 1 ] && wait || true
echo "all columns done ($AIRFOIL): $MACHS"
