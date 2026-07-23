#!/usr/bin/env bash
# Whole-table CFD sweep for a rotor section: one SU2 lm polar per Mach column,
# each at its rotor-coupled Reynolds number, then build_c81.py + Viterna assemble
# the +-180 C81 table. Generic over airfoil and blade chord.
#
#   Re(M) = M * SOUND * CHORD / NU     (fixed-RPM rotor: Re and Mach are coupled)
#
# Configure via env vars (defaults = VR12 on the UT Austin rotor):
#   AIRFOIL=vr12  CHORD=0.08  MACHS="0.1 0.2 0.3 0.4 0.5 0.6"
#   SOUND=340.3  NU=1.46e-5  AOA=-14:20:1  ITERS=10000  NP=4  TU=0.001
#
# Runs one Mach column at a time (each fits the background wall-clock budget;
# >~46 solves/job has been killed). Usage:  bash tools/run_rotor_table.sh
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
NP="${NP:-4}"
TU="${TU:-0.001}"

for M in $MACHS; do
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
    echo "=== $OUT  M$M Re$RE $REGIME  $(date) ==="
    python polar.py --airfoil "$AIRFOIL" --mach "$M" --re "$RE" --aoa "$AOA" \
      --transition lm --regime "$REGIME" --farfield "$FF" --np "$NP" \
      --iters "$ITERS" --tu "$TU" --outdir "$OUT"
    echo "=== $OUT DONE $(date) exit=$? ==="
  } &>> "$LOG"
done
