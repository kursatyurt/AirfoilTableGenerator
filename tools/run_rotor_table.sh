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
#   CONCURRENT=1  -> keep SLOTS columns running at once (SLOTS from machine.conf
#                    = cores/NP, so the box fills without oversubscribing on any
#                    machine). Warm-starts stay intact within a column. A single
#                    2D solve stops scaling past ~1-2 dozen ranks, so throughput
#                    comes from concurrent columns, not bigger NP.
#   default (CONCURRENT=0) -> one column at a time.
# NP / CORES / SLOTS default to the tuned machine.conf (run tune_np.py once);
# any of NP, SLOTS, CORES can be overridden by env.  Usage: bash tools/run_rotor_table.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source env.sh

conf() { [ -f machine.conf ] && grep -E "^$1=" machine.conf | tail -1 | cut -d= -f2 | tr -d ' '; }

AIRFOIL="${AIRFOIL:-vr12}"
CHORD="${CHORD:-0.08}"
MACHS="${MACHS:-0.1 0.2 0.3 0.4 0.5 0.6}"
SOUND="${SOUND:-340.3}"
NU="${NU:-1.46e-5}"
AOA="${AOA:--14:20:1}"
ITERS="${ITERS:-10000}"
TU="${TU:-0.001}"
CONCURRENT="${CONCURRENT:-0}"
NP="${NP:-$(conf NP)}"; NP="${NP:-4}"                 # tuned ranks/column, else 4
if [ "$CONCURRENT" = 1 ]; then
  SLOTS="${SLOTS:-$(conf SLOTS)}"; SLOTS="${SLOTS:-1}" # tuned concurrent columns
else
  SLOTS=1
fi

run_column() {  # $1 = Mach
  local M="$1" RE REGIME FF NN OUT LOG
  RE=$(python -c "print(f'{$M*$SOUND*$CHORD/$NU:.4g}')")
  # compressible everywhere (Roe + low-Mach preconditioning) so every column of the
  # C81 table is built with one consistent solver; transonic just wants a bigger farfield
  REGIME=comp
  if python -c "import sys; sys.exit(0 if $M<0.7 else 1)"; then FF=15; else FF=25; fi
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

echo "sweep: $AIRFOIL  columns[$MACHS]  np=$NP  slots=$SLOTS"
# ponytail: FIFO window (wait on the oldest), not a greedy wait -n scheduler --
# columns are similar length, and this stays portable to bash 3.2 (macOS).
pids=()
for M in $MACHS; do
  run_column "$M" &
  pids+=("$!")
  if [ "${#pids[@]}" -ge "$SLOTS" ]; then
    wait "${pids[0]}"
    pids=("${pids[@]:1}")
  fi
done
wait
echo "all columns done ($AIRFOIL): $MACHS"
