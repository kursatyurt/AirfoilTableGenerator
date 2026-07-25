#!/usr/bin/env bash
# Whole-table rotor sweep. Numerical/mesh choices are automatic unless overridden.
set -euo pipefail
cd "$(dirname "$0")/.."

usage() {
  cat <<'EOF'
Usage:
  AIRFOIL=... CHORD=... MACHS="..." SOUND=... NU=... AOA=... ITERS=... \
  TU=... TRANSITION=none|lm INC_MAX_MACH=... YPLUS=... NP=... \
  SLOTS=... OUTROOT=... STALL_DROP=... URANS_STEPS_PER_CHORD=... \
  URANS_CONVECTIVE_TIMES=... URANS_INNER_ITERS=... bash tools/run_rotor_table.sh

Re(M) is calculated as M * SOUND * CHORD / NU. Columns below INC_MAX_MACH use
INC_RANS; higher columns use compressible RANS. polar.py calibrates and caches
the farfield once per airfoil/Mach. Set optional FARFIELD=... to use a fixed
radius for a deliberate mesh-domain study; REGIME=inc|comp overrides regime selection.
EOF
}

require() { [ -n "${!1:-}" ] || { echo "missing required environment variable: $1" >&2; usage >&2; exit 2; }; }
if [ "${1:-}" = "--help" ]; then usage; exit 0; fi
for v in AIRFOIL CHORD MACHS SOUND NU AOA ITERS TU TRANSITION INC_MAX_MACH YPLUS NP SLOTS OUTROOT \
         STALL_DROP URANS_STEPS_PER_CHORD URANS_CONVECTIVE_TIMES URANS_INNER_ITERS; do require "$v"; done
if [ -n "${REGIME:-}" ]; then
  case "$REGIME" in inc|comp) ;; *) echo "REGIME must be inc or comp" >&2; exit 2;; esac
fi

source env.sh
cores=$(python -c 'import os; print(os.cpu_count() or 1)')
if (( NP * SLOTS > cores )); then
  echo "NP * SLOTS exceeds available logical cores ($NP * $SLOTS > $cores)" >&2
  exit 2
fi
mkdir -p "$OUTROOT"

run_column() {  # $1 = Mach
  local m="$1" re regime nn out log
  local -a farfield_arg=()
  re=$(python -c "print(f'{$m*$SOUND*$CHORD/$NU:.4g}')")
  if [ -n "${REGIME:-}" ]; then
    regime=$REGIME
  elif python -c "import sys; sys.exit(0 if $m < $INC_MAX_MACH else 1)"; then
    regime=inc
  else
    regime=comp
  fi
  if [ -n "${FARFIELD:-}" ]; then
    farfield_arg=(--farfield "$FARFIELD")
  fi
  nn=$(python -c "print(f'{int(round($m*100)):03d}')")
  out="$OUTROOT/${AIRFOIL}_m${nn}"
  log="$OUTROOT/${AIRFOIL}_m${nn}.log"
  {
    echo "=== $out M$m Re$re $regime $TRANSITION np$NP $(date) ==="
    python polar.py --airfoil "$AIRFOIL" --mach "$m" --re "$re" --aoa "$AOA" \
      --regime "$regime" --np "$NP" --iters "$ITERS" --yplus "$YPLUS" "${farfield_arg[@]}" \
      --transition "$TRANSITION" --tu "$TU" --outdir "$out" --stall-drop "$STALL_DROP" \
      --urans-steps-per-chord "$URANS_STEPS_PER_CHORD" \
      --urans-convective-times "$URANS_CONVECTIVE_TIMES" \
      --urans-inner-iters "$URANS_INNER_ITERS"
    echo "=== $out DONE $(date) ==="
  } &>> "$log"
}

echo "sweep: $AIRFOIL columns[$MACHS] regime=${REGIME:-auto@$INC_MAX_MACH} np=$NP slots=$SLOTS"
pids=()
for m in $MACHS; do
  run_column "$m" &
  pids+=("$!")
  if [ "${#pids[@]}" -ge "$SLOTS" ]; then
    wait "${pids[0]}"
    pids=("${pids[@]:1}")
  fi
done
wait
echo "all columns done: $AIRFOIL $MACHS"
