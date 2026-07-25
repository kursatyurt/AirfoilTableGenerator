#!/usr/bin/env bash
# Whole-table rotor sweep. Only the physical table definition is required.
set -euo pipefail
cd "$(dirname "$0")/.."

usage() {
  cat <<'EOF'
Usage:
  AIRFOIL=... MACH="..." RE=... AOA=MIN:MAX:STEP TRANSITION=none|lm \
  bash tools/run_rotor_table.sh

MACH accepts a space-separated sweep, e.g. MACH="0.1 0.3 0.5". Columns below
INC_MAX_MACH (default 0.25) use INC_RANS; higher columns use compressible RANS.
YPLUS defaults to 1 and OUTROOT defaults to runs/<airfoil>_Re<re>. polar.py
calibrates and caches the farfield once per airfoil/Mach.

Optional overrides: INC_MAX_MACH, YPLUS, OUTROOT, ITERS, TU, NP, SLOTS,
STALL_DROP, URANS_STEPS_PER_CHORD, URANS_CONVECTIVE_TIMES, URANS_INNER_ITERS,
FARFIELD, REGIME. MACHS is accepted as a backwards-compatible alias for MACH.
EOF
}

require() { [ -n "${!1:-}" ] || { echo "missing required environment variable: $1" >&2; usage >&2; exit 2; }; }
if [ "${1:-}" = "--help" ]; then usage; exit 0; fi
MACH=${MACH:-${MACHS:-}}
for v in AIRFOIL MACH RE AOA TRANSITION; do require "$v"; done
INC_MAX_MACH=${INC_MAX_MACH:-0.25}
YPLUS=${YPLUS:-1}
OUTROOT=${OUTROOT:-"runs/${AIRFOIL}_Re${RE}"}
SLOTS=${SLOTS:-1}
if [ -n "${REGIME:-}" ]; then
  case "$REGIME" in inc|comp) ;; *) echo "REGIME must be inc or comp" >&2; exit 2;; esac
fi

source env.sh
if [ -n "${NP:-}" ]; then
  cores=$(python -c 'import os; print(os.cpu_count() or 1)')
  if (( NP * SLOTS > cores )); then
    echo "NP * SLOTS exceeds available logical cores ($NP * $SLOTS > $cores)" >&2
    exit 2
  fi
fi
mkdir -p "$OUTROOT"

run_column() {  # $1 = Mach
  local m="$1" regime nn out log
  local -a farfield_arg=()
  local -a optional=()
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
  [ -n "${ITERS:-}" ] && optional+=(--iters "$ITERS")
  [ -n "${TU:-}" ] && optional+=(--tu "$TU")
  [ -n "${NP:-}" ] && optional+=(--np "$NP")
  [ -n "${STALL_DROP:-}" ] && optional+=(--stall-drop "$STALL_DROP")
  [ -n "${URANS_STEPS_PER_CHORD:-}" ] && optional+=(--urans-steps-per-chord "$URANS_STEPS_PER_CHORD")
  [ -n "${URANS_CONVECTIVE_TIMES:-}" ] && optional+=(--urans-convective-times "$URANS_CONVECTIVE_TIMES")
  [ -n "${URANS_INNER_ITERS:-}" ] && optional+=(--urans-inner-iters "$URANS_INNER_ITERS")
  nn=$(python -c "print(f'{int(round($m*100)):03d}')")
  out="$OUTROOT/${AIRFOIL}_m${nn}"
  log="$OUTROOT/${AIRFOIL}_m${nn}.log"
  {
    echo "=== $out M$m Re$RE $regime $TRANSITION $(date) ==="
    python polar.py --airfoil "$AIRFOIL" --mach "$m" --re "$RE" --aoa "$AOA" \
      --regime "$regime" --yplus "$YPLUS" "${farfield_arg[@]}" "${optional[@]}" \
      --transition "$TRANSITION" --outdir "$out"
    echo "=== $out DONE $(date) ==="
  } &>> "$log"
}

echo "sweep: $AIRFOIL columns[$MACH] Re=$RE regime=${REGIME:-auto@$INC_MAX_MACH} slots=$SLOTS"
pids=()
for m in $MACH; do
  run_column "$m" &
  pids+=("$!")
  if [ "${#pids[@]}" -ge "$SLOTS" ]; then
    wait "${pids[0]}"
    pids=("${pids[@]:1}")
  fi
done
wait
echo "all columns done: $AIRFOIL $MACH"
