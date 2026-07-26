#!/usr/bin/env bash
# Whole-table rotor sweep. Only the physical table definition is required.
set -euo pipefail
cd "$(dirname "$0")/.."

usage() {
  cat <<'EOF'
Usage:
  AIRFOIL=... MACH="..." CHORD=... AOA=MIN:MAX:STEP TRANSITION=none|lm \
  bash tools/run_rotor_table.sh

MACH accepts a space-separated sweep, e.g. MACH="0.1 0.3 0.5". Columns below
INC_MAX_MACH (default 0.25) use INC_RANS; higher columns use compressible RANS.
With CHORD, Re(M) is derived using standard sea-level SOUND=341.348 m/s and
NU=1.4607e-5 m²/s. Set RE=... instead to use one explicit Reynolds number.
YPLUS defaults to 1 and OUTROOT is derived from the supplied airfoil/chord or Re.
polar.py uses the validated Mach-based farfield formula (25c minimum; 35c at M=0.5).

Optional overrides: RE, SOUND, NU, INC_MAX_MACH, YPLUS, OUTROOT, ITERS, TU, NP, SLOTS,
STALL_DROP, URANS_STEPS_PER_CHORD, URANS_CONVECTIVE_TIMES, URANS_INNER_ITERS,
FARFIELD, REGIME, STATUS_INTERVAL. MACHS is accepted as a backwards-compatible alias for MACH.
EOF
}

require() { [ -n "${!1:-}" ] || { echo "missing required environment variable: $1" >&2; usage >&2; exit 2; }; }
if [ "${1:-}" = "--help" ]; then usage; exit 0; fi
MACH=${MACH:-${MACHS:-}}
for v in AIRFOIL MACH AOA TRANSITION; do require "$v"; done
if [ -z "${RE:-}" ] && [ -z "${CHORD:-}" ]; then
  echo "provide RE or CHORD" >&2
  usage >&2
  exit 2
fi
INC_MAX_MACH=${INC_MAX_MACH:-0.25}
YPLUS=${YPLUS:-1}
SOUND=${SOUND:-341.348}
NU=${NU:-1.4607e-5}
if [ -z "${OUTROOT:-}" ]; then
  OUTROOT="runs/${AIRFOIL}_$(if [ -n "${RE:-}" ]; then echo "Re${RE}"; else echo "c${CHORD}"; fi)"
fi
SLOTS=${SLOTS:-1}
STATUS_INTERVAL=${STATUS_INTERVAL:-15}
if [ -n "${REGIME:-}" ]; then
  case "$REGIME" in inc|comp) ;; *) echo "REGIME must be inc or comp" >&2; exit 2;; esac
fi

source env.sh
cores=$(python -c 'import os; print(os.cpu_count() or 1)')
NP=${NP:-$(python -c 'import os; from tune_np import stored_np; print(stored_np() or max(1, (os.cpu_count() or 1) // 2))')}
mkdir -p "$OUTROOT"

run_column() {  # $1 = Mach
  local m="$1" re regime nn out log
  local -a farfield_arg=()
  local -a optional=()
  if [ -n "${RE:-}" ]; then re=$RE; else re=$(python -c "print(f'{$m*$SOUND*$CHORD/$NU:.6g}')"); fi
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
    echo "=== $out M$m Re$re $regime $TRANSITION $(date) ==="
    python -u polar.py --airfoil "$AIRFOIL" --mach "$m" --re "$re" --aoa "$AOA" \
      --regime "$regime" --yplus "$YPLUS" "${farfield_arg[@]}" "${optional[@]}" \
      --transition "$TRANSITION" --outdir "$out"
    echo "=== $out DONE $(date) ==="
  } &>> "$log"
}

fmt_time() {
  local seconds=$1
  printf '%02dh:%02dm:%02ds' "$((seconds / 3600))" "$(((seconds % 3600) / 60))" "$((seconds % 60))"
}

AOA_POINTS=$(python -c "from polar import parse_aoa; print(len(parse_aoa('$AOA')))")
TOTAL_COLUMNS=$(wc -w <<<"$MACH" | tr -d ' ')
ACTIVE_SLOTS=$(( SLOTS < TOTAL_COLUMNS ? SLOTS : TOTAL_COLUMNS ))
if (( NP * ACTIVE_SLOTS > cores )); then
  ACTIVE_SLOTS=$(( cores / NP ))
fi
if (( ACTIVE_SLOTS < 1 )); then
  echo "NP=$NP exceeds available logical cores ($cores)" >&2
  exit 2
fi
table_started=$(date +%s)
pids=()
job_machs=()
job_logs=()
job_outs=()
job_starts=()
failed=0
finished=0

terminate_tree() {
  local pid=$1 child
  # mpirun and its ranks are descendants of the tracked column wrapper.  Kill
  # children first so an interrupted table never leaves orphaned solver ranks.
  for child in $(pgrep -P "$pid" 2>/dev/null || true); do
    terminate_tree "$child"
  done
  kill -TERM "$pid" 2>/dev/null || true
}

cleanup_jobs() {
  local reason=$1 pid
  trap - EXIT INT TERM HUP
  (( ${#pids[@]} )) || return 0
  echo "[$(date '+%H:%M:%S')] stopping ${#pids[@]} launched column job(s): $reason" >&2
  for pid in "${pids[@]}"; do terminate_tree "$pid"; done
  sleep 1
  for pid in "${pids[@]}"; do
    if kill -0 "$pid" 2>/dev/null; then kill -KILL "$pid" 2>/dev/null || true; fi
    wait "$pid" 2>/dev/null || true
  done
}

interrupted() {
  local signal=$1
  cleanup_jobs "received $signal"
  exit 1
}

trap 'interrupted INT' INT
trap 'interrupted TERM' TERM
trap 'interrupted HUP' HUP
trap 'cleanup_jobs "script exit"' EXIT

print_status() {
  local now elapsed i done eta current last solver
  now=$(date +%s)
  elapsed=$((now - table_started))
  echo "[$(date '+%H:%M:%S')] table: ${#pids[@]} running, $finished/$TOTAL_COLUMNS complete, elapsed $(fmt_time "$elapsed")"
  for i in "${!pids[@]}"; do
    done=$(awk '/CL=/{count++} END{print count+0}' "${job_logs[$i]}" 2>/dev/null || true)
    current=$(awk '/--- M/{line=$0} END{if (line) {sub(/^.*--- M /, "", line); sub(/ \([0-9]+ ranks\).*/, "", line); print line}}' "${job_logs[$i]}" 2>/dev/null || true)
    last=$(awk '/mesh calibration|CL=|stall detected|FAILED|URANS average/{line=$0} END{print line}' "${job_logs[$i]}" 2>/dev/null || true)
    if (( done > 0 )); then
      eta=$(fmt_time "$(( (now - ${job_starts[$i]}) * (AOA_POINTS - done) / done ))")
    else
      eta="estimating"
    fi
    printf '  M=%s  AoA %s  %s/%s complete  elapsed %s  ETA %s\n' \
      "${job_machs[$i]}" "${current:-starting}" "$done" "$AOA_POINTS" \
      "$(fmt_time "$((now - ${job_starts[$i]}))")" "$eta"
    if [ -n "$last" ]; then echo "    $last"; fi
    solver=$(latest_solver_status "${job_outs[$i]}" || true)
    if [ -n "$solver" ]; then echo "    $solver"; fi
  done
}

latest_solver_status() {
  python - "$1" <<'PY'
from pathlib import Path
import csv, sys

files = sorted(Path(sys.argv[1]).glob("history_*.csv"), key=lambda p: p.stat().st_mtime)
if not files:
    raise SystemExit
try:
    with files[-1].open(newline="") as f:
        reader = csv.DictReader(f)
        row = None
        for row in reader:
            pass
    if not row:
        raise SystemExit
    row = {k.strip().strip('"'): v.strip() for k, v in row.items() if k}
    iteration = row.get("Time_Iter") or row.get("Inner_Iter") or row.get("Outer_Iter") or "?"
    residual = next((f"{k}={float(v):.2e}" for k, v in row.items() if k.startswith("rms[") and v), "RMS=?")
    coeff = " ".join(f"{k}={float(row[k]):.5f}" for k in ("CL", "CD", "CMz") if row.get(k))
    print(f"solver {files[-1].stem}: iter={iteration} {residual} {coeff}")
except Exception:
    pass
PY
}

reap_finished() {
  local i rc
  local -a next_pids=() next_machs=() next_logs=() next_outs=() next_starts=()
  for i in "${!pids[@]}"; do
    if kill -0 "${pids[$i]}" 2>/dev/null; then
      next_pids+=("${pids[$i]}")
      next_machs+=("${job_machs[$i]}")
      next_logs+=("${job_logs[$i]}")
      next_outs+=("${job_outs[$i]}")
      next_starts+=("${job_starts[$i]}")
    else
      if wait "${pids[$i]}"; then
        finished=$((finished + 1))
        echo "[$(date '+%H:%M:%S')] M=${job_machs[$i]} completed in $(fmt_time "$(( $(date +%s) - ${job_starts[$i]} ))")"
      else
        echo "[$(date '+%H:%M:%S')] M=${job_machs[$i]} FAILED; see ${job_logs[$i]}" >&2
        failed=1
      fi
    fi
  done
  pids=("${next_pids[@]}")
  job_machs=("${next_machs[@]}")
  job_logs=("${next_logs[@]}")
  job_outs=("${next_outs[@]}")
  job_starts=("${next_starts[@]}")
}

wait_for_slot() {
  while (( ${#pids[@]} >= ACTIVE_SLOTS )); do
    print_status
    sleep "$STATUS_INTERVAL"
    reap_finished
  done
}

echo "sweep: $AIRFOIL columns[$MACH] $(if [ -n "${RE:-}" ]; then echo "Re=$RE"; else echo "chord=$CHORD"; fi) regime=${REGIME:-auto@$INC_MAX_MACH} AoAs=$AOA_POINTS"
echo "resources: $cores logical cores, NP=$NP, requested slots=$SLOTS, active columns=$ACTIVE_SLOTS, active ranks=$((NP * ACTIVE_SLOTS))"
if (( NP * ACTIVE_SLOTS < cores )); then
  echo "resources: $((cores - NP * ACTIVE_SLOTS)) cores idle; only $TOTAL_COLUMNS Mach columns can run concurrently. Raise NP only if tune_np.py proves it improves a single-column solve."
fi
for m in $MACH; do
  wait_for_slot
  run_column "$m" &
  pids+=("$!")
  nn=$(python -c "print(f'{int(round($m*100)):03d}')")
  job_machs+=("$m")
  job_logs+=("$OUTROOT/${AIRFOIL}_m${nn}.log")
  job_outs+=("$OUTROOT/${AIRFOIL}_m${nn}")
  job_starts+=("$(date +%s)")
  echo "[$(date '+%H:%M:%S')] launched M=$m; log: $OUTROOT/${AIRFOIL}_m${nn}.log"
done
while (( ${#pids[@]} )); do
  print_status
  sleep "$STATUS_INTERVAL"
  reap_finished
done
(( failed == 0 )) || exit 1
echo "all columns done: $AIRFOIL $MACH in $(fmt_time "$(( $(date +%s) - table_started ))")"
