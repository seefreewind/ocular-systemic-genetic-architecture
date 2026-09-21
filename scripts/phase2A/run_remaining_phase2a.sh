#!/bin/zsh
set -u

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PAIR_RUNS="$ROOT/results/phase2A/pair_runs"
LOG_DIR="$ROOT/logs/phase2A"
mkdir -p "$LOG_DIR"

is_valid_pair() {
  local pair="$1"
  local placo="$PAIR_RUNS/${pair}.placo.tsv.gz"
  local extreme="$PAIR_RUNS/${pair}.extreme.tsv.gz"
  local signqc="$PAIR_RUNS/${pair}.sign_qc.tsv"
  local params="$PAIR_RUNS/${pair}.params.tsv"
  [[ -s "$placo" && -s "$extreme" && -s "$signqc" && -s "$params" ]] || return 1
  gzip -t "$placo" "$extreme" >/dev/null 2>&1
}

run_one() {
  local pair="$1"
  if is_valid_pair "$pair"; then
    print -u2 "$pair already valid; skipping"
    return 0
  fi
  PHASE2A_THREADS=4 "$ROOT/scripts/phase2A/run_pair_if_needed.sh" "$pair" >"$LOG_DIR/${pair}.log" 2>&1
}

# AMD-CKD and AMD-AD were already launched by the interactive task. Wait for
# both to finish before adding new pairs, keeping the total at two pair jobs.
while ! is_valid_pair AMD-CKD || ! is_valid_pair AMD-AD; do
  sleep 60
done

remaining=(
  AMD-PD
  POAG-CAD POAG-STROKE POAG-T2D POAG-CKD POAG-AD POAG-PD
  CAT-CAD CAT-STROKE CAT-T2D CAT-CKD CAT-AD CAT-PD
)

for ((i = 1; i <= ${#remaining}; i += 2)); do
  p1="${remaining[$i]}"
  p2=""
  if (( i + 1 <= ${#remaining} )); then p2="${remaining[$((i + 1))]}"; fi
  run_one "$p1" & pid1=$!
  if [[ -n "$p2" ]]; then run_one "$p2" & pid2=$!; else pid2=""; fi
  wait "$pid1" || print -u2 "PAIR_FAILED $p1"
  if [[ -n "$pid2" ]]; then wait "$pid2" || print -u2 "PAIR_FAILED $p2"; fi
done

all_pairs=(
  AMD-CAD AMD-STROKE AMD-T2D AMD-CKD AMD-AD AMD-PD
  POAG-CAD POAG-STROKE POAG-T2D POAG-CKD POAG-AD POAG-PD
  CAT-CAD CAT-STROKE CAT-T2D CAT-CKD CAT-AD CAT-PD
)
for pair in "${all_pairs[@]}"; do
  if ! is_valid_pair "$pair"; then
    print -u2 "MISSING_OR_INVALID_PAIR $pair"
    exit 1
  fi
done

# Remove only macOS AppleDouble metadata files created in this result folder.
find "$PAIR_RUNS" -maxdepth 1 -type f -name '._*' -delete

python3 "$ROOT/scripts/phase2A/aggregate_phase2a.py" >"$LOG_DIR/aggregate_phase2a.log" 2>&1
print "PHASE2A_AGGREGATION_COMPLETE"
