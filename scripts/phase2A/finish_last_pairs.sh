#!/bin/zsh
set -u

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PAIR_RUNS="$ROOT/results/phase2A/pair_runs"
LOG_DIR="$ROOT/logs/phase2A"
mkdir -p "$LOG_DIR"

valid_pair() {
  local pair="$1"
  local p="$PAIR_RUNS/${pair}.placo.tsv.gz"
  local e="$PAIR_RUNS/${pair}.extreme.tsv.gz"
  local s="$PAIR_RUNS/${pair}.sign_qc.tsv"
  local q="$PAIR_RUNS/${pair}.params.tsv"
  [[ -s "$p" && -s "$e" && -s "$s" && -s "$q" ]] || return 1
  gzip -t "$p" "$e" >/dev/null 2>&1
}

if ! valid_pair CAT-CKD; then
  PHASE2A_THREADS=4 "$ROOT/scripts/phase2A/run_pair_if_needed.sh" CAT-CKD >"$LOG_DIR/CAT-CKD.log" 2>&1 || {
    print -u2 "CAT-CKD_FAILED"
    exit 1
  }
fi

if ! PHASE2A_THREADS=4 "$ROOT/scripts/phase2A/run_pair_if_needed.sh" CAT-AD >"$LOG_DIR/CAT-AD.log" 2>&1; then
  print -u2 "CAT-AD_FAILED"
  exit 1
fi
if ! PHASE2A_THREADS=4 "$ROOT/scripts/phase2A/run_pair_if_needed.sh" CAT-PD >"$LOG_DIR/CAT-PD.log" 2>&1; then
  print -u2 "CAT-PD_FAILED"
  exit 1
fi

for pair in AMD-CAD AMD-STROKE AMD-T2D AMD-CKD AMD-AD AMD-PD \
  POAG-CAD POAG-STROKE POAG-T2D POAG-CKD POAG-AD POAG-PD \
  CAT-CAD CAT-STROKE CAT-T2D CAT-CKD CAT-AD CAT-PD; do
  if ! valid_pair "$pair"; then
    print -u2 "MISSING_OR_INVALID_PAIR $pair"
    exit 1
  fi
done

find "$PAIR_RUNS" -maxdepth 1 -type f -name '._*' -delete
python3 "$ROOT/scripts/phase2A/aggregate_phase2a.py" >"$LOG_DIR/aggregate_phase2a.log" 2>&1
print "PHASE2A_AGGREGATION_COMPLETE"
