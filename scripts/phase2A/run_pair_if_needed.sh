#!/bin/zsh
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 PAIR" >&2
  exit 2
fi

PAIR="$1"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
INPUT="$ROOT/results/phase2A/inputs/${PAIR}.tsv.gz"
OUTDIR="$ROOT/results/phase2A/pair_runs"
PLACO="$OUTDIR/${PAIR}.placo.tsv.gz"
EXTREME="$OUTDIR/${PAIR}.extreme.tsv.gz"
SIGNQC="$OUTDIR/${PAIR}.sign_qc.tsv"
PARAMS="$OUTDIR/${PAIR}.params.tsv"

if [[ -s "$PLACO" && -s "$EXTREME" && -s "$SIGNQC" && -s "$PARAMS" ]] && gzip -t "$PLACO" "$EXTREME" >/dev/null 2>&1; then
  echo "$PAIR already has valid pair outputs; skipping"
  exit 0
fi

R_BIN="${PHASE2A_RSCRIPT:-Rscript}"
if [[ ! -x "$R_BIN" ]]; then
  echo "Rscript executable not found: $R_BIN" >&2
  exit 127
fi
exec env PHASE2A_THREADS="${PHASE2A_THREADS:-4}" "$R_BIN" "$ROOT/scripts/phase2A/run_placo_pair.R" "$PAIR" "$INPUT"
