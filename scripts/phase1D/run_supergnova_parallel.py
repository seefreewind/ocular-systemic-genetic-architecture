#!/usr/bin/env python3
"""Run independent Phase 1D SUPERGNOVA pairs with bounded parallelism.

Each pair still invokes the unchanged official pipeline through the Phase 1C
compatibility wrapper.  Two pair-level jobs are allowed concurrently to keep
the 10-core/16-GB workstation responsive while reducing wall-clock time.
"""

from __future__ import annotations

import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_supergnova_atlas import (  # noqa: E402
    LOG_DIR,
    PYTHON,
    RAW_OUT,
    RUNNER,
    SOURCE,
    REF,
    SUMSTATS,
    PAIRS,
    read_input_qc,
    run_capture,
)


def run_pair(item: Tuple[str, str, str]) -> str:
    pair, trait1, trait2 = item
    raw = RAW_OUT / f"{pair}.txt"
    overlap = RAW_OUT / f"{pair}.overlap.tsv"
    log = LOG_DIR / f"{pair}.log"
    if raw.exists() and overlap.exists() and raw.stat().st_size > 1000:
        return f"SKIP {pair} existing raw output"
    qc = read_input_qc()
    run_capture(pair, trait1, trait2, qc, raw, overlap, log)
    return f"DONE {pair}"


def main() -> int:
    RAW_OUT.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    pending = list(PAIRS)
    print(f"Parallel pair jobs: {len(pending)}; max_workers=2", flush=True)
    failures: List[str] = []
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(run_pair, item): item[0] for item in pending}
        for future in as_completed(futures):
            pair = futures[future]
            try:
                print(future.result(), flush=True)
            except Exception as exc:
                failures.append(pair)
                print(f"FAIL {pair}: {exc}", file=sys.stderr, flush=True)
    if failures:
        print("Unfinished pairs: " + ", ".join(failures), file=sys.stderr)
        return 2
    summarize = ROOT / "scripts" / "phase1D" / "summarize_supergnova_atlas.py"
    result = subprocess.run([str(PYTHON), str(summarize)], cwd=str(ROOT), text=True)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
