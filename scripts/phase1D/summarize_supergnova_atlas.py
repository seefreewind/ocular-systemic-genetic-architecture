#!/usr/bin/env python3
"""Summarize completed Phase 1D raw runs and execute the sign-flip QC.

The atlas runner contains the prespecified summary logic.  This wrapper makes
the full-run stage idempotent: existing full-pair raw outputs are read as-is,
while the custom subset runs required for sign-flip validation still execute.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_supergnova_atlas as atlas  # noqa: E402


_original_run_capture = atlas.run_capture


def skip_completed_full_run(*args, **kwargs):
    # run_capture positional layout: pair, trait1, trait2, qc, out, overlap, log
    out = kwargs.get("out")
    if out is None and len(args) >= 5:
        out = args[4]
    partition = kwargs.get("partition")
    if partition is None and out is not None:
        overlap = kwargs.get("overlap")
        if overlap is None and len(args) >= 6:
            overlap = args[5]
        if Path(out).exists() and overlap is not None and Path(overlap).exists() and Path(out).stat().st_size > 1000:
            return None
    return _original_run_capture(*args, **kwargs)


atlas.run_capture = skip_completed_full_run
raise SystemExit(atlas.main())
