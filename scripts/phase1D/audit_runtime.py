#!/usr/bin/env python3
"""Audit Phase 1D pair processes and raw-result completeness."""

from __future__ import annotations

import csv
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "results" / "phase1D" / "raw_supergnova"
LOG = ROOT / "logs" / "phase1D"
OUT = ROOT / "results" / "phase1D" / "PHASE1D_RUNTIME_AUDIT.tsv"
PAIRS = [
    "AMD-CAD", "AMD-STROKE", "AMD-T2D", "AMD-CKD", "AMD-AD", "AMD-PD",
    "POAG-CAD", "POAG-STROKE", "POAG-T2D", "POAG-CKD", "POAG-AD", "POAG-PD",
    "CAT-CAD", "CAT-STROKE", "CAT-T2D", "CAT-CKD", "CAT-AD", "CAT-PD",
]
REQUIRED = {"chr", "start", "end", "rho", "corr", "h2_1", "h2_2", "var", "p", "m"}


def active_processes():
    try:
        output = subprocess.check_output(["ps", "-axo", "pid=,command="], text=True)
    except subprocess.CalledProcessError:
        return {}
    found = {pair: [] for pair in PAIRS}
    for line in output.splitlines():
        if "run_supergnova_capture.py" not in line:
            continue
        for pair in PAIRS:
            if f"/{pair}.txt" in line or f"/{pair}.overlap.tsv" in line:
                pid = line.strip().split(None, 1)[0]
                found[pair].append(pid)
    return found


def inspect_raw(path: Path):
    if not path.exists():
        return {"exists": False}
    result = {"exists": True, "size": path.stat().st_size, "modified": datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(timespec="seconds")}
    try:
        header = set(pd.read_csv(path, sep=r"\s+", nrows=0).columns)
        result["header_ok"] = REQUIRED.issubset(header)
        frame = pd.read_csv(path, sep=r"\s+")
        result["rows"] = len(frame)
        result["chromosomes"] = int(frame["chr"].nunique()) if "chr" in frame else 0
        result["finite_rho"] = int(np.isfinite(pd.to_numeric(frame.get("rho"), errors="coerce")).sum()) if "rho" in frame else 0
        result["finite_var"] = int(np.isfinite(pd.to_numeric(frame.get("var"), errors="coerce")).sum()) if "var" in frame else 0
        result["finite_P"] = int(np.isfinite(pd.to_numeric(frame.get("p"), errors="coerce")).sum()) if "p" in frame else 0
        result["finite_corr"] = int(np.isfinite(pd.to_numeric(frame.get("corr"), errors="coerce")).sum()) if "corr" in frame else 0
        result["negative_h2_trait1"] = int((pd.to_numeric(frame.get("h2_1"), errors="coerce") < 0).sum()) if "h2_1" in frame else 0
        result["negative_h2_trait2"] = int((pd.to_numeric(frame.get("h2_2"), errors="coerce") < 0).sum()) if "h2_2" in frame else 0
        numeric_ok = all(np.isfinite(pd.to_numeric(frame[col], errors="coerce")).all() for col in ["rho", "var", "p", "h2_1", "h2_2", "m"] if col in frame)
        result["numeric_ok"] = bool(numeric_ok)
    except Exception as exc:
        result["header_ok"] = False
        result["rows"] = 0
        result["chromosomes"] = 0
        result["numeric_ok"] = False
        result["error"] = str(exc)
    return result


def main():
    active = active_processes()
    rows = []
    for pair in PAIRS:
        raw = RAW / f"{pair}.txt"
        log = LOG / f"{pair}.log"
        info = inspect_raw(raw)
        pids = active.get(pair, [])
        if pids:
            status = "RUNNING"
            action = "DO_NOT_DUPLICATE_RUNNING"
            exit_code = ""
        elif not info.get("exists"):
            status = "NOT_STARTED"
            action = "CONTINUE"
            exit_code = ""
        elif info.get("header_ok") and info.get("numeric_ok") and info.get("rows", 0) > 0 and info.get("chromosomes", 0) == 22 and log.exists() and "Computed local genetic covariance for chromosome 22" in log.read_text(errors="ignore"):
            status = "COMPLETE_VALID"
            action = "KEEP_COMPLETED_VALID"
            exit_code = "0"
        else:
            status = "PARTIAL_INVALID"
            action = "INSPECT_LOG_AND_RERUN_SPECIFIC_PAIR"
            exit_code = "UNKNOWN"
        rows.append({
            "pair": pair,
            "status": status,
            "pid_if_active": ",".join(pids),
            "raw_file": str(raw.relative_to(ROOT)),
            "raw_file_size": info.get("size", ""),
            "last_modified": info.get("modified", ""),
            "log_file": str(log.relative_to(ROOT)),
            "exit_code": exit_code,
            "action": action,
            "regions_returned": info.get("rows", ""),
            "chromosomes_observed": info.get("chromosomes", ""),
            "finite_rho": info.get("finite_rho", ""),
            "finite_var": info.get("finite_var", ""),
            "finite_P": info.get("finite_P", ""),
            "finite_corr": info.get("finite_corr", ""),
            "negative_h2_trait1": info.get("negative_h2_trait1", ""),
            "negative_h2_trait2": info.get("negative_h2_trait2", ""),
        })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    columns = ["pair", "status", "pid_if_active", "raw_file", "raw_file_size", "last_modified", "log_file", "exit_code", "action", "regions_returned", "chromosomes_observed", "finite_rho", "finite_var", "finite_P", "finite_corr", "negative_h2_trait1", "negative_h2_trait2"]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
