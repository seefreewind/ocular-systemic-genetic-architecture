#!/usr/bin/env python3
"""Run and summarize the prespecified Phase 1D SUPERGNOVA atlas.

The official SUPERGNOVA source is called through the Phase 1C compatibility
wrapper.  This script keeps the source tree untouched and writes all Phase 1D
outputs only under the Phase 1D results, logs, and figure directories.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data" / "phase1C" / "SUPERGNOVA_source"
PYTHON = ROOT / "data" / "phase1C" / "supergnova_venv" / "bin" / "python"
SUMSTATS = ROOT / "data" / "phase1D" / "sumstats"
REF = ROOT / "data" / "phase1C" / "SUPERGNOVA_reference"
RAW_OUT = ROOT / "results" / "phase1D" / "raw_supergnova"
SIGN_DIR = ROOT / "results" / "phase1D" / "signflip"
LOG_DIR = ROOT / "logs" / "phase1D"
OUT_DIR = ROOT / "results" / "phase1D"
RUNNER = ROOT / "scripts" / "phase1C" / "run_supergnova_capture.py"
INPUT_QC = OUT_DIR / "SUPERGNOVA_INPUT_QC.tsv"
GLOBAL_RG = ROOT / "results" / "phase1A" / "OCULAR_SYSTEMIC_RG.tsv"
SEED = 20260918

OCULAR = ["AMD", "POAG", "CAT"]
SYSTEMIC = ["CAD", "STROKE", "T2D", "CKD", "AD", "PD"]
PAIRS = [(f"{ocular}-{systemic}", ocular, systemic) for ocular in OCULAR for systemic in SYSTEMIC]


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def ensure_dirs() -> None:
    for path in [RAW_OUT, SIGN_DIR, LOG_DIR, OUT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def read_input_qc() -> Dict[str, Dict[str, str]]:
    with INPUT_QC.open(newline="") as handle:
        return {row["trait"]: row for row in csv.DictReader(handle, delimiter="\t")}


def read_partition_count() -> int:
    total = 0
    for chrom in range(1, 23):
        path = REF / "partition" / f"eur_chr{chrom}.bed"
        frame = pd.read_csv(path, sep=r"\s+")
        total += len(frame)
    return total


def run_capture(
    pair: str,
    trait1: str,
    trait2: str,
    qc: Dict[str, Dict[str, str]],
    out: Path,
    overlap: Path,
    log: Path,
    partition: Optional[Path] = None,
    sumstats1: Optional[Path] = None,
    sumstats2: Optional[Path] = None,
) -> None:
    partition = partition or (REF / "partition" / "eur_chr@.bed")
    sumstats1 = sumstats1 or (SUMSTATS / f"{trait1}.supergnova.sumstats.gz")
    sumstats2 = sumstats2 or (SUMSTATS / f"{trait2}.supergnova.sumstats.gz")
    command = [
        str(PYTHON), str(RUNNER),
        "--source", str(SOURCE),
        "--sumstats1", str(sumstats1),
        "--sumstats2", str(sumstats2),
        "--N1", qc[trait1]["scalar_N_used"],
        "--N2", qc[trait2]["scalar_N_used"],
        "--bfile", str(REF / "bfiles" / "eur_chr@_SNPmaf5"),
        "--partition", str(partition),
        "--out", str(out),
        "--overlap-out", str(overlap),
        "--thread", "4",
    ]
    result = subprocess.run(command, cwd=str(SOURCE), text=True, capture_output=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text((result.stdout or "") + "\n" + (result.stderr or ""))
    if result.returncode != 0:
        raise RuntimeError(f"{pair} failed with return code {result.returncode}; see {log}")


def read_overlap(path: Path) -> Dict[str, float]:
    values: Dict[str, float] = {}
    if not path.exists():
        return values
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                values[row["parameter"]] = float(row["value"])
            except (KeyError, TypeError, ValueError):
                continue
    return values


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(values) & (values >= 0) & (values <= 1)
    if not bool(valid.any()):
        return result
    p = values.loc[valid].to_numpy(dtype=float)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    result.loc[valid] = output
    return result


def region_id(row: pd.Series) -> str:
    return f"{int(row['chr'])}:{int(row['start'])}-{int(row['end'])}"


def read_raw(
    pair: str,
    trait1: str,
    trait2: str,
    raw: Path,
    overlap: Path,
    attempted: int,
) -> Tuple[pd.DataFrame, Dict[str, object], Dict[str, object]]:
    frame = pd.read_csv(raw, sep=r"\s+")
    if "p" in frame.columns:
        frame.rename(columns={"p": "P_rho"}, inplace=True)
    required = ["chr", "start", "end", "rho", "corr", "h2_1", "h2_2", "var", "P_rho", "m"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{pair} output is missing columns: {missing}")
    for column in ["chr", "start", "end", "m"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    for column in ["rho", "corr", "h2_1", "h2_2", "var", "P_rho"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["SE_rho"] = np.sqrt(frame["var"].where(frame["var"] >= 0))
    frame["Z_rho"] = frame["rho"] / frame["SE_rho"]
    frame["ocular_trait"] = trait1
    frame["systemic_trait"] = trait2
    frame["pair"] = pair
    frame["region_id"] = frame.apply(region_id, axis=1)
    frame["input_status"] = "PASS"
    valid = (
        np.isfinite(frame["rho"])
        & np.isfinite(frame["var"])
        & (frame["var"] >= 0)
        & np.isfinite(frame["P_rho"])
        & frame["P_rho"].between(0, 1)
    )
    frame["numerical_status"] = np.where(valid, "PASS", "REVIEW")
    overlap_values = read_overlap(overlap)
    finite_rho = np.isfinite(frame["rho"])
    finite_corr = np.isfinite(frame["corr"])
    finite_h1 = np.isfinite(frame["h2_1"])
    finite_h2 = np.isfinite(frame["h2_2"])
    finite_var = np.isfinite(frame["var"])
    finite_p = np.isfinite(frame["P_rho"])
    qc_row: Dict[str, object] = {
        "pair": pair,
        "ocular_trait": trait1,
        "systemic_trait": trait2,
        "regions_attempted": attempted,
        "regions_returned": len(frame),
        "regions_missing_from_output": attempted - len(frame),
        "rho_finite": int(finite_rho.sum()),
        "corr_finite": int(finite_corr.sum()),
        "h2_ocular_negative": int((frame["h2_1"] < 0).sum()),
        "h2_systemic_negative": int((frame["h2_2"] < 0).sum()),
        "h2_ocular_zero": int((frame["h2_1"] == 0).sum()),
        "h2_systemic_zero": int((frame["h2_2"] == 0).sum()),
        "variance_nonfinite": int((~finite_var).sum()),
        "variance_negative": int((frame["var"] < 0).sum()),
        "P_nonfinite": int((~finite_p).sum()),
        "P_outside_0_1": int((finite_p & ((frame["P_rho"] < 0) | (frame["P_rho"] > 1))).sum()),
        "abs_rho_gt_1": int((finite_rho & (frame["rho"].abs() > 1)).sum()),
        "low_SNP_returned_m_lt_120": int((frame["m"] < 120).sum()),
        "mean_m_returned": frame["m"].mean(),
        "median_m_returned": frame["m"].median(),
        "internal_pheno_corr": overlap_values.get("pheno_corr", np.nan),
        "internal_pheno_corr_var": overlap_values.get("pheno_corr_var", np.nan),
        "N1_used": overlap_values.get("N1", np.nan),
        "N2_used": overlap_values.get("N2", np.nan),
        "analysis_snps_before_ldscore": overlap_values.get("analysis_snps_before_ldscore", np.nan),
        "status": "PASS" if bool(valid.all()) else "REVIEW",
    }
    overlap_row: Dict[str, object] = {
        "pair": pair,
        "ocular_trait": trait1,
        "systemic_trait": trait2,
        "supergnova_internal_pheno_corr": overlap_values.get("pheno_corr", np.nan),
        "supergnova_internal_pheno_corr_var": overlap_values.get("pheno_corr_var", np.nan),
        "N1_used": overlap_values.get("N1", np.nan),
        "N2_used": overlap_values.get("N2", np.nan),
        "analysis_snps_before_ldscore": overlap_values.get("analysis_snps_before_ldscore", np.nan),
    }
    return frame, qc_row, overlap_row


def global_rg_map() -> Dict[str, Dict[str, object]]:
    if not GLOBAL_RG.exists():
        return {}
    frame = pd.read_csv(GLOBAL_RG, sep="\t")
    result: Dict[str, Dict[str, object]] = {}
    for _, row in frame.iterrows():
        t1, t2 = str(row["trait1"]), str(row["trait2"])
        key = f"{t1}-{t2}"
        reverse = f"{t2}-{t1}"
        payload = row.to_dict()
        result[key] = payload
        result[reverse] = payload
    return result


def global_class(row: Optional[Dict[str, object]]) -> Tuple[str, float]:
    if row is None:
        return "GLOBAL_RESULT_MISSING", np.nan
    rg = pd.to_numeric(pd.Series([row.get("rg")]), errors="coerce").iloc[0]
    p = pd.to_numeric(pd.Series([row.get("rg_p")]), errors="coerce").iloc[0]
    sig = str(row.get("rg_significant_fdr", "NO")) == "YES"
    if sig and np.isfinite(rg) and rg > 0:
        return "GLOBAL_POSITIVE", float(rg)
    if sig and np.isfinite(rg) and rg < 0:
        return "GLOBAL_NEGATIVE", float(rg)
    if np.isfinite(rg) and np.isfinite(p):
        return "NO_DETECTABLE_GLOBAL_RG", float(rg)
    return "GLOBAL_RESULT_MISSING", float(rg) if np.isfinite(rg) else np.nan


def local_class(frame: pd.DataFrame) -> Tuple[str, int, int, int]:
    valid = frame[frame["numerical_status"] == "PASS"]
    significant = valid[valid["q_atlas"] < 0.05]
    positive = int((significant["rho"] > 0).sum())
    negative = int((significant["rho"] < 0).sum())
    if positive and negative:
        label = "MIXED_DIRECTION"
    elif positive:
        label = "POSITIVE_ONLY"
    elif negative:
        label = "NEGATIVE_ONLY"
    else:
        label = "NO_DETECTABLE_LOCAL_SIGNAL"
    return label, positive, negative, len(significant)


def write_tsv(path: Path, rows: Sequence[Dict[str, object]], columns: Optional[Sequence[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def make_signflip_sumstats(trait: str) -> Path:
    source = SUMSTATS / f"{trait}.supergnova.sumstats.gz"
    output = SIGN_DIR / f"{trait}.signflip.sumstats.gz"
    if output.exists() and output.stat().st_mtime >= source.stat().st_mtime:
        return output
    with gzip.open(source, "rt", newline="") as src, gzip.open(output, "wt", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fields = reader.fieldnames or []
        if "Z" not in fields:
            raise ValueError(f"{trait} input has no Z column")
        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            try:
                row["Z"] = str(-float(row["Z"]))
            except (TypeError, ValueError):
                row["Z"] = row["Z"]
            writer.writerow(row)
    return output


def make_subset_partition(pair: str, regions: pd.DataFrame) -> Path:
    directory = SIGN_DIR / "partitions" / pair
    directory.mkdir(parents=True, exist_ok=True)
    for chrom in range(1, 23):
        subset = regions[regions["chr"].astype(int) == chrom]
        path = directory / f"eur_chr{chrom}.bed"
        with path.open("w") as handle:
            handle.write("chr start stop\n")
            for _, row in subset.sort_values(["start", "end"]).iterrows():
                handle.write(f"{chrom} {int(row['start'])} {int(row['end'])}\n")
    return directory / "eur_chr@.bed"


def select_signflip_regions(atlas: pd.DataFrame) -> pd.DataFrame:
    stable = atlas[
        (atlas["numerical_status"] == "PASS")
        & np.isfinite(atlas["m"])
        & (atlas["m"] >= 120)
        & np.isfinite(atlas["rho"])
        & np.isfinite(atlas["var"])
        & (atlas["var"] > 0)
        & np.isfinite(atlas["P_rho"])
    ].copy()
    rng = np.random.default_rng(SEED)
    selected: List[Dict[str, object]] = []
    keys = set()
    for pair, _, _ in PAIRS:
        choices = stable[stable["pair"] == pair]
        if len(choices) < 2:
            continue
        indices = rng.choice(choices.index.to_numpy(), size=2, replace=False)
        for index in indices:
            row = stable.loc[index]
            key = (pair, int(row["chr"]), int(row["start"]), int(row["end"]))
            keys.add(key)
            selected.append({"pair": pair, "chr": int(row["chr"]), "start": int(row["start"]), "end": int(row["end"]), "selection_source": "random_per_pair"})
    high_signal = stable.sort_values(["P_rho", "pair", "chr", "start"], kind="mergesort")
    added = 0
    for _, row in high_signal.iterrows():
        key = (str(row["pair"]), int(row["chr"]), int(row["start"]), int(row["end"]))
        if key in keys:
            continue
        keys.add(key)
        selected.append({"pair": key[0], "chr": key[1], "start": key[2], "end": key[3], "selection_source": "high_signal_atlas"})
        added += 1
        if added >= 20:
            break
    return pd.DataFrame(selected, columns=["pair", "chr", "start", "end", "selection_source"])


def run_signflip(atlas: pd.DataFrame, qc: Dict[str, Dict[str, str]]) -> pd.DataFrame:
    selected = select_signflip_regions(atlas)
    if selected.empty:
        return pd.DataFrame([{"status": "FAIL", "reason": "no stable regions available"}])
    write_tsv(SIGN_DIR / "SIGN_FLIP_SELECTION.tsv", selected.to_dict("records"), ["pair", "chr", "start", "end", "selection_source"])
    flipped_paths: Dict[str, Path] = {}
    for trait in sorted(set(SYSTEMIC)):
        flipped_paths[trait] = make_signflip_sumstats(trait)
    rows: List[Dict[str, object]] = []
    for pair, trait1, trait2 in PAIRS:
        chosen = selected[selected["pair"] == pair].copy()
        if chosen.empty:
            rows.append({"pair": pair, "status": "FAIL", "reason": "no selected regions"})
            continue
        partition = make_subset_partition(pair, chosen)
        pair_dir = SIGN_DIR / pair
        pair_dir.mkdir(parents=True, exist_ok=True)
        original_raw = pair_dir / "original.txt"
        original_overlap = pair_dir / "original.overlap.tsv"
        flipped_raw = pair_dir / "flipped.txt"
        flipped_overlap = pair_dir / "flipped.overlap.tsv"
        original_log = LOG_DIR / f"{pair}.signflip_original.log"
        flipped_log = LOG_DIR / f"{pair}.signflip_flipped.log"
        try:
            run_capture(pair, trait1, trait2, qc, original_raw, original_overlap, original_log, partition=partition)
            run_capture(pair, trait1, trait2, qc, flipped_raw, flipped_overlap, flipped_log, partition=partition, sumstats2=flipped_paths[trait2])
            original = pd.read_csv(original_raw, sep=r"\s+")
            flipped = pd.read_csv(flipped_raw, sep=r"\s+")
            original.rename(columns={"p": "P_rho"}, inplace=True)
            flipped.rename(columns={"p": "P_rho"}, inplace=True)
            merged = chosen.merge(original, on=["chr", "start", "end"], how="left", suffixes=("", "_orig"))
            merged = merged.merge(flipped, on=["chr", "start", "end"], how="left", suffixes=("_orig", "_flip"))
            for _, row in merged.iterrows():
                h1_delta = abs(float(row["h2_1_orig"]) - float(row["h2_1_flip"]))
                h2_delta = abs(float(row["h2_2_orig"]) - float(row["h2_2_flip"]))
                rho_orig = float(row["rho_orig"])
                rho_flip = float(row["rho_flip"])
                p_orig = float(row["P_rho_orig"])
                p_flip = float(row["P_rho_flip"])
                rho_sign_ok = abs(rho_orig + rho_flip) <= max(1e-8, abs(rho_orig) * 1e-6)
                abs_rho_delta = abs(abs(rho_orig) - abs(rho_flip))
                p_delta = abs(p_orig - p_flip)
                passed = (
                    h1_delta <= 1e-8
                    and h2_delta <= 1e-8
                    and rho_sign_ok
                    and abs_rho_delta <= max(1e-8, abs(rho_orig) * 1e-6)
                    and p_delta <= 1e-8
                )
                rows.append({
                    "pair": pair,
                    "chr": int(row["chr"]),
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "selection_source": row["selection_source"],
                    "h2_ocular_abs_diff": h1_delta,
                    "h2_systemic_abs_diff": h2_delta,
                    "rho_original": rho_orig,
                    "rho_flipped": rho_flip,
                    "rho_sign_reversal": "PASS" if rho_sign_ok else "FAIL",
                    "abs_rho_abs_diff": abs_rho_delta,
                    "P_original": p_orig,
                    "P_flipped": p_flip,
                    "P_abs_diff": p_delta,
                    "status": "PASS" if passed else "FAIL",
                })
            if len(merged) != len(chosen):
                rows.append({"pair": pair, "status": "FAIL", "reason": f"matched {len(merged)} of {len(chosen)} selected regions"})
        except Exception as exc:
            rows.append({"pair": pair, "status": "FAIL", "reason": str(exc)})
    result = pd.DataFrame(rows)
    columns = [
        "pair", "chr", "start", "end", "selection_source", "h2_ocular_abs_diff", "h2_systemic_abs_diff",
        "rho_original", "rho_flipped", "rho_sign_reversal", "abs_rho_abs_diff", "P_original", "P_flipped",
        "P_abs_diff", "status", "reason",
    ]
    write_tsv(OUT_DIR / "SIGN_FLIP_ATLAS_QC.tsv", result.to_dict("records"), columns)
    return result


def classify_pairs(atlas: pd.DataFrame, qc_rows: List[Dict[str, object]]) -> pd.DataFrame:
    rg_map = global_rg_map()
    output: List[Dict[str, object]] = []
    for pair, trait1, trait2 in PAIRS:
        frame = atlas[atlas["pair"] == pair]
        lclass, pos, neg, sig = local_class(frame)
        gclass, rg = global_class(rg_map.get(pair))
        if gclass == "NO_DETECTABLE_GLOBAL_RG" and lclass == "MIXED_DIRECTION":
            hidden = "HIDDEN_MIXED_LOCAL_SHARING"
        else:
            hidden = "NO"
        if lclass == "MIXED_DIRECTION":
            balance = "MIXED_LOCAL_DIRECTIONS"
        elif lclass == "POSITIVE_ONLY":
            balance = "POSITIVE_LOCAL_ONLY"
        elif lclass == "NEGATIVE_ONLY":
            balance = "NEGATIVE_LOCAL_ONLY"
        else:
            balance = "NO_DETECTABLE_LOCAL_SIGNAL"
        if gclass == "GLOBAL_POSITIVE" and lclass == "POSITIVE_ONLY":
            concordance = "GLOBAL_LOCAL_CONCORDANT_POSITIVE"
        elif gclass == "GLOBAL_NEGATIVE" and lclass == "NEGATIVE_ONLY":
            concordance = "GLOBAL_LOCAL_CONCORDANT_NEGATIVE"
        elif gclass in ["GLOBAL_POSITIVE", "GLOBAL_NEGATIVE"] and sig:
            concordance = "GLOBAL_LOCAL_DIRECTIONALLY_MIXED_OR_DISCORDANT"
        elif gclass == "NO_DETECTABLE_GLOBAL_RG" and sig:
            concordance = "LOCAL_SIGNAL_WITHOUT_DETECTABLE_GLOBAL_RG"
        else:
            concordance = "NO_LOCAL_SIGNAL"
        output.append({
            "pair": pair,
            "ocular_trait": trait1,
            "systemic_trait": trait2,
            "global_rg": rg,
            "global_class": gclass,
            "local_class": lclass,
            "local_positive_atlas_fdr_count": pos,
            "local_negative_atlas_fdr_count": neg,
            "local_total_atlas_fdr_count": sig,
            "hidden_mixed_local_sharing": hidden,
            "global_local_classification": concordance,
            "local_directional_summary": balance,
        })
    return pd.DataFrame(output)


def burden_table(atlas: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for pair, trait1, trait2 in PAIRS:
        frame = atlas[(atlas["pair"] == pair) & (atlas["numerical_status"] == "PASS")].copy()
        sig = frame[frame["q_atlas"] < 0.05]
        def metrics(data: pd.DataFrame, suffix: str) -> Dict[str, object]:
            rho = data["rho"].to_numpy(dtype=float)
            return {
                f"positive_burden{suffix}": float(np.nansum(np.maximum(rho, 0))),
                f"negative_burden{suffix}": float(np.nansum(np.abs(np.minimum(rho, 0)))),
                f"signed_net{suffix}": float(np.nansum(rho)),
                f"absolute_total{suffix}": float(np.nansum(np.abs(rho))),
                f"regions{suffix}": len(data),
            }
        row: Dict[str, object] = {"pair": pair, "ocular_trait": trait1, "systemic_trait": trait2}
        row.update(metrics(frame, ""))
        row.update(metrics(sig, "_atlas_sig"))
        absolute = float(row["absolute_total"])
        row["balance_index"] = float(row["signed_net"]) / absolute if absolute > 0 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, seed: int = SEED, n_boot: int = 10000) -> Tuple[float, float, float, float]:
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if len(x) < 3:
        return np.nan, np.nan, np.nan, np.nan
    statistic = spearmanr(x, y, nan_policy="omit")
    rho = float(statistic.statistic)
    p = float(statistic.pvalue)
    rng = np.random.default_rng(seed)
    values = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        indices = rng.integers(0, len(x), size=len(x))
        if len(np.unique(x[indices])) < 2 or len(np.unique(y[indices])) < 2:
            values[i] = np.nan
        else:
            values[i] = float(spearmanr(x[indices], y[indices]).statistic)
    finite_boot = values[np.isfinite(values)]
    low, high = (np.nan, np.nan) if len(finite_boot) == 0 else np.percentile(finite_boot, [2.5, 97.5])
    return rho, p, float(low), float(high)


def summary_correlations(classification: pd.DataFrame, burden: pd.DataFrame) -> pd.DataFrame:
    merged = classification.merge(burden, on=["pair", "ocular_trait", "systemic_trait"], how="inner")
    rg = merged["global_rg"].to_numpy(dtype=float)
    metrics = ["signed_net", "positive_burden", "negative_burden", "absolute_total", "balance_index"]
    rows: List[Dict[str, object]] = []
    for metric in metrics:
        y = merged[metric].to_numpy(dtype=float)
        rho, p, low, high = bootstrap_spearman(rg, y, seed=SEED + metrics.index(metric))
        rows.append({
            "global_measure": "global_rg",
            "local_measure": metric,
            "n_pairs": int(np.isfinite(rg).astype(int).dot(np.isfinite(y).astype(int))),
            "spearman_rho": rho,
            "spearman_p": p,
            "bootstrap_seed": SEED + metrics.index(metric),
            "bootstrap_replicates": 10000,
            "bootstrap_95ci_low": low,
            "bootstrap_95ci_high": high,
        })
    return pd.DataFrame(rows)


def atlas_decision(classification: pd.DataFrame, signflip: pd.DataFrame, atlas: pd.DataFrame) -> Tuple[str, str]:
    qc_ok = bool((signflip["status"] == "PASS").all()) if not signflip.empty and "status" in signflip.columns else False
    if not qc_ok:
        return "METHOD_UNSTABLE", "TECHNICAL_HOLD"
    valid = atlas[atlas["numerical_status"] == "PASS"]
    sig = valid[valid["q_atlas"] < 0.05]
    hidden = int((classification["hidden_mixed_local_sharing"] == "HIDDEN_MIXED_LOCAL_SHARING").sum())
    mixed = int((classification["local_class"] == "MIXED_DIRECTION").sum())
    if hidden >= 2 or mixed >= 3:
        return "STRONG_GLOBAL_LOCAL_DISCORDANCE", "WAIT_FOR_LAVA_REPLICATION"
    if hidden == 1 or mixed == 1:
        return "LIMITED_GLOBAL_LOCAL_DISCORDANCE", "WAIT_FOR_LAVA_REPLICATION"
    if len(sig) > 0:
        return "LOCAL_SHARING_WITHOUT_MIXED_DIRECTION", "WAIT_FOR_LAVA_REPLICATION"
    return "NO_ROBUST_LOCAL_SHARING", "WAIT_FOR_LAVA_REPLICATION"


def main() -> int:
    ensure_dirs()
    qc = read_input_qc()
    attempted = read_partition_count()
    raw_frames: List[pd.DataFrame] = []
    pair_qc: List[Dict[str, object]] = []
    overlap_rows: List[Dict[str, object]] = []
    failures: List[Dict[str, object]] = []
    for pair, trait1, trait2 in PAIRS:
        raw = RAW_OUT / f"{pair}.txt"
        overlap = RAW_OUT / f"{pair}.overlap.tsv"
        log = LOG_DIR / f"{pair}.log"
        print(f"[{now()}] START {pair}", flush=True)
        try:
            run_capture(pair, trait1, trait2, qc, raw, overlap, log)
            frame, qc_row, overlap_row = read_raw(pair, trait1, trait2, raw, overlap, attempted)
            raw_frames.append(frame)
            pair_qc.append(qc_row)
            overlap_rows.append(overlap_row)
            print(f"[{now()}] DONE {pair} regions={len(frame)}", flush=True)
        except Exception as exc:
            failures.append({"pair": pair, "ocular_trait": trait1, "systemic_trait": trait2, "stage": "atlas_run", "error": str(exc)})
            pair_qc.append({"pair": pair, "ocular_trait": trait1, "systemic_trait": trait2, "regions_attempted": attempted, "status": "FAIL", "error": str(exc)})
            print(f"[{now()}] FAIL {pair}: {exc}", file=sys.stderr, flush=True)
    if not raw_frames:
        write_tsv(OUT_DIR / "PHASE1D_FAILURES.tsv", failures, ["pair", "ocular_trait", "systemic_trait", "stage", "error"])
        return 2
    atlas = pd.concat(raw_frames, ignore_index=True)
    atlas["q_atlas"] = bh(atlas["P_rho"])
    atlas["q_pair"] = np.nan
    for pair, group in atlas.groupby("pair", sort=False):
        atlas.loc[group.index, "q_pair"] = bh(group["P_rho"])
    n_tests = int(atlas["q_atlas"].notna().sum())
    atlas["bonferroni_p_atlas"] = np.minimum(atlas["P_rho"] * n_tests, 1.0)
    atlas["direction"] = np.where(atlas["rho"] > 0, "POSITIVE", np.where(atlas["rho"] < 0, "NEGATIVE", "ZERO"))
    atlas_columns = [
        "ocular_trait", "systemic_trait", "pair", "region_id", "chr", "start", "end", "rho", "var", "SE_rho", "Z_rho", "P_rho",
        "q_atlas", "q_pair", "bonferroni_p_atlas", "h2_1", "h2_2", "corr", "m", "direction", "input_status", "numerical_status",
    ]
    atlas = atlas.rename(columns={"h2_1": "h2_ocular", "h2_2": "h2_systemic", "corr": "local_corr"})
    atlas_columns = ["h2_ocular" if c == "h2_1" else "h2_systemic" if c == "h2_2" else "local_corr" if c == "corr" else c for c in atlas_columns]
    write_tsv(OUT_DIR / "SUPERGNOVA_OCULAR_SYSTEMIC_ATLAS.tsv", atlas[atlas_columns].to_dict("records"), atlas_columns)
    write_tsv(OUT_DIR / "SUPERGNOVA_PAIR_QC.tsv", pair_qc)
    write_tsv(OUT_DIR / "SUPERGNOVA_OVERLAP_QC.tsv", overlap_rows)
    write_tsv(OUT_DIR / "PHASE1D_FAILURES.tsv", failures, ["pair", "ocular_trait", "systemic_trait", "stage", "error"])

    classification = classify_pairs(atlas, pair_qc)
    burden = burden_table(atlas)
    correlations = summary_correlations(classification, burden)
    classification.to_csv(OUT_DIR / "GLOBAL_LOCAL_ARCHITECTURE_CLASSIFICATION.tsv", sep="\t", index=False)
    burden.to_csv(OUT_DIR / "LOCAL_COVARIANCE_BURDEN.tsv", sep="\t", index=False)
    correlations.to_csv(OUT_DIR / "GLOBAL_LOCAL_SUMMARY_CORRELATIONS.tsv", sep="\t", index=False)

    direction_rows: List[Dict[str, object]] = []
    for pair, trait1, trait2 in PAIRS:
        frame = atlas[(atlas["pair"] == pair) & (atlas["numerical_status"] == "PASS")]
        atlas_sig = frame[frame["q_atlas"] < 0.05]
        pair_sig = frame[frame["q_pair"] < 0.05]
        direction_rows.append({
            "pair": pair,
            "ocular_trait": trait1,
            "systemic_trait": trait2,
            "regions_tested": len(frame),
            "atlas_fdr_significant_total": len(atlas_sig),
            "atlas_fdr_positive": int((atlas_sig["rho"] > 0).sum()),
            "atlas_fdr_negative": int((atlas_sig["rho"] < 0).sum()),
            "pair_fdr_significant_total": len(pair_sig),
            "pair_fdr_positive": int((pair_sig["rho"] > 0).sum()),
            "pair_fdr_negative": int((pair_sig["rho"] < 0).sum()),
            "atlas_bonferroni_significant_total": int((frame["bonferroni_p_atlas"] < 0.05).sum()),
            "minimum_P_rho": frame["P_rho"].min(),
            "minimum_q_atlas": frame["q_atlas"].min(),
        })
    write_tsv(OUT_DIR / "DIRECTIONAL_LOCAL_SIGNAL_SUMMARY.tsv", direction_rows)

    top_rows: List[Dict[str, object]] = []
    for pair, _, _ in PAIRS:
        frame = atlas[atlas["pair"] == pair].sort_values(["q_atlas", "P_rho"], kind="mergesort").head(5)
        for rank, (_, row) in enumerate(frame.iterrows(), start=1):
            top_rows.append({
                "pair": pair, "rank_within_pair": rank, "chr": int(row["chr"]), "start": int(row["start"]), "end": int(row["end"]),
                "region_id": row["region_id"], "rho": row["rho"], "SE_rho": row["SE_rho"], "P_rho": row["P_rho"], "q_atlas": row["q_atlas"],
                "direction": row["direction"], "atlas_fdr_significant": "YES" if row["q_atlas"] < 0.05 else "NO",
            })
    write_tsv(OUT_DIR / "TOP_LOCAL_REGIONS.tsv", top_rows)

    priority_rows: List[Dict[str, object]] = []
    class_map = classification.set_index("pair").to_dict("index")
    significant = atlas[(atlas["numerical_status"] == "PASS") & (atlas["q_atlas"] < 0.05)].sort_values(["q_atlas", "P_rho"], kind="mergesort")
    for rank, (_, row) in enumerate(significant.iterrows(), start=1):
        cls = class_map.get(row["pair"], {})
        if cls.get("hidden_mixed_local_sharing") == "HIDDEN_MIXED_LOCAL_SHARING":
            tier = "TIER_1_HIDDEN_MIXED_GLOBAL_NULL"
        elif cls.get("global_class") in ["GLOBAL_POSITIVE", "GLOBAL_NEGATIVE"] and cls.get("local_class") in ["POSITIVE_ONLY", "NEGATIVE_ONLY"]:
            tier = "TIER_2_GLOBAL_LOCAL_CONCORDANT"
        else:
            tier = "TIER_3_OTHER_ATLAS_SIGNAL"
        priority_rows.append({
            "priority_rank": rank, "pair": row["pair"], "ocular_trait": row["ocular_trait"], "systemic_trait": row["systemic_trait"],
            "chr": int(row["chr"]), "start": int(row["start"]), "end": int(row["end"]), "region_id": row["region_id"],
            "rho": row["rho"], "SE_rho": row["SE_rho"], "P_rho": row["P_rho"], "q_atlas": row["q_atlas"], "direction": row["direction"], "tier": tier,
        })
    write_tsv(OUT_DIR / "PHASE2_LOCUS_PRIORITY.tsv", priority_rows)

    print(f"[{now()}] START SIGN_FLIP", flush=True)
    signflip = run_signflip(atlas, qc)
    verdict, next_step = atlas_decision(classification, signflip, atlas)
    decision = [
        "PHASE1D_VERDICT\t" + verdict,
        "SIGN_FLIP_QC\t" + ("PASS" if bool((signflip.get("status") == "PASS").all()) else "FAIL"),
        f"ATLAS_VALID_TESTS\t{n_tests}",
        f"ATLAS_FDR_SIGNIFICANT_REGIONS\t{int((atlas['q_atlas'] < 0.05).sum())}",
        "NEXT_STEP\t" + next_step,
        "GENERATED_AT\t" + now(),
    ]
    (OUT_DIR / "PHASE1D_DECISION.txt").write_text("\n".join(decision) + "\n")
    print(f"[{now()}] COMPLETE verdict={verdict} next={next_step}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
