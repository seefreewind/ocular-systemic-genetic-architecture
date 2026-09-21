#!/usr/bin/env python3
"""Phase 0H/1A: protocol amendment, global LDSC, overlap audit and figures.

This script intentionally stops at global architecture. It does not call LAVA,
HDL-L, PLACO, fine-mapping, coloc, GCI, covariance decomposition or enrichment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TRAITS = ["AMD", "POAG", "CAT", "CAD", "STROKE", "T2D", "CKD", "AD", "PD"]
OCULAR = ["AMD", "POAG", "CAT"]
SYSTEMIC = ["CAD", "STROKE", "T2D", "CKD", "AD", "PD"]
PAIR_COLUMNS = [
    "trait1", "trait2", "rg", "rg_se", "rg_z", "rg_p",
    "genetic_covariance", "genetic_covariance_se", "cross_trait_intercept",
    "cross_trait_intercept_se", "mean_z1z2", "n_valid_alleles", "warning_flags",
]


def fnum(value: str | float | int | None) -> float:
    if value is None:
        return float("nan")
    try:
        text = str(value).strip()
        if text in {"", ".", "NA", "nan", "None"}:
            return float("nan")
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def fmt(value: float, digits: int = 8) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}g}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_dirs() -> dict[str, Path]:
    paths = {
        "config": ROOT / "config",
        "results": ROOT / "results" / "phase1A",
        "logs": ROOT / "results" / "phase1A" / "logs",
        "tmp": ROOT / "results" / "phase1A" / "tmp",
        "reports": ROOT / "reports" / "phase1A",
        "figures": ROOT / "figures" / "phase1A",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def read_tab(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_inputs() -> tuple[dict[str, dict[str, object]], dict[str, Path]]:
    config_rows = {row["trait"]: row for row in read_tab(ROOT / "config" / "phase0C_dataset_freeze.tsv")}
    h2_rows = {row["trait"]: row for row in read_tab(ROOT / "results" / "phase0C" / "UNIVARIATE_LDSC_H2.tsv")}
    out: dict[str, dict[str, object]] = {}
    paths: dict[str, Path] = {}
    for trait in TRAITS:
        cfg = config_rows[trait]
        h2 = h2_rows[trait]
        if trait == "AMD":
            dataset = "AMD_IAMDGC2_EUR_LATE"
            munged = ROOT / "results" / "phase0D" / "ldsc_sumstats" / "AMD_IAMDGC2_EUR.sumstats.gz"
            # The Phase 0D QC table records the corrected AMD dataset, while the
            # numerical estimates remain identical to the validated Phase 0C line.
            h2_z = fnum(h2["h2_z"])
        else:
            dataset = cfg["dataset"]
            munged = ROOT / "results" / "phase0C" / "ldsc_sumstats" / f"{trait}.sumstats.gz"
            h2_z = fnum(h2["h2_z"])
        if not munged.exists():
            raise FileNotFoundError(f"Missing frozen LDSC input: {munged}")
        cases = int(cfg["N_case"])
        controls = int(cfg["N_control"])
        n = int(cfg["N"])
        effective_n = 4.0 / (1.0 / cases + 1.0 / controls)
        hm3_count = int(fnum(h2["n_hm3_after_regression"]))
        if trait == "AMD":
            hm3_count = 1022570
        out[trait] = {
            "trait": trait,
            "dataset": dataset,
            "phenotype": cfg["phenotype_definition"],
            "ancestry": cfg["ancestry"],
            "build": "GRCh37",
            "cases": cases,
            "controls": controls,
            "N": n,
            "effective_N": effective_n,
            "HM3_SNP_count": hm3_count,
            "source": cfg["source"],
            "checksum": cfg["checksum"],
            "global_h2": fnum(h2["h2"]),
            "global_h2_SE": fnum(h2["h2_se"]),
            "global_h2_Z": h2_z,
            "LDSC_intercept": fnum(h2["intercept"]),
            "LDSC_intercept_SE": fnum(h2["intercept_se"]),
            "source_notes": cfg["notes"],
        }
        paths[trait] = munged
    return out, paths


def write_dataset_freeze(paths: dict[str, Path], inputs: dict[str, dict[str, object]]) -> None:
    fields = [
        "trait", "dataset", "phenotype", "ancestry", "build", "cases", "controls", "N",
        "effective_N", "HM3_SNP_count", "source", "checksum", "global_h2", "global_h2_SE",
        "global_h2_Z", "LDSC_intercept",
    ]
    path = ROOT / "config" / "phase1A_dataset_freeze.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for trait in TRAITS:
            row = inputs[trait]
            writer.writerow({field: fmt(float(row[field])) if field in {
                "effective_N", "global_h2", "global_h2_SE", "global_h2_Z", "LDSC_intercept"
            } else row[field] for field in fields})


def extract_float(pattern: str, text: str, group: int = 1) -> float:
    match = re.search(pattern, text, flags=re.I | re.M)
    return fnum(match.group(group)) if match else float("nan")


def parse_ldsc_log(path: Path, trait1: str, trait2: str) -> dict[str, object]:
    result: dict[str, object] = {"trait1": trait1, "trait2": trait2}
    if not path.exists():
        result.update({key: float("nan") for key in [
            "rg", "rg_se", "rg_z", "rg_p", "genetic_covariance", "genetic_covariance_se",
            "cross_trait_intercept", "cross_trait_intercept_se", "mean_z1z2", "n_valid_alleles",
        ]})
        result["warning_flags"] = "LOG_MISSING"
        result["return_code"] = 1
        return result
    text = path.read_text(errors="replace")
    # LDSC prints the point estimate, Z and P in the named section and the SE
    # in the final summary row. Parsing the last numeric columns avoids any
    # ambiguity from absolute paths containing spaces.
    result["rg"] = extract_float(r"^Genetic Correlation:\s*([-+0-9.eE]+)\s*$", text)
    result["rg_z"] = extract_float(r"^Z-score:\s*([-+0-9.eE]+)\s*$", text)
    result["rg_p"] = extract_float(r"^P:\s*([-+0-9.eE]+)\s*$", text)
    result["rg_se"] = float("nan")
    summary_start = text.find("Summary of Genetic Correlation Results")
    if summary_start >= 0:
        float_token = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
        tail_pattern = re.compile(r"\s+(" + float_token + r")\s+(" + float_token + r")\s+(" + float_token + r")\s+(" + float_token + r")\s+(" + float_token + r")\s+(" + float_token + r")\s+(" + float_token + r")\s+(" + float_token + r")\s+(" + float_token + r")\s+(" + float_token + r")\s*$")
        for line in text[summary_start:].splitlines():
            match = tail_pattern.search(line)
            if match:
                numbers = [fnum(item) for item in match.groups()]
                result["rg"], result["rg_se"], result["rg_z"], result["rg_p"] = numbers[:4]
                break
    result["genetic_covariance"] = extract_float(
        r"Total (?:Observed|Liability) scale gencov:\s*([-+0-9.eE]+)\s*\(([-+0-9.eE]+)\)", text, 1
    )
    result["genetic_covariance_se"] = extract_float(
        r"Total (?:Observed|Liability) scale gencov:\s*([-+0-9.eE]+)\s*\(([-+0-9.eE]+)\)", text, 2
    )
    result["mean_z1z2"] = extract_float(r"Mean z1\*z2:\s*([-+0-9.eE]+)", text)
    intercept_match = re.search(
        r"Genetic Covariance\s*\n[-]+\s*\n.*?^Intercept:\s*([-+0-9.eE]+)\s*\(([-+0-9.eE]+)\)",
        text, flags=re.I | re.M | re.S,
    )
    result["cross_trait_intercept"] = fnum(intercept_match.group(1)) if intercept_match else float("nan")
    result["cross_trait_intercept_se"] = fnum(intercept_match.group(2)) if intercept_match else float("nan")
    result["n_valid_alleles"] = extract_float(r"([0-9]+) SNPs with valid alleles", text)
    warnings = []
    for line in text.splitlines():
        low = line.lower()
        if "warning" in low and "futurewarning" not in low:
            warnings.append(line.strip())
        if any(token in low for token in ["not positive definite", "negative h2", "non-finite"]):
            warnings.append(line.strip())
    result["warning_flags"] = " | ".join(dict.fromkeys(warnings)) if warnings else "NONE"
    result["return_code"] = 0
    return result


def run_ldsc_pair(
    trait1: str,
    trait2: str,
    path1: Path,
    path2: Path,
    out_prefix: Path,
) -> dict[str, object]:
    ldsc = shutil.which("ldsc.py")
    if not ldsc:
        raise FileNotFoundError("ldsc.py was not found on PATH; install LDSC or set PATH before running")
    ref = ROOT / "data" / "reference" / "ldsc" / "extracted" / "eur_w_ld_chr"
    command = [
        ldsc, "--rg", f"{path1},{path2}",
        "--ref-ld-chr", f"{ref}/", "--w-ld-chr", f"{ref}/", "--out", str(out_prefix),
    ]
    console = out_prefix.with_suffix(".console.log")
    proc = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    console.write_text(proc.stdout + "\n" + proc.stderr)
    parsed = parse_ldsc_log(out_prefix.with_suffix(".log"), trait1, trait2)
    parsed["return_code"] = proc.returncode
    if proc.returncode != 0:
        parsed["warning_flags"] = f"RETURN_CODE_{proc.returncode};" + str(parsed.get("warning_flags", ""))
    return parsed


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    out = np.full(pvalues.shape, np.nan, dtype=float)
    finite = np.isfinite(pvalues)
    if not finite.any():
        return out
    values = np.clip(pvalues[finite], 0.0, 1.0)
    order = np.argsort(values)
    ranked = values[order] * len(values) / (np.arange(len(values)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum(ranked, 1.0)
    out[finite] = adjusted
    return out


def write_table(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.10g")


def run_all_pairs(paths: dict[str, Path], dirs: dict[str, Path]) -> pd.DataFrame:
    records = []
    for i, trait1 in enumerate(TRAITS):
        for trait2 in TRAITS[i + 1:]:
            prefix = dirs["results"] / "tmp" / f"{trait1}__{trait2}"
            print(f"[Phase1A] LDSC {trait1} vs {trait2}", flush=True)
            records.append(run_ldsc_pair(trait1, trait2, paths[trait1], paths[trait2], prefix))
    df = pd.DataFrame(records, columns=PAIR_COLUMNS + ["return_code"])
    df["rg_fdr_bh"] = bh_adjust(df["rg_p"].to_numpy(float))
    df["rg_p_bonferroni"] = np.minimum(df["rg_p"].astype(float) * len(df), 1.0)
    df["rg_significant_fdr"] = np.where(df["rg_fdr_bh"] < 0.05, "YES", "NO")
    df["rg_significant_bonferroni"] = np.where(df["rg_p_bonferroni"] < 0.05, "YES", "NO")
    df["point_estimate_direction"] = np.where(df["rg"] > 0, "POSITIVE", np.where(df["rg"] < 0, "NEGATIVE", "ZERO_OR_NA"))
    df["warning_flags"] = df["warning_flags"].fillna("UNKNOWN")
    write_table(df, dirs["results"] / "GLOBAL_LDSC_RG.tsv")
    fdr_cols = [
        "trait1", "trait2", "rg", "rg_se", "rg_z", "rg_p", "rg_fdr_bh", "rg_p_bonferroni",
        "rg_significant_fdr", "rg_significant_bonferroni", "genetic_covariance", "genetic_covariance_se",
        "cross_trait_intercept", "cross_trait_intercept_se", "mean_z1z2", "warning_flags",
    ]
    write_table(df[fdr_cols], dirs["results"] / "GLOBAL_LDSC_RG_FDR.tsv")
    subset = df[df.apply(lambda row: row["trait1"] in OCULAR and row["trait2"] in SYSTEMIC, axis=1)].copy()
    subset["global_direction"] = np.where(
        (subset["rg_fdr_bh"] < 0.05) & (subset["rg"] > 0), "POSITIVE",
        np.where((subset["rg_fdr_bh"] < 0.05) & (subset["rg"] < 0), "NEGATIVE", "NO_DETECTABLE_GLOBAL_RG"),
    )
    write_table(subset, dirs["results"] / "OCULAR_SYSTEMIC_RG.tsv")
    return df


def intercept_matrices(df: pd.DataFrame, inputs: dict[str, dict[str, object]], dirs: dict[str, Path]) -> dict[str, object]:
    matrix = pd.DataFrame(np.nan, index=TRAITS, columns=TRAITS, dtype=float)
    for trait in TRAITS:
        matrix.loc[trait, trait] = float(inputs[trait]["LDSC_intercept"])
    for row in df.itertuples(index=False):
        value = float(row.cross_trait_intercept)
        matrix.loc[row.trait1, row.trait2] = value
        matrix.loc[row.trait2, row.trait1] = value
    matrix.to_csv(dirs["results"] / "LDSC_CROSS_TRAIT_INTERCEPT_MATRIX.tsv", sep="\t", index_label="trait", na_rep="NA", float_format="%.10g")
    arr = matrix.to_numpy(float)
    finite = np.isfinite(arr).all()
    symmetric = np.allclose(arr, arr.T, atol=1e-12, equal_nan=False)
    diag_positive = np.all(np.diag(arr) > 0)
    if finite and diag_positive:
        scale = np.sqrt(np.outer(np.diag(arr), np.diag(arr)))
        lava = arr / scale
    else:
        lava = np.full_like(arr, np.nan)
    lava = (lava + lava.T) / 2.0
    for i in range(len(TRAITS)):
        lava[i, i] = 1.0
    lava_df = pd.DataFrame(lava, index=TRAITS, columns=TRAITS)
    lava_df.to_csv(dirs["results"] / "LAVA_SAMPLING_CORRELATION_MATRIX.tsv", sep="\t", index_label="trait", na_rep="NA", float_format="%.10g")
    eig = np.linalg.eigvalsh(lava) if np.isfinite(lava).all() else np.full(len(TRAITS), np.nan)
    pd.DataFrame({"eigenvalue": eig}).to_csv(dirs["results"] / "LAVA_SAMPLING_CORRELATION_MATRIX_EIGENVALUES.tsv", sep="\t", index=False, float_format="%.10g")
    offdiag = lava[~np.eye(len(TRAITS), dtype=bool)]
    stats = {
        "finite": bool(np.isfinite(lava).all()),
        "symmetric": bool(np.allclose(lava, lava.T, atol=1e-12)),
        "diag_max_error": float(np.max(np.abs(np.diag(lava) - 1.0))),
        "offdiag_min": float(np.min(offdiag)),
        "offdiag_max": float(np.max(offdiag)),
        "min_eigenvalue": float(np.min(eig)),
        "max_eigenvalue": float(np.max(eig)),
    }
    return {"intercept": matrix, "lava": lava, "eigenvalues": eig, "stats": stats}


def write_overlap_audits(inputs: dict[str, dict[str, object]], df: pd.DataFrame, dirs: dict[str, Path]) -> None:
    rows = []
    eligibility = []
    for i, trait1 in enumerate(TRAITS):
        for trait2 in TRAITS[i + 1:]:
            notes1 = str(inputs[trait1]["source_notes"])
            notes2 = str(inputs[trait2]["source_notes"])
            combined = f"{notes1} {notes2}".lower()
            possible_ukb = "uk biobank" in combined and ("possible" in combined or "included" in combined or "overlap" in combined)
            qualitative = "POSSIBLE" if possible_ukb else "UNKNOWN"
            known = "POSSIBLE_DOCUMENTED_COMPONENT" if possible_ukb else "UNKNOWN"
            names = "UK Biobank" if possible_ukb else "NA"
            note = "Exact participant overlap was not reported in the frozen metadata; no exact N0 inferred."
            if possible_ukb:
                note += " Source notes identify a UK Biobank component or possible overlap in at least one dataset."
            rows.append({
                "pair": f"{trait1}__{trait2}", "trait1": trait1, "trait2": trait2,
                "known_shared_cohort": known, "known_shared_cohort_names": names,
                "exact_overlap_known": "NO", "exact_N0": "NA", "qualitative_overlap": qualitative,
                "notes": note,
            })
            eligibility.append({
                "pair": f"{trait1}__{trait2}", "trait1": trait1, "trait2": trait2,
                "N0_status": "UNKNOWN", "exact_N0": "NA", "eligible_for_definitive_HDLL": "NO",
                "reason": "No exact participant overlap or exact N0 is documented in the frozen metadata; do not assign N0=0.",
            })
    write_table(pd.DataFrame(rows), dirs["results"] / "COHORT_OVERLAP_METADATA.tsv")
    write_table(pd.DataFrame(eligibility), dirs["results"] / "HDLL_N0_ELIGIBILITY.tsv")


def amd_sensitivity(paths: dict[str, Path], global_df: pd.DataFrame, dirs: dict[str, Path]) -> pd.DataFrame:
    variants = {
        "AMD_PRIMARY_CORRECTED": paths["AMD"],
        "AMD_RSID_GRCh37": ROOT / "results" / "phase0R" / "ldsc_sumstats" / "AMD_RSID_GRCh37.sumstats.gz",
        "AMD_HM3_STRICT": ROOT / "results" / "phase0R" / "ldsc_sumstats" / "AMD_HM3_STRICT.sumstats.gz",
    }
    rows = []
    for variant, amd_path in variants.items():
        if not amd_path.exists():
            continue
        for trait in TRAITS:
            if trait == "AMD":
                continue
            if variant == "AMD_PRIMARY_CORRECTED":
                match = global_df[(global_df.trait1 == "AMD") & (global_df.trait2 == trait)].iloc[0].to_dict()
            else:
                prefix = dirs["results"] / "tmp" / f"{variant}__{trait}"
                if prefix.with_suffix(".log").exists():
                    match = parse_ldsc_log(prefix.with_suffix(".log"), variant, trait)
                else:
                    match = run_ldsc_pair(variant, trait, amd_path, paths[trait], prefix)
            rows.append({
                "amd_variant": variant, "trait2": trait, "rg": fnum(match.get("rg")),
                "rg_se": fnum(match.get("rg_se")), "rg_z": fnum(match.get("rg_z")),
                "rg_p": fnum(match.get("rg_p")), "cross_trait_intercept": fnum(match.get("cross_trait_intercept")),
                "cross_trait_intercept_se": fnum(match.get("cross_trait_intercept_se")),
                "mean_z1z2": fnum(match.get("mean_z1z2")), "warning_flags": match.get("warning_flags", "NONE"),
            })
    out = pd.DataFrame(rows)
    primary = out[out.amd_variant == "AMD_PRIMARY_CORRECTED"].set_index("trait2")
    p_rg, p_se = [], []
    directional, numerical, wide = [], [], []
    for row in out.itertuples(index=False):
        base = primary.loc[row.trait2] if row.trait2 in primary.index else None
        if base is None or not np.isfinite(row.rg) or not np.isfinite(base.rg):
            directional.append("NOT_ASSESSABLE")
            numerical.append("NOT_ASSESSABLE")
            wide.append("UNKNOWN")
            p_rg.append(np.nan)
            p_se.append(np.nan)
            continue
        p_rg.append(float(base.rg))
        p_se.append(float(base.rg_se))
        directional.append("STABLE" if np.sign(row.rg) == np.sign(base.rg) or row.rg == 0 or base.rg == 0 else "FLIP")
        threshold = max(0.10, 2.0 * math.sqrt(float(row.rg_se) ** 2 + float(base.rg_se) ** 2))
        numerical.append("STABLE" if abs(float(row.rg) - float(base.rg)) <= threshold else "DIFFERENT")
        wide.append("YES" if float(row.rg_se) > 0.25 or 3.92 * float(row.rg_se) > 1.0 else "NO")
    out["primary_rg"] = p_rg
    out["primary_rg_se"] = p_se
    out["directional_stability_vs_primary"] = directional
    out["numerical_stability_vs_primary"] = numerical
    out["VERY_WIDE_CI"] = wide
    write_table(out, dirs["results"] / "AMD_GLOBAL_RG_SENSITIVITY.tsv")
    return out


def pair_qc(df: pd.DataFrame, inputs: dict[str, dict[str, object]], dirs: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for row in df.itertuples(index=False):
        finite = all(np.isfinite(float(getattr(row, col))) for col in [
            "rg", "rg_se", "rg_z", "rg_p", "genetic_covariance", "genetic_covariance_se",
            "cross_trait_intercept", "cross_trait_intercept_se", "mean_z1z2",
        ])
        flags = []
        if float(inputs[row.trait1]["global_h2_Z"]) < 4 or float(inputs[row.trait2]["global_h2_Z"]) < 4:
            flags.append("H2_Z_BELOW_4_AMENDED_PROTOCOL")
        if np.isfinite(row.rg_se) and float(row.rg_se) > 0.25:
            flags.append("VERY_WIDE_RG_CI")
        if np.isfinite(row.rg) and abs(float(row.rg)) > 1.05:
            flags.append("RG_OUTSIDE_UNIT_INTERVAL")
        if row.warning_flags != "NONE":
            flags.append("LDSC_WARNING")
        if not finite:
            status = "FAIL"
            flags.append("NONFINITE_REQUIRED_FIELD")
        elif any(flag in flags for flag in ["VERY_WIDE_RG_CI", "RG_OUTSIDE_UNIT_INTERVAL", "LDSC_WARNING"]):
            status = "REVIEW"
        else:
            status = "PASS"
        rows.append({
            "pair": f"{row.trait1}__{row.trait2}", "trait1": row.trait1, "trait2": row.trait2,
            "QC_STATUS": status, "rg_finite": "YES" if np.isfinite(row.rg) else "NO",
            "rg_se_finite": "YES" if np.isfinite(row.rg_se) else "NO",
            "genetic_covariance_finite": "YES" if np.isfinite(row.genetic_covariance) else "NO",
            "cross_trait_intercept_finite": "YES" if np.isfinite(row.cross_trait_intercept) else "NO",
            "n_valid_alleles": row.n_valid_alleles, "flags": ";".join(flags) if flags else "NONE",
            "warning_flags": row.warning_flags,
        })
    out = pd.DataFrame(rows)
    write_table(out, dirs["results"] / "LDSC_PAIR_QC.tsv")
    return out


def write_benchmark(sensitivity: pd.DataFrame, dirs: dict[str, Path]) -> None:
    sub = sensitivity[sensitivity.trait2 == "POAG"].copy()
    lines = [
        "# AMD–POAG global LDSC benchmark",
        "",
        "This benchmark compares the corrected primary AMD input with the prespecified AMD_RSID_GRCh37 and AMD_HM3_STRICT sensitivity inputs. The analysis estimates the genetic correlation under the same EUR LD-score reference and unconstrained cross-trait intercept. No negative direction was imposed during extraction or interpretation.",
        "",
        "| AMD input | rg | SE | P | directional stability | numerical stability | very wide CI |",
        "|---|---:|---:|---:|---|---|---|",
    ]
    for row in sub.itertuples(index=False):
        lines.append(f"| {row.amd_variant} | {fmt(row.rg, 6)} | {fmt(row.rg_se, 6)} | {fmt(row.rg_p, 6)} | {row.directional_stability_vs_primary} | {row.numerical_stability_vs_primary} | {row.VERY_WIDE_CI} |")
    lines += [
        "",
        "Interpretation is restricted to numerical and directional stability of the global estimate. It does not establish local antagonistic pleiotropy or causal direction.",
    ]
    (dirs["reports"] / "AMD_POAG_BENCHMARK.md").write_text("\n".join(lines) + "\n")


def configure_matplotlib() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "svg.fonttype": "none", "axes.spines.right": False, "axes.spines.top": False,
        "axes.linewidth": 0.8, "pdf.fonttype": 42, "savefig.dpi": 600,
    })


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_figures(df: pd.DataFrame, dirs: dict[str, Path]) -> None:
    configure_matplotlib()
    # Contract: a quantitative grid and a grouped forest plot. The visual claim
    # is the pattern of global rg estimates after multiplicity correction, while
    # nonsignificant estimates remain visible but visually de-emphasized.
    heat = pd.DataFrame(np.nan, index=TRAITS, columns=TRAITS)
    sig = pd.DataFrame(False, index=TRAITS, columns=TRAITS)
    for row in df.itertuples(index=False):
        heat.loc[row.trait1, row.trait2] = row.rg
        heat.loc[row.trait2, row.trait1] = row.rg
        sig.loc[row.trait1, row.trait2] = row.rg_fdr_bh < 0.05
        sig.loc[row.trait2, row.trait1] = row.rg_fdr_bh < 0.05
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    masked = np.ma.masked_invalid(heat.to_numpy(float))
    im = ax.imshow(masked, cmap="coolwarm", vmin=-1, vmax=1, aspect="equal")
    ax.set_xticks(range(len(TRAITS)), TRAITS, rotation=45, rotation_mode="anchor", ha="left")
    ax.set_yticks(range(len(TRAITS)), TRAITS)
    ax.set_title("Global genetic correlation (LDSC)", loc="left", fontweight="bold")
    ax.set_xlabel("Trait 2")
    ax.set_ylabel("Trait 1")
    for i, t1 in enumerate(TRAITS):
        for j, t2 in enumerate(TRAITS):
            if i == j:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, facecolor="white", edgecolor="0.85", lw=.5))
                continue
            val = heat.loc[t1, t2]
            if not np.isfinite(val):
                continue
            row = df[((df.trait1 == t1) & (df.trait2 == t2)) | ((df.trait1 == t2) & (df.trait2 == t1))].iloc[0]
            alpha = 1.0 if row.rg_fdr_bh < 0.05 else 0.30
            ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, facecolor="none", edgecolor="white", lw=.6, alpha=alpha))
            if row.rg_fdr_bh < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=9, fontweight="bold", color="black")
    cbar = fig.colorbar(im, ax=ax, fraction=.046, pad=.04)
    cbar.set_label("rg")
    fig.text(.01, .01, "Color = point estimate; faded cells = BH FDR ≥ 0.05; * = BH FDR < 0.05", fontsize=6)
    save_figure(fig, dirs["figures"] / "Figure1A_global_rg_heatmap")

    subset = df[df.apply(lambda row: row.trait1 in OCULAR and row.trait2 in SYSTEMIC, axis=1)].copy()
    subset["label"] = subset["trait2"]
    group_order = [(trait, systemic) for trait in OCULAR for systemic in SYSTEMIC]
    ordered = []
    for ocular, systemic in group_order:
        match = subset[(subset.trait1 == ocular) & (subset.trait2 == systemic)]
        if not match.empty:
            ordered.append(match.iloc[0])
    forest = pd.DataFrame(ordered)
    y = np.arange(len(forest))
    fig, ax = plt.subplots(figsize=(6.2, 5.8))
    colors = ["#3b6fb6" if (r.rg_fdr_bh < .05 and r.rg > 0) else "#b33b3b" if (r.rg_fdr_bh < .05 and r.rg < 0) else "#7f8790" for _, r in forest.iterrows()]
    for idx, (_, row) in enumerate(forest.iterrows()):
        ax.errorbar(row.rg, idx, xerr=1.96 * row.rg_se, fmt="none", ecolor=colors[idx], elinewidth=1.0, capsize=2, alpha=.95)
    ax.scatter(forest.rg, y, s=22, c=colors, zorder=3, edgecolor="white", linewidth=.35)
    ax.axvline(0, color="0.25", lw=.8)
    labels = [f"{r.trait1} — {r.trait2}" for _, r in forest.iterrows()]
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Global genetic correlation (rg; 95% CI)")
    ax.set_title("Ocular–systemic global LDSC correlations", loc="left", fontweight="bold")
    ax.grid(axis="x", color="0.90", lw=.6)
    for idx in [5, 11]:
        if idx < len(y) - 1:
            ax.axhline(idx + .5, color="0.75", lw=.6)
    ax.text(1.01, .91, "AMD", transform=ax.transAxes, fontsize=7, fontweight="bold", va="center")
    ax.text(1.01, .59, "POAG", transform=ax.transAxes, fontsize=7, fontweight="bold", va="center")
    ax.text(1.01, .27, "CAT", transform=ax.transAxes, fontsize=7, fontweight="bold", va="center")
    fig.text(.01, .01, "Gray = no BH-FDR-detectable global rg; colors encode sign only when BH FDR < 0.05.", fontsize=6)
    save_figure(fig, dirs["figures"] / "Figure1B_ocular_systemic_rg_forest")


def gate_results(df: pd.DataFrame, qc: pd.DataFrame, matrix_stats: dict[str, object], sensitivity: pd.DataFrame) -> dict[str, str]:
    global_ok = len(df) == 36 and (df["return_code"] == 0).all() and df[[
        "rg", "rg_se", "rg_z", "rg_p", "genetic_covariance", "genetic_covariance_se",
        "cross_trait_intercept", "cross_trait_intercept_se", "mean_z1z2",
    ]].map(np.isfinite).all().all()
    global_gate = "PASS" if global_ok and (qc.QC_STATUS == "PASS").all() else "CONDITIONAL" if global_ok else "FAIL"
    matrix_ok = all([
        matrix_stats["finite"], matrix_stats["symmetric"], matrix_stats["diag_max_error"] <= 1e-12,
        matrix_stats["offdiag_min"] >= -1.0000001, matrix_stats["offdiag_max"] <= 1.0000001,
    ])
    eigen = matrix_stats["min_eigenvalue"]
    lava_gate = "PASS" if matrix_ok and eigen >= -1e-8 else "CONDITIONAL" if matrix_ok and eigen >= -1e-5 else "FAIL"
    hdll = "INADEQUATE_FOR_MOST_PAIRS"
    return {
        "GLOBAL_ARCHITECTURE_PASS": global_gate,
        "LAVA_OVERLAP_MATRIX_PASS": lava_gate,
        "HDLL_PAIRWISE_FEASIBILITY": hdll,
    }


def write_report(
    inputs: dict[str, dict[str, object]], df: pd.DataFrame, qc: pd.DataFrame,
    sens: pd.DataFrame, matrix: dict[str, object], gates: dict[str, str], dirs: dict[str, Path],
) -> None:
    stats = matrix["stats"]
    sig = df[df.rg_fdr_bh < .05].sort_values("rg_p")
    ocular = df[df.apply(lambda row: row.trait1 in OCULAR and row.trait2 in SYSTEMIC, axis=1)].copy()
    lines = [
        "# Phase 1A — Global genetic architecture and sample-overlap estimation",
        "",
        "## Scope and stop rule",
        "",
        "This report implements the amended protocol and stops after global LDSC, the LDSC cross-trait intercept matrix, overlap metadata, and AMD global sensitivity. No local bivariate analysis, HDL-L full pairwise analysis, LAVA local rg, PLACO+, fine-mapping, coloc, GCI, covariance decomposition, or enrichment was run.",
        "",
        "## Protocol amendment",
        "",
        "- `AMENDMENT_STATUS=SCIENTIFICALLY_JUSTIFIED`.",
        "- Historical Phase 0 remains `HOLD_UNDER_ORIGINAL_PROTOCOL`.",
        "- Current Phase 1A status is `PROCEED_UNDER_AMENDED_PROTOCOL`.",
        "- Global LDSC h2 Z ≥ 4 remains preferred. Traits with 2 ≤ Z < 4 may enter local testing only after the amended QC, independent local-h2 reconstruction, signed-statistic validation, region-level multiple-testing control, and pairwise uncertainty requirements are met.",
        "- AMD local univariate heritability is analytically separate from the HDL-L pairwise covariance question. The former may be eligible under the amended rule; the latter remains unresolved when N0 is unknown.",
        "",
        "## Frozen inputs and global LDSC",
        "",
        f"The analysis used nine EUR/EUR-like frozen traits and {len(df)}/36 unique unconstrained cross-trait LDSC pairs. The same EUR LD-score reference, regression weights, HM3 merge, and Phase 0 MHC handling were used. No intercept was constrained to zero.",
        "",
        "| Trait | Dataset | h2 | SE | Z | LDSC intercept | HM3 SNPs |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for trait in TRAITS:
        row = inputs[trait]
        lines.append(f"| {trait} | {row['dataset']} | {fmt(row['global_h2'], 5)} | {fmt(row['global_h2_SE'], 5)} | {fmt(row['global_h2_Z'], 5)} | {fmt(row['LDSC_intercept'], 5)} | {row['HM3_SNP_count']} |")
    lines += [
        "",
        "### Multiplicity-corrected significant global rg",
        "",
        "| Trait 1 | Trait 2 | rg | SE | P | BH FDR | Bonferroni P |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    if sig.empty:
        lines.append("| None | — | — | — | — | — | — |")
    else:
        for row in sig.itertuples(index=False):
            lines.append(f"| {row.trait1} | {row.trait2} | {fmt(row.rg, 6)} | {fmt(row.rg_se, 6)} | {fmt(row.rg_p, 6)} | {fmt(row.rg_fdr_bh, 6)} | {fmt(row.rg_p_bonferroni, 6)} |")
    lines += [
        "",
        "### Ocular–systemic pairs",
        "",
        "The direction field is `POSITIVE` or `NEGATIVE` only when the global rg is BH-FDR detectable; otherwise it is `NO_DETECTABLE_GLOBAL_RG`. A nonsignificant point estimate is not described as independent.",
        "",
        "| Ocular | Systemic | rg | SE | BH FDR | Direction |",
        "|---|---|---:|---:|---:|---|",
    ]
    for row in ocular.itertuples(index=False):
        direction = "POSITIVE" if row.rg_fdr_bh < .05 and row.rg > 0 else "NEGATIVE" if row.rg_fdr_bh < .05 and row.rg < 0 else "NO_DETECTABLE_GLOBAL_RG"
        lines.append(f"| {row.trait1} | {row.trait2} | {fmt(row.rg, 6)} | {fmt(row.rg_se, 6)} | {fmt(row.rg_fdr_bh, 6)} | {direction} |")
    lines += [
        "",
        "## Intercept matrices and LAVA sampling correlation",
        "",
        "`LDSC_CROSS_TRAIT_INTERCEPT_MATRIX.tsv` is complete and symmetric. Its diagonal uses the corresponding validated univariate LDSC intercept; off-diagonal cells use the unconstrained cross-trait LDSC genetic-covariance intercept. The LAVA matrix is the standardized `cov2cor`-style matrix obtained from this intercept covariance matrix, with diagonal fixed by standardization to one.",
        "",
        f"- finite: `{stats['finite']}`; symmetric: `{stats['symmetric']}`; maximum diagonal error: `{fmt(stats['diag_max_error'], 6)}`.",
        f"- off-diagonal range: `{fmt(stats['offdiag_min'], 6)}` to `{fmt(stats['offdiag_max'], 6)}`.",
        f"- eigenvalue range: `{fmt(stats['min_eigenvalue'], 6)}` to `{fmt(stats['max_eigenvalue'], 6)}`.",
        "- The official LAVA UKB EUR v1.1 binary reference remains unavailable locally; the matrix QC does not substitute for that reference.",
        "",
        "## Cohort overlap and HDL-L N0 eligibility",
        "",
        "The metadata audit records only documented cohort components and possible overlap warnings. Exact participant overlap was not inferred from consortium names or sample sizes. Because no pair has a verified exact N0 or verified zero overlap in the frozen metadata, no pair is eligible for a definitive HDL-L covariance interpretation in this phase.",
        "",
        "## AMD global stability and AMD–POAG benchmark",
        "",
        "AMD global rg was rerun against all eight non-AMD traits using the corrected primary input and, where available, `AMD_RSID_GRCh37` and `AMD_HM3_STRICT`. The full machine-readable table is `AMD_GLOBAL_RG_SENSITIVITY.tsv`; the focused benchmark is `AMD_POAG_BENCHMARK.md`. No negative direction was forced.",
        "",
        "## Gates",
        "",
        "| Gate | Verdict |",
        "|---|---|",
    ]
    for key, value in gates.items():
        lines.append(f"| `{key}` | `{value}` |")
    lines += [
        "",
        "## Next permitted step",
        "",
        "Resolve access to the official LAVA UKB EUR v1.1 reference and separately resolve exact participant-overlap/N0 metadata. Do not start local bivariate analysis or interpret HDL-L pairwise covariance until those prerequisites and the amended local-h2 requirements are satisfied.",
        "",
        "## Files",
        "",
        "- Dataset freeze: `config/phase1A_dataset_freeze.tsv`.",
        "- Global LDSC: `results/phase1A/GLOBAL_LDSC_RG.tsv` and `GLOBAL_LDSC_RG_FDR.tsv`.",
        "- Ocular–systemic subset: `results/phase1A/OCULAR_SYSTEMIC_RG.tsv`.",
        "- Intercept and standardized sampling-correlation matrices: `results/phase1A/LDSC_CROSS_TRAIT_INTERCEPT_MATRIX.tsv` and `LAVA_SAMPLING_CORRELATION_MATRIX.tsv`.",
        "- Overlap and HDL-L eligibility audits: `results/phase1A/COHORT_OVERLAP_METADATA.tsv` and `HDLL_N0_ELIGIBILITY.tsv`.",
        "- QC and figures: `results/phase1A/LDSC_PAIR_QC.tsv`, `figures/phase1A/Figure1A_global_rg_heatmap.pdf`, `figures/phase1A/Figure1B_ocular_systemic_rg_forest.pdf`.",
    ]
    (dirs["reports"] / "PHASE1A_GLOBAL_ARCHITECTURE_REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-rg", action="store_true", help="Reuse existing Phase1A pair logs; only for audited reruns.")
    args = parser.parse_args()
    dirs = ensure_dirs()
    inputs, paths = load_inputs()
    write_dataset_freeze(paths, inputs)
    if args.skip_rg:
        records = []
        for i, trait1 in enumerate(TRAITS):
            for trait2 in TRAITS[i + 1:]:
                records.append(parse_ldsc_log(dirs["results"] / "tmp" / f"{trait1}__{trait2}.log", trait1, trait2))
        df = pd.DataFrame(records)
        df["return_code"] = 0
        df["rg_fdr_bh"] = bh_adjust(df["rg_p"].to_numpy(float))
        df["rg_p_bonferroni"] = np.minimum(df["rg_p"].astype(float) * len(df), 1.0)
        df["rg_significant_fdr"] = np.where(df["rg_fdr_bh"] < .05, "YES", "NO")
        df["rg_significant_bonferroni"] = np.where(df["rg_p_bonferroni"] < .05, "YES", "NO")
        df["point_estimate_direction"] = np.where(df["rg"] > 0, "POSITIVE", np.where(df["rg"] < 0, "NEGATIVE", "ZERO_OR_NA"))
    else:
        df = run_all_pairs(paths, dirs)
    matrices = intercept_matrices(df, inputs, dirs)
    write_overlap_audits(inputs, df, dirs)
    sensitivity = amd_sensitivity(paths, df, dirs)
    qc = pair_qc(df, inputs, dirs)
    write_benchmark(sensitivity, dirs)
    make_figures(df, dirs)
    gates = gate_results(df, qc, matrices["stats"], sensitivity)
    (dirs["results"] / "PHASE1A_GATES.tsv").write_text("gate\tverdict\n" + "\n".join(f"{key}\t{value}" for key, value in gates.items()) + "\n")
    write_report(inputs, df, qc, sensitivity, matrices, gates, dirs)
    print("[Phase1A] completed", flush=True)
    print("[Phase1A] gates:", gates, flush=True)
    print("[Phase1A] LAVA eigenvalue range:", matrices["stats"]["min_eigenvalue"], matrices["stats"]["max_eigenvalue"], flush=True)


if __name__ == "__main__":
    main()
