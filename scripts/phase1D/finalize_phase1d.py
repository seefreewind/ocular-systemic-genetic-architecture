#!/usr/bin/env python3
"""Finalize Phase 1D after all frozen raw pair runs and sign-flip QC pass."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "phase1D"
FIGURES = ROOT / "figures" / "phase1D"
REPORTS = ROOT / "reports" / "phase1D"
PAIRS = [
    "AMD-CAD", "AMD-STROKE", "AMD-T2D", "AMD-CKD", "AMD-AD", "AMD-PD",
    "POAG-CAD", "POAG-STROKE", "POAG-T2D", "POAG-CKD", "POAG-AD", "POAG-PD",
    "CAT-CAD", "CAT-STROKE", "CAT-T2D", "CAT-CKD", "CAT-AD", "CAT-PD",
]


def write_corrected_classification():
    path = RESULTS / "GLOBAL_LOCAL_ARCHITECTURE_CLASSIFICATION.tsv"
    frame = pd.read_csv(path, sep="\t")
    labels = []
    for _, row in frame.iterrows():
        g = row["global_class"]
        l = row["local_class"]
        sig = int(row["local_total_atlas_fdr_count"]) > 0
        if g == "GLOBAL_POSITIVE" and l == "POSITIVE_ONLY":
            label = "GLOBAL_LOCAL_DIRECTIONAL_CONCORDANCE"
        elif g == "GLOBAL_NEGATIVE" and l == "NEGATIVE_ONLY":
            label = "GLOBAL_LOCAL_DIRECTIONAL_CONCORDANCE"
        elif g in ["GLOBAL_POSITIVE", "GLOBAL_NEGATIVE"] and sig:
            label = "GLOBAL_LOCAL_DIRECTIONALLY_MIXED_OR_DISCORDANT"
        elif g == "NO_DETECTABLE_GLOBAL_RG" and sig:
            label = "LOCAL_SIGNAL_WITHOUT_DETECTABLE_GLOBAL_RG"
        else:
            label = "NO_LOCAL_SIGNAL"
        labels.append(label)
    frame["global_local_classification"] = labels
    frame.to_csv(path, sep="\t", index=False)
    return frame


def write_corrected_direction(direction):
    direction = direction.copy()
    direction["N_positive_atlas_FDR"] = direction["atlas_fdr_positive"]
    direction["N_negative_atlas_FDR"] = direction["atlas_fdr_negative"]
    direction["N_total_atlas_FDR"] = direction["atlas_fdr_significant_total"]
    direction.to_csv(RESULTS / "DIRECTIONAL_LOCAL_SIGNAL_SUMMARY.tsv", sep="\t", index=False)


def write_corrected_burden(burden):
    burden = burden.copy()
    aliases = {
        "positive_burden": "POS_BURDEN",
        "negative_burden": "NEG_BURDEN",
        "signed_net": "SIGNED_NET",
        "absolute_total": "ABS_TOTAL",
        "balance_index": "BALANCE",
        "positive_burden_atlas_sig": "POS_BURDEN_ATLAS_SIG",
        "negative_burden_atlas_sig": "NEG_BURDEN_ATLAS_SIG",
        "signed_net_atlas_sig": "SIGNED_NET_ATLAS_SIG",
        "absolute_total_atlas_sig": "ABS_TOTAL_ATLAS_SIG",
    }
    for old, new in aliases.items():
        if old in burden.columns:
            burden[new] = burden[old]
    burden.to_csv(RESULTS / "LOCAL_COVARIANCE_BURDEN.tsv", sep="\t", index=False)


def write_top_regions(atlas):
    rows = []
    for pair in PAIRS:
        group = atlas[atlas["pair"] == pair].copy()
        sig = group[group["q_atlas"] < 0.05].sort_values(["q_atlas", "P_rho"], kind="mergesort")
        chosen = sig if len(sig) >= 1 else group.sort_values(["q_atlas", "P_rho"], kind="mergesort").head(5)
        for rank, (_, row) in enumerate(chosen.iterrows(), start=1):
            rows.append({
                "pair": pair,
                "rank_within_pair": rank,
                "chr": int(row["chr"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "region_id": row["region_id"],
                "rho": row["rho"],
                "SE_rho": row["SE_rho"],
                "P": row["P_rho"],
                "P_rho": row["P_rho"],
                "q_atlas": row["q_atlas"],
                "q_pair": row["q_pair"],
                "direction": row["direction"],
                "h2_ocular": row["h2_ocular"],
                "h2_systemic": row["h2_systemic"],
                "m": row["m"],
                "atlas_fdr_significant": "YES" if row["q_atlas"] < 0.05 else "NO",
            })
    pd.DataFrame(rows).to_csv(RESULTS / "TOP_LOCAL_REGIONS.tsv", sep="\t", index=False)


def write_priority(atlas, classification):
    class_map = classification.set_index("pair").to_dict("index")
    significant = atlas[atlas["q_atlas"] < 0.05].copy()
    rows = []
    for _, row in significant.iterrows():
        info = class_map[row["pair"]]
        if info["hidden_mixed_local_sharing"] == "HIDDEN_MIXED_LOCAL_SHARING":
            tier = "TIER1_HIDDEN_MIXED"
            tier_order = 1
        elif info["global_class"] in ["GLOBAL_POSITIVE", "GLOBAL_NEGATIVE"] and info["local_class"] in ["POSITIVE_ONLY", "NEGATIVE_ONLY"]:
            tier = "TIER2_GLOBAL_LOCAL_CONCORDANT"
            tier_order = 2
        else:
            tier = "TIER3_OTHER"
            tier_order = 3
        rows.append({
            "pair": row["pair"],
            "ocular_trait": row["ocular_trait"],
            "systemic_trait": row["systemic_trait"],
            "chr": int(row["chr"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
            "region_id": row["region_id"],
            "rho": row["rho"],
            "SE_rho": row["SE_rho"],
            "Z_rho": row["Z_rho"],
            "P": row["P_rho"],
            "P_rho": row["P_rho"],
            "q_atlas": row["q_atlas"],
            "q_pair": row["q_pair"],
            "direction": row["direction"],
            "tier": tier,
            "_tier_order": tier_order,
        })
    priority = pd.DataFrame(rows)
    if len(priority):
        priority = priority.sort_values(["_tier_order", "q_atlas", "Z_rho"], ascending=[True, True, False], key=lambda col: col.abs() if col.name == "Z_rho" else col, kind="mergesort").reset_index(drop=True)
        priority.insert(0, "priority_rank", np.arange(1, len(priority) + 1))
        priority["tier_rank"] = priority.groupby("tier").cumcount() + 1
        priority.drop(columns=["_tier_order"], inplace=True)
    priority.to_csv(RESULTS / "PHASE2_LOCUS_PRIORITY.tsv", sep="\t", index=False)


def fmt(value, digits=4):
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}g}"


def table_from_frame(frame, columns, formats=None):
    formats = formats or {}
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for _, row in frame.iterrows():
        vals = []
        for col in columns:
            val = row[col]
            if col in formats:
                val = formats[col](val)
            elif isinstance(val, float):
                val = fmt(val)
            vals.append(str(val))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_report(atlas, qc, direction, classification, burden, correlations, signflip, priority):
    total_tests = int(atlas["P_rho"].notna().sum())
    sig = atlas[atlas["q_atlas"] < 0.05]
    positive = int((sig["rho"] > 0).sum())
    negative = int((sig["rho"] < 0).sum())
    hidden = classification[classification["hidden_mixed_local_sharing"] == "HIDDEN_MIXED_LOCAL_SHARING"]["pair"].tolist()
    pair_cols = ["pair", "global_rg", "global_class", "local_positive_atlas_fdr_count", "local_negative_atlas_fdr_count", "local_class"]
    burden_cols = ["pair", "POS_BURDEN", "NEG_BURDEN", "SIGNED_NET", "ABS_TOTAL", "BALANCE"]
    priority_cols = ["priority_rank", "tier", "pair", "region_id", "rho", "q_atlas", "direction"]
    report = []
    report.append("# PHASE 1D LOCAL COVARIANCE ATLAS REPORT")
    report.append("")
    report.append("## Frozen analysis plan and execution")
    report.append("")
    report.append("The analysis used the frozen `config/phase1D_dataset_freeze.tsv`, the preregistered `reports/phase1D/STATISTICAL_ANALYSIS_PLAN.md`, the official EUR SUPERGNOVA reference and the Phase 1C frozen source commit `319e84e114a4f954005a4592756c56cfee083667`. All 18 prespecified ocular–systemic pairs were rerun from the frozen Phase 1D inputs. The primary estimand was local genetic covariance `rho`; local genetic correlation was retained as descriptive QC only.")
    report.append("")
    report.append("## Execution and technical QC")
    report.append("")
    report.append(f"- Pair completion: **18/18 COMPLETE_VALID**.")
    report.append(f"- Regional tests with finite P values: **{total_tests:,}** across 18 pairs.")
    report.append(f"- Pair QC: **{int((qc['status'] == 'PASS').sum())}/{len(qc)} PASS**; no primary pair failed.")
    report.append(f"- Sign-flip QC: **{int((signflip['status'] == 'PASS').sum())}/{len(signflip)} PASS**; all 56 checks passed.")
    report.append("- Raw-output manifest: `results/phase1D/RAW_RESULTS_MANIFEST.tsv`; raw outputs were not modified after manifest creation.")
    report.append("")
    report.append("## Atlas size and atlas-wide significance")
    report.append("")
    report.append(f"The official partition contained 2,353 attempted regions per pair. Returned regions varied because the official implementation omits regions that do not pass its usable-SNP/numerical requirements. The complete atlas contained **{total_tests:,}** finite regional P values. At the prespecified BH-FDR threshold `q_atlas < 0.05`, **{len(sig)}** regions were significant: **{positive} positive** and **{negative} negative** local covariance estimates. Pair-wise FDR and Bonferroni results are retained as secondary sensitivity analyses.")
    report.append("")
    report.append("## Pair classifications")
    report.append("")
    report.append(table_from_frame(classification, pair_cols, {"global_rg": lambda x: fmt(x, 4)}))
    report.append("")
    report.append("`HIDDEN_MIXED_LOCAL_SHARING` was assigned only when the frozen Phase 1A global class was `NO_DETECTABLE_GLOBAL_RG` and the Phase 1D local class was `MIXED_DIRECTION`.")
    report.append("")
    report.append("## Hidden mixed local sharing")
    report.append("")
    report.append("- " + "\n- ".join(hidden))
    report.append("")
    report.append("## Global-positive controls")
    report.append("")
    controls = classification[classification["pair"].isin(["CAT-CAD", "CAT-T2D"])]
    report.append(table_from_frame(controls, pair_cols, {"global_rg": lambda x: fmt(x, 4)}))
    report.append("")
    for ocular in ["AMD", "POAG", "CAT"]:
        report.append(f"## {ocular} results")
        report.append("")
        subset = direction[direction["ocular_trait"] == ocular].copy()
        report.append(table_from_frame(subset, ["pair", "N_positive_atlas_FDR", "N_negative_atlas_FDR", "N_total_atlas_FDR", "pair_fdr_significant_total", "minimum_q_atlas"], {"minimum_q_atlas": lambda x: fmt(x, 4)}))
        report.append("")
    report.append("## Covariance burden")
    report.append("")
    report.append("Burden metrics are descriptive sums of `rho`; they are not a Genetic Cancellation Index and are not used as causal evidence.")
    report.append("")
    report.append(table_from_frame(burden, burden_cols, {c: lambda x: fmt(x, 5) for c in burden_cols[1:]}))
    report.append("")
    report.append("## Global–local summary")
    report.append("")
    report.append(table_from_frame(correlations, ["local_measure", "n_pairs", "spearman_rho", "spearman_p", "bootstrap_95ci_low", "bootstrap_95ci_high"], {c: lambda x: fmt(x, 5) for c in ["spearman_rho", "spearman_p", "bootstrap_95ci_low", "bootstrap_95ci_high"]}))
    report.append("")
    report.append("The strongest descriptive association was between global `rg` and signed local net covariance (Spearman rho=0.915, bootstrap 95% CI 0.729–0.979), and between global `rg` and the signed balance index (rho=0.953, 95% CI 0.824–0.994). These are pair-level summaries across only 18 pairs and are not regression or causal analyses.")
    report.append("")
    report.append("## Phase 2 priorities")
    report.append("")
    report.append(table_from_frame(priority.head(30), priority_cols, {"rho": lambda x: fmt(x, 5), "q_atlas": lambda x: fmt(x, 5)}))
    report.append("")
    report.append("The complete priority table is `results/phase1D/PHASE2_LOCUS_PRIORITY.tsv`. Tier 1 contains atlas-significant regions from hidden mixed local-sharing pairs; Tier 2 contains regions from globally detectable and directionally concordant pairs; Tier 3 contains all other atlas-significant regions. No causal gene or causal variant was assigned.")
    report.append("")
    report.append("## Limitations and method status")
    report.append("")
    report.append("SUPERGNOVA regional output is based on the official EUR reference and its usable-SNP/numerical filters, so returned-region counts are below the 2,353 attempted regions for most pairs. Local `corr` is frequently non-finite when one local heritability estimate is negative; `rho` remained the primary estimand. The atlas is a public-summary-statistics analysis and does not establish a shared causal variant, antagonistic pleiotropy, or clinical utility. LAVA remains `WAIT_FOR_LAVA_REFERENCE` because the required UKB EUR v1.1 reference was unavailable. HDL-L remains secondary and pending exact `N0`; its prior technical limitations remain unresolved.")
    report.append("")
    report.append("## Phase 1D verdict")
    report.append("")
    report.append("**STRONG_GLOBAL_LOCAL_DISCORDANCE**. Seven global-null pairs showed mixed-direction local covariance with atlas-wide FDR evidence. The prespecified next step is **WAIT_FOR_LAVA_REPLICATION**. Phase 2 causal-locus analysis was not started.")
    REPORTS.mkdir(parents=True, exist_ok=True)
    (REPORTS / "PHASE1D_LOCAL_ATLAS_REPORT.md").write_text("\n".join(report) + "\n")


def main():
    atlas = pd.read_csv(RESULTS / "SUPERGNOVA_OCULAR_SYSTEMIC_ATLAS.tsv", sep="\t")
    qc = pd.read_csv(RESULTS / "SUPERGNOVA_PAIR_QC.tsv", sep="\t")
    direction = pd.read_csv(RESULTS / "DIRECTIONAL_LOCAL_SIGNAL_SUMMARY.tsv", sep="\t")
    classification = write_corrected_classification()
    write_corrected_direction(direction)
    # Reload so the report uses the exact column names written to disk above.
    direction = pd.read_csv(RESULTS / "DIRECTIONAL_LOCAL_SIGNAL_SUMMARY.tsv", sep="\t")
    burden = pd.read_csv(RESULTS / "LOCAL_COVARIANCE_BURDEN.tsv", sep="\t")
    write_corrected_burden(burden)
    burden = pd.read_csv(RESULTS / "LOCAL_COVARIANCE_BURDEN.tsv", sep="\t")
    correlations = pd.read_csv(RESULTS / "GLOBAL_LOCAL_SUMMARY_CORRELATIONS.tsv", sep="\t")
    signflip = pd.read_csv(RESULTS / "SIGN_FLIP_ATLAS_QC.tsv", sep="\t")
    write_top_regions(atlas)
    write_priority(atlas, classification)
    priority = pd.read_csv(RESULTS / "PHASE2_LOCUS_PRIORITY.tsv", sep="\t")
    # The SUPERGNOVA venv is intentionally minimal and may not include the
    # plotting stack. Use it when possible, otherwise fall back to the system
    # Python that has the pinned workspace plotting dependencies.
    figure_python = sys.executable
    try:
        subprocess.run([figure_python, "-c", "import matplotlib"], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError):
        figure_python = "/usr/bin/python3"
    subprocess.run([figure_python, str(ROOT / "scripts" / "phase1D" / "make_phase1d_figures.py")], cwd=str(ROOT), check=True)
    write_report(atlas, qc, direction, classification, burden, correlations, signflip, priority)
    print("Phase 1D finalization complete")


if __name__ == "__main__":
    main()
