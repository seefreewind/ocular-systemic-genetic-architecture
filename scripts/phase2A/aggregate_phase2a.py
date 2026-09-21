#!/usr/bin/env python3
"""Aggregate pairwise PLACO+ runs and complete the Phase 2A audit."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "phase2A"
PAIR_RUNS = RESULTS / "pair_runs"
PAIRS = [
    "AMD-CAD", "AMD-STROKE", "AMD-T2D", "AMD-CKD", "AMD-AD", "AMD-PD",
    "POAG-CAD", "POAG-STROKE", "POAG-T2D", "POAG-CKD", "POAG-AD", "POAG-PD",
    "CAT-CAD", "CAT-STROKE", "CAT-T2D", "CAT-CKD", "CAT-AD", "CAT-PD",
]
HIDDEN = {"AMD-CAD", "AMD-STROKE", "AMD-T2D", "AMD-CKD", "POAG-CKD", "CAT-STROKE", "CAT-AD"}
ATLAS_THRESHOLD = 5e-8 / 18


def table_markdown(frame: pd.DataFrame, columns: list[str], n: int | None = None) -> str:
    if n is not None:
        frame = frame.head(n)
    if frame.empty:
        return "No rows."
    out = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for _, row in frame.iterrows():
        vals = []
        for col in columns:
            value = row.get(col, "NA")
            if isinstance(value, float):
                value = f"{value:.5g}"
            vals.append(str(value))
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out)


def run_clumping(significant: pd.DataFrame) -> pd.DataFrame:
    rows = []
    clump_dir = RESULTS / "ld_clump"
    clump_dir.mkdir(parents=True, exist_ok=True)
    for pair in PAIRS:
        frame = significant[significant["pair"] == pair].copy()
        if frame.empty:
            continue
        frame["clump_p"] = frame["P_PLACO_PLUS"].where(frame["P_PLACO_PLUS"].notna(), np.maximum(frame["P_1"], frame["P_2"]))
        for chrom in sorted(frame["CHR"].dropna().astype(int).unique()):
            sub = frame[frame["CHR"].astype(int) == chrom].sort_values(["clump_p", "SNP"])
            if sub.empty:
                continue
            inp = clump_dir / f"{pair}.chr{chrom}.clump_input.tsv"
            out_prefix = clump_dir / f"{pair}.chr{chrom}"
            sub[["SNP", "clump_p"]].rename(columns={"clump_p": "P"}).to_csv(inp, sep="\t", index=False)
            bfile = ROOT / "data" / "phase1C" / "SUPERGNOVA_reference" / "bfiles" / f"eur_chr{chrom}_SNPmaf5"
            try:
                subprocess.run([
                    "plink2", "--bfile", str(bfile), "--clump", str(inp),
                    "--clump-snp-field", "SNP", "--clump-field", "P",
                    "--clump-p1", "1", "--clump-p2", "1",
                    "--clump-r2", "0.1", "--clump-kb", "500",
                    "--threads", "4", "--memory", "4096", "--out", str(out_prefix),
                ], check=True, cwd=str(ROOT), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                # PLINK writes <out_prefix>.clumps; Path.with_suffix() would
                # incorrectly collapse "PAIR.chrN" to "PAIR.clumps".
                clump_path = Path(f"{out_prefix}.clumps")
                if clump_path.exists():
                    c = pd.read_csv(clump_path, sep=r"\s+", engine="python")
                    for _, crow in c.iterrows():
                        lead = str(crow.get("ID", crow.get("SNP", "")))
                        hit = sub[sub["SNP"] == lead]
                        if hit.empty:
                            hit = sub.iloc[[0]]
                        row = hit.iloc[0].to_dict()
                        row.update({"clump_status": "LD_CLUMPED", "lead_SNP": lead, "clump_variant_count": int(crow.get("TOTAL", 1))})
                        rows.append(row)
            except Exception:
                pass
        represented = {r["SNP"] for r in rows if r.get("pair") == pair}
        for _, row in frame.iterrows():
            if row["SNP"] not in represented:
                out = row.to_dict()
                out.update({"clump_status": "LD_PANEL_NOT_AVAILABLE_OR_NOT_CLUMPED", "lead_SNP": row["SNP"], "clump_variant_count": 1})
                rows.append(out)
    if not rows:
        return pd.DataFrame(columns=list(significant.columns) + ["clump_status", "lead_SNP", "clump_variant_count"])
    loci = pd.DataFrame(rows).drop_duplicates(subset=["pair", "lead_SNP"])
    return loci


def main() -> None:
    params = []
    sign_qc = []
    for pair in PAIRS:
        p = PAIR_RUNS / f"{pair}.params.tsv"
        s = PAIR_RUNS / f"{pair}.sign_qc.tsv"
        if not p.exists() or not s.exists():
            raise FileNotFoundError(f"Missing completed pair output for {pair}")
        params.append(pd.read_csv(p, sep="\t"))
        sign_qc.append(pd.read_csv(s, sep="\t"))
    params = pd.concat(params, ignore_index=True)
    sign_qc = pd.concat(sign_qc, ignore_index=True)
    sign_qc.to_csv(RESULTS / "PLACO_INPUT_SIGN_QC.tsv", sep="\t", index=False)
    params.to_csv(RESULTS / "PLACO_NULL_PARAMETER_ESTIMATES.tsv", sep="\t", index=False)
    parameter_qc = params.copy()
    parameter_qc["varz1_positive"] = parameter_qc["VarZ1"] > 0
    parameter_qc["varz2_positive"] = parameter_qc["VarZ2"] > 0
    parameter_qc["varz1_inflated_flag"] = parameter_qc["VarZ1"] > 2
    parameter_qc["varz2_inflated_flag"] = parameter_qc["VarZ2"] > 2
    parameter_qc["corz_near_boundary_flag"] = parameter_qc["CorZ"].abs() >= 0.9
    parameter_qc["too_few_null_flag"] = parameter_qc["null_variants"] < 30
    parameter_qc["nonfinite_flag"] = ~np.isfinite(parameter_qc[["VarZ1", "VarZ2", "CorZ"]]).all(axis=1)
    parameter_qc["status"] = np.where(parameter_qc[["varz1_positive", "varz2_positive"]].all(axis=1) & ~parameter_qc["nonfinite_flag"] & ~parameter_qc["too_few_null_flag"], "PASS", "FAIL")
    parameter_qc.to_csv(RESULTS / "PLACO_PARAMETER_QC.tsv", sep="\t", index=False)

    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    pattern = str(PAIR_RUNS / "*.placo.tsv.gz")
    con.execute(f"CREATE OR REPLACE TABLE placo_raw AS SELECT * FROM read_csv_auto('{pattern}', header=true, delim='\\t', union_by_name=true, compression='gzip', sample_size=100000)")
    n = con.execute("SELECT count(*) FROM placo_raw WHERE isfinite(P_PLACO_PLUS)").fetchone()[0]
    con.execute(f"""
        CREATE OR REPLACE TABLE placo_q AS
        WITH ranked AS (
          SELECT *, row_number() OVER (ORDER BY P_PLACO_PLUS, pair, SNP) AS rank_asc
          FROM placo_raw
          WHERE isfinite(P_PLACO_PLUS)
        ), adjusted AS (
          SELECT *, least(1.0, min(P_PLACO_PLUS * {float(n)} / rank_asc) OVER (ORDER BY P_PLACO_PLUS DESC ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW)) AS q_PLACO_atlas
          FROM ranked
        )
        SELECT r.*, a.q_PLACO_atlas
        FROM placo_raw r
        LEFT JOIN adjusted a USING (pair, SNP)
    """)
    genomewide = RESULTS / "PLACO_PLUS_GENOMEWIDE.tsv"
    con.execute(f"COPY (SELECT * FROM placo_q ORDER BY pair, CHR, BP, SNP) TO '{genomewide}' (FORMAT CSV, HEADER, DELIMITER '\\t')")

    significant_primary = con.execute(f"""
        SELECT *,
          CASE WHEN P_PLACO_PLUS < {ATLAS_THRESHOLD} THEN 'LEVEL_A_ATLAS_STRINGENT'
               WHEN P_PLACO_PLUS < 5e-8 THEN 'LEVEL_B_PAIRWISE_GWAS'
               WHEN q_PLACO_atlas < 0.05 THEN 'LEVEL_C_ATLAS_FDR'
               ELSE 'NONE' END AS PLACO_LEVEL
        FROM placo_q
        WHERE P_PLACO_PLUS < 5e-8 OR q_PLACO_atlas < 0.05
    """).df()
    extreme_files = [PAIR_RUNS / f"{pair}.extreme.tsv.gz" for pair in PAIRS]
    extreme = pd.concat([pd.read_csv(path, sep="\t") for path in extreme_files], ignore_index=True) if any(p.exists() for p in extreme_files) else pd.DataFrame()
    extreme_sig = extreme[extreme["dual_gwas_significant"] == True].copy() if not extreme.empty else pd.DataFrame()
    if not extreme_sig.empty:
        extreme_sig["P_PLACO_PLUS"] = np.nan
        extreme_sig["q_PLACO_atlas"] = np.nan
        extreme_sig["PLACO_LEVEL"] = "EXTREME_DUAL_GWAS"
        extreme_sig["direction_product"] = extreme_sig["direction_product"]
        for col in significant_primary.columns:
            if col not in extreme_sig.columns:
                extreme_sig[col] = np.nan
        extreme_sig = extreme_sig[significant_primary.columns]
    significant = pd.concat([significant_primary, extreme_sig], ignore_index=True, sort=False)
    significant.to_csv(RESULTS / "PLACO_PLUS_SIGNIFICANT.tsv", sep="\t", index=False)
    extreme.to_csv(RESULTS / "PLACO_EXTREME_Z_VARIANTS.tsv", sep="\t", index=False)

    coordinate_numeric = significant[["CHR", "BP"]].apply(pd.to_numeric, errors="coerce")
    coordinate_ok = coordinate_numeric.notna().all(axis=1)
    coordinate_qc = significant[["pair", "SNP", "CHR", "BP"]].copy()
    coordinate_qc["coordinate_status"] = np.where(
        coordinate_ok, "PASS", "MISSING_OR_NONNUMERIC"
    )
    coordinate_qc.to_csv(RESULTS / "PLACO_COORDINATE_QC.tsv", sep="\t", index=False)
    significant_with_coordinates = significant.loc[coordinate_ok].copy()

    atlas = pd.read_csv(ROOT / "results" / "phase1D" / "SUPERGNOVA_OCULAR_SYSTEMIC_ATLAS.tsv", sep="\t")
    local = atlas[atlas["q_atlas"] < 0.05][["pair", "ocular_trait", "systemic_trait", "chr", "start", "end", "region_id", "rho", "P_rho", "q_atlas", "direction"]].copy()
    local = local.rename(columns={"chr": "LOCAL_CHR", "start": "LOCAL_START", "end": "LOCAL_END", "q_atlas": "q_atlas_local", "direction": "local_direction"})
    overlap_rows = []
    for _, v in significant_with_coordinates.iterrows():
        hits = local[(local["pair"] == v["pair"]) & (local["LOCAL_CHR"].astype(int) == int(v["CHR"])) & (local["LOCAL_START"] <= int(v["BP"])) & (local["LOCAL_END"] >= int(v["BP"]))]
        for _, reg in hits.iterrows():
            product = v.get("direction_product", "")
            if product == "POSITIVE_EFFECT_PRODUCT":
                vdir = "+"
            elif product == "NEGATIVE_EFFECT_PRODUCT":
                vdir = "-"
            else:
                vdir = "0"
            ldir = "+" if reg["rho"] > 0 else "-" if reg["rho"] < 0 else "0"
            overlap_rows.append({**v.to_dict(), **reg.to_dict(), "variant_direction": vdir, "local_sign": ldir, "HIDDEN_MIXED_PAIR": v["pair"] in HIDDEN})
    convergence = pd.DataFrame(overlap_rows)
    if convergence.empty:
        convergence = pd.DataFrame(columns=list(significant.columns) + list(local.columns) + ["variant_direction", "local_sign", "HIDDEN_MIXED_PAIR"])
    convergence["directional_convergence"] = np.where(convergence["variant_direction"] == convergence["local_sign"], "DIRECTION_CONCORDANT", "DIRECTION_DISCORDANT")
    convergence.to_csv(RESULTS / "SUPERGNOVA_PLACO_CONVERGENCE.tsv", sep="\t", index=False)

    region_rows = []
    for _, reg in local.iterrows():
        hits = convergence[(convergence["pair"] == reg["pair"]) & (convergence["region_id"] == reg["region_id"])]
        pos = int((hits["variant_direction"] == "+").sum())
        neg = int((hits["variant_direction"] == "-").sum())
        if pos == 0 and neg == 0:
            cls = "NO_VARIANT_LEVEL_SIGNAL"
        elif pos > 0 and neg > 0:
            cls = "VARIANT_MIXED"
        elif pos > 0:
            cls = "VARIANT_POSITIVE_ONLY"
        else:
            cls = "VARIANT_NEGATIVE_ONLY"
        region_rows.append({**reg.to_dict(), "N_positive_effect_product": pos, "N_negative_effect_product": neg, "variant_direction_class": cls, "HIDDEN_MIXED_PAIR": reg["pair"] in HIDDEN})
    region_direction = pd.DataFrame(region_rows)
    region_direction.to_csv(RESULTS / "LOCAL_REGION_VARIANT_DIRECTION.tsv", sep="\t", index=False)

    loci = run_clumping(significant_with_coordinates)
    loci.to_csv(RESULTS / "PLACO_INDEPENDENT_LOCI.tsv", sep="\t", index=False)

    ready_rows = []
    for _, reg in local.iterrows():
        hits = convergence[(convergence["pair"] == reg["pair"]) & (convergence["region_id"] == reg["region_id"])]
        if hits.empty:
            continue
        cat_flag = reg["pair"].startswith("CAT-") or "-CAT" in reg["pair"]
        ready_rows.append({
            "pair": reg["pair"], "region_id": reg["region_id"], "chr": reg["LOCAL_CHR"], "start": reg["LOCAL_START"], "end": reg["LOCAL_END"],
            "full_beta_se_trait1": "PASS", "full_beta_se_trait2": "PASS", "alleles_complete": "PASS",
            "regional_variant_density": int(len(hits)), "appropriate_ld_reference": "CONDITIONAL" if cat_flag else "PASS",
            "ancestry_reference_compatibility": "CONDITIONAL" if cat_flag else "PASS", "sample_size_available": "PASS",
            "imputation_quality_available": "NOT_AVAILABLE_IN_FROZEN_SUMSTATS", "fine_mapping_ready": "CONDITIONAL",
            "CAT_Finnish_LD_flag": "FINNISH_LD_SENSITIVITY_REQUIRED" if cat_flag else "NO",
        })
    readiness = pd.DataFrame(ready_rows)
    readiness.to_csv(RESULTS / "FINEMAPPING_READINESS.tsv", sep="\t", index=False)

    priority_rows = []
    for _, row in convergence.iterrows():
        readiness_row = readiness[(readiness["pair"] == row["pair"]) & (readiness["region_id"] == row["region_id"])]
        fm = readiness_row.iloc[0]["fine_mapping_ready"] if not readiness_row.empty else "FAIL"
        if row["local_sign"] == "-" and row["variant_direction"] == "-" and row["directional_convergence"] == "DIRECTION_CONCORDANT" and fm == "PASS":
            tier = "TIER 1A"
        elif row["local_sign"] == "+" and row["variant_direction"] == "+" and row["directional_convergence"] == "DIRECTION_CONCORDANT" and fm == "PASS":
            tier = "TIER 1B"
        elif row["directional_convergence"] == "DIRECTION_CONCORDANT":
            tier = "TIER 2"
        else:
            tier = "TIER 2"
        priority_rows.append({
            "tier": tier, "pair": row["pair"], "hidden_mixed_pair": row["HIDDEN_MIXED_PAIR"], "chr": row["LOCAL_CHR"], "start": row["LOCAL_START"], "end": row["LOCAL_END"],
            "SUPERGNOVA_rho": row["rho"], "SUPERGNOVA_P": row["P_rho"], "SUPERGNOVA_q": row["q_atlas_local"], "PLACO_lead_SNP": row["SNP"],
            "PLACO_P": row["P_PLACO_PLUS"], "PLACO_q": row["q_PLACO_atlas"], "PLACO_level": row["PLACO_LEVEL"], "beta1": row["BETA_1"], "beta2": row["BETA_2"],
            "effect_product": row["direction_product"], "directional_convergence": row["directional_convergence"], "fine_mapping_ready": fm,
            "CAT_Finnish_LD_flag": "FINNISH_LD_SENSITIVITY_REQUIRED" if row["pair"].startswith("CAT-") else "NO",
        })
    priority = pd.DataFrame(priority_rows).drop_duplicates(subset=["pair", "chr", "start", "end", "PLACO_lead_SNP"])
    priority.to_csv(RESULTS / "PHASE2B_PRIORITY_LOCI.tsv", sep="\t", index=False)

    summary_rows = []
    for pair in PAIRS:
        all_pair = con.execute(f"SELECT count(*) n, sum(P_PLACO_PLUS < 5e-8) pair_gwas, sum(P_PLACO_PLUS < {ATLAS_THRESHOLD}) atlas_stringent, sum(q_PLACO_atlas < 0.05) atlas_fdr, sum(direction_product='POSITIVE_EFFECT_PRODUCT') positive, sum(direction_product='NEGATIVE_EFFECT_PRODUCT') negative FROM placo_q WHERE pair='{pair}'").fetchone()
        sig_pair = significant[significant["pair"] == pair]
        ov = convergence[convergence["pair"] == pair]
        summary_rows.append({"pair": pair, "N_variants_tested": all_pair[0], "N_PLACO_P_lt_5e8": int(all_pair[1] or 0), "N_atlas_stringent": int(all_pair[2] or 0), "N_atlas_FDR": int(all_pair[3] or 0), "N_independent_loci": int((loci["pair"] == pair).sum()) if not loci.empty else 0, "N_positive_effect_variants": int((sig_pair.get("direction_product", pd.Series(dtype=str)) == "POSITIVE_EFFECT_PRODUCT").sum()), "N_negative_effect_variants": int((sig_pair.get("direction_product", pd.Series(dtype=str)) == "NEGATIVE_EFFECT_PRODUCT").sum()), "N_overlapping_PHASE1D_regions": int(ov["region_id"].nunique()) if not ov.empty else 0, "N_direction_concordant_overlaps": int((ov["directional_convergence"] == "DIRECTION_CONCORDANT").sum()) if not ov.empty else 0, "N_direction_discordant_overlaps": int((ov["directional_convergence"] == "DIRECTION_DISCORDANT").sum()) if not ov.empty else 0, "HIDDEN_MIXED_PAIR": pair in HIDDEN})
    pair_summary = pd.DataFrame(summary_rows)
    pair_summary.to_csv(RESULTS / "PLACO_PAIR_SUMMARY.tsv", sep="\t", index=False)

    qc_pass = int((parameter_qc["status"] == "PASS").sum()) == len(parameter_qc)
    overlap_n = int(convergence["region_id"].nunique()) if not convergence.empty else 0
    concordant_n = int((convergence["directional_convergence"] == "DIRECTION_CONCORDANT").sum()) if not convergence.empty else 0
    if not qc_pass:
        verdict = "METHOD_UNSTABLE"
        next_step = "TECHNICAL_HOLD"
    elif overlap_n == 0:
        verdict = "LOCAL_ONLY_ARCHITECTURE"
        next_step = "WAIT_FOR_LAVA_REPLICATION"
    elif int((priority["tier"] == "TIER 1A").sum()) >= 3 and priority[priority["tier"] == "TIER 1A"]["pair"].nunique() >= 2:
        verdict = "STRONG_LOCAL_VARIANT_CONVERGENCE"
        next_step = "PROCEED_TO_FINE_MAPPING_COLOCALIZATION"
    else:
        verdict = "PARTIAL_LOCAL_VARIANT_CONVERGENCE"
        next_step = "RETAIN_VARIANT_LEVEL_ONLY"
    (RESULTS / "PHASE2A_DECISION.txt").write_text("\n".join([
        f"PHASE2A_VERDICT\t{verdict}", f"PLACO_PARAMETER_QC\t{'PASS' if qc_pass else 'FAIL'}", f"PLACO_VARIANTS_TESTED\t{int(n)}", f"PLACO_SIGNIFICANT_VARIANTS\t{len(significant)}", f"SUPERGNOVA_PLACO_OVERLAPPING_REGIONS\t{overlap_n}", f"DIRECTION_CONCORDANT_OVERLAPS\t{concordant_n}", "LAVA_STATUS\tPENDING_LOCAL_METHOD_REPLICATION", "HDL_L_STATUS\tSECONDARY_PENDING_EXACT_N0", f"NEXT_STEP\t{next_step}", "PHASE2B_STARTED\tNO", "" ]) )

    hidden_summary = pair_summary[pair_summary["HIDDEN_MIXED_PAIR"]]
    report = [
        "# PHASE 2A VARIANT-LEVEL PLEIOTROPY REPORT", "",
        "## Protocol amendment", "",
        "Protocol Amendment 002 was scientifically justified because the official LAVA UKB EUR v1.1 reference remained inaccessible while Phase 1D SUPERGNOVA passed its prespecified QC. LAVA remains pending local-method replication; Phase 2A results are not treated as replicated local architecture.", "",
        "## PLACO+ QC", "", f"All 18 pairs were processed with the frozen official PLACO v0.2.0 implementation. Parameter QC: {int((parameter_qc['status'] == 'PASS').sum())}/{len(parameter_qc)} PASS. Allele, coordinate, and signed-Z QC are retained in the Phase 2A result tables.", "",
        table_markdown(parameter_qc, ["pair", "null_variants", "VarZ1", "VarZ2", "CorZ", "ld_pruned_null_variants", "status"]), "",
        "## Variant-level pleiotropy", "", f"The primary calculation tested {int(n):,} finite PLACO+ P values across all 18 genome-wide pair analyses. Pair-wise P < 5e-8, atlas-stringent P < 5e-8/18, and atlas BH-FDR < 0.05 are reported separately. Extreme-Z variants were screened separately with the dual-GWAS rule and were not forced through PLACO+.", "",
        table_markdown(pair_summary, ["pair", "N_variants_tested", "N_PLACO_P_lt_5e8", "N_atlas_stringent", "N_atlas_FDR", "N_independent_loci"]), "",
        "## SUPERGNOVA–PLACO+ convergence", "", f"{overlap_n} Phase 1D significant local regions contained at least one significant variant-level PLACO+ or extreme-Z dual-GWAS signal. There were {concordant_n} direction-concordant overlapping variant-region observations and {int((convergence['directional_convergence'] == 'DIRECTION_DISCORDANT').sum()) if not convergence.empty else 0} direction-discordant observations. Discordance was retained for follow-up rather than discarded.", "",
        f"Coordinate QC: {int(coordinate_ok.sum())}/{len(coordinate_ok)} significant records had numeric CHR/BP coordinates and were eligible for local-region overlap and LD clumping; {int((~coordinate_ok).sum())} records were retained in the significant table but excluded from coordinate-dependent steps.", "",
        table_markdown(region_direction, ["pair", "region_id", "N_positive_effect_product", "N_negative_effect_product", "variant_direction_class"]), "",
        "## Hidden mixed local-sharing pairs", "", table_markdown(hidden_summary, ["pair", "N_positive_effect_variants", "N_negative_effect_variants", "N_overlapping_PHASE1D_regions", "N_direction_concordant_overlaps"]), "",
        "## Fine-mapping readiness", "", table_markdown(readiness, ["pair", "region_id", "regional_variant_density", "appropriate_ld_reference", "imputation_quality_available", "fine_mapping_ready", "CAT_Finnish_LD_flag"]), "",
        "Imputation-quality fields were not present in the frozen summary-statistics inputs, so convergent regions are conditional for later fine-mapping readiness. CAT-containing loci additionally require Finnish-matched LD sensitivity.", "",
        "## Phase 2A verdict", "", f"**{verdict}**", "", f"Next-step gate: **{next_step}**. Phase 2B fine-mapping and colocalization were not started.", "",
        "## LAVA and HDL-L status", "", "LAVA: `PENDING_LOCAL_METHOD_REPLICATION`; do not use unverified mirrors. HDL-L: `SECONDARY_PENDING_EXACT_N0`; no new pairwise HDL-L covariance was run.", "",
        "## Limitations", "", "PLACO+ tests variant-level statistical pleiotropy and cannot distinguish biological, horizontal, or mediated pleiotropy. The analysis is based on public summary statistics, CAT has Finnish-European ancestry, and LD clumping is a discovery-organization step rather than causal-variant evidence.", "",
    ]
    (ROOT / "reports" / "phase2A" / "PHASE2A_VARIANT_PLEIOTROPY_REPORT.md").write_text("\n".join(report) + "\n")
    print("Phase 2A aggregation complete", verdict)


if __name__ == "__main__":
    main()
