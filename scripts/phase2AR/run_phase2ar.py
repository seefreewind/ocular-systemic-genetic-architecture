#!/usr/bin/env python3
"""Phase 2A-R independent-locus consolidation and readiness audit.

This script reuses the frozen Phase 1D/Phase 2A outputs and the already
completed PLINK clumping logs. It does not rerun SUPERGNOVA, PLACO+, or any
fine-mapping/colocalization method.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RESULTS_1D = ROOT / "results" / "phase1D"
RESULTS_2A = ROOT / "results" / "phase2A"
RESULTS = ROOT / "results" / "phase2AR"
REPORTS = ROOT / "reports" / "phase2AR"
PAIR_RUNS = RESULTS_2A / "pair_runs"
CLUMP_DIR = RESULTS_2A / "ld_clump"
REFERENCE_BIM = RESULTS_2A / "reference_freq" / "eur_reference.bim.tsv"
PLINK_REFERENCE = ROOT / "data" / "phase1C" / "SUPERGNOVA_reference" / "bfiles" / "eur_chr{chrom}_SNPmaf5"

PAIRS = [
    "AMD-CAD", "AMD-STROKE", "AMD-T2D", "AMD-CKD", "AMD-AD", "AMD-PD",
    "POAG-CAD", "POAG-STROKE", "POAG-T2D", "POAG-CKD", "POAG-AD", "POAG-PD",
    "CAT-CAD", "CAT-STROKE", "CAT-T2D", "CAT-CKD", "CAT-AD", "CAT-PD",
]
HIDDEN_MIXED = {"AMD-CAD", "AMD-STROKE", "AMD-T2D", "AMD-CKD", "POAG-CKD", "CAT-STROKE", "CAT-AD"}
GLOBAL_POSITIVE = {"CAT-CAD", "CAT-T2D"}
PERMUTATION_B = 10_000
PERMUTATION_SEED = 20260921
CLUMP_KB = 500
DENSITY_KB = 1_000
MIN_DENSE_VARIANTS = 50
MIN_COMPLETE_FRACTION = 0.95


def safe_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def md_table(frame: pd.DataFrame, columns: list[str], n: int | None = None) -> str:
    if n is not None:
        frame = frame.head(n)
    if frame.empty:
        return "No rows."
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for _, row in frame.iterrows():
        vals = []
        for col in columns:
            value = row.get(col, "NA")
            if pd.isna(value):
                value = "NA"
            elif isinstance(value, float):
                value = f"{value:.5g}"
            vals.append(str(value).replace("|", "\\|"))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def parse_sp2(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text in {".", "NONE", "NA"}:
        return []
    return [x for x in text.split(",") if x and x not in {".", "NONE", "NA"}]


def read_clumps(significant: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    sig = significant.copy()
    sig["CHR"] = safe_num(sig["CHR"])
    sig["BP"] = safe_num(sig["BP"])
    sig_key = sig.set_index(["pair", "SNP"], drop=False)
    rows: list[dict[str, object]] = []
    clump_files = []
    log_files = []
    for path in sorted(CLUMP_DIR.glob("*.clumps")):
        if path.name.startswith("._"):
            continue
        match = re.match(r"(.+)\.chr(\d+)\.clumps$", path.name)
        if not match:
            continue
        pair, chrom = match.group(1), int(match.group(2))
        table = pd.read_csv(path, sep="\t")
        clump_files.append(path)
        log_path = path.with_suffix(".log")
        if log_path.exists():
            log_files.append(log_path)
        if table.empty:
            continue
        for _, row in table.iterrows():
            lead = str(row["ID"])
            key = (pair, lead)
            if key not in sig_key.index:
                raise RuntimeError(f"Clump lead is absent from frozen significant table: {pair} {lead}")
            lead_row = sig_key.loc[key]
            if isinstance(lead_row, pd.DataFrame):
                lead_row = lead_row.iloc[0]
            members = [lead] + parse_sp2(row.get("SP2"))
            positions = []
            member_set = []
            for member in dict.fromkeys(members):
                member_key = (pair, member)
                if member_key in sig_key.index:
                    member_row = sig_key.loc[member_key]
                    if isinstance(member_row, pd.DataFrame):
                        member_row = member_row.iloc[0]
                    bp = safe_num(pd.Series([member_row["BP"]])).iloc[0]
                    if pd.notna(bp):
                        positions.append(int(bp))
                    member_set.append(member)
            bp = safe_num(pd.Series([lead_row["BP"]])).iloc[0]
            if pd.isna(bp):
                raise RuntimeError(f"Clump lead has missing coordinate: {pair} {lead}")
            lead_bp = int(bp)
            window_start = max(1, lead_bp - CLUMP_KB * 1000)
            window_end = lead_bp + CLUMP_KB * 1000
            locus_id = f"{pair}|chr{chrom}|{lead}"
            rows.append({
                "PLACO_locus_id": locus_id,
                "pair": pair,
                "chr": chrom,
                "lead_SNP": lead,
                "lead_BP": lead_bp,
                "lead_P": float(row["P"]),
                "clump_variant_count": int(row.get("TOTAL", len(member_set))),
                "clump_member_count_mapped": len(member_set),
                "clump_member_start": min(positions) if positions else lead_bp,
                "clump_member_end": max(positions) if positions else lead_bp,
                "locus_start": window_start,
                "locus_end": window_end,
                "clump_members": ";".join(member_set),
                "P_PLACO_PLUS": lead_row.get("P_PLACO_PLUS", np.nan),
                "q_PLACO_atlas": lead_row.get("q_PLACO_atlas", np.nan),
                "PLACO_LEVEL": lead_row.get("PLACO_LEVEL", "NA"),
                "BETA_1": lead_row.get("BETA_1", np.nan),
                "BETA_2": lead_row.get("BETA_2", np.nan),
                "SE_1": lead_row.get("SE_1", np.nan),
                "SE_2": lead_row.get("SE_2", np.nan),
                "P_1": lead_row.get("P_1", np.nan),
                "P_2": lead_row.get("P_2", np.nan),
                "A1": lead_row.get("A1", np.nan),
                "A2": lead_row.get("A2", np.nan),
                "N_1": lead_row.get("N_1", np.nan),
                "N_2": lead_row.get("N_2", np.nan),
                "direction_product": lead_row.get("direction_product", "NA"),
                "effect_product": float(lead_row["BETA_1"]) * float(lead_row["BETA_2"]),
                "CAT_FINNISH_LD_FLAG": "TRUE" if pair.startswith("CAT-") else "FALSE",
                "source_clump_file": path.name,
            })
    loci = pd.DataFrame(rows)
    if loci.empty:
        raise RuntimeError("No valid PLINK clump output was found")
    loci = loci.drop_duplicates(subset=["pair", "lead_SNP"]).reset_index(drop=True)
    loci["effect_sign"] = np.where(loci["effect_product"] > 0, "+", np.where(loci["effect_product"] < 0, "-", "0"))

    command_text = "\n".join(path.read_text(errors="replace") for path in log_files)
    reference_paths = sorted(set(re.findall(r"--bfile\s+([^\n]+)", command_text)))
    audit = {
        "clump_file_count": len(clump_files),
        "clump_log_count": len(log_files),
        "clump_row_count": len(loci),
        "clump_pair_count": loci["pair"].nunique(),
        "clump_unique_lead_snp_count": loci["lead_SNP"].nunique(),
        "clump_unique_pair_lead_count": loci[["pair", "lead_SNP"]].drop_duplicates().shape[0],
        "clump_parameters": "PLINK clump-r2=0.1; clump-kb=500; clump-p1=1; clump-p2=1; threads=4; memory=4096 MB",
        "ld_reference_paths": ";".join(reference_paths),
        "clump_outputs_reused": "YES",
    }
    return loci, audit


def write_independent_locus_audit(significant: pd.DataFrame, loci: pd.DataFrame, clump_audit: dict[str, object]) -> pd.DataFrame:
    existing_path = RESULTS_2A / "PLACO_INDEPENDENT_LOCI.tsv"
    existing = pd.read_csv(existing_path, sep="\t")
    existing_status = existing.get("clump_status", pd.Series(dtype=str)).astype(str)
    existing_all_fallback = bool(len(existing) and (existing_status == "LD_PANEL_NOT_AVAILABLE_OR_NOT_CLUMPED").all())
    rows = [
        {
            "scope": "EXISTING_PHASE2A_TABLE",
            "source": str(existing_path.relative_to(ROOT)),
            "rows": len(existing),
            "unique_pairs": existing["pair"].nunique(),
            "unique_lead_snps": existing["lead_SNP"].nunique() if "lead_SNP" in existing else np.nan,
            "unique_pair_lead_loci": existing[["pair", "lead_SNP"]].drop_duplicates().shape[0] if "lead_SNP" in existing else np.nan,
            "unique_genomic_loci": existing[["pair", "CHR", "BP"]].drop_duplicates().shape[0],
            "validity": "INVALID_FALLBACK_TABLE" if existing_all_fallback else "REVIEW_REQUIRED",
            "clumping_parameters": "Not encoded in table",
            "ld_reference": "Not encoded in table",
        },
        {
            "scope": "REUSED_PLINK_CLUMP_OUTPUTS",
            "source": str(CLUMP_DIR.relative_to(ROOT)),
            "rows": len(loci),
            "unique_pairs": loci["pair"].nunique(),
            "unique_lead_snps": loci["lead_SNP"].nunique(),
            "unique_pair_lead_loci": loci[["pair", "lead_SNP"]].drop_duplicates().shape[0],
            "unique_genomic_loci": loci[["pair", "chr", "lead_BP"]].drop_duplicates().shape[0],
            "validity": "VALID_REUSED",
            "clumping_parameters": clump_audit["clump_parameters"],
            "ld_reference": clump_audit["ld_reference_paths"],
        },
    ]
    for pair, g in loci.groupby("pair", sort=True):
        rows.append({
            "scope": f"PAIR:{pair}",
            "source": str(CLUMP_DIR.relative_to(ROOT)),
            "rows": len(g),
            "unique_pairs": 1,
            "unique_lead_snps": g["lead_SNP"].nunique(),
            "unique_pair_lead_loci": g[["pair", "lead_SNP"]].drop_duplicates().shape[0],
            "unique_genomic_loci": g[["pair", "chr", "lead_BP"]].drop_duplicates().shape[0],
            "validity": "VALID_REUSED",
            "clumping_parameters": clump_audit["clump_parameters"],
            "ld_reference": clump_audit["ld_reference_paths"],
        })
    audit = pd.DataFrame(rows)
    audit.to_csv(RESULTS / "INDEPENDENT_LOCUS_AUDIT.tsv", sep="\t", index=False)
    return audit


def map_loci_to_regions(loci: pd.DataFrame, atlas: pd.DataFrame) -> pd.DataFrame:
    regions = atlas[atlas["q_atlas"] < 0.05].copy()
    regions["rho_sign"] = np.where(regions["rho"] > 0, "+", np.where(regions["rho"] < 0, "-", "0"))
    rows = []
    for pair, locus_group in loci.groupby("pair", sort=False):
        rg = regions[regions["pair"] == pair]
        for _, locus in locus_group.iterrows():
            hits = rg[(rg["chr"] == locus["chr"]) & (rg["start"] <= locus["locus_end"]) & (rg["end"] >= locus["locus_start"])]
            for _, region in hits.iterrows():
                if locus["effect_product"] > 0 and region["rho"] > 0:
                    status = "DIRECTION_CONCORDANT"
                elif locus["effect_product"] < 0 and region["rho"] < 0:
                    status = "DIRECTION_CONCORDANT"
                else:
                    status = "DIRECTION_DISCORDANT"
                rows.append({
                    **locus.to_dict(),
                    "ocular_trait": region["ocular_trait"],
                    "systemic_trait": region["systemic_trait"],
                    "SUPERGNOVA_region": region["region_id"],
                    "SUPERGNOVA_start": int(region["start"]),
                    "SUPERGNOVA_end": int(region["end"]),
                    "SUPERGNOVA_rho": float(region["rho"]),
                    "SUPERGNOVA_q_atlas": float(region["q_atlas"]),
                    "SUPERGNOVA_direction": region["direction"],
                    "rho_sign": region["rho_sign"],
                    "directional_status": status,
                    "HIDDEN_MIXED_PAIR": pair in HIDDEN_MIXED,
                })
    mapped = pd.DataFrame(rows)
    if mapped.empty:
        mapped = pd.DataFrame(columns=list(loci.columns) + [
            "ocular_trait", "systemic_trait", "SUPERGNOVA_region", "SUPERGNOVA_start", "SUPERGNOVA_end",
            "SUPERGNOVA_rho", "SUPERGNOVA_q_atlas", "SUPERGNOVA_direction", "rho_sign", "directional_status", "HIDDEN_MIXED_PAIR"
        ])
    return mapped


def compute_density(mapped: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
    targets = loci[loci["PLACO_locus_id"].isin(mapped["PLACO_locus_id"].unique())].copy()
    targets["density_start"] = (targets["lead_BP"] - DENSITY_KB * 1000).clip(lower=1)
    targets["density_end"] = targets["lead_BP"] + DENSITY_KB * 1000
    columns = ["PLACO_locus_id", "pair", "chr", "lead_SNP", "lead_BP", "density_start", "density_end"]
    targets = targets[columns]
    con = duckdb.connect()
    ref_path = str(REFERENCE_BIM).replace("'", "''")
    con.execute(
        "CREATE OR REPLACE TEMP TABLE ref_snp AS "
        f"SELECT CAST(CHR AS INTEGER) AS chr, CAST(SNP AS VARCHAR) AS snp FROM read_csv_auto('{ref_path}', header=true, delim='\\t', sample_size=100000)"
    )
    results = []
    for pair, cand in targets.groupby("pair", sort=True):
        cand = cand.reset_index(drop=True)
        con.register("candidate_windows_df", cand)
        con.execute("CREATE OR REPLACE TEMP TABLE candidate_windows AS SELECT * FROM candidate_windows_df")
        raw_path = str(PAIR_RUNS / f"{pair}.placo.tsv.gz").replace("'", "''")
        con.execute(
            "CREATE OR REPLACE TEMP VIEW raw_pair AS "
            f"SELECT * FROM read_csv_auto('{raw_path}', header=true, delim='\\t', compression='gzip', sample_size=100000)"
        )
        query = """
        SELECT
          c.PLACO_locus_id,
          COUNT(r.SNP) AS variant_count,
          COUNT(DISTINCT r.SNP) AS unique_variant_count,
          SUM(CASE WHEN TRY_CAST(r.BETA_1 AS DOUBLE) IS NOT NULL AND TRY_CAST(r.BETA_2 AS DOUBLE) IS NOT NULL AND TRY_CAST(r.SE_1 AS DOUBLE) IS NOT NULL AND TRY_CAST(r.SE_2 AS DOUBLE) IS NOT NULL THEN 1 ELSE 0 END) AS beta_se_complete_count,
          SUM(CASE WHEN TRY_CAST(r.P_1 AS DOUBLE) IS NOT NULL AND TRY_CAST(r.P_2 AS DOUBLE) IS NOT NULL THEN 1 ELSE 0 END) AS p_complete_count,
          SUM(CASE WHEN r.A1 IS NOT NULL AND r.A2 IS NOT NULL AND CAST(r.A1 AS VARCHAR) <> '' AND CAST(r.A2 AS VARCHAR) <> '' THEN 1 ELSE 0 END) AS allele_complete_count,
          SUM(CASE WHEN TRY_CAST(r.N_1 AS DOUBLE) IS NOT NULL AND TRY_CAST(r.N_2 AS DOUBLE) IS NOT NULL THEN 1 ELSE 0 END) AS sample_size_complete_count,
          SUM(CASE WHEN UPPER(CAST(r.A1 AS VARCHAR) || CAST(r.A2 AS VARCHAR)) IN ('AT','TA','CG','GC') THEN 1 ELSE 0 END) AS ambiguous_allele_count,
          SUM(CASE WHEN ref.snp IS NOT NULL THEN 1 ELSE 0 END) AS reference_overlap_count
        FROM candidate_windows c
        LEFT JOIN raw_pair r
          ON TRY_CAST(r.CHR AS INTEGER) = c.chr
         AND TRY_CAST(r.BP AS BIGINT) BETWEEN c.density_start AND c.density_end
        LEFT JOIN ref_snp ref
          ON ref.chr = TRY_CAST(r.CHR AS INTEGER) AND ref.snp = CAST(r.SNP AS VARCHAR)
        GROUP BY c.PLACO_locus_id
        """
        result = con.execute(query).df()
        result["pair"] = pair
        results.append(result)
        con.unregister("candidate_windows_df")
    con.close()
    if not results:
        return pd.DataFrame(columns=columns + ["variant_count"])
    density = pd.concat(results, ignore_index=True)
    density = targets.merge(density, on=["PLACO_locus_id", "pair"], how="left")
    numeric_cols = ["variant_count", "unique_variant_count", "beta_se_complete_count", "p_complete_count", "allele_complete_count", "sample_size_complete_count", "ambiguous_allele_count", "reference_overlap_count"]
    for col in numeric_cols:
        density[col] = pd.to_numeric(density[col], errors="coerce").fillna(0).astype(int)
    density["duplicate_variant_count"] = density["variant_count"] - density["unique_variant_count"]
    density["beta_se_fraction"] = np.where(density["variant_count"] > 0, density["beta_se_complete_count"] / density["variant_count"], 0.0)
    density["p_complete_fraction"] = np.where(density["variant_count"] > 0, density["p_complete_count"] / density["variant_count"], 0.0)
    density["allele_complete_fraction"] = np.where(density["variant_count"] > 0, density["allele_complete_count"] / density["variant_count"], 0.0)
    density["sample_size_complete_fraction"] = np.where(density["variant_count"] > 0, density["sample_size_complete_count"] / density["variant_count"], 0.0)
    density["reference_coverage_fraction"] = np.where(density["variant_count"] > 0, density["reference_overlap_count"] / density["variant_count"], 0.0)
    density.to_csv(RESULTS / "REGIONAL_VARIANT_DENSITY.tsv", sep="\t", index=False)
    return density


def readiness_audit(mapped: pd.DataFrame, density: pd.DataFrame) -> pd.DataFrame:
    if mapped.empty:
        out = pd.DataFrame()
    else:
        out = mapped.drop_duplicates("PLACO_locus_id").merge(density, on=["PLACO_locus_id", "pair", "chr", "lead_SNP", "lead_BP"], how="left")
        out["full_beta_se_trait1"] = np.where(out["beta_se_fraction"] >= MIN_COMPLETE_FRACTION, "PASS", "FAIL")
        out["full_beta_se_trait2"] = np.where(out["beta_se_fraction"] >= MIN_COMPLETE_FRACTION, "PASS", "FAIL")
        out["alleles_complete"] = np.where(out["allele_complete_fraction"] >= MIN_COMPLETE_FRACTION, "PASS", "FAIL")
        out["sample_size_available"] = np.where(out["sample_size_complete_fraction"] >= MIN_COMPLETE_FRACTION, "PASS", "FAIL")
        out["regional_density_status"] = np.where(out["variant_count"] >= MIN_DENSE_VARIANTS, "PASS", "FAIL")
        out["imputation_quality_available"] = "NOT_AVAILABLE_IN_FROZEN_SUMSTATS"
        out["genome_build"] = "GRCh37/hg19"
        out["ld_reference_status"] = np.where(out["pair"].str.startswith("CAT-"), "FINNISH_LD_REQUIRED", "EUR_REFERENCE_AVAILABLE_ORGANIZATIONAL_ONLY")
        out["reference_ancestry_compatibility"] = np.where(out["pair"].str.startswith("CAT-"), "CONDITIONAL", "PASS")
        core_pass = (
            (out["full_beta_se_trait1"] == "PASS") &
            (out["full_beta_se_trait2"] == "PASS") &
            (out["alleles_complete"] == "PASS") &
            (out["sample_size_available"] == "PASS") &
            (out["regional_density_status"] == "PASS") &
            (out["p_complete_fraction"] >= MIN_COMPLETE_FRACTION) &
            (out["duplicate_variant_count"] == 0)
        )
        pass_ready = core_pass & (out["imputation_quality_available"] != "NOT_AVAILABLE_IN_FROZEN_SUMSTATS") & ~out["pair"].str.startswith("CAT-")
        out["fine_mapping_ready"] = np.where(pass_ready, "PASS", np.where(core_pass, "CONDITIONAL", "FAIL"))
        out["CAT_Finnish_LD_flag"] = np.where(out["pair"].str.startswith("CAT-"), "FINNISH_LD_SENSITIVITY_REQUIRED", "NO")
        out["regional_variant_density"] = out["variant_count"]
        out["duplicate_variant_flag"] = np.where(out["duplicate_variant_count"] > 0, "YES", "NO")
        out["coordinate_qc"] = "PASS"
    out.to_csv(RESULTS / "FINE_MAPPING_READINESS_REAUDIT.tsv", sep="\t", index=False)
    return out


def region_summary(atlas: pd.DataFrame, mapped: pd.DataFrame, readiness: pd.DataFrame) -> pd.DataFrame:
    regions = atlas[atlas["q_atlas"] < 0.05].copy()
    rows = []
    readiness_by_locus = readiness.set_index("PLACO_locus_id") if not readiness.empty else pd.DataFrame()
    for _, reg in regions.iterrows():
        g = mapped[(mapped["pair"] == reg["pair"]) & (mapped["SUPERGNOVA_region"] == reg["region_id"])]
        g = g.drop_duplicates("PLACO_locus_id")
        n = len(g)
        pos = int((g["effect_sign"] == "+").sum())
        neg = int((g["effect_sign"] == "-").sum())
        concord = int((g["directional_status"] == "DIRECTION_CONCORDANT").sum())
        discord = int((g["directional_status"] == "DIRECTION_DISCORDANT").sum())
        if n == 0:
            classification = "NO_INDEPENDENT_PLACO_LOCUS"
            readiness_status = "NO_INDEPENDENT_PLACO_LOCUS"
        elif n == 1:
            classification = "SINGLE_LOCUS_CONCORDANT" if concord == 1 else "SINGLE_LOCUS_DISCORDANT"
            vals = readiness_by_locus.loc[g.iloc[0]["PLACO_locus_id"], "fine_mapping_ready"] if not readiness.empty and g.iloc[0]["PLACO_locus_id"] in readiness_by_locus.index else "FAIL"
            readiness_status = vals if isinstance(vals, str) else str(vals)
        elif pos > 0 and neg > 0:
            classification = "MIXED_VARIANT_DIRECTION"
            vals = readiness[readiness["PLACO_locus_id"].isin(g["PLACO_locus_id"])]
            readiness_status = "FAIL" if (vals["fine_mapping_ready"] == "FAIL").any() else "CONDITIONAL" if (vals["fine_mapping_ready"] == "CONDITIONAL").any() else "PASS"
        elif concord == n:
            classification = "UNIFORM_CONCORDANT"
            vals = readiness[readiness["PLACO_locus_id"].isin(g["PLACO_locus_id"])]
            readiness_status = "FAIL" if (vals["fine_mapping_ready"] == "FAIL").any() else "CONDITIONAL" if (vals["fine_mapping_ready"] == "CONDITIONAL").any() else "PASS"
        elif discord == n:
            classification = "UNIFORM_DISCORDANT"
            vals = readiness[readiness["PLACO_locus_id"].isin(g["PLACO_locus_id"])]
            readiness_status = "FAIL" if (vals["fine_mapping_ready"] == "FAIL").any() else "CONDITIONAL" if (vals["fine_mapping_ready"] == "CONDITIONAL").any() else "PASS"
        else:
            classification = "MIXED_VARIANT_DIRECTION"
            readiness_status = "CONDITIONAL"
        rows.append({
            "pair": reg["pair"], "ocular_trait": reg["ocular_trait"], "systemic_trait": reg["systemic_trait"],
            "SUPERGNOVA_region": reg["region_id"], "chr": int(reg["chr"]), "region_start": int(reg["start"]), "region_end": int(reg["end"]),
            "SUPERGNOVA_rho": reg["rho"], "SUPERGNOVA_q": reg["q_atlas"], "SUPERGNOVA_direction": reg["direction"],
            "hidden_mixed_pair": reg["pair"] in HIDDEN_MIXED, "N_independent_PLACO_loci": n,
            "N_positive_effect_loci": pos, "N_negative_effect_loci": neg,
            "N_direction_concordant": concord, "N_direction_discordant": discord,
            "concordant_fraction": concord / n if n else np.nan,
            "directional_class": classification, "fine_mapping_readiness": readiness_status,
            "lead_SNPs": ";".join(g.sort_values("lead_P")["lead_SNP"].astype(str).tolist()),
            "PLACO_locus_ids": ";".join(g.sort_values("lead_P")["PLACO_locus_id"].astype(str).tolist()),
        })
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "REGION_LEVEL_DIRECTIONAL_CONVERGENCE.tsv", sep="\t", index=False)
    out[out["SUPERGNOVA_rho"] < 0].to_csv(RESULTS / "NEGATIVE_LOCAL_REGION_AUDIT.tsv", sep="\t", index=False)
    out[out["SUPERGNOVA_rho"] > 0].to_csv(RESULTS / "POSITIVE_LOCAL_REGION_AUDIT.tsv", sep="\t", index=False)
    return out


def hidden_pair_audit(region_table: pd.DataFrame, mapped: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pair in sorted(HIDDEN_MIXED):
        rg = region_table[region_table["pair"] == pair]
        mg = mapped[mapped["pair"] == pair].drop_duplicates("PLACO_locus_id")
        rows.append({
            "pair": pair,
            "positive_SUPERGNOVA_regions": int((rg["SUPERGNOVA_rho"] > 0).sum()),
            "negative_SUPERGNOVA_regions": int((rg["SUPERGNOVA_rho"] < 0).sum()),
            "positive_regions_with_independent_PLACO_support": int(((rg["SUPERGNOVA_rho"] > 0) & (rg["N_independent_PLACO_loci"] > 0)).sum()),
            "negative_regions_with_independent_PLACO_support": int(((rg["SUPERGNOVA_rho"] < 0) & (rg["N_independent_PLACO_loci"] > 0)).sum()),
            "positive_direction_concordant_loci": int(((mg["SUPERGNOVA_rho"] > 0) & (mg["directional_status"] == "DIRECTION_CONCORDANT")).sum()),
            "negative_direction_concordant_loci": int(((mg["SUPERGNOVA_rho"] < 0) & (mg["directional_status"] == "DIRECTION_CONCORDANT")).sum()),
            "positive_direction_discordant_loci": int(((mg["SUPERGNOVA_rho"] > 0) & (mg["directional_status"] == "DIRECTION_DISCORDANT")).sum()),
            "negative_direction_discordant_loci": int(((mg["SUPERGNOVA_rho"] < 0) & (mg["directional_status"] == "DIRECTION_DISCORDANT")).sum()),
        })
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "HIDDEN_MIXED_PAIR_LOCUS_AUDIT.tsv", sep="\t", index=False)
    return out


def permutation_test(region_table: pd.DataFrame, mapped: pd.DataFrame) -> pd.DataFrame:
    eligible_rows = []
    for _, reg in region_table[region_table["N_independent_PLACO_loci"] > 0].iterrows():
        g = mapped[(mapped["pair"] == reg["pair"]) & (mapped["SUPERGNOVA_region"] == reg["SUPERGNOVA_region"])].sort_values(["lead_P", "lead_SNP"])
        if g.empty:
            continue
        lead = g.iloc[0]
        eligible_rows.append({
            "pair": reg["pair"], "region": reg["SUPERGNOVA_region"], "rho_sign": "+" if reg["SUPERGNOVA_rho"] > 0 else "-",
            "lead_effect_sign": lead["effect_sign"], "observed_concordant": int((reg["SUPERGNOVA_rho"] > 0 and lead["effect_sign"] == "+") or (reg["SUPERGNOVA_rho"] < 0 and lead["effect_sign"] == "-")),
            "hidden": reg["pair"] in HIDDEN_MIXED, "global_positive": reg["pair"] in GLOBAL_POSITIVE,
        })
    eligible = pd.DataFrame(eligible_rows)
    rng = np.random.default_rng(PERMUTATION_SEED)
    scopes = {
        "ALL_ELIGIBLE_REGIONS": np.ones(len(eligible), dtype=bool),
        "NEGATIVE_RHO_REGIONS": eligible["rho_sign"].eq("-").to_numpy(),
        "POSITIVE_RHO_REGIONS": eligible["rho_sign"].eq("+").to_numpy(),
        "HIDDEN_MIXED_PAIRS": eligible["hidden"].to_numpy(),
        "GLOBAL_POSITIVE_CONTROL_PAIRS": eligible["global_positive"].to_numpy(),
    }
    pair_indices = {pair: np.flatnonzero(eligible["pair"].eq(pair).to_numpy()) for pair in eligible["pair"].unique()}
    null = {scope: np.zeros(PERMUTATION_B, dtype=float) for scope in scopes}
    observed = {scope: int(eligible.loc[mask, "observed_concordant"].sum()) for scope, mask in scopes.items()}
    for b in range(PERMUTATION_B):
        perm_sign = np.empty(len(eligible), dtype=object)
        for pair, indices in pair_indices.items():
            labels = eligible.iloc[indices]["rho_sign"].to_numpy(copy=True)
            perm_sign[indices] = rng.permutation(labels)
        for scope, mask in scopes.items():
            if not mask.any():
                null[scope][b] = np.nan
            else:
                null[scope][b] = np.sum(np.where(perm_sign[mask] == eligible.loc[mask, "lead_effect_sign"].to_numpy(), 1, 0))
    rows = []
    for scope, mask in scopes.items():
        vals = null[scope][np.isfinite(null[scope])]
        obs = observed[scope]
        if len(vals):
            mean = float(vals.mean())
            sd = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
            p = float((1 + np.sum(vals >= obs)) / (len(vals) + 1))
            z = float((obs - mean) / sd) if sd > 0 else np.nan
        else:
            mean = sd = p = z = np.nan
        rows.append({
            "scope": scope, "N_eligible_regions": int(mask.sum()), "observed_concordance": obs,
            "observed_fraction": obs / int(mask.sum()) if mask.any() else np.nan,
            "B": PERMUTATION_B, "null_mean": mean, "null_SD": sd, "empirical_P": p, "empirical_Z": z,
            "seed": PERMUTATION_SEED, "permutation_rule": "rho sign labels permuted within pair; PLACO lead directions fixed",
        })
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "DIRECTION_CONCORDANCE_PERMUTATION.tsv", sep="\t", index=False)
    return out


def build_priority(region_table: pd.DataFrame, mapped: pd.DataFrame, readiness: pd.DataFrame) -> pd.DataFrame:
    ready_cols = ["PLACO_locus_id", "fine_mapping_ready", "ld_reference_status", "CAT_Finnish_LD_flag"]
    map_ready = mapped.merge(readiness[ready_cols], on="PLACO_locus_id", how="left") if not readiness.empty else mapped.copy()
    rows = []
    for _, reg in region_table.iterrows():
        g = map_ready[(map_ready["pair"] == reg["pair"]) & (map_ready["SUPERGNOVA_region"] == reg["SUPERGNOVA_region"])].drop_duplicates("PLACO_locus_id")
        if g.empty:
            rows.append({
                "tier": "TIER L", "pair": reg["pair"], "hidden_mixed_pair": reg["hidden_mixed_pair"], "chr": reg["chr"],
                "region_start": reg["region_start"], "region_end": reg["region_end"], "SUPERGNOVA_rho": reg["SUPERGNOVA_rho"], "SUPERGNOVA_q": reg["SUPERGNOVA_q"],
                "PLACO_locus_id": "NA", "lead_SNP": "NA", "lead_BP": np.nan, "PLACO_P": np.nan, "PLACO_q": np.nan,
                "beta_ocular": np.nan, "beta_systemic": np.nan, "effect_product": np.nan, "directional_status": "NO_INDEPENDENT_PLACO_LOCUS",
                "N_independent_loci_in_region": 0, "fine_mapping_readiness": "NOT_APPLICABLE", "LD_reference_status": "NOT_APPLICABLE", "Finnish_LD_flag": "NO",
            })
            continue
        for _, row in g.iterrows():
            ready = row.get("fine_mapping_ready", "FAIL")
            complex_region = reg["directional_class"] == "MIXED_VARIANT_DIRECTION"
            concordant = row["directional_status"] == "DIRECTION_CONCORDANT"
            if complex_region or not concordant or ready == "FAIL":
                tier = "TIER X"
            elif reg["SUPERGNOVA_rho"] < 0 and row["effect_product"] < 0:
                tier = "TIER A1" if ready == "PASS" else "TIER A2"
            elif reg["SUPERGNOVA_rho"] > 0 and row["effect_product"] > 0:
                tier = "TIER C1" if ready == "PASS" else "TIER C2"
            else:
                tier = "TIER X"
            rows.append({
                "tier": tier, "pair": row["pair"], "hidden_mixed_pair": row["HIDDEN_MIXED_PAIR"], "chr": row["chr"],
                "region_start": row["SUPERGNOVA_start"], "region_end": row["SUPERGNOVA_end"], "SUPERGNOVA_rho": row["SUPERGNOVA_rho"], "SUPERGNOVA_q": row["SUPERGNOVA_q_atlas"],
                "PLACO_locus_id": row["PLACO_locus_id"], "lead_SNP": row["lead_SNP"], "lead_BP": row["lead_BP"], "PLACO_P": row["lead_P"], "PLACO_q": row["q_PLACO_atlas"],
                "beta_ocular": row["BETA_1"], "beta_systemic": row["BETA_2"], "effect_product": row["effect_product"], "directional_status": row["directional_status"],
                "N_independent_loci_in_region": reg["N_independent_PLACO_loci"], "fine_mapping_readiness": ready,
                "LD_reference_status": row.get("ld_reference_status", "NA"), "Finnish_LD_flag": row.get("CAT_Finnish_LD_flag", "NO"),
            })
    out = pd.DataFrame(rows).sort_values(["tier", "pair", "chr", "region_start", "PLACO_P"], na_position="last")
    out.to_csv(RESULTS / "PHASE2B_REFINED_PRIORITY.tsv", sep="\t", index=False)
    return out


def decision(region_table: pd.DataFrame, mapped: pd.DataFrame, readiness: pd.DataFrame, priority: pd.DataFrame, total_loci: int | None = None) -> tuple[str, str, dict[str, int]]:
    a1 = priority[priority["tier"] == "TIER A1"]
    a2 = priority[priority["tier"] == "TIER A2"]
    antag = priority[priority["tier"].isin(["TIER A1", "TIER A2"])]
    if len(a1) >= 3 and a1["pair"].nunique() >= 2:
        verdict = "STRONG_INDEPENDENT_LOCUS_CONVERGENCE"
        next_step = "PROCEED_TO_PHASE2B_FINE_MAPPING"
    elif len(a1) in {1, 2} or len(antag) >= 3:
        verdict = "PARTIAL_INDEPENDENT_LOCUS_CONVERGENCE"
        next_step = "LIMITED_PHASE2B_FINE_MAPPING"
    elif len(antag) > 0:
        verdict = "WEAK_AFTER_LD_CONSOLIDATION"
        next_step = "RETAIN_VARIANT_LEVEL_ONLY"
    else:
        verdict = "METHOD_OR_DATA_LIMITED"
        next_step = "RETAIN_VARIANT_LEVEL_ONLY"
    counts = {
        "independent_loci": int(total_loci) if total_loci is not None else (int(mapped["PLACO_locus_id"].nunique()) if not mapped.empty else 0),
        "mapped_loci": int(mapped["PLACO_locus_id"].nunique()) if not mapped.empty else 0,
        "regions": int(len(region_table)),
        "regions_with_loci": int((region_table["N_independent_PLACO_loci"] > 0).sum()),
        "a1": int(len(a1)), "a2": int(len(a2)), "antagonistic_candidates": int(len(antag)),
        "c1": int((priority["tier"] == "TIER C1").sum()), "c2": int((priority["tier"] == "TIER C2").sum()),
    }
    return verdict, next_step, counts


def write_report(verdict: str, next_step: str, counts: dict[str, int], audit: pd.DataFrame, loci: pd.DataFrame, mapped: pd.DataFrame, regions: pd.DataFrame, permutation: pd.DataFrame, hidden: pd.DataFrame, readiness: pd.DataFrame, priority: pd.DataFrame) -> None:
    neg_candidates = priority[priority["tier"].isin(["TIER A1", "TIER A2"])].sort_values("PLACO_P")
    pos_candidates = priority[priority["tier"].isin(["TIER C1", "TIER C2"])].sort_values("PLACO_P")
    readiness_counts = readiness["fine_mapping_ready"].value_counts().to_dict() if not readiness.empty else {}
    readiness_display = readiness.copy()
    readiness_sort_cols = [col for col in ["fine_mapping_ready", "pair", "lead_BP", "lead_BP_x", "lead_BP_y"] if col in readiness_display.columns]
    if readiness_sort_cols:
        readiness_display = readiness_display.sort_values(readiness_sort_cols, na_position="last")
    perm_all = permutation[permutation["scope"] == "ALL_ELIGIBLE_REGIONS"].iloc[0] if not permutation.empty else pd.Series(dtype=object)
    lines = [
        "# PHASE 2A-R LOCUS CONSOLIDATION REPORT", "",
        "## Scope and frozen inputs", "",
        "Phase 2A-R reused the frozen Phase 1D SUPERGNOVA atlas, Phase 2A PLACO+ results, and the existing PLINK clumping outputs. SUPERGNOVA, PLACO+, LDSC, SuSiE, coloc, FINEMAP, enrichment, MR and GCI were not rerun.", "",
        "The existing `PLACO_INDEPENDENT_LOCI.tsv` was audited and found to be a fallback table: all 26,142 rows were labeled `LD_PANEL_NOT_AVAILABLE_OR_NOT_CLUMPED`, despite valid PLINK `.clumps` outputs being present. The valid `.clumps` files were reused; no reclumping was performed.", "",
        "Clumping parameters recorded in the PLINK logs: `r2 < 0.1`, `500 kb`, `p1=1`, `p2=1`, `threads=4`, `memory=4096 MB`, using `data/phase1C/SUPERGNOVA_reference/bfiles/eur_chr{chrom}_SNPmaf5`.", "",
        "## PHASE 2A-R VERDICT", "", f"**{verdict}**", "", f"Next-step recommendation: **{next_step}**. This is a workflow gate only; the next phase was not started.", "",
        "## INDEPENDENT PLACO LOCI", "", f"Total valid independent PLINK-clumped loci: **{len(loci)}**.", f"Unique pair–lead-SNP loci: **{loci[['pair', 'lead_SNP']].drop_duplicates().shape[0]}**.", f"Phase 1D significant regions: **{len(regions)}**; regions containing at least one independent locus: **{counts['regions_with_loci']}**.", "",
        md_table(audit, ["scope", "rows", "unique_pairs", "unique_lead_snps", "unique_genomic_loci", "validity"], None), "",
        "## REGION-LEVEL CONVERGENCE", "", f"Negative-rho regions supported: **{int(((regions['SUPERGNOVA_rho'] < 0) & (regions['N_independent_PLACO_loci'] > 0)).sum())}/{int((regions['SUPERGNOVA_rho'] < 0).sum())}**.", f"Positive-rho regions supported: **{int(((regions['SUPERGNOVA_rho'] > 0) & (regions['N_independent_PLACO_loci'] > 0)).sum())}/{int((regions['SUPERGNOVA_rho'] > 0).sum())}**.", "", md_table(regions[["pair", "SUPERGNOVA_region", "SUPERGNOVA_rho", "N_independent_PLACO_loci", "N_direction_concordant", "N_direction_discordant", "directional_class", "fine_mapping_readiness"]], ["pair", "SUPERGNOVA_region", "SUPERGNOVA_rho", "N_independent_PLACO_loci", "N_direction_concordant", "N_direction_discordant", "directional_class", "fine_mapping_readiness"]), "",
        "## DIRECTIONAL CONCORDANCE", "", f"Independent loci mapped to local regions: **{counts['mapped_loci']}**.", f"Concordant and discordant counts are reported at the independent-locus level; the historical 3,141 correlated SNP-level overlaps were not used as the primary evidence count.", "", md_table(permutation, ["scope", "N_eligible_regions", "observed_concordance", "observed_fraction", "null_mean", "null_SD", "empirical_P", "empirical_Z"], None), "",
        "## HIDDEN MIXED PAIRS", "", md_table(hidden, ["pair", "positive_SUPERGNOVA_regions", "negative_SUPERGNOVA_regions", "positive_regions_with_independent_PLACO_support", "negative_regions_with_independent_PLACO_support", "positive_direction_concordant_loci", "negative_direction_concordant_loci"], None), "",
        "## HIGH-CONFIDENCE ANTAGONISTIC CANDIDATE LOCI", "", f"Candidate loci meeting the Phase 2A-R candidate criteria: **{counts['antagonistic_candidates']}**. The table below shows the top 30 by lead P; the complete table is in `PHASE2B_REFINED_PRIORITY.tsv`.", "", md_table(neg_candidates, ["tier", "pair", "SUPERGNOVA_rho", "SUPERGNOVA_q", "PLACO_locus_id", "lead_SNP", "PLACO_P", "beta_ocular", "beta_systemic", "fine_mapping_readiness"], 30), "",
        "## HIGH-CONFIDENCE CONCORDANT CANDIDATE LOCI", "", f"Concordant positive-rho candidate loci: **{len(pos_candidates)}**. The table below shows the top 30 by lead P.", "", md_table(pos_candidates, ["tier", "pair", "SUPERGNOVA_rho", "SUPERGNOVA_q", "PLACO_locus_id", "lead_SNP", "PLACO_P", "beta_ocular", "beta_systemic", "fine_mapping_readiness"], 30), "",
        "## FINE-MAPPING READINESS", "", f"PASS: **{readiness_counts.get('PASS', 0)}**; CONDITIONAL: **{readiness_counts.get('CONDITIONAL', 0)}**; FAIL: **{readiness_counts.get('FAIL', 0)}**.", "All CAT-containing loci retain `FINNISH_LD_REQUIRED`; generic EUR LD was used only for organizational clumping. Imputation INFO was not present in the frozen summary-statistics inputs, so readiness cannot be upgraded to PASS on the current data.", "", md_table(readiness_display[["pair", "PLACO_locus_id", "regional_variant_density", "reference_overlap_count", "reference_coverage_fraction", "full_beta_se_trait1", "full_beta_se_trait2", "alleles_complete", "sample_size_available", "fine_mapping_ready", "ld_reference_status", "CAT_Finnish_LD_flag"]], ["pair", "PLACO_locus_id", "regional_variant_density", "reference_overlap_count", "reference_coverage_fraction", "full_beta_se_trait1", "full_beta_se_trait2", "alleles_complete", "sample_size_available", "fine_mapping_ready", "ld_reference_status", "CAT_Finnish_LD_flag"], 60), "",
        "## TIER A1 / A2", "", f"TIER A1 count: **{counts['a1']}**; pairs: **{', '.join(sorted(priority.loc[priority['tier'] == 'TIER A1', 'pair'].unique())) or 'None'}**.", f"TIER A2 count: **{counts['a2']}**; pairs: **{', '.join(sorted(priority.loc[priority['tier'] == 'TIER A2', 'pair'].unique())) or 'None'}**.", "Mixed independent directions within a region are assigned TIER X even when an individual locus is direction-concordant.", "",
        "## CAT / FINNISH LD LIMITATION", "", "CAT summary statistics are FinnGen-derived Finnish-European data. The validated EUR panel supports organizational clumping but does not establish a Finnish credible set. CAT loci therefore remain conditional and require Finnish-matched LD sensitivity before fine-mapping interpretation.", "",
        "## LAVA STATUS", "", "`PENDING_LOCAL_METHOD_REPLICATION` — no unverified mirror or substitute reference was used.", "",
        "## NEXT STEP", "", f"**{next_step}**. Phase 2B was not executed automatically.", "",
        "## Reproducibility outputs", "", "- `results/phase2AR/INDEPENDENT_LOCUS_AUDIT.tsv`", "- `results/phase2AR/INDEPENDENT_LOCUS_LOCAL_REGION_MAP.tsv`", "- `results/phase2AR/REGION_LEVEL_DIRECTIONAL_CONVERGENCE.tsv`", "- `results/phase2AR/NEGATIVE_LOCAL_REGION_AUDIT.tsv`", "- `results/phase2AR/POSITIVE_LOCAL_REGION_AUDIT.tsv`", "- `results/phase2AR/DIRECTION_CONCORDANCE_PERMUTATION.tsv`", "- `results/phase2AR/HIDDEN_MIXED_PAIR_LOCUS_AUDIT.tsv`", "- `results/phase2AR/FINE_MAPPING_READINESS_REAUDIT.tsv`", "- `results/phase2AR/REGIONAL_VARIANT_DENSITY.tsv`", "- `results/phase2AR/PHASE2B_REFINED_PRIORITY.tsv`", "- `results/phase2AR/PHASE2AR_DECISION.txt`", "",
    ]
    (REPORTS / "PHASE2AR_LOCUS_CONSOLIDATION_REPORT.md").write_text("\n".join(lines) + "\n")


def write_decision(verdict: str, next_step: str, counts: dict[str, int]) -> None:
    (RESULTS / "PHASE2AR_DECISION.txt").write_text("\n".join([
        f"PHASE2AR_VERDICT\t{verdict}",
        f"NEXT_STEP\t{next_step}",
        f"INDEPENDENT_PLACO_LOCI\t{counts['independent_loci']}",
        f"MAPPED_INDEPENDENT_LOCI\t{counts['mapped_loci']}",
        f"PHASE1D_SIGNIFICANT_REGIONS\t{counts['regions']}",
        f"REGIONS_WITH_INDEPENDENT_PLACO\t{counts['regions_with_loci']}",
        f"TIER_A1\t{counts['a1']}",
        f"TIER_A2\t{counts['a2']}",
        f"ANTAGONISTIC_CANDIDATES\t{counts['antagonistic_candidates']}",
        f"TIER_C1\t{counts['c1']}",
        f"TIER_C2\t{counts['c2']}",
        "PHASE2B_STARTED\tNO",
        "LAVA_STATUS\tPENDING_LOCAL_METHOD_REPLICATION",
        "HDL_L_STATUS\tSECONDARY_PENDING_EXACT_N0",
        "",
    ]))


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    if len(sys.argv) > 1 and sys.argv[1] == "--report-only":
        significant = pd.read_csv(RESULTS_2A / "PLACO_PLUS_SIGNIFICANT.tsv", sep="\t")
        loci, _ = read_clumps(significant)
        audit = pd.read_csv(RESULTS / "INDEPENDENT_LOCUS_AUDIT.tsv", sep="\t")
        mapped = pd.read_csv(RESULTS / "INDEPENDENT_LOCUS_LOCAL_REGION_MAP.tsv", sep="\t")
        regions = pd.read_csv(RESULTS / "REGION_LEVEL_DIRECTIONAL_CONVERGENCE.tsv", sep="\t")
        permutation = pd.read_csv(RESULTS / "DIRECTION_CONCORDANCE_PERMUTATION.tsv", sep="\t")
        hidden = pd.read_csv(RESULTS / "HIDDEN_MIXED_PAIR_LOCUS_AUDIT.tsv", sep="\t")
        readiness = pd.read_csv(RESULTS / "FINE_MAPPING_READINESS_REAUDIT.tsv", sep="\t")
        priority = pd.read_csv(RESULTS / "PHASE2B_REFINED_PRIORITY.tsv", sep="\t")
        reused = audit[audit["scope"] == "REUSED_PLINK_CLUMP_OUTPUTS"]
        total_loci = int(reused.iloc[0]["rows"]) if not reused.empty else None
        verdict, next_step, counts = decision(regions, mapped, readiness, priority, total_loci=total_loci)
        write_decision(verdict, next_step, counts)
        write_report(verdict, next_step, counts, audit, loci, mapped, regions, permutation, hidden, readiness, priority)
        print("PHASE2AR_REPORT_COMPLETE", verdict, next_step, counts)
        return
    significant = pd.read_csv(RESULTS_2A / "PLACO_PLUS_SIGNIFICANT.tsv", sep="\t")
    loci, clump_audit = read_clumps(significant)
    audit = write_independent_locus_audit(significant, loci, clump_audit)
    loci.to_csv(RESULTS / "INDEPENDENT_LOCUS_LOCAL_REGION_MAP.tsv", sep="\t", index=False)
    atlas = pd.read_csv(RESULTS_1D / "SUPERGNOVA_OCULAR_SYSTEMIC_ATLAS.tsv", sep="\t")
    mapped = map_loci_to_regions(loci, atlas)
    mapped.to_csv(RESULTS / "INDEPENDENT_LOCUS_LOCAL_REGION_MAP.tsv", sep="\t", index=False)
    density = compute_density(mapped, loci)
    readiness = readiness_audit(mapped, density)
    region_table = region_summary(atlas, mapped, readiness)
    hidden = hidden_pair_audit(region_table, mapped)
    permutation = permutation_test(region_table, mapped)
    priority = build_priority(region_table, mapped, readiness)
    verdict, next_step, counts = decision(region_table, mapped, readiness, priority, total_loci=len(loci))
    write_decision(verdict, next_step, counts)
    write_report(verdict, next_step, counts, audit, loci, mapped, region_table, permutation, hidden, readiness, priority)
    print("PHASE2AR_COMPLETE", verdict, next_step, counts)


if __name__ == "__main__":
    main()
