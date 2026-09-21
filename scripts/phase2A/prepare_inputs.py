#!/usr/bin/env python3
"""Prepare genome-wide, allele-harmonized PLACO+ inputs from frozen GWAS files."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import duckdb


ROOT = Path(__file__).resolve().parents[2]
HARM = ROOT / "data" / "harmonized"
REF_FREQ = ROOT / "results" / "phase2A" / "reference_freq" / "eur_reference.afreq.tsv"
REF_BIM = ROOT / "results" / "phase2A" / "reference_freq" / "eur_reference.bim.tsv"
OUT = ROOT / "results" / "phase2A" / "inputs"
PAIRS = [
    ("AMD", "CAD"), ("AMD", "STROKE"), ("AMD", "T2D"), ("AMD", "CKD"), ("AMD", "AD"), ("AMD", "PD"),
    ("POAG", "CAD"), ("POAG", "STROKE"), ("POAG", "T2D"), ("POAG", "CKD"), ("POAG", "AD"), ("POAG", "PD"),
    ("CAT", "CAD"), ("CAT", "STROKE"), ("CAT", "T2D"), ("CAT", "CKD"), ("CAT", "AD"), ("CAT", "PD"),
]


def q(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def columns(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return handle.readline().rstrip("\n").split("\t")


def trait_view(con: duckdb.DuckDBPyConnection, trait: str, alias: str) -> bool:
    path = HARM / f"{trait}.standardized.tsv.gz"
    cols = set(columns(path))
    eaf = "try_cast(EAF AS DOUBLE)" if "EAF" in cols else "CAST(NULL AS DOUBLE)"
    sql = f"""
    CREATE OR REPLACE TEMP VIEW {alias} AS
    SELECT SNP, A1, A2,
           try_cast(BETA AS DOUBLE) AS BETA,
           try_cast(SE AS DOUBLE) AS SE,
           try_cast(P AS DOUBLE) AS P,
           try_cast(N AS DOUBLE) AS N,
           {eaf} AS EAF,
           try_cast(CHR AS INTEGER) AS CHR,
           try_cast(BP AS BIGINT) AS BP
    FROM read_csv_auto({q(str(path))}, header=true, delim='\\t', compression='gzip',
                       sample_size=100000, nullstr='NA', ignore_errors=false)
    WHERE SNP IS NOT NULL
    QUALIFY row_number() OVER (PARTITION BY SNP ORDER BY CHR, BP) = 1
    """
    con.execute(sql)
    return "EAF" in cols


def main() -> None:
    if not REF_FREQ.exists() or not REF_BIM.exists():
        raise FileNotFoundError(f"Reference frequency/coordinate table is missing: {REF_FREQ} / {REF_BIM}")
    OUT.mkdir(parents=True, exist_ok=True)
    qc_rows: list[dict[str, object]] = []
    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    con.execute(f"""
        CREATE OR REPLACE TEMP VIEW ref_freq AS
        SELECT ID AS SNP, try_cast(CHROM AS INTEGER) AS CHR,
               try_cast(ALT_FREQS AS DOUBLE) AS ALT_FREQ,
               least(try_cast(ALT_FREQS AS DOUBLE), 1 - try_cast(ALT_FREQS AS DOUBLE)) AS REF_MAF
        FROM read_csv_auto({q(str(REF_FREQ))}, header=true, delim='\\t', sample_size=100000,
                           nullstr='NA', ignore_errors=true)
        WHERE ID IS NOT NULL
    """)
    con.execute(f"""
        CREATE OR REPLACE TEMP VIEW ref_bim AS
        SELECT SNP, try_cast(CHR AS INTEGER) AS REF_CHR, try_cast(BP AS BIGINT) AS REF_BP,
               REF AS REF_ALLELE, ALT AS ALT_ALLELE
        FROM read_csv_auto({q(str(REF_BIM))}, header=true, delim='\\t', sample_size=100000,
                           nullstr='NA', ignore_errors=true)
        WHERE SNP IS NOT NULL
        QUALIFY row_number() OVER (PARTITION BY SNP ORDER BY REF_CHR, REF_BP) = 1
    """)

    for trait1, trait2 in PAIRS:
        pair = f"{trait1}-{trait2}"
        eaf1 = trait_view(con, trait1, "t1")
        eaf2 = trait_view(con, trait2, "t2")
        con.execute("DROP VIEW IF EXISTS joined")
        con.execute("""
            CREATE TEMP VIEW joined AS
            WITH raw AS (
                SELECT t1.SNP,
                       t1.CHR AS CHR_1, t1.BP AS BP_1, t1.A1 AS A1_1, t1.A2 AS A2_1,
                       t2.CHR AS CHR_2, t2.BP AS BP_2, t2.A1 AS A1_2, t2.A2 AS A2_2,
                       t1.BETA AS BETA_1, t1.SE AS SE_1, t1.P AS P_1, t1.N AS N_1, t1.EAF AS EAF_1,
                       t2.BETA AS BETA_2_RAW, t2.SE AS SE_2, t2.P AS P_2, t2.N AS N_2, t2.EAF AS EAF_2,
                       r1.REF_MAF AS REF_MAF_1, r2.REF_MAF AS REF_MAF_2,
                       b1.REF_CHR AS REF_CHR_1, b1.REF_BP AS REF_BP_1,
                       b2.REF_CHR AS REF_CHR_2, b2.REF_BP AS REF_BP_2
                FROM t1
                INNER JOIN t2 USING (SNP)
                LEFT JOIN ref_freq r1 ON r1.SNP=t1.SNP AND r1.CHR=t1.CHR
                LEFT JOIN ref_freq r2 ON r2.SNP=t2.SNP AND r2.CHR=t2.CHR
                LEFT JOIN ref_bim b1 ON b1.SNP=t1.SNP
                LEFT JOIN ref_bim b2 ON b2.SNP=t2.SNP
            ), flags AS (
                SELECT *,
                       upper(A1_1) || upper(A2_1) IN ('AT','TA','CG','GC') AS ambiguous,
                       CHR_1 = CHR_2 AND BP_1 = BP_2 AS direct_position_match,
                       REF_CHR_1 IS NOT NULL AND REF_CHR_2 IS NOT NULL
                         AND ((CHR_1 = REF_CHR_1 AND BP_1 = REF_BP_1)
                              OR (CHR_2 = REF_CHR_2 AND BP_2 = REF_BP_2)) AS reference_position_match,
                       upper(A1_1)=upper(A1_2) AND upper(A2_1)=upper(A2_2) AS same_orientation,
                       upper(A1_1)=upper(A2_2) AND upper(A2_1)=upper(A1_2) AS swapped_orientation,
                       abs(EAF_1-EAF_2) <= 0.15 AS freq_same_ok,
                       abs(EAF_1-(1-EAF_2)) <= 0.15 AS freq_swap_ok
                FROM raw
            ), orient AS (
                SELECT *,
                       CASE
                         WHEN NOT (direct_position_match OR reference_position_match) THEN 'POSITION_MISMATCH'
                         WHEN same_orientation AND NOT ambiguous THEN 'SAME'
                         WHEN swapped_orientation AND NOT ambiguous THEN 'SWAPPED'
                         WHEN ambiguous AND same_orientation AND freq_same_ok AND NOT freq_swap_ok THEN 'SAME'
                         WHEN ambiguous AND swapped_orientation AND freq_swap_ok AND NOT freq_same_ok THEN 'SWAPPED'
                         ELSE 'UNRESOLVED'
                       END AS orientation
                FROM flags
            )
            SELECT *,
                   CASE WHEN orientation='SWAPPED' THEN -BETA_2_RAW ELSE BETA_2_RAW END AS BETA_2,
                   CASE WHEN orientation='SWAPPED' THEN -BETA_2_RAW/SE_2 ELSE BETA_2_RAW/SE_2 END AS Z_2,
                   BETA_1/SE_1 AS Z_1,
                   CASE WHEN reference_position_match THEN 'REFERENCE_RECONCILED' WHEN direct_position_match THEN 'DIRECT_MATCH' ELSE 'POSITION_MISMATCH' END AS coordinate_status,
                   least(EAF_1, 1-EAF_1) AS TRAIT_MAF_1,
                   least(EAF_2, 1-EAF_2) AS TRAIT_MAF_2,
                   coalesce(least(EAF_1, 1-EAF_1), REF_MAF_1) AS MAF_1,
                   coalesce(least(EAF_2, 1-EAF_2), REF_MAF_2) AS MAF_2,
                   CASE WHEN EAF_1 IS NOT NULL THEN 'TRAIT_EAF' WHEN REF_MAF_1 IS NOT NULL THEN 'EUR_REFERENCE_PANEL' ELSE 'MAF_FILTER_NOT_AVAILABLE' END AS MAF_SOURCE_1,
                   CASE WHEN EAF_2 IS NOT NULL THEN 'TRAIT_EAF' WHEN REF_MAF_2 IS NOT NULL THEN 'EUR_REFERENCE_PANEL' ELSE 'MAF_FILTER_NOT_AVAILABLE' END AS MAF_SOURCE_2
            FROM orient
        """)

        count_sql = """
        SELECT
          count(*) AS common_snp,
          sum(CASE WHEN direct_position_match THEN 1 ELSE 0 END) AS direct_position_match,
          sum(CASE WHEN reference_position_match THEN 1 ELSE 0 END) AS reference_position_match,
          sum(CASE WHEN direct_position_match OR reference_position_match THEN 1 ELSE 0 END) AS position_match,
          sum(CASE WHEN orientation='SAME' THEN 1 ELSE 0 END) AS same_orientation,
          sum(CASE WHEN orientation='SWAPPED' THEN 1 ELSE 0 END) AS swapped_orientation,
          sum(CASE WHEN ambiguous THEN 1 ELSE 0 END) AS strand_ambiguous,
          sum(CASE WHEN orientation='UNRESOLVED' THEN 1 ELSE 0 END) AS unresolved_orientation,
          sum(CASE WHEN orientation IN ('SAME','SWAPPED') AND isfinite(BETA_1) AND isfinite(SE_1) AND isfinite(BETA_2) AND isfinite(SE_2) AND SE_1>0 AND SE_2>0 THEN 1 ELSE 0 END) AS numeric_valid,
          sum(CASE WHEN orientation IN ('SAME','SWAPPED') AND (MAF_1 IS NULL OR MAF_2 IS NULL) THEN 1 ELSE 0 END) AS maf_unavailable,
          sum(CASE WHEN orientation IN ('SAME','SWAPPED') AND MAF_1 IS NOT NULL AND MAF_2 IS NOT NULL AND (MAF_1<0.01 OR MAF_2<0.01) THEN 1 ELSE 0 END) AS maf_below_01,
          sum(CASE WHEN orientation IN ('SAME','SWAPPED') AND MAF_1 IS NOT NULL AND MAF_2 IS NOT NULL AND MAF_1>=0.01 AND MAF_2>=0.01 AND isfinite(BETA_1) AND isfinite(SE_1) AND isfinite(BETA_2) AND isfinite(SE_2) AND SE_1>0 AND SE_2>0 THEN 1 ELSE 0 END) AS eligible
        FROM joined
        """
        stats = con.execute(count_sql).fetchone()
        names = [d[0] for d in con.description]
        qc = dict(zip(names, stats))
        qc.update({"pair": pair, "trait1": trait1, "trait2": trait2, "eaf_available_trait1": eaf1, "eaf_available_trait2": eaf2})

        out = OUT / f"{pair}.tsv.gz"
        query = """
        SELECT SNP,
               coalesce(REF_CHR_1, CHR_1) AS CHR, coalesce(REF_BP_1, BP_1) AS BP, A1_1 AS A1, A2_1 AS A2,
               BETA_1, SE_1, Z_1, P_1, EAF_1, N_1,
               BETA_2, SE_2, Z_2, P_2, EAF_2, N_2,
               MAF_1, MAF_2, MAF_SOURCE_1, MAF_SOURCE_2,
               orientation, coordinate_status
        FROM joined
        WHERE orientation IN ('SAME','SWAPPED')
          AND MAF_1 >= 0.01 AND MAF_2 >= 0.01
          AND isfinite(BETA_1) AND isfinite(SE_1) AND isfinite(BETA_2) AND isfinite(SE_2)
          AND isfinite(P_1) AND isfinite(P_2) AND SE_1 > 0 AND SE_2 > 0
        """
        con.execute(f"COPY ({query}) TO {q(str(out))} (FORMAT CSV, HEADER, DELIMITER '\\t', COMPRESSION GZIP)")
        qc_rows.append(qc)
        (OUT / f"{pair}.qc.json").write_text(json.dumps(qc, indent=2, default=str) + "\n")
        print(pair, qc["eligible"], flush=True)

    with (RESULTS := ROOT / "results" / "phase2A" / "PLACO_ALLELE_HARMONIZATION_QC.tsv").open("w", newline="") as handle:
        fields = list(qc_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(qc_rows)


if __name__ == "__main__":
    main()
