#!/usr/bin/env python3
"""Build the frozen Phase 3A locus table and prespecified robustness audits.

This script only consumes the frozen Phase 2A-R outputs.  It does not perform
new association testing, fine-mapping, colocalization, or LD inference.
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MAP = ROOT / "results/phase2AR/INDEPENDENT_LOCUS_LOCAL_REGION_MAP.tsv"
REGIONS = ROOT / "results/phase2AR/REGION_LEVEL_DIRECTIONAL_CONVERGENCE.tsv"
FREEZE = ROOT / "config/phase3A_LOCUS_FREEZE.tsv"
OUT = ROOT / "results/phase3A"


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def num(value, default=float("nan")):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def boolish(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def canonical_region(row):
    return f"{row['chr']}:{row['SUPERGNOVA_start']}-{row['SUPERGNOVA_end']}"


def main():
    loci = read_tsv(MAP)
    regions = read_tsv(REGIONS)
    assert len(loci) == 58, f"Expected 58 frozen loci, found {len(loci)}"
    supported_regions = [r for r in regions if int(float(r.get("N_independent_PLACO_loci", 0) or 0)) > 0]
    assert len(supported_regions) == 41, f"Expected 41 supported regions, found {len(supported_regions)}"

    region_lookup = {r["SUPERGNOVA_region"]: r for r in supported_regions}
    freeze_rows = []
    for row in loci:
        region = region_lookup.get(row["SUPERGNOVA_region"], {})
        rho_sign = row.get("rho_sign", "").strip()
        direction_class = (
            "NEGATIVE_DIRECTION_CONCORDANT" if rho_sign == "-" else
            "POSITIVE_DIRECTION_CONCORDANT" if rho_sign == "+" else "UNRESOLVED_DIRECTION"
        )
        freeze_rows.append({
            "pair": row["pair"],
            "hidden_mixed_pair": row.get("HIDDEN_MIXED_PAIR", "False"),
            "chr": row["chr"],
            "lead_SNP": row["lead_SNP"],
            "lead_BP": row["lead_BP"],
            "locus_start": row["locus_start"],
            "locus_end": row["locus_end"],
            "SUPERGNOVA_region": row["SUPERGNOVA_region"],
            "SUPERGNOVA_rho": row["SUPERGNOVA_rho"],
            "SUPERGNOVA_q": row["SUPERGNOVA_q_atlas"],
            "PLACO_P": row["P_PLACO_PLUS"],
            "PLACO_level": row.get("PLACO_LEVEL", ""),
            "beta_ocular": row["BETA_1"],
            "beta_systemic": row["BETA_2"],
            "effect_product": row["effect_product"],
            "direction_class": direction_class,
            "fine_mapping_status": region.get("fine_mapping_readiness", "CONDITIONAL"),
            "PLACO_locus_id": row["PLACO_locus_id"],
        })
    freeze_fields = [
        "pair", "hidden_mixed_pair", "chr", "lead_SNP", "lead_BP", "locus_start", "locus_end",
        "SUPERGNOVA_region", "SUPERGNOVA_rho", "SUPERGNOVA_q", "PLACO_P", "PLACO_level",
        "beta_ocular", "beta_systemic", "effect_product", "direction_class", "fine_mapping_status",
        "PLACO_locus_id",
    ]
    write_tsv(FREEZE, freeze_rows, freeze_fields)

    # Four prespecified nested statistical subsets.
    variant_threshold = 5e-8 / 18.0
    robust_rows = []
    for row in freeze_rows:
        q = num(row["SUPERGNOVA_q"])
        p = num(row["PLACO_P"])
        memberships = {
            "PRIMARY": True,
            "STRICT_LOCAL": q < 0.01,
            "STRICT_VARIANT": p < variant_threshold,
            "DOUBLE_STRICT": q < 0.01 and p < variant_threshold,
        }
        for subset, include in memberships.items():
            robust_rows.append({
                "subset": subset,
                "included": int(include),
                "pair": row["pair"],
                "hidden_mixed_pair": row["hidden_mixed_pair"],
                "chr": row["chr"],
                "lead_SNP": row["lead_SNP"],
                "SUPERGNOVA_region": row["SUPERGNOVA_region"],
                "rho_sign": "NEGATIVE" if "NEGATIVE" in row["direction_class"] else "POSITIVE",
                "SUPERGNOVA_q": row["SUPERGNOVA_q"],
                "PLACO_P": row["PLACO_P"],
                "local_q_threshold": "0.05" if subset == "PRIMARY" else "0.01" if subset in {"STRICT_LOCAL", "DOUBLE_STRICT"} else "0.05",
                "variant_p_threshold": "NA" if subset in {"PRIMARY", "STRICT_LOCAL"} else f"{variant_threshold:.12g}",
            })
    robust_fields = [
        "subset", "included", "pair", "hidden_mixed_pair", "chr", "lead_SNP", "SUPERGNOVA_region",
        "rho_sign", "SUPERGNOVA_q", "PLACO_P", "local_q_threshold", "variant_p_threshold",
    ]
    write_tsv(OUT / "ROBUSTNESS_SUBSETS.tsv", robust_rows, robust_fields)

    # LOCO and LOPO are leave-one-unit-out descriptive sensitivity summaries.
    def sensitivity_rows(exclude_field, values, label_field):
        rows = []
        all_pairs = {r["pair"] for r in freeze_rows}
        all_regions = {r["SUPERGNOVA_region"] for r in freeze_rows}
        all_mixed = {r["pair"] for r in freeze_rows if boolish(r["hidden_mixed_pair"])}
        for value in values:
            kept = [r for r in freeze_rows if r[exclude_field] != value]
            regs = {r["SUPERGNOVA_region"] for r in kept}
            neg = [r for r in kept if "NEGATIVE" in r["direction_class"]]
            pos = [r for r in kept if "POSITIVE" in r["direction_class"]]
            pairs = {r["pair"] for r in kept}
            mixed = {r["pair"] for r in kept if boolish(r["hidden_mixed_pair"])}
            rows.append({
                label_field: value,
                "n_loci_remaining": len(kept),
                "n_negative_loci": len(neg),
                "n_positive_loci": len(pos),
                "n_regions_remaining": len(regs),
                "n_supported_negative_regions": len({r["SUPERGNOVA_region"] for r in neg}),
                "n_supported_positive_regions": len({r["SUPERGNOVA_region"] for r in pos}),
                "n_pairs_remaining": len(pairs),
                "n_hidden_mixed_pairs_remaining": len(mixed),
                "direction_concordance_fraction": "1.0" if kept else "NA",
                "excluded_unit_note": "leave-one-chromosome-out" if label_field == "excluded_chromosome" else "leave-one-pair-out",
            })
        return rows

    chromosomes = sorted({r["chr"] for r in freeze_rows}, key=lambda x: int(x) if str(x).isdigit() else 99)
    pairs = sorted({r["pair"] for r in freeze_rows})
    audit_fields = [
        "excluded_chromosome", "n_loci_remaining", "n_negative_loci", "n_positive_loci", "n_regions_remaining",
        "n_supported_negative_regions", "n_supported_positive_regions", "n_pairs_remaining",
        "n_hidden_mixed_pairs_remaining", "direction_concordance_fraction", "excluded_unit_note",
    ]
    write_tsv(OUT / "LOCO_DIRECTION_ROBUSTNESS.tsv", sensitivity_rows("chr", chromosomes, "excluded_chromosome"), audit_fields)
    lopo_fields = [
        "excluded_pair", "n_loci_remaining", "n_negative_loci", "n_positive_loci", "n_regions_remaining",
        "n_supported_negative_regions", "n_supported_positive_regions", "n_pairs_remaining",
        "n_hidden_mixed_pairs_remaining", "direction_concordance_fraction", "excluded_unit_note",
    ]
    write_tsv(OUT / "LOPO_DIRECTION_ROBUSTNESS.tsv", sensitivity_rows("pair", pairs, "excluded_pair"), lopo_fields)

    # A base recurrence table is made before annotation and updated by the
    # annotation stage with mapped genes and functional support counts.
    recurrence = defaultdict(list)
    for row in freeze_rows:
        recurrence[row["SUPERGNOVA_region"]].append(row)
    recurrence_rows = []
    for region, members in sorted(recurrence.items()):
        pairs_here = sorted({r["pair"] for r in members})
        recurrence_rows.append({
            "genomic_region": region,
            "chr": members[0]["chr"],
            "region_start": region.split(":", 1)[1].split("-", 1)[0],
            "region_end": region.rsplit("-", 1)[1],
            "n_loci": len(members),
            "n_pairs": len(pairs_here),
            "pairs": ";".join(pairs_here),
            "n_negative_loci": sum("NEGATIVE" in r["direction_class"] for r in members),
            "n_positive_loci": sum("POSITIVE" in r["direction_class"] for r in members),
            "lead_variants": ";".join(sorted({r["lead_SNP"] for r in members})),
            "mapped_genes": "",
            "biological_support_status": "PENDING_ANNOTATION",
        })
    recurrence_fields = [
        "genomic_region", "chr", "region_start", "region_end", "n_loci", "n_pairs", "pairs",
        "n_negative_loci", "n_positive_loci", "lead_variants", "mapped_genes", "biological_support_status",
    ]
    write_tsv(OUT / "RECURRENT_CROSS_PAIR_LOCI.tsv", recurrence_rows, recurrence_fields)

    summary = [
        f"frozen_loci={len(loci)}",
        f"supported_regions={len(supported_regions)}",
        f"negative_loci={sum('NEGATIVE' in r['direction_class'] for r in freeze_rows)}",
        f"positive_loci={sum('POSITIVE' in r['direction_class'] for r in freeze_rows)}",
        f"variant_threshold={variant_threshold:.12g}",
        "input=results/phase2AR/INDEPENDENT_LOCUS_LOCAL_REGION_MAP.tsv",
        "new_association_testing=NO",
        "fine_mapping=NOT_RUN_PUBLIC_LD_ROUTE_CLOSED",
    ]
    (OUT / "CORE_AUDIT.txt").write_text("\n".join(summary) + "\n")
    print("Phase 3A core complete:", len(loci), "loci;", len(supported_regions), "supported regions")


if __name__ == "__main__":
    main()
