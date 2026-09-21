#!/usr/bin/env python3
"""Finalize the prespecified Phase 2B-R Pan-UKBB LD rescue audit."""

from __future__ import annotations

import csv
import datetime as dt
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "phase2BR"
REPORTS = ROOT / "reports" / "phase2BR"
ANALYZED = [f"P2BP{i:02d}" for i in range(2, 9)]
PAIR01 = "P2BP01"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        out = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        out.writeheader()
        out.writerows({k: "NA" if v is None else v for k, v in r.items()} for r in rows)


def as_float(value: object) -> float:
    try:
        x = float(value)
        return x if math.isfinite(x) else math.nan
    except (TypeError, ValueError):
        return math.nan


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def fmt(value: object, digits: int = 4) -> str:
    x = as_float(value)
    return "NA" if math.isnan(x) else f"{x:.{digits}f}"


def md_table(headers: list[str], rows: list[list[object]]) -> str:
    def cell(x: object) -> str:
        return str(x).replace("|", "\\|").replace("\n", " ")

    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    lines.extend("| " + " | ".join(cell(x) for x in row) + " |" for row in rows)
    return "\n".join(lines)


def main() -> None:
    access = read_tsv(RES / "PANUKBB_ACCESS_AUDIT.tsv")
    match = read_tsv(RES / "PANUKBB_VARIANT_MATCH_QC.tsv")
    allele = read_tsv(RES / "PANUKBB_ALLELE_ALIGNMENT.tsv")
    diag = read_tsv(RES / "PANUKBB_SUSIE_DIAGNOSTICS.tsv")
    comparison = read_tsv(RES / "REFERENCE_COMPARISON.tsv")
    coloc = read_tsv(RES / "COLOC_SUSIE_PANUKBB.tsv")
    priors = read_tsv(RES / "COLOC_PRIOR_SENSITIVITY.tsv")
    ser = read_tsv(RES / "SER_FALLBACK_DIAGNOSTIC.tsv")

    if not access or not match or not diag or not comparison:
        missing = [
            str(p.relative_to(ROOT))
            for p, rows in [
                (RES / "PANUKBB_ACCESS_AUDIT.tsv", access),
                (RES / "PANUKBB_VARIANT_MATCH_QC.tsv", match),
                (RES / "PANUKBB_SUSIE_DIAGNOSTICS.tsv", diag),
                (RES / "REFERENCE_COMPARISON.tsv", comparison),
            ]
            if not rows
        ]
        raise SystemExit("Phase 2B-R inputs are incomplete: " + ", ".join(missing))

    sample_rows = [r for r in access if r.get("check") == "n_samples"]
    n_samples = as_float(sample_rows[0].get("value")) if sample_rows else math.nan
    access_ok = bool(sample_rows) and n_samples == 420542 and all(
        r.get("status", "").upper() in {"PASS", "OK", "VALID"} for r in access
        if r.get("status")
    )
    match_by = {r.get("pilot_id"): r for r in match}
    analyzed_match = [match_by[p] for p in ANALYZED if p in match_by]
    coverage_pass = sum(r.get("coverage_gate") == "PASS" for r in analyzed_match)
    complete_ld = sum((RES / p / "PANUKBB_LD.tsv").exists() for p in ANALYZED)

    primary = [r for r in diag if r.get("model") == "pan_mismatch" and r.get("pilot_id") in ANALYZED]
    primary_by = {(r.get("pilot_id"), r.get("trait")): r for r in primary}
    reliable = [r for r in primary if r.get("reliability_class") == "RELIABLE" and as_bool(r.get("R_reliability_flag")) is False]
    rescued = [r for r in comparison if r.get("rescue_class") == "RESCUED"]
    partial = [r for r in comparison if r.get("rescue_class") == "PARTIALLY_RESCUED"]
    not_rescued = [r for r in comparison if r.get("rescue_class") == "NOT_RESCUED"]
    reliable_pairs = 0
    pair_rows: list[list[object]] = []
    for p in ANALYZED:
        rows = [r for r in primary if r.get("pilot_id") == p]
        both = len(rows) == 2 and all(r.get("reliability_class") == "RELIABLE" and not as_bool(r.get("R_reliability_flag")) for r in rows)
        if both:
            reliable_pairs += 1
        pair_rows.append([
            p,
            match_by.get(p, {}).get("pair", "NA"),
            match_by.get(p, {}).get("panukbb_included_variants", "NA"),
            match_by.get(p, {}).get("coverage_fraction", "NA"),
            "YES" if both else "NO",
            "; ".join(r.get("trait", "NA") + ":" + r.get("reliability_class", "NA") for r in rows) or "NA",
        ])

    strong_shared = sum(r.get("classification") == "STRONG_SHARED_SIGNAL" for r in coloc)
    moderate_shared = sum(r.get("classification") == "MODERATE_SHARED_SIGNAL" for r in coloc)
    antagonistic_upgrade = sum(r.get("upgrade_label") == "ANTAGONISTIC_UPGRADE" for r in coloc)
    concordant_upgrade = sum(r.get("upgrade_label") == "CONCORDANT_UPGRADE" for r in coloc)
    sensitivity_pass = sum(r.get("status") == "PASS" for r in priors)
    ser_n = len(ser)

    rescue_success = len(rescued) >= 8 and reliable_pairs >= 3
    if not access_ok:
        verdict = "PANUKBB_REFERENCE_ACCESS_BLOCKED"
    elif rescue_success:
        verdict = "PANUKBB_LD_RESCUE_SUCCESS"
    elif rescued or partial:
        verdict = "PARTIAL_LD_RESCUE"
    else:
        verdict = "NO_LD_RESCUE"

    if verdict == "PANUKBB_LD_RESCUE_SUCCESS" and (strong_shared or antagonistic_upgrade):
        next_step = "RESUME_LIMITED_FINE_MAPPING"
    elif verdict in {"NO_LD_RESCUE", "PANUKBB_REFERENCE_ACCESS_BLOCKED"} or reliable_pairs < 3:
        next_step = "STOP_FINE_MAPPING_PERMANENTLY_WITH_PUBLIC_LD"
    else:
        next_step = "RETAIN_LOCAL_VARIANT_LEVEL_EVIDENCE"

    now = dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")
    audit = f"""# Pan-UKBB EUR LD reference audit (Phase 2B-R)

Generated: {now}

## Scope and frozen constraints

This rescue audit re-analyzed exactly P2BP02–P2BP08 from the frozen Phase 2B-P pilot. P2BP01 was not analyzed or replaced because the original pilot had 198 included variants and failed the prespecified minimum-200-variant LD coverage gate. No additional loci, PLACO+, CAT, LAVA, enrichment, MR, GCI, or Phase 2B expansion were run.

## Official reference and access

The official resource is the Pan-UKBB EUR LD release documented at [Pan-UKBB LD documentation](https://pan.ukbb.broadinstitute.org/docs/ld/index.html), using the Hail variant index `s3://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht` and BlockMatrix `s3://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm`. The variant-index global metadata supplied the exact EUR `n_samples={n_samples:.0f}`; this value was used as B and was not hardcoded in the model.

The resource uses GRCh37 coordinates, dosage-based covariate-adjusted Pearson LD, and the published INFO>0.8/MAC>20 variant filters within a 10-Mb LD radius. Access used Hail {access[0].get('hail_version', '0.2.135') if access else '0.2.135'}, Python 3.11, OpenJDK 17, Spark/Hadoop S3A with anonymous read-only access. Only locus-specific variant-index rows and the required BlockMatrix blocks were read; the full 43.3-TB matrix was not downloaded.

## Access and regional matching

{md_table(['Metric', 'Value'], [
    ['n_samples', f'{n_samples:.0f}'],
    ['Eligible loci', f'{len(ANALYZED)}/7'],
    ['Regional coverage gates passed', f'{coverage_pass}/7'],
    ['Complete Pan LD matrices', f'{complete_ld}/7'],
    ['Allele-alignment audit records', f'{len(allele)}'],
    ['P2BP01 status', 'NOT_ANALYZED_PHASE2BR'],
])}

The detailed access and variant/allele audits are retained in `results/phase2BR/PANUKBB_ACCESS_AUDIT.tsv`, `results/phase2BR/PANUKBB_VARIANT_MATCH_QC.tsv`, and `results/phase2BR/PANUKBB_ALLELE_ALIGNMENT.tsv`.
"""
    REPORTS.mkdir(parents=True, exist_ok=True)
    (REPORTS / "PANUKBB_LD_REFERENCE_AUDIT.md").write_text(audit, encoding="utf-8")

    comparison_rows = []
    for r in comparison:
        comparison_rows.append([
            r.get("pilot_id", "NA"), r.get("trait", "NA"), r.get("old_R_reliability_flag", "NA"),
            r.get("new_R_reliability_flag", "NA"), r.get("rescue_class", "NA"),
            fmt(r.get("old_Q_art")), fmt(r.get("new_Q_art")), fmt(r.get("old_B_corrected"), 1),
            fmt(r.get("new_B_corrected"), 1), fmt(r.get("old_top_PIP")), fmt(r.get("new_top_PIP")),
            r.get("old_credible_set_sizes", "NA"), r.get("new_credible_set_sizes", "NA"),
        ])
    coloc_rows = [[
        r.get("pilot_id", "NA"), r.get("signal1", "NA"), r.get("signal2", "NA"),
        fmt(r.get("PP.H3")), fmt(r.get("PP.H4")), fmt(r.get("H4_H3")),
        r.get("classification", "NA"), r.get("status", "NA"), r.get("upgrade_label", "NA"),
    ] for r in coloc]
    ser_rows = [[r.get("pilot_id", "NA"), r.get("trait", "NA"), fmt(r.get("top_SER_PIP")), r.get("lead_variant", "NA"), r.get("SER_credible_set_width", "NA")] for r in ser]
    report = f"""# PHASE 2B-R — Pan-UKBB EUR LD rescue report

Generated: {now}

## Verdict

- Final verdict: **`{verdict}`**
- Next-step option: **`{next_step}`**
- Exact scope: 7 eligible Phase 2B-P loci (P2BP02–P2BP08), 14 trait-level fits.
- P2BP01 was not reanalyzed or replaced: its frozen Phase 2B-P status remains `LD_COVERAGE_FAIL`.

## Technical rescue results

{md_table(['Metric', 'Result', 'Prespecified interpretation'], [
    ['Pan-UKBB EUR n_samples', f'{n_samples:.0f}', 'Exact variant-index global metadata'],
    ['Regional match gate', f'{coverage_pass}/7', 'At least 80% coverage and >=200 variants where required'],
    ['Primary Pan mismatch-aware fits', f'{len(primary)}/14', 'R_finite=B, R_mismatch="eb", L=10, estimate_residual_variance=FALSE'],
    ['Reliable primary fits', f'{len(reliable)}/14', 'Converged, credible sets, purity QC, and no reliability flag'],
    ['Fits classified RESCUED', f'{len(rescued)}/14', 'Old mismatch flag TRUE; Pan primary flag FALSE'],
    ['Fits classified PARTIALLY_RESCUED', f'{len(partial)}/14', 'Reliability flag retained but diagnostics improved substantially'],
    ['Pairs with both traits reliable', f'{reliable_pairs}/7', 'Required technical success threshold: >=3/7'],
    ['Primary technical rescue gate', 'PASS' if rescue_success else 'FAIL', '>=8/14 rescued and >=3/7 pair-loci both reliable'],
])}

## Per-locus summary

{md_table(['Locus', 'Pair', 'Matched variants', 'Coverage', 'Both traits reliable', 'Primary classes'], pair_rows)}

## Trait-level old-vs-Pan comparison

The complete comparison is retained in `results/phase2BR/REFERENCE_COMPARISON.tsv`.

{md_table(['Locus', 'Trait', 'Old flag', 'Pan flag', 'Class', 'Old Q_art', 'Pan Q_art', 'Old Bcorr', 'Pan Bcorr', 'Old top PIP', 'Pan top PIP', 'Old CS sizes', 'Pan CS sizes'], comparison_rows)}

## coloc.susie and prior sensitivity

coloc.susie was run only where both Pan primary traits were reliable, converged, and passed purity QC. The primary prior was p1=p2=1e-4 and p12=5e-6; sensitivity priors were p12=1e-6 and 1e-5. The retained outputs are `results/phase2BR/COLOC_SUSIE_PANUKBB.tsv` and `results/phase2BR/COLOC_PRIOR_SENSITIVITY.tsv`.

{md_table(['Locus', 'Signal 1', 'Signal 2', 'PP.H3', 'PP.H4', 'H4/H3', 'Classification', 'Status', 'Upgrade'], coloc_rows) if coloc_rows else 'No locus passed the Pan primary reliability and purity gate; coloc.susie was not accepted for inference.'}

Prior-sensitivity rows with a successful coloc computation: **{sensitivity_pass}**. Strong shared-signal loci (H4>=0.8 and H4/H3>=3): **{strong_shared}**. Antagonistic upgrades: **{antagonistic_upgrade}**. Concordant upgrades: **{concordant_upgrade}**.

## SER fallback diagnostic

SER was used only as a diagnostic fallback when the Pan primary fit remained reliability-flagged; it was not used to upgrade the primary conclusion.

{md_table(['Locus', 'Trait', 'Top SER PIP', 'Lead', '95% CS/width'], ser_rows) if ser_rows else 'No SER fallback record was required.'}

## Boundaries and deferred analyses

- CAT-containing fine-mapping remains deferred because the required Finnish-matched LD reference was not part of this rescue.
- LAVA remains `PENDING_LOCAL_METHOD_REPLICATION`.
- HDL-L remains `SECONDARY_PENDING_EXACT_N0`.
- No Phase 2B expansion, PLACO+ rerun, CAT rerun, variant-level scan outside the frozen loci, enrichment, cell-type, MR, GCI, FINEMAP, or Phase 2B follow-on was started.

## Reproducibility files

- `scripts/phase2BR/extract_panukbb_ld.py`
- `scripts/phase2BR/run_panukbb_rescue.R`
- `results/phase2BR/PANUKBB_ACCESS_AUDIT.tsv`
- `results/phase2BR/PANUKBB_VARIANT_MATCH_QC.tsv`
- `results/phase2BR/PANUKBB_ALLELE_ALIGNMENT.tsv`
- `results/phase2BR/PANUKBB_SUSIE_DIAGNOSTICS.tsv`
- `results/phase2BR/REFERENCE_COMPARISON.tsv`
- `results/phase2BR/SER_FALLBACK_DIAGNOSTIC.tsv`
- `results/phase2BR/COLOC_SUSIE_PANUKBB.tsv`
- `results/phase2BR/COLOC_PRIOR_SENSITIVITY.tsv`
"""
    (REPORTS / "PHASE2BR_LD_RESCUE_REPORT.md").write_text(report, encoding="utf-8")

    decision = f"""PHASE 2B-R DECISION
final_verdict={verdict}
next_step={next_step}
panukbb_eur_n_samples={n_samples:.0f}
eligible_loci={len(ANALYZED)}
regional_coverage_pass={coverage_pass}/7
primary_fits={len(primary)}/14
reliable_primary_fits={len(reliable)}/14
rescued_fits={len(rescued)}/14
partially_rescued_fits={len(partial)}/14
reliable_pairs={reliable_pairs}/7
strong_shared_signals={strong_shared}
antagonistic_upgrades={antagonistic_upgrade}
concordant_upgrades={concordant_upgrade}
ser_fallback_records={ser_n}
p2bp01_status=NOT_ANALYZED_PHASE2BR
cat_status=FINNISH_LD_DEFERRED
lava_status=PENDING_LOCAL_METHOD_REPLICATION
hdl_l_status=SECONDARY_PENDING_EXACT_N0
"""
    (RES / "PHASE2BR_DECISION.txt").write_text(decision, encoding="utf-8")
    print(decision, end="")


if __name__ == "__main__":
    main()
