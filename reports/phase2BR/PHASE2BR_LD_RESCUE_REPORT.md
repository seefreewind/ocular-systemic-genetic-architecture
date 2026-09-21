# PHASE 2B-R — Pan-UKBB EUR LD rescue report

Generated: 2026-09-21T16:42:07+08:00

## Verdict

- Final verdict: **`PARTIAL_LD_RESCUE`**
- Next-step option: **`STOP_FINE_MAPPING_PERMANENTLY_WITH_PUBLIC_LD`**
- Exact scope: 7 eligible Phase 2B-P loci (P2BP02–P2BP08), 14 trait-level fits.
- P2BP01 was not reanalyzed or replaced: its frozen Phase 2B-P status remains `LD_COVERAGE_FAIL`.

## Technical rescue results

| Metric | Result | Prespecified interpretation |
|---|---|---|
| Pan-UKBB EUR n_samples | 420542 | Exact variant-index global metadata |
| Regional match gate | 7/7 | At least 80% coverage and >=200 variants where required |
| Primary Pan mismatch-aware fits | 14/14 | R_finite=B, R_mismatch="eb", L=10, estimate_residual_variance=FALSE |
| Reliable primary fits | 7/14 | Converged, credible sets, purity QC, and no reliability flag |
| Fits classified RESCUED | 13/14 | Old mismatch flag TRUE; Pan primary flag FALSE |
| Fits classified PARTIALLY_RESCUED | 1/14 | Reliability flag retained but diagnostics improved substantially |
| Pairs with both traits reliable | 1/7 | Required technical success threshold: >=3/7 |
| Primary technical rescue gate | FAIL | >=8/14 rescued and >=3/7 pair-loci both reliable |

## Per-locus summary

| Locus | Pair | Matched variants | Coverage | Both traits reliable | Primary classes |
|---|---|---|---|---|---|
| P2BP02 | AMD-T2D | 521 | 1.0 | NO | AMD:INDETERMINATE; T2D:RELIABLE |
| P2BP03 | AMD-CKD | 480 | 1.0 | NO | AMD:RELIABLE; CKD:INDETERMINATE |
| P2BP04 | AMD-STROKE | 384 | 1.0 | NO | AMD:INDETERMINATE; STROKE:RELIABLE |
| P2BP05 | POAG-CKD | 626 | 1.0 | NO | POAG:INDETERMINATE; CKD:INDETERMINATE |
| P2BP06 | AMD-T2D | 447 | 0.9977678571428571 | NO | AMD:INDETERMINATE; T2D:RELIABLE |
| P2BP07 | AMD-CAD | 283 | 1.0 | YES | AMD:RELIABLE; CAD:RELIABLE |
| P2BP08 | AMD-STROKE | 560 | 0.9982174688057041 | NO | AMD:RELIABLE; STROKE:LD_MISMATCH_LIMITED |

## Trait-level old-vs-Pan comparison

The complete comparison is retained in `results/phase2BR/REFERENCE_COMPARISON.tsv`.

| Locus | Trait | Old flag | Pan flag | Class | Old Q_art | Pan Q_art | Old Bcorr | Pan Bcorr | Old top PIP | Pan top PIP | Old CS sizes | Pan CS sizes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P2BP02 | AMD | TRUE | FALSE | RESCUED | 0.4358 | 0.0131 | 0.1 | 14.2 | 0.0247 | 0.0665 | NA | NA |
| P2BP02 | T2D | TRUE | FALSE | RESCUED | 0.4495 | 0.0146 | 8.5 | 420542.0 | 0.4333 | 0.6523 | 8 | 5 |
| P2BP03 | AMD | TRUE | FALSE | RESCUED | 0.3031 | 0.0143 | 2.3 | 420542.0 | 0.1183 | 0.1183 | 17 | 17 |
| P2BP03 | CKD | TRUE | FALSE | RESCUED | 0.4956 | 0.0230 | 0.2 | 2.5 | 0.0217 | 0.0217 | NA | NA |
| P2BP04 | AMD | TRUE | FALSE | RESCUED | 0.4306 | 0.0117 | 0.1 | 10.1 | 0.0285 | 0.0476 | NA | NA |
| P2BP04 | STROKE | TRUE | FALSE | RESCUED | 0.5550 | 0.0933 | 0.6 | 420542.0 | 0.0980 | 0.0980 | 47 | 47 |
| P2BP05 | POAG | TRUE | FALSE | RESCUED | 0.2048 | 0.0294 | 3.4 | 1332.7 | 0.0813 | 0.0813 | NA | NA |
| P2BP05 | CKD | TRUE | FALSE | RESCUED | 0.3422 | 0.0170 | 0.3 | 3.3 | 0.0254 | 0.0254 | NA | NA |
| P2BP06 | AMD | TRUE | FALSE | RESCUED | 0.2214 | 0.0053 | 6.6 | 10.3 | 0.2978 | 0.3112 | NA | NA |
| P2BP06 | T2D | TRUE | FALSE | RESCUED | 0.3284 | 0.0032 | 17.8 | 50.9 | 0.3366 | 0.2728 | 7 | 11 |
| P2BP07 | AMD | TRUE | FALSE | RESCUED | 0.1770 | 0.0107 | 10.1 | 420542.0 | 0.7780 | 0.7780 | 2 | 2 |
| P2BP07 | CAD | TRUE | FALSE | RESCUED | 0.0706 | 0.0309 | 359.5 | 420542.0 | 1.0000 | 0.5207 | 1;2;2;3;3;5;8 | 10;5 |
| P2BP08 | AMD | TRUE | FALSE | RESCUED | 0.0662 | 0.0173 | 503.0 | 420542.0 | 0.9637 | 0.0640 | 4;3;10 | 38 |
| P2BP08 | STROKE | TRUE | TRUE | PARTIALLY_RESCUED | 0.0804 | 0.1303 | 503.0 | 420542.0 | 1.0000 | 0.9671 | 1;2;2 | 1 |

## coloc.susie and prior sensitivity

coloc.susie was run only where both Pan primary traits were reliable, converged, and passed purity QC. The primary prior was p1=p2=1e-4 and p12=5e-6; sensitivity priors were p12=1e-6 and 1e-5. The retained outputs are `results/phase2BR/COLOC_SUSIE_PANUKBB.tsv` and `results/phase2BR/COLOC_PRIOR_SENSITIVITY.tsv`.

| Locus | Signal 1 | Signal 2 | PP.H3 | PP.H4 | H4/H3 | Classification | Status | Upgrade |
|---|---|---|---|---|---|---|---|---|
| P2BP02 | NA | NA | NA | NA | NA | LD_LIMITED | LD_LIMITED | NO_UPGRADE |
| P2BP03 | NA | NA | NA | NA | NA | LD_LIMITED | LD_LIMITED | NO_UPGRADE |
| P2BP04 | NA | NA | NA | NA | NA | LD_LIMITED | LD_LIMITED | NO_UPGRADE |
| P2BP05 | NA | NA | NA | NA | NA | LD_LIMITED | LD_LIMITED | NO_UPGRADE |
| P2BP06 | NA | NA | NA | NA | NA | LD_LIMITED | LD_LIMITED | NO_UPGRADE |
| P2BP07 | 1 | 1 | 1.0000 | 0.0000 | 0.0000 | DISTINCT_SIGNALS | PASS | NO_UPGRADE |
| P2BP07 | 1 | 2 | 0.0151 | 0.9756 | 64.4260 | STRONG_SHARED_SIGNAL | PASS | NO_UPGRADE |
| P2BP08 | NA | NA | NA | NA | NA | LD_LIMITED | LD_LIMITED | NO_UPGRADE |

Prior-sensitivity rows with a successful coloc computation: **6**. Strong shared-signal loci (H4>=0.8 and H4/H3>=3): **1**. Antagonistic upgrades: **0**. Concordant upgrades: **0**.

## SER fallback diagnostic

SER was used only as a diagnostic fallback when the Pan primary fit remained reliability-flagged; it was not used to upgrade the primary conclusion.

| Locus | Trait | Top SER PIP | Lead | 95% CS/width |
|---|---|---|---|---|
| P2BP08 | STROKE | 0.9671 | 12:112072424:A:G | 1 |

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
