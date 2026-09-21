# Phase 1A — Global genetic architecture and sample-overlap estimation

## Scope and stop rule

This report implements the amended protocol and stops after global LDSC, the LDSC cross-trait intercept matrix, overlap metadata, and AMD global sensitivity. No local bivariate analysis, HDL-L full pairwise analysis, LAVA local rg, PLACO+, fine-mapping, coloc, GCI, covariance decomposition, or enrichment was run.

## Protocol amendment

- `AMENDMENT_STATUS=SCIENTIFICALLY_JUSTIFIED`.
- Historical Phase 0 remains `HOLD_UNDER_ORIGINAL_PROTOCOL`.
- Current Phase 1A status is `PROCEED_UNDER_AMENDED_PROTOCOL`.
- Global LDSC h2 Z ≥ 4 remains preferred. Traits with 2 ≤ Z < 4 may enter local testing only after the amended QC, independent local-h2 reconstruction, signed-statistic validation, region-level multiple-testing control, and pairwise uncertainty requirements are met.
- AMD local univariate heritability is analytically separate from the HDL-L pairwise covariance question. The former may be eligible under the amended rule; the latter remains unresolved when N0 is unknown.

## Frozen inputs and global LDSC

The analysis used nine EUR/EUR-like frozen traits and 36/36 unique unconstrained cross-trait LDSC pairs. The same EUR LD-score reference, regression weights, HM3 merge, and Phase 0 MHC handling were used. No intercept was constrained to zero.

| Trait | Dataset | h2 | SE | Z | LDSC intercept | HM3 SNPs |
|---|---|---:|---:|---:|---:|---:|
| AMD | AMD_IAMDGC2_EUR_LATE | 0.5627 | 0.2064 | 2.7263 | 1.0827 | 1022570 |
| POAG | GCST90011766_EUR_stage1 | 0.2152 | 0.0163 | 13.202 | 1.0317 | 1183069 |
| CAT | FinnGen_R10_H7_CATARACTSENILE | 0.0555 | 0.005 | 11.1 | 1.0897 | 1160872 |
| CAD | GCST90132314_EUR_primary_discovery | 0.0747 | 0.0044 | 16.977 | 1.0192 | 1182021 |
| STROKE | GCST90104540_EUR_GIGASTROKE_AIS | 0.0283 | 0.0023 | 12.304 | 1.0719 | 1175792 |
| T2D | DIAMANTE_EUR_2022 | 0.1483 | 0.0067 | 22.134 | 0.997 | 1171068 |
| CKD | CKDGen_Wuttke2019_EA | 0.0432 | 0.0061 | 7.082 | 1.0754 | 1175015 |
| AD | GCST90027158_EUR_Bellenguez2022 | 0.0462 | 0.0045 | 10.267 | 1.0695 | 1167505 |
| PD | GCST009325_EUR_Nalls2019 | 0.0586 | 0.0054 | 10.852 | 0.9847 | 1132078 |

### Multiplicity-corrected significant global rg

| Trait 1 | Trait 2 | rg | SE | P | BH FDR | Bonferroni P |
|---|---|---:|---:|---:|---:|---:|
| CAD | T2D | 0.3944 | 0.021 | 1.7721e-78 | 6.37956e-77 | 6.37956e-77 |
| CAD | STROKE | 0.5081 | 0.0334 | 2.5141e-52 | 4.52538e-51 | 9.05076e-51 |
| STROKE | T2D | 0.3521 | 0.03 | 1.0007e-31 | 1.20084e-30 | 3.60252e-30 |
| CAT | CAD | 0.2159 | 0.0325 | 3.009e-11 | 2.7081e-10 | 1.08324e-09 |
| CAD | AD | -0.2216 | 0.0379 | 4.9093e-09 | 3.5347e-08 | 1.76735e-07 |
| CAT | T2D | 0.1459 | 0.0318 | 4.6198e-06 | 2.77188e-05 | 0.000166313 |
| AD | PD | 0.2051 | 0.0549 | 0.0002 | 0.00102857 | 0.0072 |
| CAD | CKD | 0.1344 | 0.0453 | 0.003 | 0.0135 | 0.108 |
| POAG | CAT | 0.1234 | 0.0429 | 0.004 | 0.016 | 0.144 |
| T2D | AD | -0.089 | 0.0327 | 0.0066 | 0.02376 | 0.2376 |
| STROKE | AD | -0.1396 | 0.0529 | 0.0083 | 0.0271636 | 0.2988 |
| STROKE | CKD | 0.1544 | 0.0637 | 0.0153 | 0.0459 | 0.5508 |

### Ocular–systemic pairs

The direction field is `POSITIVE` or `NEGATIVE` only when the global rg is BH-FDR detectable; otherwise it is `NO_DETECTABLE_GLOBAL_RG`. A nonsignificant point estimate is not described as independent.

| Ocular | Systemic | rg | SE | BH FDR | Direction |
|---|---|---:|---:|---:|---|
| AMD | CAD | -0.007 | 0.0383 | 0.8975 | NO_DETECTABLE_GLOBAL_RG |
| AMD | STROKE | 0.0078 | 0.0542 | 0.8975 | NO_DETECTABLE_GLOBAL_RG |
| AMD | T2D | -0.064 | 0.0374 | 0.1488 | NO_DETECTABLE_GLOBAL_RG |
| AMD | CKD | 0.0141 | 0.0571 | 0.877418 | NO_DETECTABLE_GLOBAL_RG |
| AMD | AD | -0.0368 | 0.0565 | 0.676671 | NO_DETECTABLE_GLOBAL_RG |
| AMD | PD | -0.0114 | 0.0458 | 0.877418 | NO_DETECTABLE_GLOBAL_RG |
| POAG | CAD | -0.0159 | 0.051 | 0.87689 | NO_DETECTABLE_GLOBAL_RG |
| POAG | STROKE | 0.0059 | 0.0458 | 0.8975 | NO_DETECTABLE_GLOBAL_RG |
| POAG | T2D | -0.0265 | 0.0267 | 0.4821 | NO_DETECTABLE_GLOBAL_RG |
| POAG | CKD | -0.0922 | 0.052 | 0.13752 | NO_DETECTABLE_GLOBAL_RG |
| POAG | AD | 0.0408 | 0.0457 | 0.535968 | NO_DETECTABLE_GLOBAL_RG |
| POAG | PD | -0.0336 | 0.0413 | 0.575862 | NO_DETECTABLE_GLOBAL_RG |
| CAT | CAD | 0.2159 | 0.0325 | 2.7081e-10 | POSITIVE |
| CAT | STROKE | 0.1034 | 0.0531 | 0.109482 | NO_DETECTABLE_GLOBAL_RG |
| CAT | T2D | 0.1459 | 0.0318 | 2.77188e-05 | POSITIVE |
| CAT | CKD | 0.0292 | 0.0603 | 0.780083 | NO_DETECTABLE_GLOBAL_RG |
| CAT | AD | -0.0997 | 0.0538 | 0.120316 | NO_DETECTABLE_GLOBAL_RG |
| CAT | PD | -0.0314 | 0.0496 | 0.676671 | NO_DETECTABLE_GLOBAL_RG |

## Intercept matrices and LAVA sampling correlation

`LDSC_CROSS_TRAIT_INTERCEPT_MATRIX.tsv` is complete and symmetric. Its diagonal uses the corresponding validated univariate LDSC intercept; off-diagonal cells use the unconstrained cross-trait LDSC genetic-covariance intercept. The LAVA matrix is the standardized `cov2cor`-style matrix obtained from this intercept covariance matrix, with diagonal fixed by standardization to one.

- finite: `True`; symmetric: `True`; maximum diagonal error: `0`.
- off-diagonal range: `-0.0360491` to `0.0824373`.
- eigenvalue range: `0.90892` to `1.1452`.
- The official LAVA UKB EUR v1.1 binary reference remains unavailable locally; the matrix QC does not substitute for that reference.

## Cohort overlap and HDL-L N0 eligibility

The metadata audit records only documented cohort components and possible overlap warnings. Exact participant overlap was not inferred from consortium names or sample sizes. Because no pair has a verified exact N0 or verified zero overlap in the frozen metadata, no pair is eligible for a definitive HDL-L covariance interpretation in this phase.

## AMD global stability and AMD–POAG benchmark

AMD global rg was rerun against all eight non-AMD traits using the corrected primary input and, where available, `AMD_RSID_GRCh37` and `AMD_HM3_STRICT`. The full machine-readable table is `AMD_GLOBAL_RG_SENSITIVITY.tsv`; the focused benchmark is `AMD_POAG_BENCHMARK.md`. No negative direction was forced.

## Gates

| Gate | Verdict |
|---|---|
| `GLOBAL_ARCHITECTURE_PASS` | `PASS` |
| `LAVA_OVERLAP_MATRIX_PASS` | `PASS` |
| `HDLL_PAIRWISE_FEASIBILITY` | `INADEQUATE_FOR_MOST_PAIRS` |

## Next permitted step

Resolve access to the official LAVA UKB EUR v1.1 reference and separately resolve exact participant-overlap/N0 metadata. Do not start local bivariate analysis or interpret HDL-L pairwise covariance until those prerequisites and the amended local-h2 requirements are satisfied.

## Files

- Dataset freeze: `config/phase1A_dataset_freeze.tsv`.
- Global LDSC: `results/phase1A/GLOBAL_LDSC_RG.tsv` and `GLOBAL_LDSC_RG_FDR.tsv`.
- Ocular–systemic subset: `results/phase1A/OCULAR_SYSTEMIC_RG.tsv`.
- Intercept and standardized sampling-correlation matrices: `results/phase1A/LDSC_CROSS_TRAIT_INTERCEPT_MATRIX.tsv` and `LAVA_SAMPLING_CORRELATION_MATRIX.tsv`.
- Overlap and HDL-L eligibility audits: `results/phase1A/COHORT_OVERLAP_METADATA.tsv` and `HDLL_N0_ELIGIBILITY.tsv`.
- QC and figures: `results/phase1A/LDSC_PAIR_QC.tsv`, `figures/phase1A/Figure1A_global_rg_heatmap.pdf`, `figures/phase1A/Figure1B_ocular_systemic_rg_forest.pdf`.
