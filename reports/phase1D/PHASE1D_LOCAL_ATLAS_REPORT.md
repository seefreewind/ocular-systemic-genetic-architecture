# PHASE 1D LOCAL COVARIANCE ATLAS REPORT

## Frozen analysis plan and execution

The analysis used the frozen `config/phase1D_dataset_freeze.tsv`, the preregistered `reports/phase1D/STATISTICAL_ANALYSIS_PLAN.md`, the official EUR SUPERGNOVA reference and the Phase 1C frozen source commit `319e84e114a4f954005a4592756c56cfee083667`. All 18 prespecified ocular–systemic pairs were rerun from the frozen Phase 1D inputs. The primary estimand was local genetic covariance `rho`; local genetic correlation was retained as descriptive QC only.

## Execution and technical QC

- Pair completion: **18/18 COMPLETE_VALID**.
- Regional tests with finite P values: **40,978** across 18 pairs.
- Pair QC: **18/18 PASS**; no primary pair failed.
- Sign-flip QC: **56/56 PASS**; all 56 checks passed.
- Raw-output manifest: `results/phase1D/RAW_RESULTS_MANIFEST.tsv`; raw outputs were not modified after manifest creation.

## Atlas size and atlas-wide significance

The official partition contained 2,353 attempted regions per pair. Returned regions varied because the official implementation omits regions that do not pass its usable-SNP/numerical requirements. The complete atlas contained **40,978** finite regional P values. At the prespecified BH-FDR threshold `q_atlas < 0.05`, **46** regions were significant: **19 positive** and **27 negative** local covariance estimates. Pair-wise FDR and Bonferroni results are retained as secondary sensitivity analyses.

## Pair classifications

| pair | global_rg | global_class | local_positive_atlas_fdr_count | local_negative_atlas_fdr_count | local_class |
|---|---|---|---|---|---|
| AMD-CAD | -0.007 | NO_DETECTABLE_GLOBAL_RG | 2 | 2 | MIXED_DIRECTION |
| AMD-STROKE | 0.0078 | NO_DETECTABLE_GLOBAL_RG | 1 | 1 | MIXED_DIRECTION |
| AMD-T2D | -0.064 | NO_DETECTABLE_GLOBAL_RG | 3 | 2 | MIXED_DIRECTION |
| AMD-CKD | 0.0141 | NO_DETECTABLE_GLOBAL_RG | 1 | 1 | MIXED_DIRECTION |
| AMD-AD | -0.0368 | NO_DETECTABLE_GLOBAL_RG | 0 | 0 | NO_DETECTABLE_LOCAL_SIGNAL |
| AMD-PD | -0.0114 | NO_DETECTABLE_GLOBAL_RG | 0 | 2 | NEGATIVE_ONLY |
| POAG-CAD | -0.0159 | NO_DETECTABLE_GLOBAL_RG | 0 | 4 | NEGATIVE_ONLY |
| POAG-STROKE | 0.0059 | NO_DETECTABLE_GLOBAL_RG | 0 | 1 | NEGATIVE_ONLY |
| POAG-T2D | -0.0265 | NO_DETECTABLE_GLOBAL_RG | 0 | 3 | NEGATIVE_ONLY |
| POAG-CKD | -0.0922 | NO_DETECTABLE_GLOBAL_RG | 1 | 2 | MIXED_DIRECTION |
| POAG-AD | 0.0408 | NO_DETECTABLE_GLOBAL_RG | 0 | 0 | NO_DETECTABLE_LOCAL_SIGNAL |
| POAG-PD | -0.0336 | NO_DETECTABLE_GLOBAL_RG | 0 | 0 | NO_DETECTABLE_LOCAL_SIGNAL |
| CAT-CAD | 0.2159 | GLOBAL_POSITIVE | 4 | 0 | POSITIVE_ONLY |
| CAT-STROKE | 0.1034 | NO_DETECTABLE_GLOBAL_RG | 1 | 1 | MIXED_DIRECTION |
| CAT-T2D | 0.1459 | GLOBAL_POSITIVE | 4 | 3 | MIXED_DIRECTION |
| CAT-CKD | 0.0292 | NO_DETECTABLE_GLOBAL_RG | 0 | 1 | NEGATIVE_ONLY |
| CAT-AD | -0.0997 | NO_DETECTABLE_GLOBAL_RG | 2 | 3 | MIXED_DIRECTION |
| CAT-PD | -0.0314 | NO_DETECTABLE_GLOBAL_RG | 0 | 1 | NEGATIVE_ONLY |

`HIDDEN_MIXED_LOCAL_SHARING` was assigned only when the frozen Phase 1A global class was `NO_DETECTABLE_GLOBAL_RG` and the Phase 1D local class was `MIXED_DIRECTION`.

## Hidden mixed local sharing

- AMD-CAD
- AMD-STROKE
- AMD-T2D
- AMD-CKD
- POAG-CKD
- CAT-STROKE
- CAT-AD

## Global-positive controls

| pair | global_rg | global_class | local_positive_atlas_fdr_count | local_negative_atlas_fdr_count | local_class |
|---|---|---|---|---|---|
| CAT-CAD | 0.2159 | GLOBAL_POSITIVE | 4 | 0 | POSITIVE_ONLY |
| CAT-T2D | 0.1459 | GLOBAL_POSITIVE | 4 | 3 | MIXED_DIRECTION |

## AMD results

| pair | N_positive_atlas_FDR | N_negative_atlas_FDR | N_total_atlas_FDR | pair_fdr_significant_total | minimum_q_atlas |
|---|---|---|---|---|---|
| AMD-CAD | 2 | 2 | 4 | 4 | 3.944e-06 |
| AMD-STROKE | 1 | 1 | 2 | 2 | 0.004444 |
| AMD-T2D | 3 | 2 | 5 | 6 | 0.002827 |
| AMD-CKD | 1 | 1 | 2 | 3 | 0.01097 |
| AMD-AD | 0 | 0 | 0 | 0 | 0.1443 |
| AMD-PD | 0 | 2 | 2 | 2 | 0.02044 |

## POAG results

| pair | N_positive_atlas_FDR | N_negative_atlas_FDR | N_total_atlas_FDR | pair_fdr_significant_total | minimum_q_atlas |
|---|---|---|---|---|---|
| POAG-CAD | 0 | 4 | 4 | 4 | 2.445e-09 |
| POAG-STROKE | 0 | 1 | 1 | 1 | 0.001743 |
| POAG-T2D | 0 | 3 | 3 | 3 | 0.02055 |
| POAG-CKD | 1 | 2 | 3 | 3 | 0.03817 |
| POAG-AD | 0 | 0 | 0 | 0 | 0.09792 |
| POAG-PD | 0 | 0 | 0 | 0 | 0.3196 |

## CAT results

| pair | N_positive_atlas_FDR | N_negative_atlas_FDR | N_total_atlas_FDR | pair_fdr_significant_total | minimum_q_atlas |
|---|---|---|---|---|---|
| CAT-CAD | 4 | 0 | 4 | 9 | 0.000499 |
| CAT-STROKE | 1 | 1 | 2 | 3 | 0.01963 |
| CAT-T2D | 4 | 3 | 7 | 7 | 2.059e-05 |
| CAT-CKD | 0 | 1 | 1 | 1 | 0.002379 |
| CAT-AD | 2 | 3 | 5 | 5 | 0.0002752 |
| CAT-PD | 0 | 1 | 1 | 1 | 3.276e-06 |

## Covariance burden

Burden metrics are descriptive sums of `rho`; they are not a Genetic Cancellation Index and are not used as causal evidence.

| pair | POS_BURDEN | NEG_BURDEN | SIGNED_NET | ABS_TOTAL | BALANCE |
|---|---|---|---|---|---|
| AMD-CAD | 0.05135 | 0.051194 | 0.00015576 | 0.10254 | 0.0015189 |
| AMD-STROKE | 0.065192 | 0.062421 | 0.0027714 | 0.12761 | 0.021717 |
| AMD-T2D | 0.069759 | 0.084494 | -0.014735 | 0.15425 | -0.095522 |
| AMD-CKD | 0.082253 | 0.074903 | 0.0073499 | 0.15716 | 0.046768 |
| AMD-AD | 0.05691 | 0.06392 | -0.0070098 | 0.12083 | -0.058013 |
| AMD-PD | 0.081525 | 0.081549 | -2.402e-05 | 0.16307 | -0.0001473 |
| POAG-CAD | 0.037314 | 0.042131 | -0.0048162 | 0.079445 | -0.060623 |
| POAG-STROKE | 0.049545 | 0.046602 | 0.0029433 | 0.096146 | 0.030613 |
| POAG-T2D | 0.055719 | 0.057768 | -0.0020488 | 0.11349 | -0.018053 |
| POAG-CKD | 0.052443 | 0.066678 | -0.014234 | 0.11912 | -0.11949 |
| POAG-AD | 0.047585 | 0.041496 | 0.0060884 | 0.089081 | 0.068347 |
| POAG-PD | 0.059405 | 0.060623 | -0.0012175 | 0.12003 | -0.010144 |
| CAT-CAD | 0.027919 | 0.014322 | 0.013597 | 0.042242 | 0.32189 |
| CAT-STROKE | 0.027266 | 0.022812 | 0.004454 | 0.050078 | 0.08894 |
| CAT-T2D | 0.037226 | 0.02482 | 0.012406 | 0.062045 | 0.19996 |
| CAT-CKD | 0.031078 | 0.030009 | 0.0010696 | 0.061087 | 0.01751 |
| CAT-AD | 0.022105 | 0.025784 | -0.003679 | 0.04789 | -0.076822 |
| CAT-PD | 0.029652 | 0.034254 | -0.0046022 | 0.063906 | -0.072015 |

## Global–local summary

| local_measure | n_pairs | spearman_rho | spearman_p | bootstrap_95ci_low | bootstrap_95ci_high |
|---|---|---|---|---|---|
| signed_net | 18 | 0.91538 | 1.0113e-07 | 0.72944 | 0.97912 |
| positive_burden | 18 | -0.2549 | 0.30735 | -0.74167 | 0.32014 |
| negative_burden | 18 | -0.51909 | 0.027279 | -0.87989 | 0.057601 |
| absolute_total | 18 | -0.35397 | 0.14955 | -0.80418 | 0.22269 |
| balance_index | 18 | 0.95253 | 1.117e-09 | 0.82445 | 0.9937 |

The strongest descriptive association was between global `rg` and signed local net covariance (Spearman rho=0.915, bootstrap 95% CI 0.729–0.979), and between global `rg` and the signed balance index (rho=0.953, 95% CI 0.824–0.994). These are pair-level summaries across only 18 pairs and are not regression or causal analyses.

## Phase 2 priorities

| priority_rank | tier | pair | region_id | rho | q_atlas | direction |
|---|---|---|---|---|---|---|
| 1 | TIER1_HIDDEN_MIXED | AMD-CAD | 5:120963248-121481183 | -0.00035629 | 3.9437e-06 | NEGATIVE |
| 2 | TIER1_HIDDEN_MIXED | CAT-AD | 1:206706419-208158787 | -0.00028484 | 0.00027524 | NEGATIVE |
| 3 | TIER1_HIDDEN_MIXED | AMD-CAD | 12:110112932-113027651 | 0.0006149 | 0.0010234 | POSITIVE |
| 4 | TIER1_HIDDEN_MIXED | AMD-CAD | 7:98152191-100555361 | -0.00039314 | 0.0019213 | NEGATIVE |
| 5 | TIER1_HIDDEN_MIXED | AMD-T2D | 4:103388441-104802530 | 0.0004845 | 0.0028273 | POSITIVE |
| 6 | TIER1_HIDDEN_MIXED | AMD-CAD | 16:74730819-75517115 | 0.00065262 | 0.0044438 | POSITIVE |
| 7 | TIER1_HIDDEN_MIXED | AMD-STROKE | 17:3426133-4321755 | -0.00051527 | 0.0044438 | NEGATIVE |
| 8 | TIER1_HIDDEN_MIXED | CAT-AD | 11:59620206-61870732 | -0.00026002 | 0.0062178 | NEGATIVE |
| 9 | TIER1_HIDDEN_MIXED | AMD-CKD | 7:104158491-105425027 | -0.00048831 | 0.010972 | NEGATIVE |
| 10 | TIER1_HIDDEN_MIXED | AMD-CKD | 12:110112932-113027651 | 0.00046886 | 0.019634 | POSITIVE |
| 11 | TIER1_HIDDEN_MIXED | CAT-STROKE | 17:17299591-18446092 | 0.00017872 | 0.019634 | POSITIVE |
| 12 | TIER1_HIDDEN_MIXED | AMD-STROKE | 12:110112932-113027651 | 0.00065601 | 0.019634 | POSITIVE |
| 13 | TIER1_HIDDEN_MIXED | CAT-STROKE | 10:12799286-13718384 | -0.00024438 | 0.020554 | NEGATIVE |
| 14 | TIER1_HIDDEN_MIXED | AMD-T2D | 17:3426133-4321755 | -0.00068679 | 0.023307 | NEGATIVE |
| 15 | TIER1_HIDDEN_MIXED | AMD-T2D | 5:73508509-75240469 | 0.00059996 | 0.02643 | POSITIVE |
| 16 | TIER1_HIDDEN_MIXED | AMD-T2D | 7:104158491-105425027 | 0.00059621 | 0.026812 | POSITIVE |
| 17 | TIER1_HIDDEN_MIXED | AMD-T2D | 2:43309247-44048346 | -0.00067835 | 0.031473 | NEGATIVE |
| 18 | TIER1_HIDDEN_MIXED | POAG-CKD | 7:77102731-77953248 | -0.00029619 | 0.038168 | NEGATIVE |
| 19 | TIER1_HIDDEN_MIXED | CAT-AD | 7:98152191-100555361 | -0.00014233 | 0.039202 | NEGATIVE |
| 20 | TIER1_HIDDEN_MIXED | POAG-CKD | 4:142808503-145021868 | -0.00032524 | 0.047253 | NEGATIVE |
| 21 | TIER1_HIDDEN_MIXED | POAG-CKD | 1:154547292-156334116 | 0.00033024 | 0.047253 | POSITIVE |
| 22 | TIER1_HIDDEN_MIXED | CAT-AD | 2:26894589-28597305 | 0.00011555 | 0.047253 | POSITIVE |
| 23 | TIER1_HIDDEN_MIXED | CAT-AD | 2:202830132-204814896 | 0.00012766 | 0.047253 | POSITIVE |
| 24 | TIER2_GLOBAL_LOCAL_CONCORDANT | CAT-CAD | 11:110702636-111981853 | 0.00013289 | 0.00049896 | POSITIVE |
| 25 | TIER2_GLOBAL_LOCAL_CONCORDANT | CAT-CAD | 1:175385076-176416712 | 0.00010414 | 0.011008 | POSITIVE |
| 26 | TIER2_GLOBAL_LOCAL_CONCORDANT | CAT-CAD | 16:74730819-75517115 | 0.00011296 | 0.02371 | POSITIVE |
| 27 | TIER2_GLOBAL_LOCAL_CONCORDANT | CAT-CAD | 7:98152191-100555361 | 0.00010821 | 0.031473 | POSITIVE |
| 28 | TIER3_OTHER | POAG-CAD | 17:1481897-2613925 | -0.00060075 | 2.4452e-09 | NEGATIVE |
| 29 | TIER3_OTHER | CAT-PD | 17:14860784-16319512 | -0.00041032 | 3.2758e-06 | NEGATIVE |
| 30 | TIER3_OTHER | CAT-T2D | 4:5971634-7135912 | 0.00039729 | 2.0594e-05 | POSITIVE |

The complete priority table is `results/phase1D/PHASE2_LOCUS_PRIORITY.tsv`. Tier 1 contains atlas-significant regions from hidden mixed local-sharing pairs; Tier 2 contains regions from globally detectable and directionally concordant pairs; Tier 3 contains all other atlas-significant regions. No causal gene or causal variant was assigned.

## Limitations and method status

SUPERGNOVA regional output is based on the official EUR reference and its usable-SNP/numerical filters, so returned-region counts are below the 2,353 attempted regions for most pairs. Local `corr` is frequently non-finite when one local heritability estimate is negative; `rho` remained the primary estimand. The atlas is a public-summary-statistics analysis and does not establish a shared causal variant, antagonistic pleiotropy, or clinical utility. LAVA remains `WAIT_FOR_LAVA_REFERENCE` because the required UKB EUR v1.1 reference was unavailable. HDL-L remains secondary and pending exact `N0`; its prior technical limitations remain unresolved.

## Phase 1D verdict

**STRONG_GLOBAL_LOCAL_DISCORDANCE**. Seven global-null pairs showed mixed-direction local covariance with atlas-wide FDR evidence. The prespecified next step is **WAIT_FOR_LAVA_REPLICATION**. Phase 2 causal-locus analysis was not started.
