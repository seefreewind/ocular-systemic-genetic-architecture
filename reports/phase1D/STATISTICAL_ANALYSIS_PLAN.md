# Phase 1D statistical analysis plan

## Scope

Phase 1D is the prespecified discovery atlas of local genetic covariance between three ocular traits (AMD, POAG, CAT) and six systemic traits (CAD, STROKE, T2D, CKD, AD, PD), comprising 18 ocular–systemic pairs. The seven Phase 1C pairs will be rerun from the frozen inputs so that the final atlas is homogeneous.

The primary estimand is SUPERGNOVA local genetic covariance, `rho`. Local correlation is descriptive because Phase 1C showed frequent negative local h² estimates and non-finite local correlation values while rho remained estimable. Positive and negative rho are reported as `POSITIVE_LOCAL_COVARIANCE` and `NEGATIVE_LOCAL_COVARIANCE`. They are not called concordant or antagonistic pleiotropy in this phase.

## Frozen software, reference, and inputs

- SUPERGNOVA source commit: `319e84e114a4f954005a4592756c56cfee083667`.
- Official SUPERGNOVA EUR bfiles and chromosome partitions from the Phase 1C reference audit.
- Exact Phase 1A/Phase 1C frozen GWAS inputs listed in `config/phase1D_dataset_freeze.tsv`.
- The official Phase 1C runtime and compatibility wrapper are retained; source code is not updated during this phase.

## Multiplicity hierarchy

The primary testing family is every local covariance test across all 18 ocular–systemic pairs. After all 18 runs are complete, BH-FDR is applied once across the complete pair-by-region family. This adjusted value is named `q_atlas`, and `q_atlas < 0.05` is the primary discovery criterion.

BH-FDR within each pair is also calculated as `q_pair` and is a secondary descriptive sensitivity analysis. Bonferroni correction across the complete atlas is reported as a conservative sensitivity analysis. The multiplicity hierarchy is frozen before the remaining runs and will not be changed after results are observed.

## Required regional fields

Each returned region will retain trait labels, pair, chromosome, start, end, region identifier, rho, variance, SE, Z, P, local h² for both traits, local corr, and the returned SNP count `m`, together with input and numerical status fields.

## QC and sign-flip validation

For every pair, the analysis will report attempted and returned regions, finite rho, variance, P and local corr, negative local h² counts, failed rows, and low-SNP rows. Rho, variance and P must be finite and P must lie in [0,1] for the numerical gate to pass.

The sign-flip spot check will select two stable regions per pair using a fixed seed and 20 additional high-signal regions after the atlas run. The systemic-trait Z will be sign-flipped. A check passes only when both local h² estimates and P are invariant, rho reverses sign, and absolute rho is invariant. Every discrepancy will be investigated.

## Prespecified summaries

For each pair, atlas-wide and pair-wise positive/negative significant-region counts will be reported. Pair-level local classes are `MIXED_DIRECTION`, `POSITIVE_ONLY`, `NEGATIVE_ONLY`, or `NO_DETECTABLE_LOCAL_SIGNAL`, based only on atlas-wide FDR results.

Global LDSC classes use Phase 1A BH-FDR only: `GLOBAL_POSITIVE`, `GLOBAL_NEGATIVE`, or `NO_DETECTABLE_GLOBAL_RG`. `HIDDEN_MIXED_LOCAL_SHARING` is assigned only to global-null pairs with atlas-wide `MIXED_DIRECTION` local results. Global–local concordance and discordance are reported descriptively.

Covariance burden summaries (positive burden, negative burden, signed net, and absolute total) are descriptive and are not named a cancellation index. Their association with global rg is summarized with Spearman correlation and a bootstrap confidence interval across the 18 pairs; no regression is performed.

## Exclusions and stop rule

GCI, PLACO+, fine-mapping, coloc, functional enrichment, cell-type analysis, MR, causal inference, and any expansion beyond the prespecified atlas are excluded. The phase stops after the 18-pair atlas, QC, figures, priority set, and decision are complete. LAVA remains replication-pending because the official UKB EUR reference is unavailable. HDL-L remains non-primary because exact N0 is unresolved.
