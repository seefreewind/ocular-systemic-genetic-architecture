# Protocol Amendment 002 — Phase 2A PLACO+ before LAVA replication

Date: 2026-09-19 (Asia/Shanghai)

## Scientific basis

The original sequence required LAVA replication before Phase 2. The official LAVA documentation still identifies the UK Biobank EUR v1.1 binary LD reference as the recommended European-ancestry reference, with download links hosted through the SURF distribution endpoint. The required reference has remained inaccessible from this project after repeated validated access attempts, including the official endpoint check at the start of Phase 2A.

Phase 1C and Phase 1D SUPERGNOVA analyses used the frozen official source/reference, passed numerical and sign-flip validation, and produced atlas-wide FDR-controlled local covariance signals across all 18 prespecified ocular–systemic pairs. This provides a justified basis for a distinct variant-level statistical pleiotropy audit while the local-method replication remains unavailable.

## Amendment

Phase 2A genome-wide PLACO+ testing may proceed before LAVA replication. The analysis will use the unchanged Phase 1A/Phase 1D GWAS inputs and the official PLACO v0.2.0 implementation. PLACO+ results will be interpreted as variant-level statistical pleiotropy evidence only.

LAVA remains mandatory as a later replication step and is recorded as:

`PENDING_LOCAL_METHOD_REPLICATION`

Phase 2A findings must not be described as replicated local architecture until LAVA replication is completed with an official or maintainer-approved UKB EUR v1.1 reference.

## Prespecified safeguards

- Run PLACO+ genome-wide for all 18 primary pairs; do not restrict the calculation to the 46 Phase 1D regions.
- Keep the frozen Phase 1A/Phase 1D datasets, genome build, ancestry labels, sample-size convention, and effect-allele convention unchanged.
- Use signed `Z = BETA / SE`; do not reconstruct unsigned Z scores from P values.
- Harmonize alleles using SNP identity plus CHR/BP and allele checks; remove unresolved strand-ambiguous variants.
- Apply the primary common-variant rule `MAF >= 0.01` when reliable trait or compatible reference-panel frequency is available. Record unavailable-frequency exclusions explicitly.
- Estimate PLACO+ null parameters at `p.threshold = 1e-4` using the official `var.placo()` and `cor.pearson()` functions.
- Exclude variants with `Z1^2 > 80` or `Z2^2 > 80` from the primary PLACO+ calculation and apply the prespecified dual-GWAS-significance screen separately.
- Report pair-wise genome-wide P < 5e-8, atlas-stringent P < 5e-8/18, and atlas BH-FDR < 0.05 separately.
- Flag all CAT-containing LD-clumping results as `FINNISH_LD_SENSITIVITY_REQUIRED`.
- Stop after Phase 2A. Do not start fine-mapping, coloc, enrichment, cell-type analysis, MR, GCI, or other Phase 2B analyses automatically.

## Amendment status

`AMENDMENT_STATUS = SCIENTIFICALLY_JUSTIFIED`

This amendment changes the order of the statistical audit only. It does not lower significance thresholds, change the primary pair family, replace frozen datasets, or convert statistical pleiotropy into a causal or mechanistic claim.
