# PLACO+ software audit

Audit date: 2026-09-19 (Asia/Shanghai)

## Official source

- Repository: https://github.com/RayDebashree/PLACO
- Frozen source file: `data/phase2A/PLACO_source/PLACO_v0.2.0.R`
- Git commit: `3ba3cae1d323ad117fb4540e620bcefa79f70663`
- Source SHA256: `fb684a8ed88f27dd138f5c2e8613b092904e84f364030d036058f17db95b7124`
- Official function used: `placo.plus()`
- Official parameter helpers used: `var.placo()` and `cor.pearson()`
- Statistical source functions: not modified

## Runtime

- R: `R version 4.4.3 (2025-02-28)`
- `data.table`: `1.18.4`
- `MASS`: `7.3-65`
- `mvtnorm`: `1.4-2`
- `stats`: `4.4.3`

The official PLACO source is sourced directly from the frozen local copy. No replacement implementation is used for the primary PLACO+ P values.

## Reference and input status

- Summary-statistics genome build: GRCh37.
- Primary ancestry: European or predominantly European; CAT is Finnish-European.
- Available EUR LD panel for discovery organization and sensitivity: `data/phase1C/SUPERGNOVA_reference/bfiles/eur_chr{1..22}_SNPmaf5`.
- CAT-containing clumping results will carry `FINNISH_LD_SENSITIVITY_REQUIRED`.
- LAVA UKB EUR v1.1 reference: unavailable at the Phase 2A start check; status remains `PENDING_LOCAL_METHOD_REPLICATION`.

## Reproducibility files

The Phase 2A preparation and execution scripts are under `scripts/phase2A/`. Logs are written to `logs/phase2A/`, and all result tables are written under `results/phase2A/`.
