# Pan-UKBB EUR LD reference audit (Phase 2B-R)

Generated: 2026-09-21T16:42:07+08:00

## Scope and frozen constraints

This rescue audit re-analyzed exactly P2BP02–P2BP08 from the frozen Phase 2B-P pilot. P2BP01 was not analyzed or replaced because the original pilot had 198 included variants and failed the prespecified minimum-200-variant LD coverage gate. No additional loci, PLACO+, CAT, LAVA, enrichment, MR, GCI, or Phase 2B expansion were run.

## Official reference and access

The official resource is the Pan-UKBB EUR LD release documented at [Pan-UKBB LD documentation](https://pan.ukbb.broadinstitute.org/docs/ld/index.html), using the Hail variant index `s3://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht` and BlockMatrix `s3://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.bm`. The variant-index global metadata supplied the exact EUR `n_samples=420542`; this value was used as B and was not hardcoded in the model.

The resource uses GRCh37 coordinates, dosage-based covariate-adjusted Pearson LD, and the published INFO>0.8/MAC>20 variant filters within a 10-Mb LD radius. Access used Hail 0.2.135, Python 3.11, OpenJDK 17, Spark/Hadoop S3A with anonymous read-only access. Only locus-specific variant-index rows and the required BlockMatrix blocks were read; the full 43.3-TB matrix was not downloaded.

## Access and regional matching

| Metric | Value |
|---|---|
| n_samples | 420542 |
| Eligible loci | 7/7 |
| Regional coverage gates passed | 7/7 |
| Complete Pan LD matrices | 7/7 |
| Allele-alignment audit records | 3303 |
| P2BP01 status | NOT_ANALYZED_PHASE2BR |

The detailed access and variant/allele audits are retained in `results/phase2BR/PANUKBB_ACCESS_AUDIT.tsv`, `results/phase2BR/PANUKBB_VARIANT_MATCH_QC.tsv`, and `results/phase2BR/PANUKBB_ALLELE_ALIGNMENT.tsv`.
