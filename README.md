# Ocular–systemic genetic architecture atlas

This public release contains the reproducibility code, prespecified configuration, derived summary tables, audit reports and figures for a cross-trait genetic-architecture atlas spanning three ocular phenotypes and six systemic phenotypes.

## Authors and correspondence

- Da Lin¹ — Department of Ophthalmology, The Second Affiliated Hospital of Wenzhou Medical University
- Ying Chen² — Wenzhou Medical University
- Yue Liu² — Wenzhou Medical University
- Yu Zhang¹ — Department of Ophthalmology, The Second Affiliated Hospital of Wenzhou Medical University

¹ Department of Ophthalmology, The Second Affiliated Hospital of Wenzhou Medical University, No. 109 Xueyuan West Road, Lucheng District, Wenzhou, Zhejiang Province, China<br>
² Wenzhou Medical University, Wenzhou, Zhejiang Province, China

Corresponding author: Yu Zhang — zhangyu1@wzhealth.com<br>
ORCID: [0000-0001-8579-3692](https://orcid.org/0000-0001-8579-3692)

Funding: none.<br>
Competing interests: none declared.

The preferred citation metadata are provided in [`CITATION.cff`](CITATION.cff). The Zenodo archive is versioned separately from the GitHub repository. This release is archived as [10.5281/zenodo.22874563](https://doi.org/10.5281/zenodo.22874563); the concept DOI for all versions is [10.5281/zenodo.22874303](https://doi.org/10.5281/zenodo.22874303).

The release is organized around the completed analysis layers:

- genome-wide LDSC pairwise genetic correlation;
- atlas-wide SUPERGNOVA local covariance and directional classification;
- LD-clumped PLACO+ variant-level follow-up;
- locus consolidation and directional permutation checks;
- limited public-LD fine-mapping readiness assessment;
- functional annotation and prespecified robustness summaries.

## What is included

- `scripts/` — portable analysis scripts for the released analysis layers;
- `config/` — dataset-freeze files, thresholds and locus-selection manifests;
- `results/` — derived tables, QC outputs and decision records;
- `reports/` — analysis and audit reports;
- `figures/` — final derived figures.

## What is intentionally excluded

Raw GWAS summary statistics, consortium-restricted files, LD reference panels, large genome-wide intermediate tables, local machine environments, credentials, logs and submission-specific materials are not redistributed. The full PLACO+ genome-wide table is therefore excluded because of size and source-data considerations. Source access URLs, checksums and dataset identifiers are retained in the configuration manifests; users must obtain each source dataset under its original access and citation terms.

The released tables are derived statistical outputs. They do not establish causal variants, colocalization, biological antagonism or participant-level clinical comorbidity.

## Reproduction notes

Run scripts from the repository root and provide the required source data and reference panels through portable relative paths. External tools such as LDSC, R, PLACO, PLINK, Hail and the public Pan-UKBB resources must be installed separately. Environment-specific executable paths were removed from the released scripts; use `PATH` or explicit environment variables such as `PHASE2A_RSCRIPT`.

The released numerical outputs are frozen results from the completed analysis. Regenerating them requires the exact source datasets, reference resources, software versions and parameter settings documented in `config/` and `reports/`.

## Scope boundary

This repository is a reproducibility release for the genetic-architecture analysis. It does not include journal-specific manuscript files or an assertion of independent replication. PLACO+ and SUPERGNOVA use the same underlying GWAS summary statistics, so their agreement is interpreted as cross-scale internal consistency.
