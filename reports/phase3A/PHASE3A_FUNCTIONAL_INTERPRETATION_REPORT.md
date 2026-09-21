# Phase 3A — Robust Locus-Level Functional Interpretation of Mixed-Direction Ocular–Systemic Genetic Sharing

Date: 2026-09-21 (Asia/Shanghai)

## Executive conclusion

Phase 3A supports a robust **directional local genetic-sharing architecture** but does not establish functional convergence between negative-direction and positive-direction loci. The prespecified final verdict is:

`DIRECTIONAL_ARCHITECTURE_WITHOUT_FUNCTIONAL_CONVERGENCE`

The next step is:

`RETAIN_STATISTICAL_ARCHITECTURE_ONLY`

The 58 frozen independent loci remain valid descriptive evidence within the Phase 1D regions. Their positional annotations, xQTL lookup records, and recurrent-region structure can guide later biological work. They must not be used to call variants or genes causal, to reopen public-LD fine-mapping, or to label negative-direction loci as confirmed antagonistic pleiotropy.

## 1. Frozen evidence and scope

The project-wide frozen evidence comprises 18 ocular–systemic pairs and 40,978 finite local tests. Phase 1D identified 46 atlas-wide FDR-significant SUPERGNOVA regions, including 19 positive-rho and 27 negative-rho regions. Seven pair-level global-null pairs showed mixed-direction local genetic sharing: AMD–CAD, AMD–STROKE, AMD–T2D, AMD–CKD, POAG–CKD, CAT–STROKE, and CAT–AD.

Phase 2A-R established 2,096 valid independent PLACO loci. Fifty-eight loci map to 41 supported Phase 1D pair-region entries: 24 of 27 negative-rho entries and 17 of 19 positive-rho entries. The lead independent locus was direction-concordant in all 41 supported entries; the prespecified 10,000-permutation audit gave empirical P=0.0001 and Z=4.35. The 58 loci comprise 32 negative-direction and 26 positive-direction loci.

The public-LD fine-mapping route was closed by Phase 2B-R. Pan-UKBB EUR rescued 13/14 trait fits, but only 7/14 were fully reliable and only 1/7 pair-loci had both traits reliable. The technical gate failed. Accordingly, Phase 3A did not run new fine-mapping, coloc, PLACO, SUPERGNOVA, MR, GCI, or causal mediation analyses.

## 2. Prespecified analysis

The primary evidence units were the 41 supported pair-region entries and the 58 independent loci. Correlated PLACO overlaps and raw significant PLACO records were not treated as independent loci.

Four nested robustness subsets were retained:

| Subset | Definition | Included loci |
|---|---|---:|
| PRIMARY | q_atlas<0.05 and independent support | 58 |
| STRICT_LOCAL | q_atlas<0.01 and independent support | 26 |
| STRICT_VARIANT | PLACO P<5e-8/18 and q_atlas<0.05 | 12 |
| DOUBLE_STRICT | q_atlas<0.01 and PLACO P<5e-8/18 | 10 |

Leave-one-chromosome-out (LOCO) and leave-one-pair-out (LOPO) audits were descriptive sensitivity analyses. They were not used to select a preferred subset.

## 3. Directional robustness

LOCO retained 47–57 loci and 23–27 supported regions across chromosome exclusions. LOPO retained 48–57 loci and 24–28 supported regions across pair exclusions. Every exclusion retained a direction-concordance fraction of 1.0. This indicates that the observed lead-locus directional pattern is not attributable to one chromosome or one pair under the frozen mapping.

The robustness result supports the statistical phrase “mixed-direction local genetic sharing.” It does not identify a molecular mechanism and does not convert negative-direction local sharing into confirmed antagonistic pleiotropy.

## 4. Variant and region annotation

### VEP lead-variant annotation

Lead variants were queried against the Ensembl GRCh37 REST VEP endpoint using release 116 and the annotation date 2026-09-21. The output contains 632 transcript, regulatory, or intergenic annotation rows for 56 unique lead variant IDs shared by 58 locus records. There were 608 transcript annotation rows, 10 regulatory annotation rows, and 14 intergenic annotation rows. Fifty-seven gene symbols occurred in the VEP output. The dominant consequence category was intron_variant; the annotation was used for positional interpretation only.

### Protein-coding region inventory

The 41 supported pair-region entries generated 927 inventory records. These comprised 2 Level-1 transcript-disrupting/coding or splice records, 808 Level-2 protein-coding gene-body overlap records, and 117 Level-3 nearest protein-coding-gene records. The inventory contained 530 unique gene symbols across the queried regions. The mapping is positional and descriptive; “mapped gene” does not mean “causal gene.”

### xQTL lookup

The GTEx v8 API was queried for the frozen lead variants and relevant tissues, including brain, arterial, heart, kidney, liver, adipose, blood, nerve, muscle, spleen, pancreas, small-intestine, and skin tissues. The output contains 1,159 raw returned lookup records involving 47 unique lead variants, 167 gene symbols, and 30 tissues. Four lead variants had no GTEx v8 internal variant record and five had no record in the filtered relevant-tissue scope. EyeGEx/retina was recorded as an explicit optional-source status and was not silently interpreted as a negative result.

These xQTL records are lookup support only. No colocalization, mediation, or causal gene inference was performed.

## 5. Functional background and pathway analysis

The primary background was constructed from the 2,096 valid independent PLACO lead loci. Ensembl VEP produced 6,838 locus–transcript rows representing 984 unique protein-coding genes; 705 of the 2,096 lead loci had no protein-coding VEP transcript and were retained as unmapped loci. The background was not replaced by all human genes.

Direction-specific gene sets were formed from Level-1 to Level-3 mappings. A one-locus-one-gene sensitivity set selected the highest-priority available mapping deterministically. The predefined interpretive domains were RETINAL/OCULAR, VASCULAR/ENDOTHELIAL, IMMUNE/INFLAMMATORY, LIPID/METABOLIC, ECM, NEURONAL/NEURODEGENERATIVE, RENAL, and OXIDATIVE/MITOCHONDRIAL.

Primary GO Biological Process enrichment used the PANTHER over-representation service with a custom reference list derived from the full 2,096-locus background. No reported direction-specific term passed FDR<0.05 in the one-locus-one-gene sensitivity. Negative-direction one-locus-one-gene analysis returned 6,562 GO-BP records, with the best reported FDR approximately 0.986; after excluding the three most recurrent regions, the best reported FDR was approximately 0.857. The positive-direction custom-background requests did not return an interpretable result field for the supplied symbol set, so they were recorded as an API/symbol-coverage limitation rather than as evidence of biological absence. Reactome results were retained as descriptive default-background lookups and yielded no FDR<0.05 result; they were not used as the primary enrichment conclusion.

Excluding the single most recurrent region or the three most recurrent regions did not reveal a direction-specific FDR-significant pathway. Thus, the available functional results do not establish convergent biology separating negative-direction and positive-direction loci.

## 6. Recurrent cross-pair loci

The 58 loci occupy 28 unique genomic intervals after coordinate collapsing. Nine intervals recur across more than one pair, with a maximum of three pairs per interval. The most recurrent interval is 12:110112932–113027651, which contains eight frozen loci across AMD–CAD, AMD–CKD, and AMD–STROKE. The next recurrent intervals include 16:74730819–75517115, 17:17299591–18446092, and 7:104158491–105425027, each represented across three pairs. These recurrence summaries support cross-pair prioritization, not shared-causal interpretation.

## 7. Claim boundary

The strongest supported claim is that global and local summaries can show mixed-direction local genetic sharing, with independent PLACO loci maintaining lead-locus directional concordance in the validated regions and under LOCO/LOPO audits. Functional annotation provides positional and lookup context but does not establish causal variants, causal genes, mechanisms, or clinical utility.

The phrase “negative-direction local sharing” is retained for rho<0 regions with opposite PLACO effect directions. The phrase “positive-direction local sharing” is retained for rho>0 regions with same-direction effects. “Mixed-direction local genetic sharing” refers to pair-level coexistence of both directions. “Antagonistic pleiotropy” is not used as a confirmed result.

## 8. Decision and next step

The Phase 3A decision is recorded in `results/phase3A/PHASE3A_DECISION.txt`:

`DIRECTIONAL_ARCHITECTURE_WITHOUT_FUNCTIONAL_CONVERGENCE`

The next step is `RETAIN_STATISTICAL_ARCHITECTURE_ONLY`. No Phase 2 downstream analysis is started automatically. LAVA remains pending local-method replication, HDL-L remains pending exact N0, and the public-LD fine-mapping stop remains in force.

## 9. Reproducibility and source endpoints

- Frozen locus source: `results/phase2AR/INDEPENDENT_LOCUS_LOCAL_REGION_MAP.tsv`
- Frozen region source: `results/phase2AR/REGION_LEVEL_DIRECTIONAL_CONVERGENCE.tsv`
- Ensembl GRCh37 REST: https://grch37.rest.ensembl.org/
- GTEx API specification: https://github.com/broadinstitute/gtex-public/blob/master/gtex_api/gtex_api.yml
- PANTHER over-representation service: https://pantherdb.org/
- Reactome Analysis Service: https://reactome.org/AnalysisService

