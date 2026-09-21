# Phase 1D 当前情况

更新时间：2026-09-19 05:05（Asia/Shanghai）

## 总体状态

Phase 1D 已完成，当前不再有运行中的 SUPERGNOVA 任务。

- 预设 ocular–systemic pairs：18/18 `COMPLETE_VALID`
- 有限 regional P 值：40,978
- Atlas-wide BH-FDR 显著区域：46
  - 正向 local covariance：19
  - 负向 local covariance：27
- Pair QC：18/18 `PASS`
- Sign-flip QC：56/56 `PASS`
- Raw-output manifest：已创建并核对，18/18 SHA256 一致
- Phase 1D 自动化：已暂停，避免重复运行

## 主要判定

Phase 1D 判定为：`STRONG_GLOBAL_LOCAL_DISCORDANCE`。

7 个 Phase 1A global-null pairs 显示 mixed-direction local covariance，并达到 atlas-wide FDR 证据：

- AMD-CAD
- AMD-STROKE
- AMD-T2D
- AMD-CKD
- POAG-CKD
- CAT-STROKE
- CAT-AD

这些结果提示部分 pair 可能存在被 global summary 遗漏的局部方向异质性。该结果属于 summary-statistics 层面的局部协方差证据，不能单独证明共享因果变异、拮抗性多效性或临床效用。

## Global-positive controls

- CAT-CAD：global `rg = 0.2159`；4 个 positive atlas-FDR regions，0 个 negative regions；`POSITIVE_ONLY`
- CAT-T2D：global `rg = 0.1459`；4 个 positive、3 个 negative atlas-FDR regions；`MIXED_DIRECTION`

## 方向性和全局–局部汇总

描述性 pair-level Spearman 分析结果：

- global `rg` vs signed local net covariance：rho = 0.9154，P = 1.01 × 10⁻⁷，bootstrap 95% CI 0.7294–0.9791
- global `rg` vs signed balance index：rho = 0.9525，P = 1.12 × 10⁻⁹，bootstrap 95% CI 0.8244–0.9937
- global `rg` vs negative burden：rho = −0.5191，P = 0.0273；bootstrap 95% CI −0.8799–0.0576

这些相关性基于 18 个 pair 的描述性汇总，不是回归分析，也不构成因果证据。

## 技术和方法状态

- 使用冻结的 `config/phase1D_dataset_freeze.tsv`。
- 使用预注册分析计划 `reports/phase1D/STATISTICAL_ANALYSIS_PLAN.md`。
- 使用 Phase 1C 冻结的官方 SUPERGNOVA 源码/reference，commit：`319e84e114a4f954005a4592756c56cfee083667`。
- 保留 `thread=4`、sample-size convention 和 sign convention。
- `rho` 是主要局部协方差估计；local `corr` 仅作为描述性 QC，因为部分区域的 local heritability 为负时 `corr` 不可计算。
- LAVA：`WAIT_FOR_LAVA_REFERENCE`；所需 UKB EUR v1.1 reference 仍不可用。
- HDL-L：`SECONDARY_PENDING_EXACT_N0`；不对 pairwise covariance 做最终生物学解释。

## 当前未执行的分析

以下分析均未启动：

- Phase 2 causal-locus analysis
- PLACO+
- variant-level scan
- fine-mapping、SuSiE、coloc
- functional enrichment、cell-type analysis
- MR
- GCI validation

## 已生成文件

- 最终报告：`reports/phase1D/PHASE1D_LOCAL_ATLAS_REPORT.md`
- Phase 1D 决策：`results/phase1D/PHASE1D_DECISION.txt`
- Raw manifest：`results/phase1D/RAW_RESULTS_MANIFEST.tsv`
- Regional atlas：`results/phase1D/SUPERGNOVA_OCULAR_SYSTEMIC_ATLAS.tsv`
- Pair QC：`results/phase1D/SUPERGNOVA_PAIR_QC.tsv`
- Sign-flip QC：`results/phase1D/SIGN_FLIP_ATLAS_QC.tsv`
- Directional summary：`results/phase1D/DIRECTIONAL_LOCAL_SIGNAL_SUMMARY.tsv`
- Global–local classification：`results/phase1D/GLOBAL_LOCAL_ARCHITECTURE_CLASSIFICATION.tsv`
- Local covariance burden：`results/phase1D/LOCAL_COVARIANCE_BURDEN.tsv`
- Spearman bootstrap summary：`results/phase1D/GLOBAL_LOCAL_SUMMARY_CORRELATIONS.tsv`
- Top local regions：`results/phase1D/TOP_LOCAL_REGIONS.tsv`
- Phase 2 priority table：`results/phase1D/PHASE2_LOCUS_PRIORITY.tsv`
- Figure 2A–2C：`figures/phase1D/`
- 项目状态：`CURRENT_STATUS.md`

## 下一步

下一步是等待可用的官方 LAVA UKB EUR v1.1 reference，并按预设计划进行 LAVA replication。获得该 reference 前，不启动 Phase 2 或其他未解锁的下游分析。
