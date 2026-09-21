# Phase 1D 当前情况

更新时间：2026-09-19 04:58（Asia/Shanghai）

## 当前状态：已完成

- 18/18 个预设 ocular–systemic pair 均为 `COMPLETE_VALID`。
- 已冻结 `results/phase1D/RAW_RESULTS_MANIFEST.tsv`；manifest 创建后未修改 raw SUPERGNOVA 输出。
- 完成 40,978 个有限 regional P 值的 atlas-wide BH-FDR 分析，得到 46 个显著区域：19 个正向、27 个负向。
- Pair QC：18/18 `PASS`；sign-flip QC：56/56 `PASS`。
- Phase 1D 判定：`STRONG_GLOBAL_LOCAL_DISCORDANCE`。
- Hidden mixed local sharing：AMD-CAD、AMD-STROKE、AMD-T2D、AMD-CKD、POAG-CKD、CAT-STROKE、CAT-AD。
- 下一步：`WAIT_FOR_LAVA_REPLICATION`；LAVA 仍等待 UKB EUR v1.1 reference。
- HDL-L 保持 `SECONDARY_PENDING_EXACT_N0`，不作最终 pairwise covariance 解释。
- 未启动 Phase 2、PLACO+、variant-level scan、fine-mapping、SuSiE、coloc、functional enrichment、cell-type analysis、MR 或 GCI validation。

最终报告：`reports/phase1D/PHASE1D_LOCAL_ATLAS_REPORT.md`。

## 已完成

- 已建立并冻结 `config/phase1D_dataset_freeze.tsv`，与 Phase 1A 冻结清单逐字一致。
- 已预注册 `reports/phase1D/STATISTICAL_ANALYSIS_PLAN.md`。
- 9 个 trait 的 Phase 1D SUPERGNOVA 输入已生成，输入质控表为 `results/phase1D/SUPERGNOVA_INPUT_QC.tsv`。
- 已完成 18 对 Phase 1D 全量重跑原始输出：`results/phase1D/raw_supergnova/`。
- Phase 1D atlas、并行调度器、汇总、sign-flip QC、最终化脚本和三幅图的生成脚本位于 `scripts/phase1D/`。

## 运行状态

- `run_supergnova_parallel.py` 与 SUPERGNOVA pair 子进程已退出；没有需要继续运行或恢复的 pair。
- 使用官方 SUPERGNOVA 源码与 Phase 1C 固定 commit/reference，保留 `thread=4`、冻结数据集、统计分析计划、sample-size convention 与 sign convention。
- 日志目录：`logs/phase1D/`；原始结果目录：`results/phase1D/raw_supergnova/`。

## 已生成的 Phase 1D 输出

- 全 18 对 regional atlas、pair QC、sign-flip QC；
- atlas-wide BH-FDR、pair-wise FDR 和 Bonferroni 敏感性结果；
- global–local architecture classification、方向性计数、局部协方差 burden 及 Spearman bootstrap 汇总；
- Figure 2A–2C、`TOP_LOCAL_REGIONS.tsv`、`PHASE2_LOCUS_PRIORITY.tsv`；
- `PHASE1D_DECISION.txt`、最终 Phase 1D 报告和项目 `CURRENT_STATUS.md`。

Phase 1D 已完成。后续仅在获得可用 LAVA reference 并按预设计划进行 LAVA replication 后再考虑下一步；本轮不启动 Phase 2 因果分析。
