# CIMA peak–gene 因果连接：可视化（Visualization）

> 目标：用 Claude_project 新跑出的 MR 和验证结果，复现论文画图目录
> `/scratch/user/tfcui2025/CIMA/Results/RunMR/Materials_runMR/20261006_paper`（下称 **paper 目录**，只读）中的图。
> 画图思路、版式、配色都沿用 paper 目录；只换输入数据和输出位置。
> 运行环境、文件安全规则见 `Claude_project/CLAUDE.md`；分析流程见 `Scripts/Pipeline.md`。

## 1. 文件位置

| 内容 | 位置 |
|---|---|
| 画图代码 | `Scripts/Visualization/Plot<N>_*.R`（重的另有 `*_slurm.sh`） |
| 共用设置 | `Scripts/Visualization/00_common.R`（路径、细胞类型、方法名/配色、主题、`read_mr()`） |
| Example pairs | `Scripts/Visualization/01_example_pairs.R`（24×3 个 pair、4 个类别、`prep_examples()`、`check_examples()`；Fig1 检查、Fig2 画图共用这一份） |
| SLURM 日志 | `Scripts/Visualization/logs/` |
| 图 / 表 | `Results/Figures/` |
| 临时文件（PLINK LD 等） | `tmp/` |

## 2. 共同约定

- **输入**：
  - MR：`Results/MR/<CT>_AllPairs_results.tsv.gz`（已删 NA、`ivw_sigma2 <= 10`、gene 内 BH）。画图时**不再**做 sigma2 过滤。
  - 验证：`Results/Validation/Valid<k>_<名>/`（见 `Pipeline.md` 第 3 节）。
- **细胞类型**：`Bn_TCL1A`、`CD4_Tn_CCR7`、`CD8_Tn_CCR7`。
- **方法显示名、顺序和配色**（所有图统一，与 paper 目录相同）：

  | 显示名 | 数据中的 key | 颜色 |
  |---|---|---|
  | Causal SMR | `smr` | `#F0A14A` |
  | Causal IVW | `ivw` | `#E64B35` |
  | Causal Egger | `egger` | `#3C5488` |
  | Causal Joint | `cct` | `#8491B4` |
  | Marginal Multiomics | `wnn` | `#4DBBD5` |
  | Marginal Pseudobulk | `pseudo` | `#00A087` |

  另有 `COL_CONC = #D5DDE5`、`COL_DISC = #7E6148`、`COL_INK = #2F2F2F`。
- **FDR 阈值**：沿用 paper 目录——Fig1/Fig2 用 0.1；Fig3 依赖 Valid1（0.1）；Fig5、Fig7 依赖 Valid2/Valid5（0.05）。
- **运行**：轻的脚本在登录节点 `bash -lc 'module load ...; Rscript Scripts/Visualization/PlotN_*.R'`；Fig2 用 `sbatch`。

## 3. 各图

| 图 | 本项目脚本 | paper 目录原脚本 | 状态 |
|---|---|---|---|
| Fig1 MR 总体结果 | `Plot1_MRresults.R` | `Plot1_MRresults.sh` | ✅ |
| Fig2 Example pairs（v1 + v2） | `Plot2_examples_v1.R` / `_v2.R` | `Source_visual/Plot2_scripts.R` / `_v2.R` | ✅ |
| Fig3 超几何富集 | `Plot3_Hyper.R` | `Plot3_val_Hyper.sh` | ✅ |
| Fig4 ROC AUC | `Plot4_AUC.R` | `Plot4_val_AUC.sh` | ✅ |
| Fig5 S-LDSC | `Plot5_SLDSC.R` | `Plot5_val_SLDSC.sh` Part 2 | ✅（只做第 1 页） |
| Fig6 Motif | `Plot6_Motif.R` | `Plot6_val_Motif.sh` | ✅（只出 `Motif_all.tsv` + 复制 xlsx） |
| Fig7 CIMA GRN | `Plot7_GRN.R` | `Plot7_val_GRN.sh` 第 6 节 | ✅ |
| Fig8 位点 link 图 | `Plot8_Link.R` | `Plot8_val_Link.sh` | ✅（只做 GWAS 模式；依赖 Valid6_GWAS） |

### 3.1 Fig1 — MR 总体结果 ✅
- **输入**：`Results/MR/<CT>_AllPairs_results.tsv.gz`
- **脚本**：`Scripts/Visualization/Plot1_MRresults.R`
- **方法**：FDR < 0.1 判显著；每个方法的显著 pair 数柱图；每个 CT 一张 UpSet（ComplexHeatmap，前 25 个组合）；六方法全显著 pair 中 Causal SMR/IVW/Egger 与 Marginal 的效应方向一致性。
- **输出**：`Results/Figures/Fig1_MR_overview.pdf`，5 页，11×5.2 in，cairo_pdf。
- **Example pairs 检查**：脚本开头用 `check_examples()` 在新表和旧表（`MR_results_formal/*.rds`，只读）下分别检查 72 个 example pair 是否满足自己的类别，结果写 `tmp/example_pairs_check.tsv`。
  2026-10-09 结果：72 个全部在新表中、新旧 FDR 下都满足类别 → **沿用 paper 的 example pairs，不重挑**。

### 3.2 Fig2 — Example pairs ✅
- **输入**：MR 结果；eQTL、caQTL（全量 + LD clump）；PLINK + `data/ref_chr/hg38`；pseudobulk；WNN（`1_pbmc_multiomics.reprocessed.rda`、`ATAC_RNA_WNN.rda`）。
- **方法**：对每个 example pair 重跑 QTL 对齐 → 去重 → PLINK LD → GLS IVW/Egger + SMR，拿到 SNP 层数据（p/FDR/beta 一律取 MR 结果表）。
  - v1：每页 4 个 panel（caQTL/eQTL Manhattan、MR 散点、WNN 相关、pseudobulk 相关），21×7 in。
  - v2：左 2/3 为上述 2×2 面板，右 1/3 为位点 link 图（只突出当前 pair）。
- **脚本**：`Scripts/Visualization/Plot2_examples_v1.R`、`Plot2_examples_v2.R`（由 paper 的 `Plot2_scripts.R` / `_v2.R` 复制，作图代码不动），`Plot2_examples_slurm.sh`（array 0–2，每个 CT 一个 task，依次跑 v1 和 v2；4 核、96G、3 h；日志 `Scripts/Visualization/logs/`）。
  单独运行：`Rscript Scripts/Visualization/Plot2_examples_v1.R --celltype <CT|all> --n <N|0>`（`--n 0` = 全部 example，否则文件名带 `_topN`）。
- **与 paper 的差异**：
  - 读 `Results/MR/<CT>_AllPairs_results.tsv.gz`，不再做 sigma2 过滤；LD 参考用 `data/ref_chr/hg38`；example pairs 改为 source `01_example_pairs.R`。
  - 新 MR 表没有 `egger_int` 列（paper 的 Egger 线截距取自旧表）→ 改用重跑得到的 Egger 截距（`alpha`）；日志打印重跑 Egger 斜率与表中 `egger_beta` 的差，作为对账。
  - PLINK 临时文件、awk 过滤用的 key 文件放 `tmp/plink_fig2/ld_fig2_<pid>/`；SLURM 中 `TMPDIR` 设为 `tmp/plink_fig2/Rtmp_<job>_<task>/`，结束删除。
- **输出**：`Results/Figures/Fig2_examples_<CT>.pdf`（v1）和 `Fig2_examples_<CT>_v2.pdf`（v2），每个 pair 一页（同一 pair 属于多个类别时只画一页）。
- **2026-10-09 运行**：job 20080188，每个 CT 约 9 分钟、内存约 4.8G；作图 pair 数 Bn 22 / CD4 23 / CD8 24，全部 status = ok；重跑的 IVW 斜率、SMR lead SNP 与 MR 表完全一致，重跑 Egger 斜率与表中 `egger_beta` 差 ≤ 5e-15。

### 3.3 Fig3 — 超几何富集 ✅
- **输入**：`Results/Validation/Valid1_Hyper/<CT>_hyper_long.tsv`（method 为 `SMR/IVW/Egger/CCT/WNN/Pseudobulk`，database 为 `enhancerAtlas/FANTOM5/4DGenome/GTEx/HiChIP/PCHiC`，与旧表一致）。
- **脚本**：`Scripts/Visualization/Plot3_Hyper.R`
- **方法**：只画 `N > 0, K > 0, n_sig > 0, N > K` 的行；每个 gene 的超几何 p 值画直方图（bin 0.05），facet = 6 方法 × 6 金标准，按方法上色；每个格子把 N/K/n_sig/k 求和后算一个 global hypergeometric P 标在右上角。日志打印每格 gene 数和 global P。
- **输出**：`Results/Figures/Fig3_Hyper.pdf`，每个 CT 一页，13×10.5 in。

### 3.4 Fig4 — ROC AUC ✅
- **输入**：`Results/Validation/Valid4_AUC/<CT>_AUC_long.tsv.gz`（`gold_standard` 为 `enhancerAtlas/FANTOM5/...`，`method` 为小写 key）。
- **脚本**：`Scripts/Visualization/Plot4_AUC.R`
- **方法**：按 CT × 金标准 × 方法汇总 mean/median/sd；小提琴 + 箱线图，facet = 6 个金标准，0.5 虚线为随机水平。金标准显示名沿用 paper（EnhancerAtlas、FANTOM、4DGenome、GTEx cisQTL、HiChIP、PCHiC）。
- **输出**：`Results/Figures/Fig4_AUC.pdf`，3 页，13×8.5 in。

### 3.5 Fig5 — S-LDSC ✅
- **输入**：`Results/Validation/Valid2_SLDSC/SLDSC_all.tsv`（`cell_type`、`set`、9 个短名 trait；只有 ext500）。
- **脚本**：`Scripts/Visualization/Plot5_SLDSC.R`
- **方法**：enrichment ± s.e. 柱图（ggh4x `facet_grid2`，行 = CT（顺序 CD4 / Bn / CD8，沿用 paper），列 = 9 个 trait：AS, CD, IBD, RA, SLE, MS, T1D, Asthma, PBC），6 个方法 + All peaks（`pool`）。日志打印 enrichment 表（`enrich_p < 0.05` 标 *）。
- **输出**：`Results/Figures/Fig5_SLDSC.pdf`，**只做第 1 页**（paper 的第 2、3 页需要 MR-All 集合，本项目没有跑），18×7.5 in。

### 3.6 Fig6 — Motif ✅
- **输入**：`Results/Validation/Valid3_Motif/<CT>_motif_enrichment.tsv.gz`
- **脚本**：`Scripts/Visualization/Plot6_Motif.R`
- **方法**：只把 3 个 CT 的 motif 富集结果合并成 `Motif_all.tsv`（列与 paper 的 `Motif_details/Motif_all.tsv` 相同，method 名仍为 `IVW/MR_Egger/SMR/CCT/WNN/Pseudo`）；不画 logo、不做 16 个 motif 的表、不生成新的 xlsx。
  日志里另做两项对照（只打印）：与 paper 旧 `Motif_all.tsv` 逐行比较 p 值和显著数；xlsx 中 16 个目标 motif 在新旧结果下是否显著（p.adjust < 0.05 且 FE > 1）。
- **2026-10-09 对照结果**：6 个 CT×方法组合（如 Bn 的 IVW/SMR/CCT）与旧结果完全相同，其余前景 peak 数差 1–27 个，p 值随之变化（每组 2–12 个 motif 的显著性翻转）；原因应是旧 Motif 代码的前景用的是旧 FDR 且保留了 sigma2 为 NA 的 pair，前景变了，按 GC/可及性匹配的背景也跟着变（未逐一核实）。
  16 个目标 motif 中 14 个结论不变；**CD4 JUND（MA0492.1）在 IVW 下不再显著**（SMR 下仍显著）；**CD8 FOSL2::JUND（MA1145.1）在 MR_Egger 下不再显著**（新结果中 6 个方法都不显著）；Bn ELF2 在 MR_Egger 下新变为显著。复制过来的 xlsx 仍是旧数字，用时需注意这两行。
- **输出**：`Results/Figures/Motif_all.tsv`；另把 paper 目录的 `Materials/Fig6_Motif_table.xlsx` 原样复制到 `Results/Figures/`（该 xlsx 是本地 Claude Code 结合外部知识手工整理的，数字来自旧结果）。

### 3.7 Fig7 — CIMA GRN ✅
- **输入**：`Results/Validation/Valid5_GRN/GRN_sig_counts_by_method.tsv`（统计量在 Valid5 中已算好：FDR < 0.05，2×2 Fisher）。
- **脚本**：`Scripts/Visualization/Plot7_GRN.R`
- **方法**：上图 = 每个方法显著且与 GRN 重叠的 pair 数（行 = CT / CT_HQ，列 = 细胞类型）；下图 = OR（GRN vs 非 GRN），CT 与 CT_HQ 并排。
- **输出**：`Results/Figures/Fig7_GRN.pdf`，1 页，14×11 in。

### 3.8 Fig8 — 位点 link 图 ✅
- **输入**：`Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv`（替代 paper 的 `Valid8_matchGWAS/G1_Focus_MRonly_GWAS.tsv`）；MR 结果；QTL；PLINK + `data/ref_chr/hg38`；GWAS sumstats（9 个）；EnsDb.Hsapiens.v86。
- **只做 GWAS 模式**：一页 = 一个 link × 一个 trait；输出 `Results/Figures/Fig8_Link_GWAS.pdf`。不做 example 模式（`Fig8_Link.pdf`）。
- **脚本**：`Scripts/Visualization/Plot8_Link.R`（由 paper 的 `Plot8_val_Link.sh` 复制，作图代码不动），`Plot8_Link_slurm.sh`（4 核、96G、2 h；日志 `Scripts/Visualization/logs/`）。
- **与 paper 的差异**：GWAS pairs 读 Valid6 输出；`TRAIT_FILES` 改为项目的 9 个 GWAS（标签 AS/CD/IBD/RA/SLE/MS/T1D/Asthma/PBC）；MR 读新表、不做 sigma2 过滤；LD 参考 `data/ref_chr/hg38`；example pairs 改为 source `01_example_pairs.R`；PLINK / 临时文件放 `tmp/plink_fig8/`，SLURM 中 `TMPDIR` 也设在项目内。
- **版式**（与 paper 相同）：左栏 = 该 trait 的 GWAS / caQTL（全量）/ eQTL 三条 −log10 P 轨（灰 = 全部 SNP，蓝 = MR instruments，金三角 = SMR lead，紫菱形 = peak 内 GWAS SNP，GWAS 轨虚线 = 5e-8）；右栏 = Gene / Peaks / 6 个方法的 links（FDR < 0.05，颜色 = 方向，线宽 = −log10 FDR，三个 CT 共用标尺）；example peak 处紫色竖带；22×11 in，无图例。
- **2026-10-09 运行**：job 20080684，约 10 分钟、内存约 2.7G；27 个 pair、31 个 SNP-pair 行，共 32 页（CT × gene × trait），无作图失败；MR SNP 对账全部一致（SNP 数、lead SNP：Bn 7/7、CD4 10/10、CD8 10/10）；GWAS SNP 在 QTL 中找到 Bn 8/8、CD4 8/12、CD8 11/11（CD4 缺的 4 个旧版同样缺：旧版 11/15；找不到时 GWAS 点按 Valid6 表中的 P 补画）。

## 4. Update History
- 2026-10-09：建立本文件和 `Scripts/Visualization/`；写好 `00_common.R`、`Plot3_Hyper.R`、`Plot7_GRN.R` 并运行，输出 `Results/Figures/Fig3_Hyper.pdf`、`Fig7_GRN.pdf`。确定：Fig2 出 v1 和 v2；Fig5 只做第 1 页；Fig6 只出 `Motif_all.tsv` 并复制旧 xlsx；Fig8 暂不做。
- 2026-10-09：新增 `01_example_pairs.R`（example pairs 只保留一份）；写好并运行 `Plot1_MRresults.R`、`Plot4_AUC.R`、`Plot5_SLDSC.R`，输出 Fig1（5 页）、Fig4（3 页）、Fig5（第 1 页）；72 个 example pair 在新 FDR 下全部仍满足类别，沿用不重挑。
- 2026-10-09：写好并运行 `Plot6_Motif.R`：输出 `Results/Figures/Motif_all.tsv`，复制 `Fig6_Motif_table.xlsx`；与旧结果对照，16 个目标 motif 中 CD4 JUND（IVW）、CD8 FOSL2::JUND（MR_Egger）不再显著。
- 2026-10-09：移植 Fig2 v1 / v2（`Plot2_examples_v1.R`、`_v2.R`、`Plot2_examples_slurm.sh`），Egger 线截距改用重跑结果；job 20080188 跑完 3 个 CT，输出 6 个 PDF。
- 2026-10-09：移植 Fig8（`Plot8_Link.R`、`Plot8_Link_slurm.sh`，只做 GWAS 模式，读 Valid6_GWAS，9 个 GWAS）；job 20080684 输出 `Results/Figures/Fig8_Link_GWAS.pdf`（32 页）。至此 Fig1–Fig8 全部完成（Fig5 只有第 1 页、Fig6 只有表）。
