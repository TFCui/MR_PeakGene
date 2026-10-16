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
| Example pairs 表 | `Scripts/Visualization/Example_pairs.tsv`：`pair_set = original_example`（paper 的 24×3 个，`ord/type/category`）/ `Fig2_selected`（Fig2 实际画的 pair，`page/row/type`） |
| Example pairs 检查 | `Scripts/Visualization/01_example_pairs.R`（`check_examples()`，Fig1 用） |
| SLURM 日志 | `Scripts/Visualization/logs/` |
| 图 / 表 | 当前版 `Results/Figures/Figures_20261019/`；上一版 `Results/Figures/Figures_20261012/`（只作对比） |
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
- **FDR 阈值**：Fig1/Fig2 用 0.1；Fig8 link 用 0.05（Valid6 选 pair 用 MR FDR < 0.001）；Fig3 依赖 Valid1（0.1）；Fig5 依赖 Valid2（MR-Only，0.1）；Fig7 依赖 Valid5（0.05）。
- **运行**：轻的脚本在登录节点 `bash -lc 'module load ...; Rscript Scripts/Visualization/PlotN_*.R'`；Fig2、Fig8 需要 PLINK，用 `sbatch`。
- **检查图**：登录节点没有 pdftoppm，用 `module load GCCcore/13.3.0 Ghostscript/10.03.1` 后 `gs -q -dNOPAUSE -dBATCH -sDEVICE=png16m -r80 -sOutputFile=tmp/x_p%d.png <pdf>` 转 PNG 查看。

## 3. 各图

| 图 | 本项目脚本 | paper 目录原脚本 | 状态 |
|---|---|---|---|
| Fig1 MR 总体结果 | `Plot1_MRresults.R` | `Plot1_MRresults.sh` | ✅ |
| Fig2 Example pairs | `Plot2_examples_v2.R`（`_v1.R` 不再维护） | `Source_visual/Plot2_scripts_v2.R` | ✅ |
| Fig3 超几何富集 | `Plot3_Hyper.R` | `Plot3_val_Hyper.sh` | ✅（global test 表，不再画图） |
| Fig4 ROC AUC | `Plot4_AUC.R` | `Plot4_val_AUC.sh` | ✅ |
| Fig5 S-LDSC | `Plot5_SLDSC.R` | `Plot5_val_SLDSC.sh` Part 2 | ✅（只做第 1 页） |
| Fig6 Motif | `Plot6_Motif.R` | `Plot6_val_Motif.sh` | ✅（只出 `Motif_all.tsv` + 复制 xlsx） |
| Fig7 CIMA GRN | `Plot7_GRN.R` | `Plot7_val_GRN.sh` 第 6 节 | ✅ |
| Fig8 Link + GWAS | `Plot8_Link.R` | `Plot8_val_Link.sh`（已重新设计） | ✅（依赖 Valid6_GWAS） |

Fig1、Fig2、Fig3（改为表）、Fig8 为 2026-10-13 新版（`Figures_20261019/`）；Fig5 已于 2026-10-13 重跑（`Figures_20261019/`）；Fig4、Fig7 本轮未重跑，现有图在 `Figures_20261012/`，重跑后写到 `Figures_20261019/`。

### 3.1 Fig1 — MR 总体结果 ✅
- **脚本**：`Plot1_MRresults.R`（登录节点）；输入 `Results/MR/<CT>_AllPairs_results.tsv.gz`，FDR < 0.1 判显著。
- **输出**：`Fig1_MR_overview.pdf`，4 页，11×5.2 in。
  - 第 1–3 页：每个 CT 一张 UpSet（前 25 个组合）；右侧条形图标题 "Significant pairs"（`upset_right_annotation` 把 "Set size" 写死，生成后改 `@anno_list[[1]]@label`）。各方法显著数只打印到日志。
  - 第 4 页：六方法全显著的 pair，每个 CT 一个环形图——扇区 = Causal IVW 与 Marginal Multiomics / Pseudobulk 的方向关系（都一致 / 只与 Multiomics 不一致 / 只与 Pseudobulk 不一致 / 都不一致），中心 = n 和 SMR/IVW/Egger 内部方向一致比例。
- **Example pairs 检查**：`check_examples()` 检查 72 个 original example 在新旧 FDR 下是否满足类别，写 `tmp/example_pairs_check.tsv`（2026-10-09：全部满足）。

### 3.2 Fig2 — Example pairs ✅
- **脚本**：`Plot2_examples_v2.R` + `Plot2_examples_slurm.sh`（单个作业，约 7 分钟）。
- **选 pair**：`Example_pairs.tsv` 中 `pair_set = Fig2_selected` 的行（`page`、`row`、`type`）。第 1 页 = 4 个指定 pair；第 2–6 页 = 每页 Type 1–4 各一个。改 pair 只需改这张表再重跑。
- **方法**：对每个 pair 重跑 QTL 对齐 → 去重 → PLINK LD → GLS IVW/Egger + SMR 取 SNP 层数据；p/FDR/beta 取 MR 表（新 MR 表没有 Egger 截距，用重跑的截距画 Egger 线）。日志对账 IVW/Egger 斜率和 SMR lead SNP。
- **输出**：`Fig2_examples.pdf`，6 页，22×20 in；每页 4 行，每行一个 pair（标题 = gene | peak (CT)，副标题 = Type、peak–TSS 距离、MR SNP 数），行内 4 个 panel：
  - a：caQTL / eQTL −log10 P（peak/TSS 中点 ±500 kb；红三角 = MR SNP，金三角 = SMR lead；caQTL 轨不画 P < 1e-8 且不在 MR 中的 variant；角标 = SNPs in MR）；
  - b：MR 散点 + Causal IVW / Egger 线（alpha 0.75、线宽 1.6；不画 SMR 线，保留 lead SNP 标注）；
  - c：Marginal Multiomics (across metacells / cell types)，WNN metacell Spearman；
  - d：Marginal Pseudobulk (across individuals)，个体层 Spearman。

### 3.3 Fig3 — 超几何富集（global test 表）✅
- **输入**：`Results/Validation/Valid1_Hyper/<CT>_hyper_long.tsv`（method 为 `SMR/IVW/Egger/CCT/WNN/Pseudobulk`，database 为 `enhancerAtlas/FANTOM5/4DGenome/GTEx/HiChIP/PCHiC`，与旧表一致）。
- **脚本**：`Scripts/Visualization/Plot3_Hyper.R`（登录节点，几秒）。
- **方法**：不再画 gene 级 p 值直方图。只用 `N > 0, K > 0, n_sig > 0, N > K` 的 (gene, 方法, 金标准) 行；对每个 CT × 方法 × 金标准，把各 gene 的 N/K/n_sig/k 分别求和，得到 ΣN、ΣK、Σn_sig、Σk，再算两个量：
  - **Global P** = `phyper(Σk − 1, ΣK, ΣN − ΣK, Σn_sig, lower.tail = FALSE)`，即对合并计数做一次超几何检验（上尾）。用 `log.p = TRUE` 算 log10 P，极小的 P 不会下溢（旧图写 "< 1e-300"，现在给出具体指数，如 Bn SMR/HiChIP = 1.0e-416）。
  - **Fold enrichment（pooled）** = `(Σk / Σn_sig) / (ΣK / ΣN)`，代码里是 `(sum_k / sum_nsig) / (sum_K / sum_N)`。
- **符号含义**（每个 gene 内，来自 Valid1_Hyper）：N = 该 gene 在 All_GenePeak 中的全部候选 peak 数（背景）；K = 其中与该 gene 金标准区域重叠的 peak 数；n_sig = 该方法在该 gene 的显著 peak 数（FDR < 0.1）；k = 其中与金标准重叠的 peak 数。
- **fold 的含义**：方法选出的显著 peak 落在金标准里的比例（命中率 Σk/Σn_sig），除以候选 peak 落在金标准里的比例（背景率 ΣK/ΣN）。1 = 与随机相同，越大富集越强。Global P 检验的正是这两个比例的差距：fold 是效应大小，P 只回答"是否偶然"。
- **读表注意**：
  - 合并的是计数，不是把各 gene 的 fold 取平均，所以候选 peak 多的 gene 权重更大；
  - 背景是每个 gene 自己的候选 peak 集，不是全基因组 peak；
  - n_sig = 0（该方法在此 gene 无显著 peak）、K = 0（金标准无覆盖）或 N = K 的 gene 不参与；
  - 显著 peak 越多的方法，global P 越容易极小，但 fold 往往更低（如 Bn 的 Marginal Multiomics：P 最小、fold 最低），所以 P 和 fold 要一起看，不同方法的 fold 也不要脱离各自的显著 peak 数来比较；
  - fold 没有置信区间（未做 bootstrap）。
- **输出**（`Figures_20261019/`）：
  - `Fig3_Hyper_table.xlsx`：每个 CT 一个 sheet；上块 Global P（文本科学计数，行 = 6 方法，列 = 6 金标准），下块 pooled fold（每列最高 fold 的方法高亮），底部有脚注；
  - `Fig3_Hyper_global.tsv`：长表，108 行（3 CT × 6 方法 × 6 金标准），列 `n_gene, sum_N, sum_K, sum_nsig, sum_k, fold, log10_gp, gp, p_txt`（`gp` 在 < 1e-300 时下溢为 0，以 `log10_gp` 为准）。
  - 旧版直方图 `Figures_20261012/Fig3_Hyper.pdf` 保留作对比。

### 3.4 Fig4 — ROC AUC ✅
- **输入**：`Results/Validation/Valid4_AUC/<CT>_AUC_long.tsv.gz`（`gold_standard` 为 `enhancerAtlas/FANTOM5/...`，`method` 为小写 key）。
- **脚本**：`Scripts/Visualization/Plot4_AUC.R`
- **方法**：按 CT × 金标准 × 方法汇总 mean/median/sd；小提琴 + 箱线图，facet = 6 个金标准，0.5 虚线为随机水平。金标准显示名沿用 paper（EnhancerAtlas、FANTOM、4DGenome、GTEx cisQTL、HiChIP、PCHiC）。
- **输出**：`Fig4_AUC.pdf`，3 页，13×8.5 in。

### 3.5 Fig5 — S-LDSC ✅
- **输入**：`Results/Validation/Valid2_SLDSC/SLDSC_all.tsv`（`cell_type`、`set`、9 个短名 trait；只有 ext500）。
- **脚本**：`Scripts/Visualization/Plot5_SLDSC.R`
- **方法**：enrichment 柱图，不画 s.e.，enrichment raw P < 0.1 的柱子上方标粗体黑色 `*`（ggh4x `facet_grid2`，行 = CT（顺序 CD4 / Bn / CD8，沿用 paper），列 = 9 个 trait：AS, CD, IBD, RA, SLE, MS, T1D, Asthma, PBC）。每个面板 2 根柱：MR-Only（`MROnly`，浅天蓝 `#6FA8DC`）vs All peaks（`pool`，`COL_CONC`）；MR-Only = IVW/Egger/SMR 任一 FDR < 0.1 且 WNN、Pseudo 都不显著的 link 的 peak（2026-10-14 起，不再画 6 个方法）。日志打印 enrichment 表（`enrich_p < 0.1` 标 *）。
- **第 2 页**：Combined（3 个 CT 的 MR-Only 并集 vs 3 个 CT 的全部候选 peak 并集），1 行 × 9 个 trait，样式同第 1 页（`aspect.ratio = 2`）。
- **输出**：`Figures_20261019/Fig5_SLDSC.pdf`，2 页，12×7.5 in（旧版 6 方法、FDR 0.05、带 s.e. 的图留在 `Figures_20261012/`）。

### 3.6 Fig6 — Motif ✅
- **输入**：`Results/Validation/Valid3_Motif/<CT>_motif_enrichment.tsv.gz`
- **脚本**：`Scripts/Visualization/Plot6_Motif.R`
- **方法**：只把 3 个 CT 的 motif 富集结果合并成 `Motif_all.tsv`（列与 paper 的 `Motif_details/Motif_all.tsv` 相同，method 名仍为 `IVW/MR_Egger/SMR/CCT/WNN/Pseudo`）；不画 logo、不做 16 个 motif 的表、不生成新的 xlsx。
  日志里另做两项对照（只打印）：与 paper 旧 `Motif_all.tsv` 逐行比较 p 值和显著数；xlsx 中 16 个目标 motif 在新旧结果下是否显著（p.adjust < 0.05 且 FE > 1）。
- **2026-10-09 对照结果**：6 个 CT×方法组合（如 Bn 的 IVW/SMR/CCT）与旧结果完全相同，其余前景 peak 数差 1–27 个，p 值随之变化（每组 2–12 个 motif 的显著性翻转）；原因应是旧 Motif 代码的前景用的是旧 FDR 且保留了 sigma2 为 NA 的 pair，前景变了，按 GC/可及性匹配的背景也跟着变（未逐一核实）。
  16 个目标 motif 中 14 个结论不变；**CD4 JUND（MA0492.1）在 IVW 下不再显著**（SMR 下仍显著）；**CD8 FOSL2::JUND（MA1145.1）在 MR_Egger 下不再显著**（新结果中 6 个方法都不显著）；Bn ELF2 在 MR_Egger 下新变为显著。复制过来的 xlsx 仍是旧数字，用时需注意这两行。
- **输出**：`Motif_all.tsv`；另把 paper 目录的 `Materials/Fig6_Motif_table.xlsx` 原样复制到图目录（该 xlsx 是本地 Claude Code 结合外部知识手工整理的，数字来自旧结果）。

### 3.7 Fig7 — CIMA GRN ✅
- **输入**：`Results/Validation/Valid5_GRN/GRN_sig_counts_by_method.tsv`（统计量在 Valid5 中已算好：FDR < 0.05，2×2 Fisher）。
- **脚本**：`Scripts/Visualization/Plot7_GRN.R`
- **方法**：上图 = 每个方法显著且与 GRN 重叠的 pair 数（行 = CT / CT_HQ，列 = 细胞类型）；下图 = OR（GRN vs 非 GRN），CT 与 CT_HQ 并排。
- **输出**：`Fig7_GRN.pdf`，1 页，14×11 in。

### 3.8 Fig8 — Link + GWAS ✅
- **脚本**：`Plot8_Link.R` + `Plot8_Link_slurm.sh`（约 3 分钟）。
- **输入**：`Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv`（MR-only 强 link 且 peak 内有 GWAS P < 5e-8 的 SNP；按文件顺序给 pair 编号）；MR 结果；QTL + PLINK（只为取 MR SNP）；EnsDb.Hsapiens.v86。
- **输出**：`Fig8_Link_GWAS.pdf`，一页一个 pair（27 页），14×9 in。只画这一条 peak–gene link，不画其他 peak。
  - 轨道（上 → 下）：Gene / SNPs / Peak / 5 个方法（SMR、IVW、Egger、Multiomics、Pseudobulk；不画 Causal Joint）各一行 link / x 轴。
  - Link：FDR < 0.05 = 实线，颜色 = 方向，线宽 ∝ −log10 FDR（所有页共用标尺）；不显著 = 灰色虚线；行标签写 FDR。
  - 分段 x 轴：peak 中点 ± 1 kb 放大到 35% 宽（灰底），两侧按 bp 线性压缩，刻度写真实位置；peak–TSS < 2 kb 时不压缩。参数 `ZOOM_HALF`、`ZOOM_FRAC`。
  - SNPs 轨只画放大区内：蓝竖线 = MR instrument（重跑，与 `ivw_n_snp` / SMR lead 对账），金三角 = SMR lead，紫菱形 + rsID = GWAS SNP；副标题写 MR SNP 在窗口内/总数和 GWAS SNP 的 trait、P。
  - 部分 pair（如 DR1）MR SNP 都在 ±1 kb 外，SNPs 轨只有 GWAS SNP。

## 4. Update History
- 2026-10-09：建立本文件和 `Scripts/Visualization/`；写好 `00_common.R`、`Plot3_Hyper.R`、`Plot7_GRN.R` 并运行，输出 `Results/Figures/Fig3_Hyper.pdf`、`Fig7_GRN.pdf`。确定：Fig2 出 v1 和 v2；Fig5 只做第 1 页；Fig6 只出 `Motif_all.tsv` 并复制旧 xlsx；Fig8 暂不做。
- 2026-10-09：新增 `01_example_pairs.R`（example pairs 只保留一份）；写好并运行 `Plot1_MRresults.R`、`Plot4_AUC.R`、`Plot5_SLDSC.R`，输出 Fig1（5 页）、Fig4（3 页）、Fig5（第 1 页）；72 个 example pair 在新 FDR 下全部仍满足类别，沿用不重挑。
- 2026-10-09：写好并运行 `Plot6_Motif.R`：输出 `Results/Figures/Motif_all.tsv`，复制 `Fig6_Motif_table.xlsx`；与旧结果对照，16 个目标 motif 中 CD4 JUND（IVW）、CD8 FOSL2::JUND（MR_Egger）不再显著。
- 2026-10-09：移植 Fig2 v1 / v2（`Plot2_examples_v1.R`、`_v2.R`、`Plot2_examples_slurm.sh`），Egger 线截距改用重跑结果；job 20080188 跑完 3 个 CT，输出 6 个 PDF。
- 2026-10-09：移植 Fig8（`Plot8_Link.R`、`Plot8_Link_slurm.sh`，只做 GWAS 模式，读 Valid6_GWAS，9 个 GWAS）；job 20080684 输出 `Results/Figures/Fig8_Link_GWAS.pdf`（32 页）。至此 Fig1–Fig8 全部完成（Fig5 只有第 1 页、Fig6 只有表）。
- 2026-10-13：Fig1（去计数柱图，方向一致性改为环形图，UpSet 右侧改 "Significant pairs"）、Fig2（每页 4 个指定 pair、去掉位点图）、Fig8（只画单条 link、分段 x 轴、去掉 GWAS/QTL 轨和 Causal Joint）重画，输出 `Results/Figures/Figures_20261019/`；pair 表合并为 `Example_pairs.tsv`。
- 2026-10-13：Fig3 去掉直方图，改为 global test 表（Global P + pooled fold enrichment，P 用 `log.p` 计算不下溢），输出 `Figures_20261019/Fig3_Hyper_table.xlsx` 和 `Fig3_Hyper_global.tsv`；fold 定义写入 3.3 节。
- 2026-10-13：Valid2_SLDSC 改为 FDR < 0.1 的显著 peak 重跑；Fig5 去掉 s.e. 误差棒，raw enrichment P < 0.1 标粗体 `*`，输出到 `Figures_20261019/`。
- 2026-10-14：Fig5 改为 MR-Only vs All peaks（每面板 2 根柱，MR-Only 浅天蓝），不再画 6 个方法；12×7.5 in。
- 2026-10-14：Fig5 加第 2 页 Combined（Bn + CD4 + CD8 合并的 MR-Only vs All peaks）；画图部分抽成 `page_dt()` / `plot_page()`。
