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
| 图 / 表 | 论文图表新一版 `Results/Figures/Figures_20261102/`（逐图更新中，尚未更新的图仍以 `Figures_20261027/` 为准）；`Figures_20261027/`（含 Figma 手工拼的 `Figure_1.pdf`）；Fig5–Fig8 等旧图在 `Figures_20261019/`；更早一版 `Figures_20261012/`（只作对比） |
| 临时文件（PLINK LD 等） | `tmp/` |

## 2. 共同约定

- **输入**：
  - MR：`Results/MR/<CT>_AllPairs_results.tsv.gz`（已删 NA、`ivw_sigma2 <= 10`、gene 内 BH）。画图时**不再**做 sigma2 过滤。
  - 验证：`Results/Validation/Valid<k>_<名>/`（见 `Pipeline.md` 第 3 节）。
- **细胞类型**：`Bn_TCL1A`、`CD4_Tn_CCR7`、`CD8_Tn_CCR7`；图中标签统一写 **Bn / CD4 T cell / CD8 T cell**（`00_common.R` 的 `CELL_LAB`，规则见 `CLAUDE.md` 第 1 节；2026-11-02 起逐图改用）。
- **输出格式**（2026-11-02 起）：除 `Figure_1B_C` 只出 SVG（`svglite`，供 Figma 拼 Fig 1）外，所有图**只输出 PDF**，不再出 SVG。
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
| Fig10 S-LDSC + Link 合成 | `Plot10_Combined.R` | —（Fig5 第 2 页 + Fig8 的代码合并） | ✅（初版，版面已定，pair 可换） |

Fig1、Fig2、Fig3（改为表）、Fig8 为 2026-10-13 新版（`Figures_20261019/`）；Fig5 已于 2026-10-13 重跑（`Figures_20261019/`）；Fig4、Fig7 本轮未重跑，现有图在 `Figures_20261012/`，重跑后写到 `Figures_20261019/`。

### 3.1 Fig1 — MR 总体结果 ✅
- **脚本**：`Plot1_MRresults.R` + `Plot1_MRresults_slurm.sh`（约 1 分钟）；输入 `Results/MR/<CT>_AllPairs_results.tsv.gz`，FDR < 0.1 判显著（图中写作 p.adj < 0.1）。
- **输出**（`Figures_20261102/`）：`Figure_1B_C.svg`（CD4 T cell，只出 SVG，`svglite`，9×3.6 in；Figma 拼 `Figure_1` 用，脚本默认不重写，设 `PLOT1_WRITE_1BC=1`（slurm 脚本已设）才重出）+ `Figure_S3.pdf`（只出 PDF，Bn 在上、CD8 T cell 在下，同一页，9×7.2 in，补充用；每个细胞类型一个 panel 字母 A / B，画在每行左上角）。每个 CT 一行，左 UpSet（列宽 0.60）、右环形图（0.40）；两侧标题（细胞类型标签 / “Effect direction”）由 `plot_row()` 用 `grid.text` 画在同一行，环形图标题下两行灰字说明（“Significant in all six methods (p.adj < 0.1)” / “Slices: sign of Causal IVW vs. Marginal Multiomics / Pseudobulk”）：
  - 左：UpSet（前 15 个组合；柱顶数字 45° 倾斜；上方 y 轴名 “Size”（竖排）；右侧条形图标题 "Significant pairs"，两个名字都靠改 `@anno_list[[1]]@label` 实现）；字号：方法名 9、柱顶数字 7.5、右侧数字 8。各方法显著数只打印到日志。
  - 右：六方法全显著的 pair 的环形图——扇区 = Causal IVW 与 Marginal Multiomics / Pseudobulk 的方向关系（都一致 / 只与 Multiomics 不一致 / 只与 Pseudobulk 不一致 / 都不一致），≥3% 扇区标百分比，中心写 n 和 SMR / IVW / Egger 方向一致比例，底部单列图例（2026-11-02 恢复 10-25 之前的写法）；`coord_polar` 四周固定留约 10% 空白，用负的 `plot.margin` / `legend.box.spacing` 抵消。
- **Figure_1A**：手工示意图（`Figures_20261019/Fig0_Schematics.pdf` 的修改版），`Figures_20261027/Figure_1A.{pdf,svg}`。2026-11-02 版：用户更新的 `Figures_20261102/Fig0_Schematics.pdf`（639×315 pt），SVG 用 poppler 的 `pdftocairo -svg` 转出 `Fig0_Schematics.svg`（`module load GCC/13.2.0 poppler/24.04.0`；文字转为字形路径，不可直接编辑），供 Figma 拼接。
- **Figure_1**：论文 Fig 1 = 在 Figma 中把 `Fig0_Schematics.svg`（1A）与 `Figure_1B_C.svg`（CD4 T cell）拼接排版的结果，`Figures_20261102/Figure_1.pdf`（手工，非脚本生成，图内已有 A / B / C 标签）。Figma 导出含带 SMask 的透明位图，部分阅读器显示黑底，论文用 Ghostscript 拍平版 `Figure_1_flat.pdf`（`sbatch tmp/Fig1_flatten/flatten_slurm.sh <in.pdf> <out.pdf>`：PDF 1.3、整页 720 ppi、Flate 无损）。
- **Example pairs 检查**：`check_examples()` 检查 72 个 original example 在新旧 FDR 下是否满足类别，写 `tmp/example_pairs_check.tsv`（2026-10-09：全部满足）。

### 3.2 Fig2 — Example pairs ✅
- **脚本**：`Plot2_examples_v2.R` + `Plot2_examples_slurm.sh`（单个作业，约 7 分钟）。
- **选 pair**：`Example_pairs.tsv` 中 `pair_set = --set` 的行（`page`、`row`、`type`，默认 `Figure_2`）。`Figure_2` = 1 页：SIDT2 (Bn)、RPL11、ERAP2、ARPC1B (CD4)；`Fig2_selected` = 上一版 6 页。改 pair 只需改这张表再重跑。
- **方法**：对每个 pair 重跑 QTL 对齐 → 去重 → PLINK LD → GLS IVW/Egger + SMR 取 SNP 层数据；p/FDR/beta 取 MR 表（新 MR 表没有 Egger 截距，用重跑的截距画 Egger 线）。日志对账 IVW/Egger 斜率和 SMR lead SNP。
- **输出**：`Figures_20261102/Figure_2.pdf`（只出 PDF），1 页，22×20 in；每页 4 行（行标签 A–D，行内 panel 不标字母），每行一个 pair（标题 = gene | peak (细胞类型，`CELL_LAB`)，无副标题；分类型的 pair 集才写 “Type k: …”），列宽 1.2 : 1 : 1 : 1，行内 4 个 panel：
  - a：caQTL / eQTL −log10 P（peak/TSS 中点 ±500 kb；红三角 = MR SNP，金三角 = SMR lead；caQTL 轨不画 P < 1e-8 且不在 MR 中的 variant；两轨间距收紧；角标 = 红三角 + “N SNPs”）；
  - b：MR 散点 + Causal IVW / Egger 线（alpha 0.75、线宽 1.6；不画 SMR 线，保留 lead SNP 标注）；
  - c：Marginal Multiomics (across metacells / cell types)，WNN metacell Spearman；
  - d：Marginal Pseudobulk (across individuals)，个体层 Spearman。c / d 的线都是 `lm` 回归实线 + 置信带（c 灰色，d 深绿）。

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
- **输出**（`Figures_20261027/`；2026-10-20 前在 `Figures_20261019/`，文件名 `Fig3_Hyper_table.xlsx` / `Fig3_Hyper_global.tsv`）：
  - `Table_1_2.xlsx`（论文 S1 Table = Global P，S2 Table = pooled fold；`Figures_20261102/` 中的副本改名 `Table_S1_S2.xlsx` / `Table_S1_S2_global.tsv`）：3 个 CT 合并在**单个 sheet**（2026-10-23 起）；上块 Global P（文本科学计数，行 = 3 CT × 6 方法，列 = 6 金标准），下块 pooled fold（每个 CT 内每列最高 fold 的方法高亮），底部有脚注；
  - `Table_1_2_global.tsv`：长表，108 行（3 CT × 6 方法 × 6 金标准），列 `n_gene, sum_N, sum_K, sum_nsig, sum_k, fold, log10_gp, gp, p_txt`（`gp` 在 < 1e-300 时下溢为 0，以 `log10_gp` 为准）。
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
- **MR 特异 motif 表**：Valid3 步骤 2 现在直接输出 `Results/Validation/Valid3_Motif/Motif_MR_specific.tsv`（47 个）和文献注释 `Motif_MR_specific_literature.tsv`（见 `Pipeline.md` 3.3）；上面的旧 xlsx（16 个 motif，旧数字）只作对照。

### 3.7 Fig7 — CIMA GRN ✅（2026-10-29 起只保留 OR，并入论文 Figure 3C；S5 Fig 已删除）
- **输入**：`Results/Validation/Valid5_GRN/GRN_sig_counts_by_method.tsv`、`GRN_summary.tsv`（统计量在 Valid5 中已算好：p.adj（FDR）< 0.05，2×2 Fisher）。
- **脚本**：旧版 `Plot7_GRN.R`（输出 `Fig7_GRN.pdf`，14×11 in，不再维护）；论文图用 `Plot9_GO_GRN.R`（与 GO 合并，见 3.10；原 `Plot_Sup2_GRN.R` 已删除）。
- **方法**：OR（GRN vs 非 GRN），每个细胞类型一个小面板（Bn / CD4 T cell / CD8 T cell），CT 与 CT_HQ 并排（浅蓝 / 亮蓝 `#A6CEE3` / `#1F78B4`，ColorBrewer Paired），虚线 OR = 1。日志打印 eRegulon 数、与 MR 重叠的 GRN pair 数、计数表、OR 表和最大 Fisher P（1.2e-31，所有比较 P < 1e-30）。原 S5 的 panel A（重叠 pair 数）不再画。

### 3.8 Fig8 — Link + GWAS ✅
- **脚本**：`Plot8_Link.R` + `Plot8_Link_slurm.sh`（约 3 分钟）。
- **输入**：`Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv`（MR-only 强 link 且 peak 内有 GWAS P < 5e-8 的 SNP；按文件顺序给 pair 编号）；MR 结果；QTL + PLINK（只为取 MR SNP）；EnsDb.Hsapiens.v86。
- **输出**：`Fig8_Link_GWAS.pdf`，一页一个 pair（27 页），14×9 in。只画这一条 peak–gene link，不画其他 peak。
  - 轨道（上 → 下）：Gene / SNPs / Peak / 5 个方法（SMR、IVW、Egger、Multiomics、Pseudobulk；不画 Causal Joint）各一行 link / x 轴。
  - Link：FDR < 0.05 = 实线，颜色 = 方向，线宽 ∝ −log10 FDR（所有页共用标尺）；不显著 = 灰色虚线；行标签写 FDR。
  - 分段 x 轴：peak 中点 ± 1 kb 放大到 35% 宽（灰底），两侧按 bp 线性压缩，刻度写真实位置；peak–TSS < 2 kb 时不压缩。参数 `ZOOM_HALF`、`ZOOM_FRAC`。
  - SNPs 轨只画放大区内：蓝竖线 = MR instrument（重跑，与 `ivw_n_snp` / SMR lead 对账），金三角 = SMR lead，紫菱形 + rsID = GWAS SNP；副标题写 MR SNP 在窗口内/总数和 GWAS SNP 的 trait、P。
  - 部分 pair（如 DR1）MR SNP 都在 ±1 kb 外，SNPs 轨只有 GWAS SNP。

### 3.9 Fig10 — S-LDSC（Combined）+ Link/GWAS 合成图 ✅
- **脚本**：`Plot10_Combined.R` + `Plot10_Combined_slurm.sh`（约 3–5 分钟，需要 PLINK + EnsDb，用 sbatch）。由 `Plot5_SLDSC.R`（Combined 页）和 `Plot8_Link.R`（单 pair 位点图）的代码合并拷贝而成，**Fig5、Fig8 及其脚本不动**。
- **输出**（2026-10-29 起）：`Figures_20261102/Figure_4.pdf`（只出 PDF；一页，panel B 1 行 × 5 个 pair，24×11.5 in；页高 = `PAGE_H_A` 4.5 + `PAGE_H_ROW` 7.0 × 行数）；不再输出 `Figure_S4`。行名列宽 0.3，pair 内容缩放 0.99、轨道左右边距缩小，例子之间更紧凑。panel B 的行名（Gene / SNPs / Peak / 5 个方法）只在最左边一列写一次。
- **panel a**：S-LDSC Combined（3 个 CT 并集）的 MR-Only vs All peaks，2026-10-29 起只画 7 个 GWAS（去掉 Asthma、PBC，`SL_TRAITS`）；单个面板，9 个 GWAS 共用同一个 y 轴，GWAS 名写在 x 轴，每个 GWAS 的两根柱紧贴（2026-10-23 起；原为 9 个独立 y 轴的小面板）；不放标题和解释文字，只留图例、y 轴标题和 enrichment P < 0.1 的 `*`。MR-Only 柱为亮蓝 `#1F78B4`（2026-10-29 由紫红改，与 Fig 3C 同色）。数据 `Valid2_SLDSC/SLDSC_all.tsv`。
- **panel b**：1 × 5 并排（2026-10-29：去掉每个 pair 顶部的 Gene / Peak 文字；显著 link 和星号统一红色，不再按效应方向分红 / 蓝；底部细胞类型用 `CELL_LAB`），内容同 Fig8（Gene / SNPs / Peak / 5 个方法 link / 分段 x 轴），无边框，相邻 pair 之间一条浅灰细竖线（`SEP_COL`）；每个 pair 顶部 `Gene: XXX  Peak: chr:start-end`（hg38，取 `peak_id`），底部 `细胞类型短名 | GWAS trait`（trait 去重）；左侧行名字号 `LABEL_SIZE`；x 轴标题缩成 `chrN (Mb)`，`N_SNP_TARGET = 6`，`BASE_SIZE = 12`。
- **选 pair**：脚本里的 `PAGES`（`c(cell_type, gene)`，来自 `GWAS_MRonly_overlap.tsv`）。当前（2026-10-29）一页 5 个：Bn RGL2 / CD4 ERBB3 / CD4 POLR2E / CD4 UBAC2 / CD8 GNAS（chr20:59254001-59254502，IBD、CD；原在 Figure_S4）。
- **注意**：
  - 先开 PDF 设备再构建图会让设备自动开出一张空白第 1 页（patchworkGrob 量布局时在设备上画东西）；第一页 `print(newpage = FALSE)` 解决。
  - 设环境变量 `FIG10_PNG_DIR` 时改输出 PNG（`pg1.png`…，每页一张），用于检查版面，不写 Figures；预览脚本在 `tmp/Plot10_test/`。
  - peak 离 TSS 很近的 pair（如 POLR2E）link 弧线窄成一条线，是 Fig8 就有的特点。

### 3.10 Figure 3 — GO 富集（DAVID）+ GRN odds ratio ✅
- **脚本**：`Plot9_GO_GRN.R`（登录节点，几秒；2026-10-29 由 `Plot9_GO.R` 与 `Plot_Sup2_GRN.R` 合并）。输入 3 个 CT 合并的 MR-Only gene（`Valid7_GO/DAVID_input/MROnly_genes_all.txt`，4,509 个）在 DAVID 网站的结果 `DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv`（只含 GO BP / MF）。
- **输出**：`Figures_20261102/Figure_3.pdf`（只出 PDF），13×11 in（上下行高 1 : 1.1）；上行 A = GO BP、B = GO MF 并排（图例收在右侧），下行 C = GRN odds ratio 通栏（见 3.7；`free()` 使 C 不与上行轴对齐）；GO 部分：dotplot，每个 panel 按 DAVID FDR 取前 10 个 term（x = GeneRatio，大小 = Count，颜色 = DAVID FDR，图例写作 p.adj，< 1e-10 截为最深色）。

### 3.11 Figure S1 — One-sample MR QC ✅
- **脚本**：`Plot_Sup1_OneSampleMR.R`（登录节点）；输入 `Results/QC/1_OneSampleMR/<CT>_pair_corr.tsv.gz`、`<CT>_global_rho.tsv`。
- **输出**：`Figures_20261102/Figure_S1.pdf`（只出 PDF），1 页，9.5×4.8 in；2 个 panel：A = pair 内相关 ρ̂，B = 相关性检验的 P 值；x 轴 = pair 内共享 null SNP 数分组，每组 3 根并排箱子 = Bn / CD4 T cell / CD8 T cell（`CELL_LAB`；ColorBrewer Dark2：紫 `#7570B3` / 绿 `#1B9E77` / 橙 `#D95F02`，箱体半透明、描边同色加深）；图例写各 CT 的全局 ρ 和 P；图内不放 n 和虚线说明（pair 数、SNP 数、虚线含义写在论文 caption）。日志打印各 CT 的 pair 数 / SNP 数 / 全局 ρ、P 和每组中位数。

### 3.12 Table 4 — 三个细胞类型的数据概览 ✅（论文 S4 Table）
- **脚本**：`Table4_DataOverview.R`（`Table4_DataOverview_slurm.sh`，约 25 分钟，主要是 QTL 汇总统计的全文件扫描）；输入：`libsize_<CT>.csv`、pseudobulk 矩阵、eQTL / caQTL 汇总统计的 `sample_sizes` 列、`caQTL_sig/<CT>_Significant.txt.gz`、`Scripts/Model/logs/MR_AllPairs_20074030_{0,1,2}.out`（漏斗）、`Results/MR/<CT>_AllPairs_results.tsv.gz`（核对和 `ivw_n_snp`）。
- **输出**：`Figures_20261027/Table_4.xlsx`（论文 S3 Table，2026-10-29 前为 S4；`Figures_20261102/` 中的副本改名 `Table_S3.xlsx` / `Table_S3.tsv`，xlsx 标题行改为 “S3 Table …”；单 sheet，4 块：单细胞 / cis-QTL / 候选 pair 与 MR 分析集 / 边际方法样本）+ `Table_4.tsv`（长表）。漏斗末行的 pair / 基因 / peak 数与 MR 结果表用 `stopifnot` 核对。
- **口径**：caQTL peak 数 = 名义 P < 1e-8（选工具变量的阈值），不是 CIMA 的 FDR caPeak；没有 eGene 数。论文 tex 由 `Manuscript` 侧 `tmp/Manuscript_tables/make_data_table.py` 从 `Table_4.tsv` 生成。

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
- 2026-10-16：新增 Fig10（`Plot10_Combined.R`）：Fig5 Combined 柱状图（panel a）+ Fig8 单 pair 位点图 1×4（panel b），3 页，24×11 in，输出 `Figures_20261019/Fig10_SLDSC_Link.pdf`；MR-Only 柱改紫红色；修复第 1 页空白。
- 2026-10-19：论文图整理到 `Results/Figures/Figures_20261026/`，每张图 PDF + SVG（`svglite`，文字可编辑），图中 FDR 一律写作 p.adj：Figure_1A（复制 Fig0）、Figure_1B_C（每 CT 一页 UpSet + 环形图）、Figure_2（4 个 pair，行标签 A–D）、Figure_3（DAVID GO BP / MF dotplot，3 CT 合并 gene list）、Figure_4（原 Fig10，panel B 行名只写一次、弧线标星号）、Figure_S1（OneSampleMR 3 CT 合并为 1 页）、Table_1（复制 Fig3 xlsx）。Table_2（MR 特异 motif + logo + 文献注释）由 `Scripts/Validation/Valid3_Motif_3_table.R` 生成，见 `Pipeline.md` 3.3。
- 2026-10-20：Figure_S2（QC2 flip 前后对比，`Scripts/QC/QC2_2_FlipPoints_plot.R`，只出 PDF）；Fig6 说明加上 Valid3 新的 MR 特异 motif 表。
