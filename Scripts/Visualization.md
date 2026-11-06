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

脚本编号 `Plot<N>_` 是本项目内部的固定编号，按新增顺序往后排（下一个是 `Plot16_`），**不跟论文图号对应**：论文图会不断修改、合并，图号会变，改号时不改脚本名，只在本表和各小节标题中写“论文当前 … Fig”。输出文件名（`Figure_S1.pdf` 等）跟当前论文图号走。

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
| Fig11 One-sample MR QC（论文当前 S1 Fig） | `Plot11_OneSampleMR.R`（原 `Plot_Sup1_OneSampleMR.R`） | — | ✅ |
| Fig12 跨细胞类型 Venn（论文当前 S4 Fig） | `Plot12_Venn.R`（原 `Plot_Sup4_Venn.R`） | — | ✅ |
| Fig13 MR-Only 观测 vs 背景（论文当前 S5 Fig） | `Plot13_MROnly_ObsExp.R` | — | ✅ |
| Fig14 5 个方法 effect 两两相关与密度矩阵（论文暂定 S6 Fig） | `Plot14_MethodCorr.R` | `Server6_wgen1_plot.sh`（`pairs()` + `smoothScatter`） | ✅ |
| Fig15 全基因组 Manhattan / Miami 图（论文暂定 S7 Fig） | `Plot15_Manhattan.R` | — | ✅ |
| Table 5 GWAS 性状汇总（论文暂定 S4 Table；表，不是 Plot） | `Table5_GWASTraits.R` | — | ✅ |

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

### 3.7 Fig7 — CIMA GRN ✅（2026-10-29 起只保留 OR，并入论文 Figure 3C；2026-11-04 起 Figure 3E 只画 CD4 T cell，Bn / CD8 T cell 移到 Figure S8；S5 Fig 已删除）
- **输入**：`Results/Validation/Valid5_GRN/GRN_sig_counts_by_method.tsv`、`GRN_summary.tsv`（统计量在 Valid5 中已算好：p.adj（FDR）< 0.05，2×2 Fisher）。
- **脚本**：旧版 `Plot7_GRN.R`（输出 `Fig7_GRN.pdf`，14×11 in，不再维护）；论文图用 `Plot9_GO_GRN.R`（与 GO 合并，见 3.10；原 `Plot_Sup2_GRN.R` 已删除）。
- **方法**：OR（GRN vs 非 GRN），每个细胞类型一个小面板（Bn / CD4 T cell / CD8 T cell），CT 与 CT_HQ 并排（2026-11-04 起为 NPG 藏蓝系浅 / 深 `#A9B4D6` / `#3C5488`，原为 ColorBrewer Paired 蓝），虚线 OR = 1。日志打印 eRegulon 数、与 MR 重叠的 GRN pair 数、计数表、OR 表和最大 Fisher P（1.2e-31，所有比较 P < 1e-30）。原 S5 的 panel A（重叠 pair 数）不再画。

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
- **2026-11-06 微调**：panel b 的 pair 改为 RGL2 / ERBB3 / POLR2E / ERAP1 (Bn, CD) / GNAS（UBAC2 只有 Asthma，去掉；读 `GWAS_MRonly_overlap.tsv` 时过滤 `trait != "Asthma"`）；配色改 Set3（柱 MR-Only `#8dd3c7` / All peaks `#d9d9d9`，Gene `#80b1d3`，MR SNP `#a29cc9`，GWAS SNP `#bc80bd`，Peak `#fdb462`，红色 link 不变）；`BASE_SIZE` 12→16、A 的 `BS` 12→18、tag 18→32、行名列宽 0.3→0.4、`TICK_SEP` 0.14→0.27、`PAGE_H_A` 4.5→5.0；panel a 图例放图内左上角。
- **panel a**：S-LDSC Combined（3 个 CT 并集）的 MR-Only vs All peaks，2026-10-29 起只画 7 个 GWAS（去掉 Asthma、PBC，`SL_TRAITS`）；单个面板，9 个 GWAS 共用同一个 y 轴，GWAS 名写在 x 轴，每个 GWAS 的两根柱紧贴（2026-10-23 起；原为 9 个独立 y 轴的小面板）；不放标题和解释文字，只留图例、y 轴标题和 enrichment P < 0.1 的 `*`。MR-Only 柱为亮蓝 `#1F78B4`（2026-10-29 由紫红改，与 Fig 3C 同色）。数据 `Valid2_SLDSC/SLDSC_all.tsv`。
- **panel b**：1 × 5 并排（2026-10-29：去掉每个 pair 顶部的 Gene / Peak 文字；显著 link 和星号统一红色，不再按效应方向分红 / 蓝；底部细胞类型用 `CELL_LAB`），内容同 Fig8（Gene / SNPs / Peak / 5 个方法 link / 分段 x 轴），无边框，相邻 pair 之间一条浅灰细竖线（`SEP_COL`）；每个 pair 顶部 `Gene: XXX  Peak: chr:start-end`（hg38，取 `peak_id`），底部 `细胞类型短名 | GWAS trait`（trait 去重）；左侧行名字号 `LABEL_SIZE`；x 轴标题缩成 `chrN (Mb)`，`N_SNP_TARGET = 6`，`BASE_SIZE = 12`。
- **选 pair**：脚本里的 `PAGES`（`c(cell_type, gene)`，来自 `GWAS_MRonly_overlap.tsv`）。当前（2026-10-29）一页 5 个：Bn RGL2 / CD4 ERBB3 / CD4 POLR2E / CD4 UBAC2 / CD8 GNAS（chr20:59254001-59254502，IBD、CD；原在 Figure_S4）。
- **注意**：
  - 先开 PDF 设备再构建图会让设备自动开出一张空白第 1 页（patchworkGrob 量布局时在设备上画东西）；第一页 `print(newpage = FALSE)` 解决。
  - 设环境变量 `FIG10_PNG_DIR` 时改输出 PNG（`pg1.png`…，每页一张），用于检查版面，不写 Figures；预览脚本在 `tmp/Plot10_test/`。
  - peak 离 TSS 很近的 pair（如 POLR2E）link 弧线窄成一条线，是 Fig8 就有的特点。

### 3.10 Figure 3 — Miami 图 + GO 富集（DAVID）+ HiChIP 校准 + GRN odds ratio ✅（另出 Figure S8）
- **脚本**：`Plot9_GO_GRN.R`（登录节点，几秒；2026-10-29 由 `Plot9_GO.R` 与 `Plot_Sup2_GRN.R` 合并）。GO 输入为 3 个 CT 合并的 MR-Only gene（`Valid7_GO/DAVID_input/MROnly_genes_all.txt`，4,509 个）在 DAVID 网站的结果 `DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv`（只含 GO BP / MF）。两个 panel 由其他脚本构建、存 rds 后读入：A 来自 `Plot15_Manhattan.R`（`tmp/Plot15_Manhattan/miami_CD4.rds`，见 3.16），D 来自 `Plot13_MROnly_ObsExp.R`（`tmp/Plot13_MROnly/hichip_fig3.rds`，见 3.14）。改这两个 panel 后先跑对应脚本，再跑 Plot9。
- **输出**：`Figures_20261102/Figure_3.pdf`（只出 PDF），13×14 in。版式（2026-11-04 第 2 轮）：A 通栏，下面是 2×2 网格（左列 B、C，右列 D、E；A : 网格高度 1.05 : 2，左 : 右宽度 0.85 : 1）。网格是扁平的 `wrap_plots(…, ncol = 2)`，同一行两图的绘图区上下对齐；网格按行排、字母按列走，用自定义 tag 顺序。tag 画在图左上角（`plot.tag.position = c(0, 1)`），不额外占一行，以减少 panel 之间的空白。
  - 各 panel 不写 title / subtitle，统一用浅薄荷色小标签条（`00_common.R` 的 `theme_strip_mint()`）写简洁标题，括号内写细胞类型 / 数据库：A “Genome-wide peak–gene associations (CD4 T cell)”、B “GO biological process”、C “GO molecular function”、D “MR-only link overlap (HiChIP)”、E “CIMA GRN overlap (CD4 T cell)”。B、C、D、E 都用完整加粗边框（`TH_BOX`），字号 base 11，tag 18。
  - A：CD4 T cell 的 Miami 图（见 3.16），用 `wrap_elements` 作为一个 panel；方法名用论文写法（`lab_paper()`：Causal joint / Marginal multiomics）。
  - B = GO BP、C = GO MF（上下排）：dotplot，每个类别按 DAVID FDR 取前 10 个 term（x = GeneRatio，大小 = Count，颜色 = DAVID FDR，图例写作 p.adj，< 1e-10 截为最深色）；term 名按 24 字符折行；不写 gene 数。两个 panel 共用标尺，图例只放在 C 内右下角（p.adj 色条与 Count 图例顶端对齐、同高）。
  - D = MR-Only link 与 HiChIP 的分箱校准图（S5 panel E 的同一口径；`panel(fig3 = TRUE, cart = TRUE)`：边框 / 刻度 / 分箱线条加粗，去掉 0.05 刻度，右下角加粗写 `Fold = 3.28 (P = 6.4e-620)`，用 `coord_cartesian`（x、y 同范围）而不是固定纵横比，以便与 E 在网格中对齐，因此对角线不一定是 45°）；E = GRN odds ratio，只画 CD4 T cell（见 3.7；方法名两行水平、论文写法，图例在图内右上角）。
- **Figure S8**：`Figures_20261102/Figure_S8.pdf`（10×5.2 in），Bn 和 CD8 T cell 的 GRN odds ratio（2026-11-04 从 Figure 3 移出，口径、配色与 E 相同；标签条只写细胞类型）。
- **配色**（2026-11-04，薄荷色板 `PAL_MINT`，定义在 `00_common.R`：`#EBF5E6` `#C4E9CA` `#9ACDA2` `#6FBD8F` `#47A49D` `#336E99` `#2E4E77`）：A causal 钢蓝 `#336E99`、marginal 薄荷绿 `#6FBD8F`、非显著点浅冷灰交替、例子 pair 柔和珊瑚 `#E07A5F`（全图唯一暖色）；B/C p.adj 浅薄荷 `#C4E9CA` → 深蓝 `#2E4E77`（点带深蓝细描边）；D 分箱深蓝 `#2E4E77`、gene 点 `#9ACDA2`、富集侧底色 `#F5FAF3`、对角线灰；E / S8 cell-type eRegulons 浅薄荷 `#9ACDA2` / high-quality eRegulons 亮青蓝 `#3492B2`；轴线 / 边框统一 `#374151`。

### 3.11 Fig11 — One-sample MR QC ✅（论文当前 S1 Fig）
- **脚本**：`Plot11_OneSampleMR.R`（原 `Plot_Sup1_OneSampleMR.R`；登录节点）；输入 `Results/QC/1_OneSampleMR/<CT>_pair_corr.tsv.gz`、`<CT>_global_rho.tsv`。
- **输出**：`Figures_20261102/Figure_S1.pdf`（只出 PDF），1 页，9.5×4.8 in；2 个 panel：A = pair 内相关 ρ̂，B = 相关性检验的 P 值；x 轴 = pair 内共享 null SNP 数分组，每组 3 根并排箱子 = Bn / CD4 T cell / CD8 T cell（`CELL_LAB`；ColorBrewer Dark2：紫 `#7570B3` / 绿 `#1B9E77` / 橙 `#D95F02`，箱体半透明、描边同色加深）；图例写各 CT 的全局 ρ 和 P；图内不放 n 和虚线说明（pair 数、SNP 数、虚线含义写在论文 caption）。日志打印各 CT 的 pair 数 / SNP 数 / 全局 ρ、P 和每组中位数。

### 3.12 Table 4 — 三个细胞类型的数据概览 ✅（论文 S4 Table）
- **脚本**：`Table4_DataOverview.R`（`Table4_DataOverview_slurm.sh`，约 25 分钟，主要是 QTL 汇总统计的全文件扫描）；输入：`libsize_<CT>.csv`、pseudobulk 矩阵、eQTL / caQTL 汇总统计的 `sample_sizes` 列、`caQTL_sig/<CT>_Significant.txt.gz`、`Scripts/Model/logs/MR_AllPairs_20074030_{0,1,2}.out`（漏斗）、`Results/MR/<CT>_AllPairs_results.tsv.gz`（核对和 `ivw_n_snp`）。
- **输出**：`Figures_20261027/Table_4.xlsx`（论文 S3 Table，2026-10-29 前为 S4；`Figures_20261102/` 中的副本改名 `Table_S3.xlsx` / `Table_S3.tsv`，xlsx 标题行改为 “S3 Table …”；单 sheet，4 块：单细胞 / cis-QTL / 候选 pair 与 MR 分析集 / 边际方法样本）+ `Table_4.tsv`（长表）。漏斗末行的 pair / 基因 / peak 数与 MR 结果表用 `stopifnot` 核对。
- **口径**：caQTL peak 数 = 名义 P < 1e-8（选工具变量的阈值），不是 CIMA 的 FDR caPeak；没有 eGene 数。论文 tex 由 `Manuscript` 侧 `tmp/Manuscript_tables/make_data_table.py` 从 `Table_4.tsv` 生成。

### 3.13 Fig12 — 各方法显著 link 的跨细胞类型 Venn 图 ✅（论文当前 S4 Fig）
- **脚本**：`Plot12_Venn.R`（原 `Plot_Sup4_Venn.R`；登录节点，约 30 秒）；输入 `Results/MR/<CT>_AllPairs_results.tsv.gz`。
- **口径**（2026-11-04 起）：只用 3 个 CT 都检验过的 14,890 个候选 pair（各 CT 检验的 pair：Bn 40,087 / CD4 T cell 85,527 / CD8 T cell 50,931），避免“另一个 CT 没检验”造成的假 unique。每个方法在某 CT 中 p.adj（gene 内 BH，各 CT 自己的 MR 表）< 0.1 且属于共有 pair 的 link 作为该 CT 的 call set，用 `peak_id + gene_id` 匹配。区域内写 count 和（占三者并集的 %）；圆外标签写 `共有 pair 中的显著数 / 该 CT 全部显著数`（下一行为占比），例如 IVW Bn `1,562 / 3,977 (39.3%)`。
- **输出**：`Figures_20261102/Figure_S4.pdf`（只出 PDF，10×8 in）+ `Figure_S4_counts.tsv`（长表：`region_<Bn CD4 CD8 位>` 的 n / 并集 / %（共有 pair 口径），`common_<CT>` 的 n_sig_common / n_sig_all / %，`unique_<CT>` 的 n / n_sig_common / %）。
- **版式**：2×3 panel（A–F，按 `METHODS` 顺序），每个 panel 一个固定大小（不按比例）的三圆 Venn；标题用方法配色（2026-11-04 起标题下不再写并集数，并集只写入 tsv）；圆用低饱和的清新配色（2026-11-04 起：Bn 雾蓝 `#9DBEDF` / CD4 T cell 鼠尾草绿 `#A9D3A2` / CD8 T cell 杏色 `#F3C29B`，alpha 0.45，边用同色加深；与 Figure S1 的 Dark2 不同），圆内数字加粗、字号 3.5；圆外标签（Bn 左上、CD4 T cell 右上、CD8 T cell 正下，两行，用该 CT 的加深色）；细胞类型名称靠页底共用图例（Bn / CD4 T cell / CD8 T cell）；无页顶说明（2026-11-04 起）。
- **写 caption 时注意**：圆内只统计共有 pair，所以 unique 是真正“检验过但只在本 CT 显著”；圆外占比说明各 CT 的显著 link 有多少落在共有 pair 中（18.8–42.0%）。WNN 用的是外部 PBMC metacell，不分 CT，共有 pair 上三个 CT 的 WNN 原始 p 相同，所以 Marginal Multiomics 三者共享 91.5%（6,816 / 7,451）；其余方法三者共享 8.8–22.1%。

### 3.14 Fig13 — MR-Only link 的金标准命中：观测 vs 背景（分箱校准图）✅（论文当前 S5 Fig）
- **脚本**：`Plot13_MROnly_ObsExp.R`（登录节点，几秒）；输入 `Results/Validation/Valid1_Hyper/MROnly_gene_overlap_combined.tsv`（`Scripts/Validation/Valid1_Hyper_MROnly.R` 生成，口径见 `Pipeline.md` 3.1 “MR-Only 合并扩展”）。
- **口径**：每个 gene 一行（3 个 CT 合并；HiChIP / PCHiC 用 3 个 CT loop 的并集）；x = K / N（全部候选 peak 的背景命中率），y = k / n_sig（MR-Only peak 的命中率）；只用 K > 0、N > K 的 gene（同 Fig3 / Table S1–S2）。
- **为什么不画示例那样的逐 gene 2D 密度**：MR-Only 在 gene 层很稀疏，55% 的 gene 只有 1 个 MR-Only peak（n_sig 中位 = 1），K/N 中位只有 0.01–0.08，所以单个 gene 的 y 几乎只能取 0 或 1。逐 gene 画会得到 y = 0 的一条浓带（58–94% 的点在对角线以下），这和整体 1.9–3.3 倍的富集正好相反；y 是离散值，2D KDE 也不成立。2026-11-04 与用户讨论后改为分箱校准图。
- **版式**：2×3 panel（A–F：GTEx、EnhancerAtlas、FANTOM5、4DGenome、HiChIP、PCHiC），无页顶标题和说明。每个库按 K/N **等数量**分 8 箱（不用等距：K/N 严重右偏，等距时约 70% 的 gene 挤在第 1 箱、右侧几箱只有个位数 gene 甚至为空；2026-11-04 与用户确认），每箱画 (ΣK/ΣN, Σk/Σn_sig)，带 95% CI（按 gene 的 bootstrap，见下）；小点 = 落在 panel 范围内的 gene（y 很高的如 y = 1 不画，y = 0 的 gene 在底部形成一条横带）；虚线 y = x，对角线以上极浅底色。每个 panel 的 x、y 用同一范围（按数据放大），两轴用带符号的平方根变换（同一变换，`coord_equal`，对角线 45°；两端各留 3% 空白）。轴名只在最下一行 / 最左一列显示。**2026-11-04 起与 Figure 3D 同一风格**（`panel(fig3 = TRUE)`）：薄荷色板（分箱深蓝 `#2E4E77`、点 size 2.9、误差棒 0.7、连线 0.75；gene 点 `#9ACDA2` alpha 0.35；富集侧 `#F5FAF3`；对角线灰），浅薄荷标签条写库名（代替 panel 标题），边框 / 刻度加粗；不再写副标题，右下角对角线以下用加粗文字（白底）写 `Fold = 2.17` / `(P = 2.8e-73)` 两行（pooled fold 与合并计数的超几何 P，口径同 Table S1/S2；Figure 3D 写成一行、字号更大）；gene 数只写 summary tsv。范围大的 panel（lim > 0.3，即 HiChIP）去掉 0.05 刻度。tag 画在左上角不占行。命中 gene 的观测 vs 期望 %（期望 = 各 gene 的 1 − P(k = 0)，超几何精确计算后取平均）只写 summary tsv。
- **CI**（2026-11-04 起）：每箱内按 gene 有放回重抽（bootstrap，B = 2000，`set.seed(20261104)`，百分位法），每次重算 Σk/Σn_sig；这样同一 gene 内 peak 的相关性和 gene 之间的差异都计入。原来的精确二项（Clopper–Pearson）CI 把每个 peak 当独立试验，仍写在 `Figure_S5_bins.tsv`（`ci_*_binom`）作对照；两者宽度比 0.75–1.31，差别很小。bootstrap CI 下限 ≤ 背景比的箱为 9 / 48（精确二项为 8 / 48），都在命中数少的低背景箱。
- **输出**：`Figures_20261102/Figure_S5.pdf`（只出 PDF，10×7.2 in）+ `Figure_S5_summary.tsv`（每个库：n、n_sig = 1 的比例、命中 gene 的观测 / 期望 %、对角线以上 %、ΣN / ΣK / Σn_sig / Σk、fold、log10 P）+ `Figure_S5_bins.tsv`（每个库 8 箱：背景比、观测比、CI）。另存 `tmp/Plot13_MROnly/hichip_fig3.rds`：HiChIP panel 的 Figure 3 版（`panel(fig3 = TRUE, strip = …)`：薄荷色板、标签条、线条加粗、`coord_cartesian`、无副标题数字），供正文 Figure 3D 使用（见 3.10）；S5 不受影响。
- **结果（2026-11-04）**：pooled fold：GTEx 2.17、EnhancerAtlas 2.35、FANTOM5 2.84、4DGenome 2.87、HiChIP 3.28、PCHiC 1.89（log10 P 在 −28 到 −619 之间）；48 个箱的点估计全部在对角线以上；CI 下限跨过对角线的箱都在背景比最低、命中数少的位置（精确二项 8 个，bootstrap 9 个）。

### 3.15 Fig14 — 5 个方法 effect 的两两相关性与密度矩阵 ✅（论文暂定 S6 Fig）
- **脚本**：`Plot14_MethodCorr.R`（登录节点，约 1.5 分钟）；输入 `Results/MR/<CT>_AllPairs_results.tsv.gz`。思路来自旧 `Server6_wgen1_plot.sh` 的 `pairs()` 图（旧图实际画的是 signed −log10 P，截断 ±30）。
- **口径**：5 个方法（Causal SMR / IVW / Egger，Marginal Multiomics / Pseudobulk；CCT 没有 effect，不画），每个 CT 用全部候选 pair（Bn 40,087 / CD4 T cell 85,527 / CD8 T cell 50,931）。两种刻度（每页 panel A / B）：
  - **A = effect**：SMR `smr_b`、IVW `ivw_beta`、Egger `egger_beta`、WNN `wnn_rho`、Pseudobulk `pseudo_rho`。每个方法一个对称轴范围 = 3 个 CT 中 |x| 的 99% 分位的最大值（向上取到 1 / 1.5 / 2 / 2.5 …）：SMR、IVW ±1.5，Egger ±2.5，WNN ±0.75，Pseudobulk ±0.35。
  - **B = signed −log10 P**：`sign(effect) × −log10(P)`。所有方法共用 ±11（合并后 |x| 的 99% 分位）。
  - 超出范围的 pair 不进密度（直接舍去，不压到边上，避免边缘条带），显示比例为 98.1–99.8%（写入 tsv）。上三角 r 是全部 pair 的 Spearman 相关，不截断。
- **版式**（2026-11-04 起）：一个 PDF，每个 CT 一页（Bn / CD4 T cell / CD8 T cell），8.5×17.6 in；页顶只写细胞类型（不写 pair 数）；每页两个 5×5 矩阵上下排（补充材料页面不够宽）：A Effect estimate 在上，B Signed −log10(P) 在下，panel 名放在矩阵上方单独一行。
  - 下三角：`MASS::kde2d` 二维密度（160×160，带宽 1.5 × `bandwidth.nrd`），填色为 density^0.25（同 smoothScatter），蓝色系 白 → `#DCEAF5` → `#7FB2D9` → `#2F6FA8`；x = 0、y = 0 用灰点线；y = x 用柔和橙红虚线 `#E08A6B`（A 只在两轴都是 MR beta 时画，B 每格都画）。
  - 上三角：`r = 0.87` 加粗，字号随 |r| 变大，背景按 r 用低饱和发散色（浅红 ← 白 → 浅蓝）。
  - 对角线：方法名（两行，方法配色）+ 该方法自身数值（A 为 effect，B 为 signed −log10 P）在全部 pair 中的分布曲线（`density()`，同色浅填充；高度缩放到格子下半部，只看形状，纵轴无读数）；横轴用本方法的范围，与同列密度图共用刻度。
  - 刻度只在最下一行 / 最左一列，不放两端刻度（避免相邻格子的刻度字相撞）；左上角（SMR 对角线）格子不画纵轴刻度（第一行没有密度图，刻度易被误读为曲线读数）。
- **输出**：`Figures_20261102/Figure_S6.pdf`（3 页，只出 PDF；原单独的 `Figure_S7.pdf` 已合并进来并删除）+ `Figure_S6_cor.tsv`（`type = cor`：刻度 × CT × 方法对的 r、n；`type = range`：每个方法的轴范围、显示比例）。
- **结果（2026-11-04）**：effect 的 Spearman r 在 3 个 CT 中很接近：SMR–IVW 0.86–0.87，IVW–Egger 0.55–0.56，SMR–Egger 0.39–0.40；MR 与 Pseudobulk 0.16–0.30；WNN 与其余方法都 ≤ 0.07。B（signed −log10 P）的 r：SMR–IVW 0.85，IVW–Egger 0.70–0.71，其余与 effect 相近。

### 3.16 Fig15 — 全基因组 Manhattan / Miami 图 ✅（论文暂定 S7 Fig + Figure 3A）
- **脚本**：`Plot15_Manhattan.R`（登录节点，约 1 分钟）；输入 `Results/MR/<CT>_AllPairs_results.tsv.gz`、`Example_pairs.tsv`（`pair_set = Figure_2`）。
- **口径**：6 个方法，每个 CT 用全部候选 pair。x = peak 中心，按 hg38 chr1–22 长度累加；y = −log10 P，线性等距（2026-11-04 起不再分段压缩），每条轨道用自己的范围（各方法量级差很多，例如 CD4 T cell 中 Joint 最大 189，Multiomics 31）。阈值线 = Bonferroni（0.05 / 该 CT pair 数，−log10 为 5.9 / 6.2 / 6.0）；FDR 是 gene 内 BH，没有统一阈值线。
- **配色**：不按方法分色。过线的点 causal 用亮橙 `#FF9F1C`，marginal 用亮绿 `#3BB273`；其余按奇偶染色体交替浅灰 / 中灰；点用 `ggrastr` 栅格化（300 dpi），点较大。
- **标注**：过线 pair 中每个 gene 取最强 pair，按 P 贪心挑选，同一染色体 5 Mb 内只标一个（Miami 每半边 10 个，第 2–4 页每轨道 8 个）。Figure 2 的例子 pair 在其所在 CT 的页面上每条轨道都标（不显著也标在实际位置）：洋红菱形 + 浅洋红底色的粗体标签；与 top 基因标签在同一个 `geom_label_repel` 图层，互不遮挡。
- **版式**：`Figure_S7.pdf` 3 页（8.5 × 11 in），每个 CT 一页，6 条轨道（SMR / IVW / Egger / Joint / Multiomics / Pseudobulk）；Bn 页标 SIDT2，CD4 T cell 页标 ERAP2 / ARPC1B / RPL11，CD8 T cell 页没有例子。轨道名写在 y 轴标题（方法名 + −log10 P）；页标题只写细胞类型；图中不写 pair 数等数字。CD4 T cell 的 Miami 图（上半 Causal joint 向上，下半 Marginal multiomics 镜像向下：画 y = −(−log10 P)，刻度写绝对值；下半高度为上半的 0.75（2026-11-04 第 2 轮由 0.55 加高，放得下竖排标题），两半之间无空白，中间一条染色体编号；例子 ERAP2 / ARPC1B / RPL11 两半都标）原为 S7 第 1 页，2026-11-04 起移入正文 Figure 3A：本脚本只构建并存为 `tmp/Plot15_Manhattan/miami_CD4.rds`，由 `Plot9_GO_GRN.R` 拼图（见 3.10）。Miami 的配色和标签条通过 `manhattan_track()` 的 `cols` / `strip` 参数传入（S7 用默认配色，不受影响）。两条轨道的 y 轴标题不用 ggplot 的 y 标题，改为左侧单独一列文字格（两个标题在同一竖线上）；**两条轨道都必须是 facet 图**（下半用 `strip = ""`，facet 但不显示标签条）：若只有上半是 facet 图，patchwork 加宽共用的 y 轴列时，非 facet 图的 y 轴线会画在该列最左侧，与上半错开。
- **输出**：`Figures_20261102/Figure_S7.pdf`（只出 PDF）+ `Figure_S7_data.tsv`（一个文件；label 部分用 `panel` 列区分 `Figure_3A` / `S7_p1`–`S7_p3`：`type = count` 为每个 CT × 方法的 pair 数、Bonferroni、过线 pair / gene 数、FDR < 0.1 数、最大 −log10 P；`type = label` 为标注的基因：page、CT、方法、rank、gene、peak、−log10 P、FDR、is_example）。
- **结果（2026-11-04）**：过 Bonferroni 的 pair 数（SMR / IVW / Egger / Joint / Multiomics / Pseudobulk）：Bn 680 / 526 / 347 / 749 / 3,650 / 403；CD4 T cell 1,976 / 1,606 / 1,038 / 2,167 / 5,946 / 1,630；CD8 T cell 1,069 / 885 / 560 / 1,179 / 4,432 / 836。例子在 Joint 中很强（SIDT2 171、ERAP2 133、ARPC1B 47、RPL11 15），在 Multiomics 中弱（8.0 / 3.2 / 0.7 / 3.6）；Pseudobulk 中 SIDT2 72.7、ERAP2 43.7，ARPC1B、RPL11 均为 0.1。Multiomics 过线 pair 多但都不很强（最大 26–31），top 基因（NOSIP、PARVB、CD7、ETS1）与 MR 的（ACYP2、LCLAT1、ERAP2、TMEM8A）基本不重合。

### 3.17 Table 5 — S-LDSC / Fig 4B 的 GWAS 性状汇总 ✅（论文暂定 S4 Table）
- **脚本**：`Table5_GWASTraits.R`（登录节点，约 1 分钟）；输入：`selected_sumstats/*.sumstats.gz`（7 个 LDSC 输入文件）、`GWAS_data/onk1k_GWAS/{AS,CD}.log`（cases / controls）、`Valid2_SLDSC/SLDSC_all.tsv`（Combined pool 的 h2 和 SE）、`Valid2_SLDSC/h2/Combined/Combined_pool_<trait>.log`（SNP 数核对）、`Valid6_GWAS/GWAS_MRonly_overlap.tsv` + `GWAS_pair_counts.tsv`（Fig 4B 的 link 统计和核对）、`Plot10_Combined.R` 的 `SL_TRAITS` / `PAGES`（读出展示的性状和 5 个例子）。
- **口径**：只列 Fig 4A 展示的 7 个性状（AS、CD、IBD、RA、SLE、MS、T1D；Asthma、PBC 虽跑了 S-LDSC 但不展示、不入表）。N 和 SNP 数只写 LDSC 输入文件里的数（不引用论文的样本量）：AS、CD 的 sumstats 有约 3 万行 Z / N 缺失，LDSC 读入时丢弃，SNP 数按完整行统计并用 `stopifnot` 与 S-LDSC 日志的 “Read summary statistics for X SNPs” 核对；T1D 的 N 随 SNP 变化，写中位数和范围。来源和 PMID 是人工整理（PubMed 核对）：AS = Pan-UK Biobank，CD = IIBDGC 网站，MS = IMSGC 2011（PMID 21833088，ImmunoBase 受保护数据），IBD / RA / SLE / T1D = Zenodo 10515792。Fig 4B 的 link 统计仍在脚本里算并核对，但 2026-11-05 起不输出（论文 S4 Table 去掉了 “Fig 4B: MR-only links” 列）；表只留来源、N、SNP 数、h²、是否在 Fig 4A。
- **输出**：`Figures_20261102/Table_S4.xlsx`（单 sheet，标题 “S4 Table …”，含脚注）+ `Table_S4.tsv`。
- **结果（2026-11-05）**：h2（obs）：AS 0.0052、CD 0.0649、IBD 0.3202、RA 0.0132、SLE 0.4567、MS 0.0800、T1D 0.0397。

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
- 2026-11-03：新增 Figure S4（`Plot_Sup4_Venn.R`，同日改名 `Plot12_Venn.R`）：6 个方法各一个 Bn / CD4 T cell / CD8 T cell 三圆 Venn 图（count + 占并集 %，圆外标签写 n 和 % unique），见 3.13。
- 2026-11-03：画图脚本统一按 `Plot<N>_` 顺序编号，与论文图号脱钩：`Plot_Sup1_OneSampleMR.R` → `Plot11_OneSampleMR.R`，`Plot_Sup4_Venn.R` → `Plot12_Venn.R`；第 3 节表格和小节标题写明“论文当前 … Fig”。
- 2026-11-03：新增 Fig13（`Plot13_MROnly_ObsExp.R`，论文当前 S5 Fig）：MR-Only link 的金标准命中率观测 vs 背景（3 CT 合并，6 个库），见 3.14。
- 2026-11-04：Fig13 由逐 gene 2D 密度改为分箱校准图（8 箱，精确 95% CI，平方根坐标；原因是 MR-Only 在 gene 层太稀疏，见 3.14），新增 `Figure_S5_bins.tsv`。
- 2026-11-04：Fig12（S4）去掉页顶说明和圆外标签，panel 标题下加并集数，页底加细胞类型图例；Fig13（S5）去掉页顶标题和说明，panel 副标题改为一行 `genes | Fold (P)`，灰点和分箱点加深加粗，坐标两端留白；分箱保持等数量。
- 2026-11-04：Fig12（S4）圆内数字加粗加大；细胞类型改为低饱和配色（雾蓝 / 鼠尾草绿 / 杏色），让黑字更清楚。
- 2026-11-04：Fig13（S5）CI 改为按 gene 的 bootstrap（B = 2000，百分位法；精确二项 CI 留在 bins tsv 作对照）；去掉 “enriched” 文字；分箱点、误差棒、连线调细；配色改为同一蓝色系（分箱深蓝 `#1F4E79`、gene 淡蓝 `#8FB3D3`、对角线柔和红 `#C44E52`）。
- 2026-11-04：Fig12（S4）改为只统计 3 个 CT 都检验过的 14,890 个候选 pair；圆外加 `共有 pair 中的显著数 / 全部显著数 (占比)` 标签。
- 2026-11-04：Fig12（S4）去掉 panel 标题下的并集数；S4、S5 已放入论文的 `supplementary.tex`（caption 见 `Manuscript/notes/REFERENCE.md` 第 1 节）。
- 2026-11-04：新增 Fig14（`Plot14_MethodCorr.R`，论文暂定 S6 / S7 Fig）：5 个方法两两的 effect（S6）和 signed −log10 P（S7）密度 + Spearman r 矩阵，每个 CT 一页，见 3.15。
- 2026-11-04：Fig14 的 S6 / S7 合并为一个 `Figure_S6.pdf`：每个 CT 一页，A effect 在上、B signed −log10 P 在下；去掉 pair 数副标题和左上角格子的纵轴刻度。
- 2026-11-04：新增 Fig15（`Plot15_Manhattan.R`，论文暂定 S7 Fig）：MR 方法的全基因组 Manhattan 图，第 1 页 Causal Joint × 3 CT（高亮 Figure 2 例子），第 2–4 页每个 CT 4 个方法，见 3.16。
- 2026-11-04：Fig15（S7）第 1 页只画 Causal Joint 的 CD4 T cell（8.5 × 4.6 in），例子 pair 改用洋红色高亮；各页散点放大。
- 2026-11-04：Fig15（S7）图中去掉 pair 数和过线数（只保留在 tsv）；两个 tsv 合并为 `Figure_S7_data.tsv`。
- 2026-11-04：Fig15（S7）y 轴改为线性（每条轨道自己的范围）；第 1 页改为 CD4 T cell 的 Miami 图（Joint 向上 / Multiomics 向下）；第 2–4 页每页 6 个方法；配色改为 causal 亮橙、marginal 亮绿；Figure 2 例子在所在 CT 页每条轨道都标。
- 2026-11-04：Fig15（S7）第 1 页 Miami 图压缩 y 轴：页高 4.6 in，下半高度为上半的 0.55，两半之间去掉空白。
- 2026-11-04：Figure 3 加入 panel A = CD4 T cell 的 Miami 图（`Plot15` 存 rds，`Plot9_GO_GRN.R` 读入），原 A–C 顺延为 B–D；GO 配色改为暖色（深红橙 → 浅杏），GRN 改为 NPG 藏蓝系；Fig15（S7）去掉 Miami 页，只留 3 页。
- 2026-11-04：Figure 3 改版为三行：A Miami 图（加高）/ B、C GO（压缩行高）/ D HiChIP 校准图（`Plot13` 存暖色版 rds，S5 不变）+ E GRN odds ratio（只画 CD4 T cell）；Bn / CD8 T cell 的 GRN 结果移到新的 Figure S8（`Plot9_GO_GRN.R` 一并输出）。
- 2026-11-04：Figure 3 微调（见 3.10）：版式改为 A 通栏 + 2×2 网格（左列 B GO BP / C GO MF，右列 D HiChIP / E GRN，同行绘图区对齐）；去掉所有 title / subtitle，统一用浅薄荷色标签条写简洁标题（括号内细胞类型 / 数据库）；整体换为薄荷色板 `PAL_MINT`（S8 随 E 换色，S5 / S7 不变）；修正 Miami 两半 y 轴错位（两条轨道都做成 facet 图，y 标题放在左侧单独一列），方法名改为论文写法（Causal joint / Marginal multiomics）。
- 2026-11-04：`Scripts/Visualization.md` 在 21:52 被一次外部 “Workspace Edit” 覆盖成上午的旧版本；已用 Claude Code 的 file-history 备份（当日 15:58 版）加上之后的修改恢复，旧文件留在 `tmp/Visualization.md.stale_overwritten_20261104_2152`。
- 2026-11-04：Figure S5 改为 Figure 3D 的风格（薄荷色板、库名标签条、加粗边框、右下角加粗 Fold / P，去掉副标题和 gene 数）；Figure 3D 右下角加 Fold / P。Figure 3、S5、S7、S8 已更新到 `Manuscript/`（见 `Manuscript/notes/REFERENCE.md`）。
- 2026-11-05：新增 `Scripts/Visualization/Table5_GWASTraits.R`：S-LDSC（Fig 4A）和 Fig 4B 用到的 7 个 GWAS 性状汇总表，输出 `Figures_20261102/Table_S4.xlsx` / `Table_S4.tsv`（`Scripts/Visualization.md` 3.17）；尚未写入 `Manuscript/`。
