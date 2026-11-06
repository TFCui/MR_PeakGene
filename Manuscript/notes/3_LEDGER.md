# 3_LEDGER：台账（数字、图表、投稿、同步、更新记录）

润色时需要数字就查第 1 节，其余需要时再读。路径相对项目根 `Claude_project/`。正文里的数字以本节和 `Results/` 为准；发现不一致只报告，不改。

## 1. 已核对的数字

2026-11-01 已对照结果表或代码核对（R1–R5 证据卡、A1 审计），2026-11-04 补了新图的数字。阈值都指 $P_{\mathrm{adj}}$（每个 gene 内部、每个方法分别做 BH）。

### 1.1 数据与候选 pair（Methods、R1）

| 项 | Bn | CD4 | CD8 | 来源 |
|---|---|---|---|---|
| 候选 pair（最终分析集） | 40,087 | 85,527 | 50,931 | `Results/MR/<CT>_AllPairs_results.tsv.gz` |
| 任一方法 $P_{\mathrm{adj}}<0.1$ | 22,766 | 50,704 | 31,480 | 同上 |
| pseudobulk 个体数（两个模态都通过 QC） | 335 | 372 | 346 | `Data/.../QC_pseudobulk/` |
| 工具变量数 J ≤ 5 的 pair | 11.0% | 10.2% | 12.7% | `ivw_n_snp` |
| J ≤ 10 的 pair（全部 / MR 任一显著的 link） | 26.5% / 15.0% | 25.9% / 14.2% | 28.7% / 16.3% | 同上 |

- 每个 pair 的 J：中位数 22（IQR 约 10–40）。3 个细胞类型都检验过的 pair：14,890（S4 Fig）。
- CIMA：428 名健康中国成年人；全部细胞的 count 矩阵为 RNA 6,484,974 个细胞 × 36,326 个基因，ATAC 3,762,242 个细胞 × 338,036 个 peak；pseudobulk 个体过滤：细胞数 ≥ 100 且 ≤ median + 3·MAD。
- 工具变量：caQTL 名义 $P<10^{-8}$；逐 pair clumping，$r^2$ 0.99，窗口 20,000 kb；至少 3 个共享变异。
- LD 参考：1000 Genomes GRCh38（3,202 人）中的 EAS 子集，n = 585；PLINK v1.90b7.2。GLS 截断：保留 99.9% 方差、至少 5 个方向、相对特征值阈值 $10^{-10}$；只保留 $\hat\phi_{\mathrm{IVW}} \le 10$ 的 pair。
- marginal multiomics：外部 PBMC multiome 的 80 个 WNN metacell（`wsnn_res.15`，不重新聚类），在 338,036 个 CIMA peak 上重新计数；3 个细胞类型共用一个矩阵。

### 1.2 R1：发现数、重叠和方向（Fig 1、S3 Fig）

$P_{\mathrm{adj}}<0.1$ 的 link 数（IVW / Egger / SMR / joint / multiomics / pseudobulk）：

| | IVW | Egger | SMR | joint | multiomics | pseudobulk |
|---|---|---|---|---|---|---|
| Bn | 3,977 | 3,361 | 3,534 | 4,406 | 18,179 | 4,924 |
| CD4 | 9,569 | 7,265 | 8,627 | 10,597 | 37,181 | 17,181 |
| CD8 | 5,722 | 4,538 | 5,220 | 6,451 | 24,114 | 9,618 |

- 六个方法都显著：469 / 1,259 / 759；三个 MR 方法同号：469/469、1,258/1,259、758/759。
- CD4 的四个扇区（IVW 与两个 marginal 的符号）：都一致 918（72.9%）/ 只与 multiomics 不一致 304（24.1%）/ 只与 pseudobulk 不一致 17（1.4%）/ 都不一致 20（1.6%）。与 multiomics 不一致合计 25.7%，与 pseudobulk 3.0%；Bn 24.2% vs 1.0%，CD8 24.8% vs 2.1%。

### 1.3 R2：四个例子（Fig 2；`Scripts/Visualization/Example_pairs.tsv`）

| 例子 | peak | IVW β（$P_{\mathrm{adj}}$） | Egger / SMR β | multiomics ρ（$P_{\mathrm{adj}}$） | pseudobulk ρ（$P_{\mathrm{adj}}$） |
|---|---|---|---|---|---|
| *SIDT2*（Bn） | chr11:117185715-117186216 | 1.024（1.6e-171），194 个 SNP | 1.026 / 1.013 | 0.586 | 0.792（2.0e-72） |
| *RPL11*（CD4） | chr1:23691494-23691995 | 1.127（6.0e-11） | — | 0.400（3.9e-4） | 0.010（0.92） |
| *ERAP2*（CD4） | chr5:96885894-96886395 | −1.525（4.2e-133） | −1.510 / −1.470 | +0.375（P 6.1e-4，$P_{\mathrm{adj}}$ 1.5e-3） | −0.641（9.9e-44） |
| *ARPC1B*（CD4） | chr7:99376204-99376705 | 0.740（1.8e-47） | 0.719 / 0.784 | −0.149（0.37） | 0.016（0.79） |

SIDT2 的 P 值排序：IVW < joint < Egger < SMR。结果表里没有 Egger 截距。

### 1.4 R3：GO、参考数据集、GRN（Fig 3、S1 / S2 Table、S5 Fig、S8 Fig）

- **GO（DAVID，2026-10-19）**：MR-only 基因 4,509 个（3 个细胞类型合并；`Results/Validation/Valid7_GO/DAVID_input/MROnly_genes_all.txt`），有 BP 注释 4,031、MF 注释 4,185；背景是 DAVID 默认的全基因组。DNA-templated transcription：689 个基因，fold 1.36，FDR 5.18e-18。protein binding fold 1.16，FDR 2.81e-85；RNA binding fold 1.41，FDR 3.58e-15；zinc ion binding fold 1.33，FDR 3.58e-15。
- **超几何（S1 / S2 Table）**：所有方法 × 参考 × 细胞类型都富集，P ≤ 3.6e-17，fold 1.68–6.18。causal SMR 在 17 / 18 个组合中 fold 最大（例外：CD4 PCHi-C，IVW 2.61 vs SMR 2.55）；marginal multiomics 在 18 个组合中都最小。Bn：FANTOM5 SMR 6.18 vs multiomics 3.56，HiChIP 4.46 vs 3.08。
- **MR-only 与参考（S5 Fig、Fig 3D）**：约 55% 的 gene 只有 1 个 MR-only peak；分箱 95% CI 来自按 gene 的 bootstrap（B = 2,000）。数值见 FIG `Figure_S5_summary.tsv`、`Figure_S5_bins.tsv`。
- **GRN（Fig 3E、S8 Fig）**：CIMA SCENIC+ 1.0.1，237 个 high-quality eRegulon。被检验的 pair：cell-type eRegulon 864 / 2,520 / 1,364，high-quality 817 / 458 / 382；HQ eRegulon 个数 Bn 29、CD4 2、CD8 8。阈值 $P_{\mathrm{adj}}<0.05$，Fisher 检验，所有 P < 1e-30。OR（SMR / IVW / joint / Egger / pseudobulk / multiomics）：

| | SMR | IVW | joint | Egger | pseudobulk | multiomics |
|---|---|---|---|---|---|---|
| Bn cell-type | 9.60 | 7.33 | 7.27 | 5.11 | 5.09 | 2.30 |
| CD4 cell-type | 7.21 | 5.83 | 6.13 | 4.86 | 3.15 | 2.77 |
| CD8 cell-type | 7.48 | 5.61 | 6.01 | 4.23 | 3.35 | 3.57 |
| Bn HQ | 9.52 | 7.18 | 7.13 | 5.04 | 5.12 | 2.30 |
| CD4 HQ | 10.75 | 8.28 | 9.19 | 6.70 | 4.79 | 5.24 |
| CD8 HQ | 11.60 | 7.83 | 8.74 | 5.46 | 4.74 | 6.51 |

Egger 的例外：Bn HQ 5.04 < pseudobulk 5.12；CD8 HQ 5.46 < multiomics 6.51。

### 1.5 R4：S-LDSC 与 GWAS 重叠（Fig 4）

- **S-LDSC**：MR-only peak 3,777 个（hg19；hg38 下 3,818 个）；对照 = 全部候选 peak 并集，hg19 下 316,798 个。覆盖的 SNP 比例：MR-only 0.25%，全部 peak 15.1%。baselineLD v2.2（97 个注释），两个注释分别拟合；EUR LD 参考。
  enrichment（MR-only vs 全部 peak）：AS 17.7 vs 3.05，CD 20.0 vs 6.36，IBD 18.0 vs 5.12，RA 10.1 vs 5.73，SLE 14.0 vs 6.49，MS 50.8 vs 26.0，T1D 19.7 vs 4.64。MR-only 的 P：CD 0.084，IBD 0.050，MS 0.069，T1D 0.060（AS 0.30，RA 0.25，SLE 0.25）。全部 peak 的 P 除 AS（0.20）外都 < 1e-7。
- **GWAS 重叠**（strong MR-only link + peak 内有 GWAS P < 5e-8 的 SNP；Fig 4B 中红色 link 的阈值是 0.05，`Plot10_Combined.R:40`）：
  - *ERBB3*（CD4，chr12:56041553-56042054）：rs705704，T1D P 1.5e-62；IVW $P_{\mathrm{adj}}$ 1.4e-4
  - *RGL2*（Bn，chr6:33313910-33314411）：rs3130099 / rs3106189，T1D P 7.5e-54，RA P 7.1e-15
  - *POLR2E*（CD4，chr19:1095225-1095726）：rs1046040 / rs11084882，CD P 1.5e-8
  - ~~*UBAC2*（CD4，chr13:99375239-99375740）：Asthma P 2.8e-9~~（2026-11-06 去掉：只有 Asthma，不在 7 个展示的 GWAS 内）
  - *ERAP1*（Bn，chr5:96896259-96896760，2026-11-06 新加）：rs2549785，CD P 4.8e-11；SMR / IVW / Egger FDR 2.6e-5 / 8.0e-6 / 4.8e-4，WNN / Pseudobulk 不显著（0.19 / 0.96）
  - *GNAS*（CD8，chr20:59254001-59254502）：rs259956，IBD P 2.6e-9，CD P 8.4e-9

### 1.6 R5：motif（Table 1）

- JASPAR2020 CORE vertebrates 非冗余集，746 个 motif；背景是 50,000 个按 GC 和可及性匹配的 peak（universe：Bn 316,117 个）。
- MR-specific motif 47 个（Bn 23、CD4 14、CD8 10）。Table 1 的 6 个都在 Egger 中富集，其中 FOSL2::JUND、MNT、POU3F1、TEF 只由 Egger 选出。
- MNT（CD4）：Egger $P_{\mathrm{adj}}$ 2.9e-3，fold 1.34（5.07% vs 3.79%）。FOSL2::JUND（Bn）：1.1e-3，fold 1.57。ATF3（Bn，最强）：1.5e-4，fold 1.75；在 multiomics 中 $P_{\mathrm{adj}}$ 0.072，在 pseudobulk 中 0.056。CENPB（CD4）：3.6e-3，fold 1.55。POU3F1（CD8）：6.3e-3，fold 1.38。TEF（CD8）：0.010，fold 1.74。

### 1.7 QC（S1 Fig、S2 Fig）

- **样本重叠诊断（S1 Fig，`Results/QC/1_OneSampleMR/`）**：null SNP 是两侧 P > 0.05 的变异；全局 ρ：Bn 0.0035、CD4 0.0060、CD8 0.0053；P 分别为 2.3e-48、1.8e-262、5.6e-134。观测数 17,074,938 / 33,263,695 / 21,537,822（同一 SNP 出现在多个 pair 中，观测不独立）。pair 数 547,452 / 1,051,953 / 686,829。CD4 逐 pair：平均 ρ 0.004–0.007；ρ 的 SD 从 0.54（3–5 个 SNP）降到 0.15（> 50 个 SNP）。
- **等位方向扫描（S2 Fig，`Results/QC/2_FlipPoints/`）**：共扫描 176,545 个 pair；不 flip 时 $P_{\mathrm{adj}}<0.05$、flip 后 ≥ 0.1 的 pair 有 700 / 1,763 / 1,056 个（其中 SNP ≥ 50 的 187 / 507 / 305 个）；例：*PIGX*（Bn，136 个 SNP，P 1.3e-11 → 0.175）。

### 1.8 软件版本（2026-11-01 读自服务器当前安装，不一定是运行时版本）

R 4.4.2；PLINK v1.90b7.2；Signac 1.14.0；JASPAR2020 0.99.10；TFBSTools 1.44.0；motifmatchr 1.28.0；BSgenome.Hsapiens.UCSC.hg38 1.4.5；rtracklayer 1.66.0；LDSC v1.0.1；Seurat 5.1.0；data.table 1.16.2；GenomicRanges 1.58.0；EnsDb.Hsapiens.v86 2.99.0；ieugwasr 1.1.0；TwoSampleMR 0.7.1；scanpy 1.10.4（pseudobulk 聚合）；SCENIC+ 1.0.1（CIMA）。

## 2. 图表台账

**FIG** = `Results/Figures/Figures_20261102/`。hash 取 SHA-256 的前 12 位，刷新命令：`sha256sum Manuscript/figures/*.pdf`。2026-11-01 核对过：论文里的 7 个 PDF 与 FIG 中的同名文件逐字节一致。图只出 PDF（`Figure_1B_C` 只出 SVG，供 Figma 拼图用）；图中的细胞类型标签用 `Scripts/Visualization/00_common.R` 的 `CELL_LAB`。

| 条目 | 论文内文件 / label | 来源与脚本 | hash | 内容与 caption 备注 |
|---|---|---|---|---|
| Fig 1 | `figures/Figure_1.pdf` / `fig:overview` | FIG `Figure_1_flat.pdf`。用户在 Figma 里把 `Fig0_Schematics.svg`（A）和 `Figure_1B_C.svg`（B–C，`Plot1_MRresults.R`）拼在一起，再用 Ghostscript 拍平（`tmp/Fig1_flatten/flatten_slurm.sh`，PDF 1.3，720 ppi） | `8c0d08b7fabf` | A：研究框架；B：CD4 的 UpSet（前 15 个组合）；C：1,259 个 link 的环形图，只标了 73% 和 24%。拍平后图内文字是位图（Figma 导出的 PDF 中带 SMask 的透明层在部分阅读器里显示为黑底）。Figma 原文件 hash `81411a569c5d` |
| Fig 2 | `figures/Figure_2.pdf` / `fig:examples` | `Plot2_examples_v2.R`；pair 列表 `Scripts/Visualization/Example_pairs.tsv`（`Figure_2`），每个 pair 都重新跑了 MR | `ae5b94e5c822` | 4 行：A SIDT2（Bn），B RPL11、C ERAP2、D ARPC1B（都是 CD4）。每行 4 列：QTL 轨道 / MR 散点 / WNN 相关 / pseudobulk 相关。caption 没写 A–D 对应的基因；受 6.5 in 宽度限制，图内文字偏小 |
| Fig 3 | `figures/Figure_3.pdf` / `fig:go` | `Plot9_GO_GRN.R` 拼图（13×14 in）：A 由 `Plot15_Manhattan.R`、D 由 `Plot13_MROnly_ObsExp.R` 存 rds 后读入；B/C 用 DAVID 结果；E 用 `Results/Validation/Valid5_GRN/` | `ae7df7ad7d96` | 2026-11-04 改版：A CD4 T cell Miami 图（上 causal joint / 下 marginal multiomics，Bonferroni 0.05/85,527，例子 RPL11 / ERAP2 / ARPC1B）；B GO BP、C GO MF（各前 10 个 term，4,031 / 4,185 个有注释基因）；D MR-only link 与 HiChIP 的校准图（同 S5 E，右下 Fold = 3.28）；E GRN OR 只画 CD4 T cells（2,520 / 458 个 pair），Bn / CD8 在 S8。各 panel 用浅薄荷标签条，薄荷色板。tex 中图宽 = `\linewidth`（不缩小）；caption 只在图页写标题和 A，B–E 用 `\ContinuedFloat` + `\caption*{Fig 3 (continued). …}` 接到下一页（2026-11-04 作者要求） |
| Fig 4 | `figures/Figure_4.pdf` / `fig:disease` | `Plot10_Combined.R`（`SL_TRAITS` 控制展示哪些 GWAS）；A 用 `Valid2_SLDSC/SLDSC_all.tsv`，B 用 `Valid6_GWAS/GWAS_MRonly_overlap.tsv` | `e5f1b6fa1da5` | A：7 个 GWAS 的 MR-only 与全部 peak 对比，`*` 表示 P < 0.1，不画 s.e.；图例在图内左上角；B：RGL2、ERBB3、POLR2E、ERAP1、GNAS（2026-11-06 UBAC2 → ERAP1，作图时去掉 Asthma 行），红色 link 的阈值是 0.05。2026-11-06 配色改 Set3（MR-Only 薄荷青 #8dd3c7 / All peaks 灰 / Gene 蓝 #80b1d3 / MR SNP 淡紫 / GWAS SNP 紫 #bc80bd / Peak 橙 #fdb462），字号约 ×1.3–1.5，A / B 字母 32 pt |
| Table 1 | `tables/motif.tex` + `figures/motif_logos/*.png` / `tab:motif` | `Valid3_Motif/Motif_MR_specific.tsv` + `Motif_table_main.tsv`（解释和 PMID 是人工整理的，只有这一个来源；2026-11-06 起 ATF3、FOSL2::JUND、MNT、TEF 四行的解释和 PMID 由脚本里的 `BRIEF_OVERRIDE` 覆盖为人类研究）→ `tmp/Manuscript_tables/make_motif_table.py`；完整 xlsx 在 FIG `Table_1.xlsx`（`537c47ae7b9f`） | — | 6 行，PMID 用 `\url`，不进 bib。PLOS 不允许表格里放图片，投稿前要处理 logo |
| S1 Fig | `Figure_S1.pdf` / `S1_Fig`、`fig:supp-qc` | `Plot11_OneSampleMR.R`；上游 `Scripts/QC/QC1_*` | `c0aa5df5b020` | null SNP 的 z 相关（A）和 P 值（B） |
| S2 Fig | `Figure_S2.pdf` / `S2_Fig`、`fig:supp-flip` | `Scripts/QC/QC2_2_FlipPoints_plot.R`；数值在 `Results/QC/2_FlipPoints/FigureS2_estimates.tsv` | `21cbfe159813` | A：flip 前后差别小的例子（SIDT2 / RPL11 / LCLAT1）；B：差别大的例子（PIGX / PDE4DIP / KLRC3） |
| S3 Fig | `Figure_S3.pdf` / `S3_Fig`、`fig:supp-bn-cd8` | `Plot1_MRresults.R` | `c8ba628007de` | Bn 和 CD8 的 UpSet + 环形图；caption 里写死了 "Fig 1B–C" |
| S4 Fig | `Figure_S4.pdf` / `S4_Fig`、`fig:supp-venn` | `Plot12_Venn.R`；计数在 FIG `Figure_S4_counts.tsv` | `52cae9934182` | 6 个方法各一个三圆 Venn，只统计 3 个细胞类型都检验过的 14,890 个 pair；圆内为数量（占并集 %），圆外为共有 pair 中的显著数 / 该细胞类型全部显著数（%）；caption 说明 WNN 不分细胞类型 |
| S5 Fig | `Figure_S5.pdf` / `S5_Fig`、`fig:supp-mronly-ref` | `Valid1_Hyper_MROnly.R` → `Plot13_MROnly_ObsExp.R`；数值在 FIG `Figure_S5_summary.tsv`、`Figure_S5_bins.tsv` | `2120a8356acf` | MR-only link（3 个细胞类型合并）的金标准命中率 vs 背景 K/N：按 K/N 等数量 8 组，95% CI 为按 gene 的 bootstrap。2026-11-04 起为 Fig 3D 的风格：库名标签条，右下角 pooled fold 和合并超几何 P（不再有副标题和 gene 数）；E 即 Fig 3D。caption 中引用表格用 `\ref{tab:supp-hyper-*}`（补充文档里没有 `S1_Table` 这类 label） |
| S6.1–S6.3 Fig | `Figure_S6.pdf`（3 页：Bn / CD4 T cells / CD8 T cells）/ `S6_1_Fig`–`S6_3_Fig`、`fig:supp-method-cor`、`-cd4`、`-cd8` | `Plot14_MethodCorr.R`；r、轴范围和显示比例在 FIG `Figure_S6_cor.tsv` | `4d80d1e15919` | 5 个有 effect 的方法两两的密度矩阵 + Spearman r，每页 A effect、B signed −log10 P 上下排。三页各为一个图（`page=` 选项），编号 S6.1 / S6.2 / S6.3：`\SixSub{n}` 在 figure 环境内局部改 `\thefigure` 和 `\theHfigure`，S6.2、S6.3 先 `\addtocounter{figure}{-1}`，计数器保持 6，之后的图不受影响。三张图同样大小（`height=0.86\textheight`）；完整说明放在三张图之前单独一页的 “Legend for S6.1–S6.3 Figs”（前后 `\clearpage`），每张图的 caption 只写标题、pair 数、A / B 和 “See the shared legend”；S6.3 后加 `\FloatBarrier`，防止 S1 Table 插到 legend 页和 S6.1 之间；轴范围（±1.5 / 2.5 / 0.75 / 0.35；±11）和舍去比例 0.3–1.9% 来自 tsv |
| S7.1–S7.3 Fig | `Figure_S7.pdf`（3 页：Bn / CD4 T cells / CD8 T cells）/ `S7_1_Fig`–`S7_3_Fig`、`fig:supp-manhattan`、`-cd4`、`-cd8` | `Plot15_Manhattan.R`；计数和标注基因在 FIG `Figure_S7_data.tsv` | `f8c0daf7459e` | 6 个方法的 Manhattan 图，每页 6 条轨道；Bonferroni 线；过线点 causal 橙 / marginal 绿；每轨道最多 8 个基因标签；Fig 2 例子洋红菱形。与 S6 同样的 `\SevenSub` 编号方式和共用 legend 页（2026-11-04 新增） |
| S8 Fig | `Figure_S8.pdf` / `S8_Fig`、`fig:supp-grn` | `Plot9_GO_GRN.R`；`Results/Validation/Valid5_GRN/` | `9418a7a48923` | Bn 和 CD8 T cells 的 GRN OR（口径同 Fig 3E；cell-type 864 / 1,364、high-quality 817 / 382 个 pair；2026-11-04 从 Fig 3 移出） |
| S1 / S2 Table | `tables/supplementary_hyper_{p,fold}.tex` / `S1_Table`、`S2_Table` | `Valid1_Hyper.R` → `Plot3_Hyper.R` → FIG `Table_S1_S2.xlsx`（`721d5095aeb8`），从 xlsx 逐格转写 | — | 18 行 × 6 个参考：Global P / pooled fold；方法名与图一致 |
| S3 Table | `tables/supplementary_data_overview.tex` / `S3_Table` | `Table4_DataOverview.R` → `Figures_20261027/Table_4.tsv` → `tmp/Manuscript_tables/make_data_table.py`（不手抄） | — | 数据概览；表注中的 CIMA 全队列数字待作者核对 |
| S4 Table | `tables/supplementary_gwas.tex` / `S4_Table`（`tab:supp-gwas`） | `Table5_GWASTraits.R` → `Figures_20261102/Table_S4.tsv` → `tmp/Manuscript_tables/make_gwas_table.py`（不手抄）；7 个 GWAS 的来源、N（LDSC 输入文件）、SNP 数、h²（2026-11-05 去掉了 “Fig 4B: MR-only links” 列）；PMID 经 PubMed 核对 | — | AS = Pan-UK Biobank，CD = IIBDGC 网站，MS = IMSGC 2011（ImmunoBase 受保护数据），其余 4 个来自 Zenodo 10515792；Asthma、PBC 不入表 |

- 补充条目放在 `supplementary.tex`。正文用 `\nameref{S1_Fig}` 引用，这类 label 定义在 `sections/supporting_info.tex`；补充 PDF 里的 caption 另用一套 `fig:supp-*` label。新增或重排补充条目时，这三处要一起改。
- `figures/Figure_1A.pdf`、`Figure_1B_C.pdf` 已经没有任何地方引用，暂时保留。
- 图的分析编号（`Scripts/Visualization.md` 用的 Fig1–Fig10）≠ 论文图号，对应关系以上表“来源”一列为准。

## 3. 投稿（PLOS Genetics，Research Article）

- **模板**：`Manuscript/template_reference/` 里的 PLOS LaTeX v3.8；`plos2025.bst` 不改。双倍行距，连续行号，最多三级标题，不用脚注；正文写 "Fig 1"，补充条目写 S1 Fig / S1 Table。
- **长度**：全标题 ≤ 200 字符，短标题 ≤ 70 字符（`\ShortTitle` 已定义但没有使用）。Abstract < 300 词，不含引用。Author Summary 150–200 词，非技术性，第一人称。
- **章节顺序**：Abstract → Author Summary → Introduction → Results → Discussion → Methods → Acknowledgments → References → Supporting Information。
- **投稿版改造**（草稿期不做）：合成单个 `.tex`，内联各 section 和参考文献；`\draftfiguresfalse`；去掉全部 TODO 和占位；补充条目只留 caption；出版用图另外准备 TIFF / EPS（不接受 LaTeX 生成的矢量 EPS）；表格不能有 `\multirow` 和图片。
- **必须提交**：STROBE-MR 核对表（标页码，作为补充材料，**现在还没有**）；Methods 里的 AI 使用披露（工具、用法、怎么评估输出、影响了哪些内容；素材见第 5 节）；Data Availability Statement（不接受 "data not shown"）。
- **投稿系统字段**（都还是 TODO；不能推断或编造，funding 不写进 Acknowledgments）：Data Availability（accession / DOI / 访问限制）；CRediT 作者贡献；Funding（资助方、编号、受助人、资助方的作用）；Competing Interests；作者信息和 ORCID（作者顺序和单位需要确认）；至少 4 位建议审稿人；投稿时选的 section。
- **参考文献**：Vancouver 编号；前 6 位作者后加 et al.；目前有些条目以 "others" 结尾，要补全；软件、资源和 GWAS 都要引用。
- 官方指南：https://journals.plos.org/plosgenetics/s/submission-guidelines 、https://journals.plos.org/plosgenetics/s/latex

## 4. Overleaf 同步

`push` 会改动用户的 Overleaf 项目，**只在用户明确要求时执行**。

```
Overleaf 项目 https://www.overleaf.com/project/6ae4b2115cc5cee7da165450
   │ git：https://git.overleaf.com/6ae4b2115cc5cee7da165450（用户名 git，密码 = token）
   ▼
tmp/overleaf_draft/my_paper/   Overleaf 的 git 克隆（不要在这里直接改文件）
   ▲ rsync，只同步白名单：main.tex supplementary.tex references.bib plos2025.bst sections/ figures/ tables/
   ▼
Manuscript/                    服务器上的工作目录
```

- **命令**（在项目根目录运行）：`bash tmp/overleaf_sync.sh status`（dry-run，不改任何东西）/ `pull` / `push`（会先自动 pull，没有变化时提示 nothing to sync）。文件按内容比较（rsync `-c`）。默认不同步删除；要同步删除时，先 `DELETE=1 … status` 预览，再 `DELETE=1 … push`。
- **token** 在 `tmp/.overleaf_token`（权限 600），由 `tmp/overleaf_askpass.sh` 通过 `GIT_ASKPASS` 提供。**不要读出、打印或提交 token。** 更换：Overleaf → Account Settings → Git Integration 重新生成，然后 `printf '%s' '新token' > tmp/.overleaf_token; chmod 600 tmp/.overleaf_token`。
- **pull 之前**：先 `git status --short Manuscript`。如果服务器上有还没 commit 的改动，先和用户商量，不要直接覆盖。pull 之后用 `git diff` 汇报改动，再编译。
- **冲突**：不 force。fast-forward 失败或 push 被拒时，把克隆的 `git log --oneline --all` 和 `git status` 报给用户，一起决定以哪边为准。两边改了同一个文件时：先把 `Manuscript/` 里的那个文件备份到 `tmp/`，pull，再手工合并。认证失败就换 token。
- Overleaf 一次只编译一个根文档；看补充材料要在 Menu → Main document 里切换到 `supplementary.tex`。
- `tmp/` 不进 git，所以下面保存了脚本全文。**从零重建**：(1) 拿 token，`umask 077; printf '%s' '<token>' > tmp/.overleaf_token`；(2) 新建 `tmp/overleaf_askpass.sh`（`chmod 700`），内容如下；(3) `GIT_ASKPASS=$PWD/tmp/overleaf_askpass.sh GIT_TERMINAL_PROMPT=0 git clone https://git.overleaf.com/6ae4b2115cc5cee7da165450 tmp/overleaf_draft/my_paper`，然后在克隆里 `git config user.name "Tengfei Cui"`、`git config user.email "tfcui2025@tamu.edu"`；(4) 把下面的脚本存为 `tmp/overleaf_sync.sh`（`chmod +x`），先 `status`，再 `pull`。脚本改了要同步更新这里。

```bash
#!/bin/bash
case "$1" in
  Username*) echo git ;;
  Password*) cat "$(dirname "$(readlink -f "$0")")/.overleaf_token" ;;
esac
```

`tmp/overleaf_sync.sh` 全文：

```bash
#!/bin/bash
# Overleaf <-> Manuscript/ 同步
#   bash tmp/overleaf_sync.sh status   # 对比两边（dry-run，不改任何东西）
#   bash tmp/overleaf_sync.sh pull     # Overleaf -> Manuscript/
#   bash tmp/overleaf_sync.sh push     # Manuscript/ -> Overleaf（先自动 pull 检查）
# 加 DELETE=1 才会同步删除（先用 status 看 dry-run）。
set -euo pipefail

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project
MS=$PROJ/Manuscript
CLONE=$PROJ/tmp/overleaf_draft/my_paper
ITEMS=(main.tex supplementary.tex references.bib plos2025.bst sections figures tables)

export GIT_ASKPASS=$PROJ/tmp/overleaf_askpass.sh
export GIT_TERMINAL_PROMPT=0

RSYNC_OPTS=(-rcv --exclude='*.before' --exclude='.DS_Store')
[ "${DELETE:-0}" = 1 ] && RSYNC_OPTS+=(--delete)

# rsync 白名单：每一项 src -> dst
sync_items() {  # $1=src root  $2=dst root  $3...=额外 rsync 参数
  local src=$1 dst=$2; shift 2
  for it in "${ITEMS[@]}"; do
    if [ -d "$src/$it" ]; then
      rsync "${RSYNC_OPTS[@]}" "$@" "$src/$it/" "$dst/$it/"
    else
      rsync "${RSYNC_OPTS[@]}" "$@" "$src/$it" "$dst/$it"
    fi
  done
}

git_pull_clone() {
  git -C "$CLONE" diff --quiet && git -C "$CLONE" diff --cached --quiet \
    || { echo "克隆目录有未提交改动，请先处理：$CLONE"; exit 1; }
  git -C "$CLONE" pull --ff-only
}

case "${1:-}" in
  status)
    git -C "$CLONE" fetch origin
    echo "== Overleaf 上有、本地克隆没有的提交："; git -C "$CLONE" log --oneline HEAD..origin/main
    echo "== Manuscript/ -> 克隆 的差异（dry-run，即 push 会改的内容）："
    sync_items "$MS" "$CLONE" -n -i | grep -E '^[<>c*][fLDS]' || echo "(无)"
    ;;
  pull)
    git_pull_clone
    sync_items "$CLONE" "$MS" -i | grep -E '^[<>c*][fLDS]' || echo "(Manuscript/ 已是最新)"
    ;;
  push)
    git_pull_clone
    sync_items "$MS" "$CLONE"  > /dev/null
    if [ -z "$(git -C "$CLONE" status --porcelain)" ]; then
      echo "nothing to sync"; exit 0
    fi
    git -C "$CLONE" status --short
    git -C "$CLONE" add -A
    git -C "$CLONE" commit -q -m "Sync from server $(date +%F\ %H:%M)"
    git -C "$CLONE" push origin HEAD:main
    echo "pushed: $(git -C "$CLONE" log --oneline -1)"
    ;;
  *)
    sed -n '2,6p' "$0"; exit 1 ;;
esac
```

## 5. 更新记录（只追加，一行一件事：日期 | 工具和模型 | 改了什么 | 怎么核对）

所有改动都记在这里：每轮润色（W2 / W3）、Figures 更新、格式 / 模板 / 补充材料版面改动、commit 和 Overleaf push。

AI 使用披露的素材就来自这里。2026-11-01 之前的条目是压缩过的里程碑，细节见 git 历史。

- 2026-10-15 | CodeX | 初稿 `Manuscript/`（PLOS 模板 + 旧图）；可从提交 7cf0b1a 恢复
- 2026-10-20 | Overleaf 导入 | `Manuscript_MR_v2` 替换初稿；新增 Data Preprocessing
- 2026-10-20 | CodeX | 按代码核对后写入 SMR、pseudobulk、WNN 三个 Methods 小节
- 2026-10-20 | Claude Code | 补充图表拆成独立的 `supplementary.tex`，正文用 `\nameref` 引用
- 2026-10-23~25 | Claude Code | Fig 1–4、S1–S3 重画替换；超几何表和 motif 表与图统一；Fig 1 拍平透明度
- 2026-10-29 | Overleaf 导入 + Claude Code | 导入用户在 Overleaf 上的修改；细胞类型名称全文统一；`Figures_20261102` 一轮（GRN 并入 Fig 3C、GNAS 并入 Fig 4B、Fig 4A 只展示 7 个 GWAS、motif 表进正文作为 Table 1、数据概览改为 S3 Table）
- 2026-11-01 | Claude Code | 建立 Overleaf git 同步，首次 pull；删除 `CODEX_METHODS_NOTES.md`
- 2026-11-01 | Claude Code (claude-opus-5-5) | A1 审计 Methods 与代码的对照（只读）；Results 证据卡 R1–R5（只读，数字从 `Results/MR` 重算）
- 2026-11-01 | Claude Code (claude-opus-5-5) | 笔记精简：17 个 md 合并为 `polish/PLAYBOOK.md`、`polish/REVISION_PLAN.md`、`notes/REFERENCE.md` 和 `AGENTS.md`；旧文件打包存为 `tmp/manuscript_notes_backup_20261101.tar.gz` 后删除 | 脚本核对 A1/P/X 编号全部保留在 REVISION_PLAN 中；sync 脚本与 `tmp/overleaf_sync.sh` 用 diff 比对一致
- 2026-11-01 | Claude Code (claude-opus-5-5) | P2：REVISION_PLAN 第 3 节（Abstract 逐句标注 A0 + 规划 A1–A8）和第 5 节（Discussion D0–D8，局限性 7 项全覆盖）；未改 .tex | S1 Fig 全局 ρ 读自 `Results/QC/1_OneSampleMR/*_global_rho.tsv`；工具变量数分布由 `Results/MR/*_AllPairs_results.tsv.gz` 的 `ivw_n_snp` 计算；其余数字引用第 4、6 节已 ✓ 的条目
- 2026-11-01 | Claude Code (claude-opus-5-5) | P3（Methods + Results + Discussion 清理）+ P4（34 条现稿段意）：Methods 按 M0 重排，新建 Study design / Filtering / GRN / GWAS overlap / Allele orientation 小节，推导和截断细节移到新文件 `sections/s1_text.tex`（S1 Text，由 `supplementary.tex` 引入，补充文档加了 `cite` 和参考文献表）；`\q` 改为 $P_{\mathrm{adj}}$；Results 按 R1–R5 修正并新增 R3-P3 超几何段；Discussion 删提纲、加 D1–D8 注释骨架；Abstract 未动 | invariants 差异逐项核对（均为计划内的移动和修正）；`build_ms.sh` 两个文档编译通过，main.log 无 undefined；supplementary 手动补跑 bibtex 后无 undefined（`build_ms.sh` 尚未包含 supplementary 的 bibtex）
- 2026-11-02 | Claude Code (claude-opus-5-5) | 阶段 2 规则改为可循环（PLAYBOOK 第 3、7 节，REVISION_PLAN 第 0 节；去掉“[阶段 2]”标志，现稿段意带日期）；P5 润色 Methods 的 Data sources、Study design and terminology（按作者段意把列表改成第二段文字）、Data Preprocessing、Mendelian randomization、Causal IVW、Causal Egger、Causal SMR | `ms_check.py invariants --base 69915a8`（润色前的 `git stash create` 快照）：数字、引用、label、公式无差异，只有 M0b 改写造成 1 处 $P$ / $P_{\mathrm{adj}}$ 记号增减；`build_ms.sh` 通过，main.log 无 undefined
- 2026-11-02 | Claude Code (claude-opus-5-5) | P5 按作者意见改 M2（删 Final analysis set）、M3（一般 MR → MR-Egger 式 → α₀ = 0 即 IVW → 本研究）、M4（行内公式、Eckart–Young、推断段，目标函数移到 S1）和 S1 Text（新增 NOME、IVW 一致性、InSIDE、相关变异，转述 Burgess & Thompson 2017 附录 A.2 / A.4；新增随机效应 / overdispersion 小节） | `ms_check.py invariants --base 69915a8`：dp 只少 3 个 pair 数和 2 个 `\nameref`，methods 只少 1 处重复引用，公式差异均为计划内的移动 / 合并 / 新增（s1_text 不在快照中）；分散度实现与 `MR_AllPairs.R:98-126` 核对；`build_ms.sh` 通过，supplementary 手动补跑 bibtex 后两个 log 都没有 undefined
- 2026-11-02 | Claude Code (claude-sonnet-5-5) | Table 1（`tables/motif.tex`）最后一列改为 "Brief interpretation"，去掉 PubMed 链接，句末改 `\cite{key}`；`references.bib` 新增 6 条（nus2017atf3、ubieta2017fra2、dezfouli2006mnt、parveen1995cenpb、wu2023pou3f1、inukai2005tef）；`tmp/Manuscript_tables/make_motif_table.py` 同步；results.tex 的 `[CITE: PMID …]` 未动 | 6 条元数据（作者、期刊、年份、卷页、DOI）用 PubMed 记录逐条核对；编译无 undefined citation
- 2026-11-02 | Claude Code (claude-sonnet-5-5) | 标题错别字 ehancer → enhancer（`main.tex:86`）；Figure_1 重新拍平并更新 `figures/Figure_1.pdf`，Figure_2 去掉 "Example peak-gene pairs" 标题并更新 `figures/Figure_2.pdf` | 编译无 undefined；`ms_check.py invariants` 只显示 Table 1 的 6 个新 `\cite` 和被去掉的 PMID 链接
- 2026-11-04 | Claude Code (claude-opus-5-5) | 补充材料新增 S4 Fig（各方法显著 link 的跨细胞类型 Venn）和 S5 Fig（MR-only link 与 6 个参考的观测 vs 背景），写入 caption；`sections/supporting_info.tex` 加 `S4_Fig`、`S5_Fig` 条目；正文尚未引用这两张图 | caption 中的数字（14,890；40,087 / 85,527 / 50,931；约 55% 的 gene 只有 1 个 MR-only peak；B = 2,000）与 `Plot12` / `Plot13` 输出核对；`build_ms.sh` 通过，两个 log 都没有 undefined；figures 中的 PDF 与 FIG 同名文件 hash 一致
- 2026-11-04 | Claude Code (claude-opus-5-5) | 补充材料新增 S6 Fig（5 个方法两两的 effect / signed −log10 P 密度矩阵和 Spearman r，3 页），写入 caption；`sections/supporting_info.tex` 加 `S6_Fig`；补充材料开头改为 S1–S6 Figs；正文尚未引用 | caption 中的数字（40,087 / 85,527 / 50,931；轴范围；0.3–1.9%）与 `Figure_S6_cor.tsv` 核对；`build_ms.sh` 通过，两个 log 都没有 undefined；S6 在 `supplementary.pdf` 第 13–15 页，已转 PNG 查看；figures 中的 PDF 与 FIG 同名文件 hash 一致
- 2026-11-04 | Claude Code (claude-opus-5-5) | S6 Fig 拆为 S6.1 / S6.2 / S6.3 Fig（Bn / CD4 T cells / CD8 T cells 各一个编号和 caption，替代 `\ContinuedFloat`）；`supporting_info.tex` 的 `S6_Fig` 改为 `S6_1_Fig`–`S6_3_Fig`；补充材料开头改为 “S1–S5 Figs, S6.1–S6.3 Figs” | `build_ms.sh` 通过，两个 log 都没有 undefined，也没有重复 hyperref 锚点；图目录列出 S6.1–S6.3（第 13–15 页），已转 PNG 查看
- 2026-11-04 | Claude Code (claude-opus-5-5) | S6.1–S6.3 Fig 改为同样大小；共用的完整说明移到三张图之前单独一页（Legend for S6.1–S6.3 Figs），每张图的 caption 改为短 caption；S6.3 后加 `\FloatBarrier` | `build_ms.sh` 通过，两个 log 都没有 undefined；页序为 legend 第 13 页、S6.1–S6.3 第 14–16 页、S1–S3 Table 第 17–19 页（共 20 页），已转 PNG 查看
- 2026-11-04 | Claude Code (claude-opus-5-5) | `supplementary.tex` 首页去掉论文标题和 “S1 Text, … accompany the main manuscript” 一句；“Supporting” 改为 “Supplementary”（页眉、首页标题、目录名）；新增 S1 Text 目录 “Supplementary text”（`\tableofcontents`；`\SoneTOC` 让 S1 Text 的 `\subsection*` 也写入目录，`s1_text.tex` 不改）。正文 `sections/supporting_info.tex` 的 “Supporting information” 是 PLOS 要求的标题，未改 | `build_ms.sh` 通过，两个 log 都没有 undefined；首页转 PNG 查看：依次为 Supplementary text（S1 Text 及 8 个小节和页码）、Supplementary figures、Supplementary tables
- 2026-11-04 | Claude Code (claude-opus-5-5) | Figure 3 换新版（A Miami / B GO BP / C GO MF / D HiChIP 校准 / E GRN CD4），重写 caption，正文和 Methods 中 Fig 3 的字母改为新版（GO → B/C，GRN → E + S8 Fig）；Figure S5 换为 Fig 3D 风格并改 caption（点的颜色、Fold / P 位置）；新增 S7.1–S7.3 Fig（Manhattan，共用 legend 页）和 S8 Fig（Bn / CD8 GRN），`supporting_info.tex` 加 4 个 label；R3 开头留 TODO（Fig 3A、3D、S7 尚无正文引用）；Fig 3 图高限为 0.53\textheight 让 caption 同页 | `build_ms.sh` 6 步都返回 0，两个 log 无 undefined reference；caption 中的数字逐个对照 `Figure_S7_data.tsv`、`Figure_S5_summary.tsv`、`Valid5_GRN/GRN_summary.tsv` 和原 caption；渲染 main 第 9 页、supplementary 第 18–21 页目检 |
- 2026-11-04 | Claude Code (claude-opus-5-5) | Fig 3 恢复为全文宽（`width=\linewidth`，取消 0.53\textheight 的缩小）；caption 拆成两部分：图页写标题 + A，B–E 放到下一页的 `\ContinuedFloat` 续 caption（“Fig 3 (continued).”），文字不变 | `build_ms.sh` 6 步返回 0；main.log 无 undefined reference、无 “Float too large”；渲染第 9–10 页目检 |
- 2026-11-04 | Claude Code (claude-opus-5-5) | 管理文件重组：`polish/PLAYBOOK.md`、`polish/REVISION_PLAN.md`、`notes/REFERENCE.md` 改为 `notes/1_WORKFLOW.md`（分工、W1–W3 工作卡流程、润色规范、scMultiMap 写法对照表）、`notes/2_BOARD.md`（统一表达、全文地图、全文问题、小节工作区）、`notes/3_LEDGER.md`（新增第 1 节已核对数字；原第 4 节工作记录改为第 5 节更新记录）；旧文件按作者意见不备份、直接删除；`AGENTS.md`、项目 `CLAUDE.md`、.tex 中指向旧文件的注释和 `\TODO` 指针同步改名；`tools/ms_check.py invariants` 新增 `--snapshot DIR` | 迁移逐条核对（遗留问题、已定、✓ 数字均有去处）；`--snapshot` 用改动一个数字的快照测试，能报出差异；`build_ms.sh` 通过
- 2026-11-05 | Claude Code (claude-sonnet-5-5) | 新增 S4 Table（GWAS 性状汇总）：`tables/supplementary_gwas.tex` 加入 `supplementary.tex`（S3 Table 之后），`supporting_info.tex` 加 S4 Table 条目；`references.bib` 新增 5 条（delange2017ibd、saevarsdottir2022ra、bentham2015sle、chiou2021t1d、sawcer2011ms，PMID 28067908 / 35470158 / 26502338 / 34012112 / 21833088 均经 PubMed 核对）；Methods S-LDSC 段：4 处 `[CITE: …]` 改为 `\cite`，AS / CD / MS 来源句替换原 `\TODO`，并引用 `\nameref{S4_Table}`（`Methods` 其余 `\TODO` 未动，包括 “9 or 7 GWAS”） | `build_ms.sh` 通过 + `bibtex supplementary` 重跑；两个 log 无 undefined；`ms_check invariants --snapshot tmp/ms_rounds/20261105_s4table` 只显示 5 个新 cite、`S4_Table` 的 label / nameref 和被替换的 `[CITE]` 里的年份；S4 Table 为 supplementary.pdf 第 25 页，已转 PNG 查看
- 2026-11-05 | Claude Code (claude-sonnet-5-5) | S4 Table 去掉最后一列 “Fig 4B: MR-only links” 及对应表注，列宽重排；`Table_S4.xlsx` / `.tsv` 同步（`Table5_GWASTraits.R` 不再输出这些列） | `supplementary.tex` 重编两遍无 undefined，S4 Table 仍在第 25 页，已转 PNG 查看
- 2026-11-06 | Claude Code (claude-sonnet-5-5) | Fig 4 微调：B 去掉 UBAC2（只有 Asthma）换成 ERAP1（Bn，CD）、ERBB3 标签去掉 Asthma；A 图例放入图内、字号和 A / B 字母放大；柱子 / Gene / SNPs / Peak 改 Set3 配色（红色 link 不变）；`results.tex` 的 UBAC2 句改为 ERAP1，caption 的 “dark-blue ticks / gold bar” 改为 “lavender ticks / orange bar”，Methods 写 “seven GWAS” 并删除对应 `\TODO` | `Plot10_Combined.R` 预览 PNG 逐轮检查；S-LDSC 数值与 1.5 节一致；`build_ms.sh` 6 步返回 0，supplementary 无 undefined，main 只有 motif 表中原有的 4 条 undefined 引用；图 hash 与 FIG 一致 `e5f1b6fa1da5`
- 2026-11-06 | Claude Code (claude-opus-5-5) | Table 1：ATF3、FOSL2::JUND、MNT、TEF 四行的小鼠文献换成人类研究（zhang2018atf3、wu2024jund、vermeer2008mnt、rowsey2019tef；bib 删掉 4 条旧条目），表头和表注的方法名改为 Causal IVW / Causal Egger / Causal SMR / Causal Joint / Marginal Multiomics / Marginal Pseudobulk；`results.tex` R5-P2 中 MNT、FOSL2::JUND 的两处 `[CITE: PMID]`（小鼠研究）改为新的 `\cite`，句子随之改写；`make_motif_table.py` 加 `BRIEF_OVERRIDE` 同步 | 4 篇用 PubMed 记录核对作者 / 期刊 / 卷页 / DOI；脚本重跑后与快照 diff 只有文字和方法名；编译无 undefined；`ms_check.py invariants --snapshot tmp/ms_rounds/20261106_table1` 只报引用和注释文字
