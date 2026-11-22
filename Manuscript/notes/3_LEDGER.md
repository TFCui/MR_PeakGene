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
- **2026-11-22 重算更正**（`Results/MR/*_AllPairs_results.tsv.gz`，按 IVW 符号）：与 multiomics / pseudobulk 不一致 Bn 24.1% / 1.1%（113、5 / 469），CD4 25.7% / 2.9%（324、37 / 1,259），CD8 24.8% / 2.1%（188、16 / 759）；上一行 3.0%、24.2%、1.0% 是扇区四舍五入后相加的结果，正文用本行。
- MR-only link（IVW / Egger / SMR 任一 $P_{\mathrm{adj}}<0.1$，两个 marginal 都 ≥ 0.1）：Bn 2,262 / CD4 4,749 / CD8 2,886。CD4 UpSet 中只含 causal 方法的最大组合 870 = SMR + IVW + joint（四个 causal 都显著、marginal 都不显著的为 680）。
- S4 Fig（`Figure_S4_counts.tsv`，14,890 个共有 pair，region_111 / 并集）：SMR 22.1%、IVW 14.9%、Egger 8.8%、joint 16.4%、multiomics 91.5%、pseudobulk 12.5%。
- S6 Fig（`Figure_S6_cor.tsv`，全部 pair 的 effect Spearman r，Bn / CD4 / CD8）：IVW–pseudobulk 0.30 / 0.24 / 0.27；IVW–multiomics 0.034 / 0.042 / 0.041；multiomics–pseudobulk 0.038 / 0.065 / 0.048；IVW–SMR 0.86 / 0.87 / 0.87；IVW–Egger 0.56 / 0.56 / 0.55。signed −log10P：IVW–pseudobulk 0.31 / 0.26 / 0.28，IVW–multiomics 0.04。
- Fig 3A（CD4）：causal joint 最小 P 8.2e-190；multiomics 最小 P 8.4e-32（−log10 31.1）。

### 1.3 R2：四个例子（Fig 2；`Scripts/Visualization/Example_pairs.tsv`）

| 例子 | peak | IVW β（$P_{\mathrm{adj}}$） | Egger / SMR β | multiomics ρ（$P_{\mathrm{adj}}$） | pseudobulk ρ（$P_{\mathrm{adj}}$） |
|---|---|---|---|---|---|
| *SIDT2*（Bn） | chr11:117185715-117186216 | 1.024（1.6e-171），194 个 SNP | 1.026 / 1.013 | 0.586 | 0.792（2.0e-72） |
| *RPL11*（CD4） | chr1:23691494-23691995 | 1.127（6.0e-11） | — | 0.400（3.9e-4） | 0.010（0.92） |
| *ERAP2*（CD4） | chr5:96885894-96886395 | −1.525（4.2e-133） | −1.510 / −1.470 | +0.375（P 6.1e-4，$P_{\mathrm{adj}}$ 1.5e-3） | −0.641（9.9e-44） |
| *ARPC1B*（CD4） | chr7:99376204-99376705 | 0.740（1.8e-47） | 0.719 / 0.784 | −0.149（0.37） | 0.016（0.79） |

SIDT2 的 P 值排序：IVW < joint < Egger < SMR。结果表里没有 Egger 截距。

### 1.4 R3：GO、参考数据集、GRN（Fig 3、S1 / S2 Table、S5 Fig、S8 Fig）

- **GO（DAVID，2026-10-19）**：MR-only 基因 4,509 个（3 个细胞类型合并；`Results/Validation/Valid7_GO/DAVID_input/MROnly_genes_all.txt`），有 BP 注释 4,031、MF 注释 4,185；背景是 DAVID 默认的全基因组。DNA-templated transcription：689 个基因，fold 1.36，FDR 5.18e-18。（2026-11-22 核对 `DAVIDChartReport_A1BG_2026-10-19.csv`：正文的 5.20e-18 是 Benjamini 列，即 BH；正文其余 term 的值也都是 Benjamini 列。*SRSF7*、*DDX5* 在 RNA splicing / RNA binding / mRNA binding 三个 term 中，*UBE2N* 在 protein ubiquitination 和 ubiquitin-protein transferase activity 中。）protein binding fold 1.16，FDR 2.81e-85；RNA binding fold 1.41，FDR 3.58e-15；zinc ion binding fold 1.33，FDR 3.58e-15。
- **超几何（S1 / S2 Table）**：所有方法 × 参考 × 细胞类型都富集，P ≤ 3.6e-17，fold 1.68–6.18。causal SMR 在 17 / 18 个组合中 fold 最大（例外：CD4 PCHi-C，IVW 2.61 vs SMR 2.55）；marginal multiomics 在 18 个组合中都最小。Bn：FANTOM5 SMR 6.18 vs multiomics 3.56，HiChIP 4.46 vs 3.08。
- **MR-only 与参考（S5 Fig、Fig 3D）**：约 55% 的 gene 只有 1 个 MR-only peak；分箱 95% CI 来自按 gene 的 bootstrap（B = 2,000）。数值见 FIG `Figure_S5_summary.tsv`、`Figure_S5_bins.tsv`。 pooled fold：GTEx 2.17、EnhancerAtlas 2.35、FANTOM5 2.84、4DGenome 2.87、HiChIP 3.28、PCHi-C 1.89；HiChIP 8 组 k/n 全部高于 K/N（第 1 组 0.087 vs 0.014，第 8 组 0.584 vs 0.264）。
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
  - 总数（2026-11-11 核对，`GWAS_example_candidates.tsv` 中 `in_old_strict`）：7 个 GWAS 下共 24 个 link（Bn 7 / CD4 7 / CD8 10），18 个基因、12 个 peak，性状 CD / IBD / RA / T1D；6 个例子都在其中。*DR1* 在 3 个 CT、*GNAS* 在 CD4 / CD8 都是 strong MR-only link。
  - R4-P2 草稿另用到的数字（来源 `tmp/Valid6_GWAS/annot/GWAS_example_{candidates,external_snp,external_peak}.tsv`、`GWAS_Fig4B_examples.tsv`）：peak 中心到 TSS 距离 ERAP1 88 kb / DR1 166 kb / ERBB3 35 kb / GNAS 415 kb / TMBIM1 5.3 kb / THBS3 18.6 kb；lead SNP 对该 peak 的 caQTL P 5.2e-27 / 4.1e-30 / 4.5e-100 / 9.8e-28 / 6.0e-10 / 2.5e-22；eQTL P DR1 1.9e-33、TMBIM1 1.3e-10；与 MR SNP 最大 r²（EAS）1.00 / 1.00 / 0.99 / 本身 / 1.00 / 本身；Open Targets target H4 0.65 / 0.98 / 0.99 / 0.995 / 0.98 / 0.99；FAVOR rs13031757 CADD 11.0；RegulomeDB rs259956 1b；GeneHancer TMBIM1 排 1/15。
  - *ERBB3*（CD4，chr12:56041553-56042054）：rs705704，T1D P 1.5e-62；IVW $P_{\mathrm{adj}}$ 1.4e-4
  - ~~*RGL2*（Bn，chr6:33313910-33314411）：rs3130099 / rs3106189，T1D P 7.5e-54，RA P 7.1e-15~~（2026-11-10 去掉）
  - ~~*POLR2E*（CD4，chr19:1095225-1095726）：rs1046040 / rs11084882，CD P 1.5e-8~~（2026-11-10 去掉）
  - ~~*UBAC2*（CD4，chr13:99375239-99375740）：Asthma P 2.8e-9~~（2026-11-06 去掉：只有 Asthma，不在 7 个展示的 GWAS 内）
  - *ERAP1*（Bn，chr5:96896259-96896760，2026-11-06 新加）：rs2549785，CD P 4.8e-11；SMR / IVW / Egger FDR 2.6e-5 / 8.0e-6 / 4.8e-4，WNN / Pseudobulk 不显著（0.19 / 0.96）
  - *GNAS*（CD8，chr20:59254001-59254502）：rs259956，IBD P 2.6e-9，CD P 8.4e-9
  - 2026-11-10 新加（Fig 4B 第 2 版，选例依据见 `Results/Validation/Valid6_GWAS/GWAS_Fig4B_examples.tsv`；FDR 依次为 SMR / IVW / Egger / Joint，WNN / Pseudobulk）：
    - *DR1*（CD4，chr1:93179440-93179941）：rs1335679，CD P 3.9e-10；2.4e-16 / 3.8e-17 / 5.4e-8 / 9.8e-17；0.94 / 0.45
    - *TMBIM1*（CD8，chr2:218287070-218287571）：rs13031757，CD P 2.0e-8；4.0e-5 / 1.9e-4 / 0.53 / 1.0e-4；0.17 / 0.094
    - *THBS3*（CD8，chr1:155227359-155227860）：rs2990245，IBD P 7.2e-9；1.3e-4 / 2.8e-4 / 0.023 / 2.7e-4；0.067 / 0.45
  - 2026-11-11 正文只留 ERAP1 / ERBB3 / THBS3；IVW β（`Results/MR/*_AllPairs_results.tsv.gz`，2026-11-22）ERAP1 −0.359（Egger −0.327、SMR −0.434）、ERBB3 +0.242（0.203 / 0.255）、THBS3 −0.303（−0.225 / −0.506）；同一 peak 的邻近基因：ERAP1 peak → *ERAP2* β +2.28（IVW $P_{\mathrm{adj}}$ 1.6e-38），ERBB3 peak → *RPS26* β −1.23（9.6e-46）、*SUOX* +0.39（2.3e-9）；外部证据数值：HiChIP（Chandra 2021 Supp Table 3，hg19 5 kb bin liftOver 后含 peak）ERAP1 naive B raw count 44、Q 8.4e-3（bin chr5:96230000-96235000 hg19），THBS3 naive CD8 raw count 1,300、Q 1.9e-48（bin chr1:155195000-155200000 hg19）；GTEx（`gtex.rda` 的 `meta.pval`）rs2990245 对 THBS3 P 5.2e-33；RegulomeDB 三个 lead 都是 1f、probability 0.55；FAVOR aPC-TF / aPC-EpiActive：rs705704 29.1 / 18.1，rs2990245 21.4 / 13.0，rs2549785 16.2 / 18.3（PHRED 标度：20 = 前 1%，30 = 前 0.1%）

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
| Fig 3 | `figures/Figure_3.pdf` / `fig:go` | `Plot9_GO_GRN.R` 拼图（13×15 in，2026-11-06 起）：A 由 `Plot15_Manhattan.R`、D 由 `Plot13_MROnly_ObsExp.R` 存 rds 后读入；B/C 用 DAVID 结果；E 用 `Results/Validation/Valid5_GRN/` | `485013d70aa7` | **2026-11-06 第 3 轮（当前）**：字母 A–D。A 三条同向 Manhattan（Causal joint 薄荷绿 / Marginal multiomics 青绿 / Marginal pseudobulk 青蓝，无标题条，基因名加粗放大）；B = GO BP（左）+ MF（右），共用字母；C = HiChIP 校准图（L 形箭头坐标轴，轴名 Background ratio (K / N) / MR-Only ratio (k / n)）；D = GRN OR 横向条形（箭头坐标轴，cell-type #AADCA9 / high-quality #58B8D1）。配色 `PAL_FRESH`。正文 / Methods / 补充材料中的 Fig 3 字母引用已同步（MF 写作 “B, right”）。以下为更早版本描述：**2026-11-06 改版**：A 改为 CD4 T cell 的三条同向 Manhattan（Causal joint / Marginal multiomics / Marginal pseudobulk，每条最多 15 个基因名 + 例子）；B GO BP | C GO MF 同一行；D HiChIP | E GRN 同一行、等宽对齐；caption A 加 “marginal pseudobulk”。以下为 2026-11-04 版描述：A CD4 T cell Miami 图（上 causal joint / 下 marginal multiomics，Bonferroni 0.05/85,527，例子 RPL11 / ERAP2 / ARPC1B）；B GO BP、C GO MF（各前 10 个 term，4,031 / 4,185 个有注释基因）；D MR-only link 与 HiChIP 的校准图（同 S5 E，右下 Fold = 3.28）；E GRN OR 只画 CD4 T cells（2,520 / 458 个 pair），Bn / CD8 在 S8。各 panel 用浅薄荷标签条，薄荷色板。tex 中图宽 = `\linewidth`（不缩小）；caption 只在图页写标题和 A，B–E 用 `\ContinuedFloat` + `\caption*{Fig 3 (continued). …}` 接到下一页（2026-11-04 作者要求） |
| Fig 4 | `figures/Figure_4.pdf` / `fig:disease` | `Plot10_Combined.R`（`SL_TRAITS` 控制展示哪些 GWAS）；A 用 `Valid2_SLDSC/SLDSC_all.tsv`，B 用 `tmp/Valid6_GWAS/annot/GWAS_example_candidates_snps.tsv`（2026-11-10 起；6 个 pair 都在 strong 集 `GWAS_MRonly_overlap.tsv` 内） | `f21ddcdb962e` | A：7 个 GWAS 的 MR-only 与全部 peak 对比，`*` 表示 P < 0.1，不画 s.e.；图例在图内左上角；B：2026-11-11 起 1 行 × 3（16×14 in，行高 9 in）：Bn ERAP1 / CD4 ERBB3 / CD8 THBS3（FIG `Figure_4.pdf`；11-10 的 2 行 × 3 试看版 `Figure_4_CC1.pdf` 已删除），红色 link 的阈值是 0.05。2026-11-06 配色改 Set3（MR-Only 薄荷青 #8dd3c7 / All peaks 灰 / Gene 蓝 #80b1d3 / MR SNP 淡紫 / GWAS SNP 紫 #bc80bd / Peak 橙 #fdb462），字号约 ×1.3–1.5，A / B 字母 32 pt |
| Table 1 | `tables/motif.tex` + `figures/motif_logos/*.png` / `tab:motif` | `Valid3_Motif/Motif_MR_specific.tsv` + `Motif_table_main.tsv`（解释和 PMID 是人工整理的，只有这一个来源；2026-11-06 起 ATF3、FOSL2::JUND、MNT、TEF 四行的解释和 PMID 由脚本里的 `BRIEF_OVERRIDE` 覆盖为人类研究）→ `tmp/Manuscript_tables/make_motif_table.py`；完整 xlsx 在 FIG `Table_1.xlsx`（`537c47ae7b9f`） | — | 6 行，PMID 用 `\url`，不进 bib。PLOS 不允许表格里放图片，投稿前要处理 logo |
| Fig A in S1 Text | TikZ，写在 `sections/s1_text.tex` “Marginal association and causal effect” 小节内（非浮动体，不占 S 图编号） | 依据作者推导图 `paper_reference/image.png` 重画；符号改为 γ（Z→X）/ ψ_X / ψ_Y，与 S1 Text 一致 | — | 因果图 Z→X→Y，U→X、U→Y；图注为 minipage 小字 |
| S1 Fig | `Figure_S1.pdf` / `S1_Fig`、`fig:supp-qc` | `Plot11_OneSampleMR.R`；上游 `Scripts/QC/QC1_*` | `c0aa5df5b020` | null SNP 的 z 相关（A）和 P 值（B） |
| S2 Fig | `Figure_S2.pdf` / `S2_Fig`、`fig:supp-flip` | `Scripts/QC/QC2_2_FlipPoints_plot.R`；数值在 `Results/QC/2_FlipPoints/FigureS2_estimates.tsv` | `21cbfe159813` | A：flip 前后差别小的例子（SIDT2 / RPL11 / LCLAT1）；B：差别大的例子（PIGX / PDE4DIP / KLRC3） |
| S3 Fig | `Figure_S3.pdf` / `S3_Fig`、`fig:supp-bn-cd8` | `Plot1_MRresults.R` | `c8ba628007de` | Bn 和 CD8 的 UpSet + 环形图；caption 里写死了 "Fig 1B–C" |
| S4 Fig | `Figure_S4.pdf` / `S4_Fig`、`fig:supp-venn` | `Plot12_Venn.R`；计数在 FIG `Figure_S4_counts.tsv` | `52cae9934182` | 6 个方法各一个三圆 Venn，只统计 3 个细胞类型都检验过的 14,890 个 pair；圆内为数量（占并集 %），圆外为共有 pair 中的显著数 / 该细胞类型全部显著数（%）；caption 说明 WNN 不分细胞类型 |
| S5 Fig | `Figure_S5.pdf` / `S5_Fig`、`fig:supp-mronly-ref` | `Valid1_Hyper_MROnly.R` → `Plot13_MROnly_ObsExp.R`；数值在 FIG `Figure_S5_summary.tsv`、`Figure_S5_bins.tsv` | `a3de0f3043c4` | MR-only link（3 个细胞类型合并）的金标准命中率 vs 背景 K/N：按 K/N 等数量 8 组，95% CI 为按 gene 的 bootstrap。2026-11-04 起为 Fig 3D 的风格：库名标签条，右下角 pooled fold 和合并超几何 P（不再有副标题和 gene 数）；E 即 Fig 3D。caption 中引用表格用 `\ref{tab:supp-hyper-*}`（补充文档里没有 `S1_Table` 这类 label） |
| S6 Fig | `Figure_S6.pdf`（3 页 = panel A Bn / B CD4 T cells / C CD8 T cells）/ `S6_Fig`、`fig:supp-method-cor` | `Plot14_MethodCorr.R`；r、轴范围和显示比例在 FIG `Figure_S6_cor.tsv` | `48801737085b` | 5 个有 effect 的方法两两的密度矩阵 + Spearman r；每页（每个 panel）上 effect、下 signed −log10 P 两个矩阵（页内不标字母，只写小标题），panel 字母写在页标题。2026-11-06 起为一张图（不再用 S6.1–S6.3）：图前单独一页 “Legend for S6 Fig”，第 1 页带编号 caption（A），第 2、3 页 `\ContinuedFloat` + `\caption*{S6 Fig (continued). B / C …}`，后接 `\FloatBarrier`；轴范围（±1.5 / 2.5 / 0.75 / 0.35；±11）和舍去比例 0.3–1.9% 来自 tsv |
| S7 Fig | `Figure_S7.pdf`（3 页 = panel A Bn / B CD4 T cells / C CD8 T cells）/ `S7_Fig`、`fig:supp-manhattan` | `Plot15_Manhattan.R`；计数和标注基因在 FIG `Figure_S7_data.tsv` | `72f0ca82a903` | 6 个方法的 Manhattan 图，每页 6 条轨道；Bonferroni 线；过线点 causal 橙 / marginal 绿；每轨道最多 8 个基因标签；Fig 2 例子洋红菱形。与 S6 相同的单图 + 续页 caption + 共用 legend 页（2026-11-06 起，不再用 S7.1–S7.3） |
| S8 Fig | `Figure_S8.pdf` / `S8_Fig`、`fig:supp-grn` | `Plot9_GO_GRN.R`；`Results/Validation/Valid5_GRN/` | `8541635dc9d0` | Bn 和 CD8 T cells 的 GRN OR（口径同 Fig 3E；cell-type 864 / 1,364、high-quality 817 / 382 个 pair；2026-11-04 从 Fig 3 移出） |
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
- 2026-11-06 | Claude Code (claude-opus-5-5) | Overleaf push（81c4dc5，含上次同步后积累的全部稿件改动）；GitHub 提交 496efc8 并 push（.gitignore 加 __pycache__/） | push 输出确认；git status 干净
- 2026-11-06 | Claude Code (claude-opus-5-5) | Overleaf pull（7591385）：作者改写 `results.tex` 中 Fig 1–3 的 caption（未完成），服务器同步；未改文字 | pull 前本地与克隆一致（status 无差异）；`build_ms.sh` 返回 0，两份文档无 undefined；figures 与 FIG 同名 PDF 逐字节一致
- 2026-11-06 | Claude Code (claude-opus-5-5) | Fig 3 改版（`Plot15_Manhattan.R` + `Plot9_GO_GRN.R`）：A 改为三条同向 Manhattan（加 Marginal pseudobulk、每条 15 个基因名），C / D 互换，HiChIP 与 GRN 同一行等宽对齐；caption A 加 “marginal pseudobulk”，第 26 行 `\TODO` 中 “Miami plot” → “Manhattan plots”；`Figure_S7.pdf` / `Figure_S8.pdf` 重新生成（内容不变）一并复制 | PNG 预览逐轮检查（标题不重叠、标签不出界、D/E 对齐）；`build_ms.sh` 返回 0，无 undefined，Fig 3 与 caption 在第 9 页；图 hash 与 FIG 一致（Fig 3 `c1663b4aa269`）
- 2026-11-06 | Claude Code (claude-opus-5-5) | Fig 3 第 2–3 轮定稿进论文：字母 A–E → A–D（GO BP / MF 合为 B，HiChIP = C，GRN = D）；caption 改为 “B, … (left) … (right)”、C、D；`results.tex` 第 26、28、30 行、`methods.tex` 第 329 行、`supplementary.tex` S5 / S8 caption 中的 Fig 3 字母同步；复制 Figure_3 / S5 / S7 / S8（S8 high-quality 颜色改 #58B8D1，S5 / S7 内容不变） | `build_ms.sh` 返回 0，两份文档无 undefined；转 PNG 查看第 9–10 页；4 个 PDF hash 与 FIG 一致
- 2026-11-06 | Claude Code (claude-opus-5-5) | S6.1–S6.3 / S7.1–S7.3 Fig 合并为 S6 Fig / S7 Fig，三页作为 panel A / B / C（Bn / CD4 T cells / CD8 T cells）：`Plot14` / `Plot15` 页标题加字母，S6 页内矩阵去掉 A / B 字母；`supplementary.tex` 删 `\SixSub` / `\SevenSub`，改为第 1 页编号 caption + 第 2–3 页 `\ContinuedFloat` 续页 caption，legend 中 “A / B” 改写为 top / bottom（effect estimates / signed −log10 P）；`supporting_info.tex` 改为 `S6_Fig`、`S7_Fig` 两个条目；`results.tex` 第 26 行 TODO 改引 `\nameref{S7_Fig}` | 全文无 S6.x / S7.x 残留；`build_ms.sh` 返回 0，无 undefined；补充 PDF 第 13–20 页转 PNG 查看；hash 与 FIG 一致
- 2026-11-06 | Claude Code (claude-opus-5-5) | Overleaf push（7591385 → 10f2aa3）：本日 Fig 3 定稿、S6 / S7 合并的 tex（`results.tex`、`methods.tex`、`supporting_info.tex`、`supplementary.tex`）和 Figure_3 / S5–S8 PDF；作者接着在 Overleaf 润色文字 | push 前 status 显示 Overleaf 无新提交、差异只有这 9 个文件；push 后 status 无差异
- 2026-11-09 | Claude Code (claude-opus-5-5) | Overleaf pull（10f2aa3 → 7d3fbfc）：作者在 Overleaf 改写 `results.tex`（Results 第 1 节改为 MR 框架 + 模型公式，第 2–3 节重写，Fig 3 / Fig 4 caption）、`introduction.tex`（加一句占位）、`tables/motif.tex` 脚注，`references.bib` 新增 7 条（su2025scmultimap、calderon2019landscape、sherman2022david、ogrady2022splicing、zisman2025splicing、fang2001cblb、yamamoto2006ubc13） | pull 前 status 显示 Manuscript/ 与克隆无差异；`build_ms.sh` 返回 0，无 undefined reference / citation
- 2026-11-09 | Claude Code (claude-opus-5-5) | Methods 补写外部数据来源：M1 Data sources 写明 10x granulocyte-sorted 10k PBMC multiome，使用 TRIPOD（jiang2022tripod）处理的版本；M8 补 QC 后 7,790 细胞、CITE-seq 标签转移（162,000 细胞 / 228 抗体）、14 个细胞类型、sctransform / TF-IDF、WNN 80 个 metacell、求和 + size factor；M10 重写 4 个参考的证据类型（FANTOM5 eRNA、4DGenome 染色质互作、EnhancerAtlas 2.0 enhancer 注释 + 靶基因预测、GTEx whole-blood cis-eQTL），删 2 个 TODO（M8 只留 Signac 版本）；`references.bib` 新增 7 条（jiang2022tripod、hao2021seurat、andersson2014atlas、teng2015fourdgenome、wang2019hacer、gao2020enhanceratlas、gtex2017genetic）；10x 页面、sctransform、Signac 留 `[CITE]` | 事实来自 TRIPOD 全文（PubMed PMC9509445）和 GitHub 脚本，文献元数据经 PubMed 核对；invariants（快照 `tmp/ms_rounds/20261109_methods_data`）只多出新数字和新 cite；`build_ms.sh` 返回 0，无 undefined，7 条都在 main.bbl
- 2026-11-09 | Claude Code (claude-opus-5-5) | Methods M10 Hypergeometric enrichment 重写：原第 2 段改为 5 段（Overlap counts：背景 = 工具变量筛选前的 1 Mb 候选 peak、≥1 bp 重叠、纳入条件；Pooled hypergeometric test：phyper log.p、fold、不校正；Links of each method：108 个检验 → S1/S2 Table；MR-only 3 CT 合并：N / n_sig / HiChIP·PCHi-C 取并集 → Fig 3C、S5 Fig；Stratification：按 K/N 分 8 组、按 gene bootstrap 2,000 次百分位 CI）；留 `[CITE: GenomicRanges]` | 细节逐条对照 `Valid1_Hyper.R`、`Valid1_Hyper_MROnly.R`、`Plot3_Hyper.R`、`Plot13_MROnly_ObsExp.R`；invariants（快照 `tmp/ms_rounds/20261109_methods_hyper`）只多出 1 Mb / 1 bp / 108 / 2,000 / 95% 和已有 label 的引用；`build_ms.sh` 返回 0，无 undefined
- 2026-11-09 | Claude Code (claude-opus-5-5) | 10x 数据集引用：作者从 10x 官网复制页面信息，新增 `tenx2020pbmc10k`（@misc，2020-09-09 发布，Cell Ranger ARC，CC BY 4.0；url 暂用 datasets 首页），Data sources 替换 `[CITE]` 并补 “25-year-old healthy female donor” | build_ms.sh 返回 0，无 undefined；bbl 中条目正常
- 2026-11-09 | Claude Code (claude-opus-5-5) | M10 PCHi-C / HiChIP 引用：作者提供文献，新增 `javierre2016pchic`（Cell 2016）、`chandra2021pieqtl`（Nat Genet 2021），替换 2 处 `[CITE]`；句中补 “17 primary human hematopoietic cell types” / “five common immune cell types”（来自两篇摘要），HiChIP 改为 promoter 与 H3K27ac 区域之间的接触 | PubMed 核对元数据；build_ms.sh 返回 0，无 undefined；两条都在 main.bbl
- 2026-11-09 | Claude Code (claude-opus-5-5) | 记号统一：显著 peak 数 `n_{\mathrm{sig}}` 改为 `n`（gene 级 `n_g`），涉及 `methods.tex` M10、`supplementary.tex` S5 Fig caption、`tables/supplementary_hyper_{p,fold}.tex` 表注；与 Fig 3C 一致。S5 Fig 图内 y 轴仍为 `k / n_sig`（Plot13 脚本） | invariants（快照 `tmp/ms_rounds/20261109_n_notation`）只有记号变化；`build_ms.sh` 返回 0，无 undefined
- 2026-11-09 | Claude Code (claude-opus-5-5) | Overleaf push（7d3fbfc → e81840c）：`references.bib`（新增 10 条：TRIPOD、Seurat v4、FANTOM5、4DGenome、HACER、EnhancerAtlas 2.0、GTEx 2017、10x 数据集、PCHi-C、HiChIP）、`methods.tex`（M1 / M8 / M10 外部数据与验证资源，另一会话的 M10 重写）、`supplementary.tex` 和 `tables/supplementary_hyper_{p,fold}.tex`（另一会话的 n_sig → n 记号统一） | push 前 status 显示 Overleaf 无新提交；build_ms.sh 返回 0；push 后 status 无差异
- 2026-11-10 | Claude Code (claude-sonnet-5-5) | Overleaf pull（e81840c → 625d815，作者 11-10 在 Overleaf 的改动）：只改 `sections/results.tex`——小节标题改为 “Biological evidence for MR-inferred enhancer--gene links” 和 “Disease heritability enrichment of MR-supported peak–gene links”；S-LDSC 段落被作者重写（删除 `\TODO{R4-P1}` 和 OneK1K 句，新段仍待润色：“six … T1D” 列了 7 个 GWAS、“Fig 4A” 应为 `\ref{fig:disease}A`、有拼写错误）；pull 前 Manuscript/ 与上次同步一致，无覆盖；`build_ms.sh` 编译通过 | `overleaf_sync.sh status` 无差异；`git diff` 只涉及 results.tex
- 2026-11-10 | Claude Code (claude-opus-5-5) | Fig 4 更新为 FIG `Figure_4_CC1.pdf`（B 改为 2 行 × 3：Bn ERAP1 / CD4 DR1 / CD4 ERBB3；CD8 GNAS / CD8 TMBIM1 / CD8 THBS3，去掉 RGL2、POLR2E）；`results.tex` R4-P2 的 “Additional examples” 句改为新例子；1.5 节和 FIG 台账同步；caption 与 Methods 不变（6 个 pair 都是 strong MR-only link） | 新例子的 rsID、性状、FDR 与 `GWAS_MRonly_overlap.tsv` 核对；`build_ms.sh` 返回 0，两个 log 无 undefined；`figures/Figure_4.pdf` 与 FIG `Figure_4_CC1.pdf` hash 一致 `f21ddcdb962e`；改前文件快照在 `tmp/ms_rounds/fig4_20261110/`
- 2026-11-11 | Claude Code (claude-opus-5-5) | R4-P2 重写为 8 段草稿（仿 scMultiMap PICALM / INPP5D 段：导语 + Fig 4B 6 个例子各一段 + 收尾统一交代多候选基因）；新增 4 个已有 `\cite`（chandra2021pieqtl、javierre2016pchic、gao2020enhanceratlas、gtex2017genetic）和约 20 处 `[CITE]`；1.5 节补 24 个 link 的计数和草稿数字；2_BOARD R4 记待确认句和 M14 待办 | Overleaf status 无新提交；快照 `tmp/ms_rounds/r4p2_20261111/`；`build_ms.sh` 返回 0，main.log 无 undefined；invariants 只有 results.tex 的计划内新增（删掉的 “1” 来自旧句 “type 1 diabetes”）；第 5–7 页转 PNG 查看
- 2026-11-11 | Claude Code (claude-opus-5-5) | R4-P2 第 2 版：按 scMultiMap PICALM / INPP5D 段逐个例子重写（叙事式，每段以 “This implies that … may contribute to …” 收尾），导语改为 scMultiMap 开头的写法；新增文献 Albert 2010、Hsu 2015、Kiselev 2024、Zhang 2025（PMID 来自已核对文献表）；按 FAVOR 原始返回修正 rs1335679 位于 *TMED5* 内含子（距 *CCDC18* TSS 293 bp）、rs2549785 位于 lncRNA 内含子；TMBIM1 的变异到 TSS 距离用 5.4 kb（SNP 层，peak 中心为 5.3 kb） | `build_ms.sh` 返回 0，main.log 无 undefined；invariants（对 `r4p2_20261111`）只有计划内新增，另少了 `\q<0.001`（strong MR-only 定义改为只引 Methods）
- 2026-11-11 | Claude Code (claude-opus-5-5) | R4 第 3 版（作者决定）：Fig 4B 改为 1 行 3 个 pair（Bn ERAP1 / CD4 ERBB3 / CD8 THBS3，行高 9 in，FIG `Figure_4.pdf`，删除试看版 CC1 / Codex1）；`results.tex` 删去 DR1 / GNAS / TMBIM1 三段和收尾段；三段补 peak hg38 坐标和外部证据数值（HiChIP raw count + Q、RegulomeDB 1f + probability、FAVOR aPC、GTEx P = 5.2e-33，原稿 “P < 0.05” 为少报）；P1 润色（seven GWAS、`\ref{fig:disease}A`、删 OneK1K）；`[CITE]` 全部换成 `\cite`，`references.bib` 新增 12 条（PMID / DOI 经 PubMed 核对） | HiChIP 读自 Chandra 2021 Supp Table 3，GTEx 读自 `gtex.rda`；`build_ms.sh` 返回 0，main.log 无 undefined，Fig 4 的 “Float too large” 警告消失；`figures/Figure_4.pdf` 与 FIG hash 一致 `a7d50daf1418`；快照 `tmp/ms_rounds/20261111_fig4b/`；Overleaf status 无新提交，push 被 auto mode 拦下，待作者执行
- 2026-11-18 | Claude Code (claude-opus-5-5) | Overleaf pull（1ba44fc → 0664cfc，作者在 Overleaf 的大改）：`methods.tex` 删去 “Study design and terminology” 小节（`sec:terminology`，类别名 / $\q$ / MR-only 的主定义随之消失），MR / Causal IVW / Causal Egger 改写（公式改为 `equation*`、记号去掉 hat），DAVID 段改为 “exports only terms with $P<0.1$” 并加运行日期；`data_preprocessing.tex` 去掉粗体小标题、合并段落、删 A1-09 TODO；`introduction.tex` 新写草稿（5 段 + 占位句）；`results.tex` R4 标题改为 “Causal enhancer–gene links provide insights into disease mechanisms”，S-LDSC 首句和 Fig 4B 三个例子段落改写；`s1_text.tex`、`supplementary.tex` 重写；`tables/supplementary_{data_overview,gwas,hyper_p,hyper_fold}.tex` 表注改写；`references.bib` +19 行 | pull 前 status 显示 Manuscript/ 与克隆无差异，无覆盖；无删除 / 新增文件；`build_ms.sh` 全部返回 0（main 32 页、supplementary 25 页）；main.log 有 2 处 undefined：methods.tex M11 / M13 仍 `\nameref{sec:terminology}`；M11 留有 `[CITE: DAVID]`
- 2026-11-18 | Claude Code (claude-fable-5-2) | Introduction 改写（`sections/introduction.tex`，作者草稿 571 词 → 6 段约 1,090 词）：与 scMultiMap / TRIPOD / FusioMR 的 Introduction 雷同的句子全部重写；P1 拆成背景和现有方法两段；补写最后一段（应用与结果，定性，无数字）；作者的主旨句改为文件顶部注释；删去文件末尾空行。`references.bib` 新增 12 条（取自三篇论文的参考文献表）：maurano2012regulatory、gasperini2020enhancers、stuart2021signac、ma2020shareseq、sakaue2024scent、daveysmith2003mendelian、lawlor2008mendelian、burgess2013summarized、richardson2020twmr、gleason2021mrmtrobin、lin2024cismrcml、kang2026fusiomr（PMID 记在不打印的 `pmid` 字段）；留 1 处 `[CITE]`（跨个体的 peak–gene 相关） | 12 条文献的标题 / 作者 / 卷期页 / DOI 逐条经 PubMed 核对；与三篇论文原文的共有 5-gram 由 81 / 65 / 111 降到 11 / 10 / 4（剩下的是专有名词），共有 8-gram 由 69 / 40 / 79 降到 1 / 3 / 0；invariants（快照 `tmp/ms_rounds/20261118_intro`）只有 Introduction 的新 `\cite`、428、`P`、`^+`；`build_ms.sh` 返回 0（main 34 页），12 条都在 main.bbl，无 undefined citation（`sec:terminology` 2 处是 pull 进来的旧问题）；结果句与 `results.tex` 和第 1 节数字核对
- 2026-11-18 | Claude Code (claude-fable-5-2) | Abstract 按全文重写（`sections/abstract.tex`，旧稿 6 句 137 词 → 10 句 292 词）：替换强于证据的 4 句，去掉 428 individuals，加入三个细胞类型、两种对照方法、link 数量和方向、eRegulon OR、6 个外部参考、MR-only peak 的 HiChIP 3.28 倍、7 个疾病的遗传力点估计和 GWAS 位点例子；`2_BOARD` 2.2 和 Abstract 小节同步 | 每句对照 `results.tex` 和第 1 节（link 数表：三个细胞类型里每个 MR 检验都少于两个 marginal 方法；同号 469/469、1,258/1,259、758/759；GRN OR 表；18 个组合全部富集；S-LDSC 7 个疾病点估计都更高）；词数 292（< 300），无 `\cite` / `\ref`；invariants（快照 `tmp/ms_rounds/20261118_abstract`）只有 +3.28、−428、`^+`；`build_ms.sh` 返回 0（main 35 页），无新增 undefined
- 2026-11-18 | Claude Code (claude-opus-5-5) | 补引用：统计全文缺引用 24 处；`references.bib` 新增 17 条（ArchR、scMEGA、Gallagher 2018、Eckart–Young 1936、sctransform、SCENIC+、JASPAR2020、Gazal 2017 baselineLD、Gazal Zenodo、rtracklayer、PLINK 1.9、1000G 高覆盖、Burgess 2016 样本重叠、Bowden 2016 NOME、Kamal 2023 GRaNIE、Liu 2015 IBD、Karczewski 2025 Pan-UKB）；introduction / results / methods / data_preprocessing / discussion / s1_text / S4 Table 的占位符换成 `\cite`，methods 两处 “(Supp file)” 改为 `\nameref{S1_Text}`；`make_gwas_table.py` 的 AS / CD 来源同步加 cite；快照 `../tmp/ms_rounds/cite_20261118/`，bib 备份 `../tmp/references.bib.before_cite_20261118` | PubMed 元数据核对 PMID / DOI（Eckart–Young、Zenodo 不在 PubMed）；main + supplementary 编译（supplementary 手动跑 bibtex）无 undefined citation；grep 无占位符残留 |
- 2026-11-18 | Claude Code (claude-opus-5-5) | 补引第 2 轮：bib 加 zhang2026ra（Crossref 核对）、taylorweiner2019tensorqtl（PubMed 核对，CIMA 第 43 号文献），删 fang2001cblb；Intro P2 加 Zhang 2026；burgess2016overlap（D5、overlap 小节）、bowden2016nome（NOME）、lawlor2008mendelian（IV 假设）进正文；Results R4 首句 yazar 引用改为 `\nameref{S4_Table}`，Methods S-LDSC 加 “AS 和 CD 为 Yazar et al. 所用版本” 一句；Methods Data sources 写明 CIMA 用 TensorQTL 1.0.10；Methods 新加 “Terminology” 小节（`sec:terminology`，3 句取自被删原文）；main / supplementary 页脚去掉日期（`\lfoot{}`）；`tmp/build_ms.sh` 为 supplementary 加 bibtex 和第 3 遍 pdflatex；快照 `../tmp/ms_rounds/cite2_20261118/` | build 全部返回 0，两份 log 无 undefined citation / reference；PDF 中不再出现日期 |
- 2026-11-18 | Claude Code (claude-fable-5-2) | Introduction 新增一段讲 marginal association 的问题（原第 2 段末 4 句拆出并改写，6 段 → 7 段，约 1,230 词）：跨细胞类型引 Xie 2023（r = 0.87、38 个神经元类型、> 23,000 个其他染色体上的 enhancer、8 vs 38 分组）和 TRIPOD；跨个体引 scMultiMap、GRaNIE Box 1，加一句作者关于 pseudobulk 的推论；`references.bib` 新增 `xie2023robust`（PMID 37492103） | 作者给的三条支持逐条对照原文：Xie 2023 和 Kamal 2023 读 PMC 全文（PMC10363915、PMC10258561），scMultiMap 读 PDF；Xie 的基因名 *Stim2* 未能从 PMC 正文确认，正文写 “one gene”；`zhang2026ra` 在 PubMed 检索不到，待作者核对；invariants（快照 `tmp/ms_rounds/20261118_intro_marginal`）只多出 23,000 / 38 / 8 / r=0.87 和两个 `\cite`；与两篇原文引句无 6 词以上相同片段；`build_ms.sh` 返回 0，`xie2023robust` 在 main.bbl，无 undefined
- 2026-11-18 | Claude Code (claude-fable-5-2) | Overleaf push（0664cfc → 4da17fe，作者要求）：11 个文件——`sections/abstract.tex`（重写）、`sections/introduction.tex`（改写 + marginal association 一段）、`references.bib`，以及服务器上已有改动的 `main.tex`、`supplementary.tex`、`sections/{data_preprocessing,discussion,methods,results,s1_text}.tex`、`tables/supplementary_gwas.tex` | push 前 status 显示 Overleaf 无新提交；`build_ms.sh` 返回 0，main / supplementary 的 log 都没有 undefined；push 后 status 无差异
- 2026-11-20 | Claude Code (claude-opus-5-5) | Overleaf pull 87165bf 后：删 Intro 重复段（"We first applied the MR framework…"）；Intro 作者–年份引用和 "(reference)" 占位改为数字编号，新增 bib `dorans2025pgboost`、`holmes2014bmi`（PubMed 核对），pseudobulk 混杂句补引 scMultiMap；Methods 3 处 "X et al." 改为研究名（TRIPOD / OneK1K）；修 Intro / Abstract / Discussion 拼写和格式；BOARD 1.1 写入引用格式 | 编译通过、无 undefined citation、正文 0 处 "et al."
- 2026-11-20 | Claude Code (claude-opus-5-5) | S1 Text 新增首节 “Marginal association and causal effect”（作者推导：β̃ = β + ψ_Xψ_Y Var(U)/Var(X)，并接工具变量比值 = β）和 Fig A（TikZ 重画）；Intro P2 第 2 句加理论说明并引用 “Fig A in S1 Text”；`supplementary.tex` 加 tikz；supporting_info 的 S1 Text 说明更新；Overleaf 7e98393（只删重复段）未 pull，仅快进本地克隆 | 编译通过、无 undefined，supplementary 第 2 页目视检查
- 2026-11-20 | Claude Code (claude-opus-5-5) | 结构调整：Results motif 小节（R5）并入 R3，放在 eRegulon 段后（首句改为衔接句）；Discussion 的 allele orientation 段并入 Methods “Allele orientation”（与原操作描述合成一段），Discussion 局限改为 First / Second；Results 4 个、Methods 6 个小节标题更新（label 不变）；BOARD 1.1 改为标题用 enhancer–gene、2.2 地图更新 | 编译通过、无 undefined；`ms_check.py invariants --base HEAD` 本轮只有移动（Discussion 的 cite / S2_Fig 移到 Methods）
- 2026-11-20 | Claude Code (claude-opus-5-5) | Overleaf push（7e98393 → 2436a76，作者要求）：10 个文件——`supplementary.tex`、`references.bib`、`sections/{abstract,data_preprocessing,discussion,introduction,methods,results,s1_text,supporting_info}.tex`（本日三轮：引用统一与拼写、S1 Text 推导和 Fig A、结构调整与标题） | push 前 status 无 Overleaf 新提交；push 后 status 两边一致
- 2026-11-21 | Claude Code (claude-opus-5-5) | Sample overlap + NOME 统一论述（按 `../tmp/sample_overlap_revision_brief.md`）：S1 Text 删独立 NOME 小节，新增 “Sample overlap and measurement error”（ρ_ε → ρ_j → Cov(γ̂_j, ε_j) → bias ≈ (κβ̃ − β)/F̄），MR model 加样本集合和残差相关，全篇统一符号（Z → G、ξ、ψ_Y U_i + ξ_Yi、估计值一律带 hat）；Methods 该节重写并改名 “Error correlation due to sample overlap”；Discussion 两段局限合并；S1 Fig 图注补全局 ρ̂；S3 Table 脚注改写 | `build_ms.sh` 返回 0，无 undefined；invariants（对 `20261121_overlap`）只有计划内变化
- 2026-11-22 | Claude Code (claude-opus-5-5) | Results R1 扩写为 7 段（IV 理由、IVW / Egger 性质、NOME 与 overlap、SMR / joint、marginal 方法思路、术语）；Methods M3–M8 补 S1 Text 的性质，Causal SMR 按 Zhu 2016 Online Methods 重写（SMR P 值用 χ²₁，已对照 `Scripts/Model/MR_AllPairs.R`），修 `indivisual` 等笔误；随后 commit + GitHub push + Overleaf push（作者要求） | `build_ms.sh` 返回 0，无 undefined；invariants（对 `20261121_R1`、`20261122_methods`）只有计划内变化
- 2026-11-22 | Claude Code (claude-opus-5-5) | W3 收卡（sample overlap + NOME、R1 扩写）：2_BOARD 删两张工作卡，长期决定写入 R1 / Discussion / M3–M8 / M17 / S1 Text 的“已定”，过时的遗留问题删除；更新 2.2（R1、Discussion、M3–M8、M17、S1 Text 小节表）和 2.3；S1 Text 中 “described by Burgess et al.” 改为不写作者名（1.1 规则） | `build_ms.sh` 返回 0
- 2026-11-22 | Claude Code (claude-opus-5-5) | Results R2–R4 审读 + 润色（计划经作者两轮批注）：R2 补各方法 link 数、MR-only 数、S4 / S6 Fig 结果、方向比例（重算更正见 1.2），拆出例子段与机制段（回扣 S1 Text Sec 1 的混杂 / 掩盖）；R3 按 panel 顺序重排（3A 描述 → GO 重写 → S1 / S2 Table + 3C + S5 → 3D → motif），删 Egger 等位方向句和重复的 MR-only 定义；R4 S-LDSC 改为点估计并补覆盖比例，三个 GWAS 例子平行化并填入 β；Fig 1–4 图注重写（Fig 3D 阈值改为 0.05，Fig 1 去掉 “Median”），motif 表注背景按正文改为 50,000 个匹配 peak | 快照 `../tmp/ms_rounds/20261122_Results/`；invariants 只有计划内差异（新数字均对照 Results 或本台账，无 cite / label 变动）；placeholders 0；`build_ms.sh` 返回 0，无 undefined
- 2026-11-22 | Claude Code (claude-opus-5-5) | Discussion 润色（只改现有 3 段）：D1 改为 estimate causal effects、补核心发现；D2 与 Results 对齐（GO 主题改为转录调控 / RNA 加工 / 泛素化，fold 写 generally larger，补 nominal S-LDSC），删 “strong biological support”、“therapeutic intervention”；局限段 “true causal effect” 改为 nonzero，拆 which + 分号长句；Results 三句收尾句作者确认 | 快照 `../tmp/ms_rounds/20261122_Discussion/`；invariants 无差异（无数字、cite 变动）；`build_ms.sh` 返回 0，无 undefined
- 2026-11-22 | Claude Code (claude-opus-5-5) | 全文时态与称谓统一：时态 4 处（Intro 末段 may differ → differed；R4 ERBB3 scores were → are；Methods multiomics had → has；S1 Sec 8 was estimated → is estimated）；称谓：正文 enhancer–gene linkages ×3 → links / regulation，Abstract / Intro 的 enhancer–gene pairs → peak–gene pairs / links，causal links → MR-only links / the MR estimates（标题按作者意见不动）；顺带修 Methods GRN 的 Fig 3E → 3D，以及 Intro 两处语法 | 快照 `../tmp/ms_rounds/20261122_Tense_Terms/`；invariants 无差异；grep linkage 只剩标题和 LD；`build_ms.sh` 返回 0，无 undefined
- 2026-11-22 | Claude Code (claude-opus-5-5) | 作者要求：git commit 22ae473 + GitHub push（2188981 → 22ae473）；Overleaf push（a6e09d7 → 54e135e） | push 前 Overleaf 无新提交
