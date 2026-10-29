# Codex 论文科学审阅记录

本文件按论文结构保存 Codex 阅读 Methods、检查公式、讨论科学逻辑和拟定段落结构的记录。尚未讨论的章节留空，不代表 LaTeX 中没有内容。Claude Code 负责 LaTeX 整理、格式、交叉引用和编译，并维护 `Notes_Claude_Code.md`。

## Title Page

## Abstract

## Author Summary

## Introduction

## Results

### A common analysis set enables comparison of six linking methods

### Shared links can have different MR and marginal effect directions

### Selected links are enriched for regulatory-reference support

### GO summaries characterize biological processes and molecular functions

### MR-only regions show exploratory disease-heritability enrichment

### GWAS-overlapping peaks nominate disease-associated regulatory examples

## Discussion

## Materials and Methods

### Data sources

### Data Preprocessing

### Candidate peak–gene pairs

### Instrument selection and LD preparation

### MR Framework

### Inverse-variance weighted method

### MR-Egger regression

### GLS-IVW with correlated variants

### GLS MR-Egger with correlated variants

### Summary-data-based Mendelian randomization

2026-10-20：按用户要求将下述计划落实到 `sections/methods.tex` 的 SMR 小节，替换原占位说明与 TODO；包含原理简述、ratio、p-derived z、SE、卡方检验，以及协方差近似和未实施 HEIDI 的简短说明。

验证：运行现有 `../tmp/build_ms.sh`，三轮 pdfLaTeX 和 BibTeX 均成功；无未定义引用，PDF 为 32 页。仍有原先 supplementary 第 77–111 行附近的 14.22636 pt overfull，与本次 SMR 修改无关。未修改 `Notes_Claude_Code.md`。

**依据**：本地 `../paper_reference/` 中 Zhu et al. (2016) 的 Online Methods，尤其 delta-method 方差（Eq. 4）和 SMR 检验（Eq. 5）；引用键 `zhu2016smr`。实现已对照 `../Scripts/Model/MR_AllPairs.R` 的 `run_smr_pair()`。

**段落安排**：写两个简短段落。SMR 是已有方法，只概括原理，不重复完整推导。

1. **原理**：从同一工具 SNP 对 exposure 和 outcome 的汇总关联构造 ratio estimate，以 delta-method 近似方差构造 SMR 检验。将原文的 expression→trait 映射为本项目的 accessibility→expression，沿用前文的 gamma/Gamma/beta 记号。给出 ratio 和核心检验公式，不展开 2SLS 或 HEIDI 推导。
2. **实现**：从最终保留的共享 SNP 中选择 caQTL nominal p 最小者；用两侧 nominal p 值反推非负 z 幅度，设 p 值下限为 1e-300；使用自由度为 1 的卡方分布上尾计算 SMR p 值。效应方向来自已对齐的 slope ratio。标准误可紧凑说明，BH 和 CCT 留在相应小节。

**公式备忘**（不要求全部作为独立公式放入正文；下标 * 表示 lead SNP）：

```latex
\hat\beta_{\mathrm{SMR}}=\hat\Gamma_*/\hat\gamma_*.
\tilde p_X=\max(p_X,10^{-300}),\qquad
\tilde p_Y=\max(p_Y,10^{-300}).
z_X=\Phi^{-1}(1-\tilde p_X/2),\qquad
z_Y=\Phi^{-1}(1-\tilde p_Y/2).
\widehat{\mathrm{SE}}(\hat\beta_{\mathrm{SMR}})
=|\hat\beta_{\mathrm{SMR}}|\sqrt{z_X^{-2}+z_Y^{-2}}.
T_{\mathrm{SMR}}=\frac{z_X^2z_Y^2}{z_X^2+z_Y^2},\qquad
p_{\mathrm{SMR}}=\Pr(\chi_1^2\ge T_{\mathrm{SMR}}).
```

**表述要点**：代码用 `qnorm(p/2, lower.tail=FALSE)` 取得 z，SMR 标准误不直接代入输入的 `slope_se`；lead SNP 的 p 值并列时按 caQTL 的绝对 slope/SE 排序。原文一般方差包含两侧估计协方差，本实现忽略该项，不能据此称 CIMA 两侧 QTL 来自独立样本。简短说明未做 HEIDI，避免把 SMR 显著写成已区分 mediation、horizontal pleiotropy 和 linkage；不展开长篇讨论。

### Cauchy combination

### Pseudobulk correlations

2026-10-20：已按用户要求将合并小节拆分并写入 `sections/methods.tex`。Pseudobulk 按“观测单位与数据构建 → 标准化 → 相关性与检验”组织，包含 CPM/log1p、平均秩 Spearman、双侧 t 近似及未调整协变量的说明。同时更新 `data_preprocessing.tex` 的聚合顺序和细胞数 QC 定义，加入相关小节交叉引用。

- **数据与观测单位**：CIMA 内按 cell type × individual 分别汇总 RNA 和 ATAC；相关性在同一细胞类型的匹配个体间计算，不把细胞当重复样本，也不要求 RNA 与 ATAC 是同一细胞。最终匹配个体数 Bn/CD4/CD8 = 335/372/346。
- **上游证据**：`../Scripts_reference/Server1_ExploreCIMA1_scCount.sh` 的 Python 部分按个体、细胞类型累加输入 count 矩阵，随后才标准化。每个模态分别按该细胞类型内每个个体的细胞数筛选：至少 100，至多 median + 3 × MAD；MAD 使用 1.4826 缩放，再取两模态通过个体的交集。全特征汇总 counts 的总和作为各模态 library size。已核对三个 CT 的矩阵列名、libsize sample 和 valid_samples 完全一致，细胞数下限符合上述规则；未重新从原始 h5ad 验证全部生成过程。
- **标准化**：当前 `MR_AllPairs.R` 按 sample ID 匹配两模态及 library size。对 peak p、gene g、individual i 定义 A_ip = log[1 + 10^6 C^ATAC_ip / L^ATAC_i]，E_ig = log[1 + 10^6 C^RNA_ig / L^RNA_i]。分母为对应模态的全部输入特征，不是候选 pair 子集。log 是自然对数；log1p 不改变 CPM 的秩。
- **统计量**：rho^PB_pg = Corr(rank_i(A_ip), rank_i(E_ig))；ties 采用平均秩。两侧检验 H0: rho_S=0，用 `cor.test(..., method="spearman", exact=FALSE)`；本地 R 源码已核对，其 p 值采用 t = rho sqrt[(n−2)/(1−rho²)] 的 t_(n−2) 近似，不称为精确检验或置换检验。不做 covariate residualization，不按个体细胞数加权。
- **写作衔接**：聚合和样本 QC 的完整细节可放 Data Preprocessing，本小节简述并引用，避免两处重复。现有 preprocessing 的“aggregating normalized signals”应在正式写作时协调为先汇总 counts 后标准化。共同 pair 过滤、gene 内 BH 和 q<0.1 引用统一小节。若有缺失，cor.test 使用完整观测；当前输出 n 是向量长度，不能未核实就称为所有 pair 的完整观测数。

### WNN metacell correlations

2026-10-20：已写入独立 WNN metacell correlations 小节，按“外部 multiome 数据和既有分组 → feature 对齐及信号汇总 → 跨 metacell 相关性与检验”组织。原始数据引用、QC、WNN 构建参数与 RNA 聚合/标准化保留一个明确 TODO，未补造细节。未修改 Claude Code 维护的 `MANUSCRIPT_CONTEXT.md`。

本次两个小节及 preprocessing 联动修改后，三轮 pdfLaTeX 和 BibTeX 均通过，PDF 为 34 页，无未定义引用；日志仅保留原 supplementary 的 14.22636 pt overfull，没有新增公式溢出。

- **数据与观测单位**：外部 PBMC multiome 的 80 个 WNN metacells；三个 CIMA 分析共用这组矩阵，未分别只取 Bn/CD4/CD8。WNN 用于定义多模态细胞邻域/分组，最终 link score 是跨 metacell 的 Spearman，而不是 WNN 邻居权重。
- **已找到的上游证据**：`../Scripts_reference/1_Read_Data.R` 读取 `ATAC_RNA_WNN.rda`，将已有 `wsnn_res.15` 赋给 `pbmc$wnn`，沿用已有 `wnn.rna`；从外部 fragments 在 CIMA hg38 peak 区域重新计数，先排除与其他输入 peak 重叠的区域，再按既有分组汇总 ATAC counts 并做每 metacell 的 CPM。输出 `1_pbmc_multiomics.reprocessed.rda` 正是当前模型读取的路径对应文件。该脚本重算 TF-IDF/LSI 和 WNN neighbors，但没有重新聚类来替换既有 metacell 标签。
- **尚不能补写成事实的构建细节**：原始外部数据的正式来源、细胞 QC、原始 WNN 聚类步骤、80 个组的生成及 `wnn.rna` 的聚合/标准化过程仍需查其原始生成脚本。不能仅凭 `wsnn_res.15` 推断全部参数，也不能把重处理脚本的新 neighbors 写成生成现有 80 个组的步骤。当前未载入 1.4 GB RDA 做全对象审计。
- **统计量**：令 U_mp、V_mg 为已提供的 peak/RNA metacell 信号，当前代码分别 log1p，再计算 rho^WNN_pg = Corr(rank_m(log(1+U_mp)), rank_m(log(1+V_mg)))，m 遍历 80 个 metacells。不要在核实前把 RNA 矩阵也定义成 CPM。使用与 pseudobulk 相同的两侧近似检验，无 covariate residualization、无细胞类型调整、不按 metacell 细胞数加权。可引用前一小节的统计定义，避免重复推导。
- **解释与校正**：该相关性包含跨细胞类型/状态差异，不能称为各 CIMA 细胞类型内独立估计的相关性。同一 pair 的原始 WNN rho/p 在三个分析中共用；gene 内 BH 的测试集合可能不同，因此 q 值可不同。80 个数据驱动的 metacells 不能写成 80 个独立供体；名义 p 值未额外校正分组构建的不确定性或潜在组间依赖。

### Filtering and multiple testing

### Regulatory-reference enrichment

### Gene Ontology summaries

### Disease heritability and GWAS overlap

### Null-SNP diagnostic

### Allele-orientation sensitivity of MR-Egger

### Motif enrichment

## Acknowledgments

## References

## Supporting Information
