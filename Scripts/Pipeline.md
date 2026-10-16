# CIMA peak–gene 因果连接：分析流程（Pipeline）

> 运行环境、数据路径、文件安全规则见 `Claude_project/CLAUDE.md`，这里不再重复。
> 旧代码和旧做法见 `Scripts_reference/`（`00_README_Pipeline_Summary.md`）；本文件写的是**当前采用的做法**，以本文件为准。
> 路径约定：`BASE=/scratch/user/tfcui2025/CIMA`，`PROJ=$BASE/Claude_project`，
> 旧结果 `WD=$BASE/Results/RunMR/Materials_runMR/20260430_LinkPeakGene_more/MRresult_20260911`。

## 0. 目标

用 CIMA（OMIX009786）细胞类型特异的 caQTL（peak 开放度）和 eQTL，以遗传变异为工具变量做 MR，推断 peak → gene 的因果调控连接；
再和只看相关性的 marginal 方法比较，用外部证据（Valid1–Valid5，见第 3 节）验证。

- 细胞类型：`Bn_TCL1A`、`CD4_Tn_CCR7`、`CD8_Tn_CCR7`
- 6 个方法：

| 类别 | 方法 | 列前缀 |
|---|---|---|
| Causal | GLS-IVW | `ivw_` |
| Causal | GLS-MR Egger | `egger_` |
| Causal | SMR | `smr_` |
| Causal | Cauchy 合并检验（SMR + IVW + Egger） | `cct_` |
| Marginal | WNN metacell Spearman | `wnn_` |
| Marginal | 个体层 pseudobulk Spearman | `pseudo_` |

## 1. 文件位置

| 内容 | 位置 |
|---|---|
| MR 代码 | `Scripts/Model/MR_AllPairs.R`、`Scripts/Model/MR_AllPairs_slurm.sh` |
| 验证代码 | `Scripts/Validation/` |
| QC 代码 | `Scripts/QC/`（日志 `Scripts/QC/logs/`；结果 `Results/QC/`；中间文件 `tmp/QC_OneSampleMR/`） |
| SLURM 日志 | 和脚本放一起：`Scripts/Model/logs/`、`Scripts/Validation/logs/` |
| MR 结果 | `Results/MR/<CT>_AllPairs_results.tsv.gz` |
| 验证结果 | `Results/Validation/Valid<k>_<名>/` |
| 中间文件 | `tmp/MR_chunks/<TAG>/`（chunk）、`tmp/plink/<TAG>_chunkXXX/`（PLINK 临时文件，用完即删）、`tmp/Valid2_SLDSC/`（S-LDSC 的 annot / ldscore） |
| 一次性 / 检查代码 | `tmp/`（例如对账脚本 `tmp/MR_check_vs_old.R`） |

`TAG` = 细胞类型名；冒烟测试时为 `<CT>_TEST<N>`，结果与正式结果分开。

---

## 2. MR Pipeline

代码移植自 `Scripts_reference/Server5_wgen6_MR9.sh`，核心统计函数（`dedup_keep`、`make_plink_id`、`run_plink_ld`、`run_gls_mr`、
`run_smr_pair`、`cauchy_combine_matrix`、`run_one_pair`、`spearman_pairs`）保持不变；只改了输出路径和 2.5 节的合并/过滤步骤。

### 2.1 候选 pair
1. `All_GenePeak_<CT>.rds` 与 clumped caQTL 按 (peak, gene) 取交集，得到 (peak, gene, variant)。peak 格式由 `chr-start-end` 转成 `chr:start-end`。
2. 只保留在该 gene 的 eQTL 中也出现的 variant；每个 pair 共享 variant 数 `n_variant ≥ 3`（`MIN_SNP`）。
3. gene 和 peak 必须同时出现在 pseudobulk 矩阵和 WNN 矩阵中（Seurat 里 peak 写作 `chr-start-end`）。
4. fork 之前把 caQTL / eQTL 表缩到只含这些 pair 的行并建 key，子进程共享内存。

### 2.2 单个 pair 的 Causal 分析
1. **对齐**：按 `variant_id` 合并 caQTL 和 eQTL；两边 A1/A2 必须完全一致，不一致的 SNP 直接丢弃（不做 harmonize）。
2. **近似去重**：caQTL slope 和 eQTL slope 都相差 ≤ 1e-3 的 SNP 视为重复，只保留第一个。去重后仍需 ≥ 3 个 SNP。
3. **LD**：
   - variant 转成 `CHR:BP:REF:ALT`（X→23），按 peak 所在染色体选 `data/ref_chr/hg38/chr*`（1000G EAS，585 人）。
   - `plink --extract --keep-allele-order --r square --write-snplist`，得到带符号的 r。
   - 行列顺序以 PLINK 输出为准，据此重排 QTL 表；含 NaN（单态）的行列删除，删后仍需 ≥ 3 个 SNP。
4. **Flip**：让 beta_x 全为正，beta_y 和 R 同步变号。
5. **GLS-IVW / GLS-Egger**：
   - Σ = diag(se_y) · R · diag(se_y)，特征分解。
   - 保留方向：累计方差 < 99.9% 的方向或前 5 个，且特征值 > 1e-10 × λ₁；其余方向权重置 0。
   - 对称白化 W = V·diag(λ^{-1/2})·Vᵀ，`lm(W·by ~ W·X − 1)`；Egger 的 X 多一个截距列。
   - df = n − p（IVW p=1，Egger p=2），σ² = Q / df；SE × max(1, σ)（只修正 over-dispersion）；p 值用 t(df)。
6. **SMR**：
   - top SNP = 最终 SNP 集合中 caQTL nominal p 最小者（|Z| 打破平局）。
   - z 由 p 反推（p 下限 1e-300），不用 slope_se。
   - b = beta_y / beta_x；χ² = z_x² z_y² / (z_x² + z_y²)，df = 1。

每个 pair 记录一个 status（`ok`、`partial`、`no_qtl`、`few_snp_qtl`、`plink_fail`、`few_snp_ld`、`few_snp_nan`、`eigen_fail`、`gls_fail`、`error: …`），
只打印分布到日志，不保存。

### 2.3 Marginal 分析
- **Pseudobulk**：该 CT 自己的个体层 counts，取 RNA 与 ATAC 列名交集的个体；按 libsize 文件做 log1p(count / lib × 1e6)；
  每个 pair `cor.test(..., method = "spearman", exact = FALSE)`，记录 n、rho、p。
- **WNN**：外部 10x PBMC 多组学的 WNN metacell（80 个），log1p 后在 metacell 层算 Spearman。三个 CT 共用同一份矩阵。
- 注意：
  - 两者都**不用任何 covariate**（年龄、性别、批次、PC 都不校正，也不先回归残差）。
  - Spearman 只看秩：pseudobulk 中真正起作用的只有 CPM 标准化，WNN 中 log1p 不影响结果。
  - WNN 数据不是 CIMA、不分细胞类型，相关性里含跨细胞类型差异。
  - pseudobulk 未按每个个体的细胞数过滤；QTL summary statistics 的 covariate 校正由 OMIX 发布方完成，不在本流程控制范围内。

### 2.4 并行与断点续跑
- pair 轮流分配到 144 个 chunk，`mclapply` 并行（核数 = `SLURM_CPUS_PER_TASK`）。
- 每个 chunk 先写 `.tmp` 再改名为 `chunk_XXX.rds`；已存在的 chunk 跳过，重投同一作业即可续跑。全部完成并合并后删除 chunk 目录。
- 环境变量 `TEST_N > 0` 时只跑前 N 个 pair（冒烟测试）。

### 2.5 合并、过滤、校正（按顺序）
1. 把 IVW、Egger、SMR、WNN、Pseudobulk 合并为一行一个 pair。
2. **删除 NA**：5 个基础方法中任一 p 值为 NA 的 pair 整行删除。
3. **删除拟合不好的 pair**：只保留 `ivw_sigma2 <= 10`。
4. **CCT**：对 (smr_p, ivw_p, egger_p) 等权 Cauchy 合并，p < 1e-15 或 > 1 − 1e-15 时用尾部近似防止 tan 溢出。
5. **FDR**：6 个 p 值分别在每个 gene 内部做 BH，得到 `*_fdr`（在过滤后的 pair 上做）。

日志中打印每一步前后的 pair / gene / peak 数，以及各方法 p<0.05、FDR<0.05、FDR<0.1 的显著数。
下游所有验证模块直接使用这张表，不再自己做 sigma2 / NA 过滤或重新做 BH。

### 2.6 输出 `Results/MR/<CT>_AllPairs_results.tsv.gz`（32 列，不含 Egger 截距和 status）

tab 分隔（旧结果虽然叫 `.tsv.gz`，实际是逗号分隔；`fread` 两种都能自动识别）。

| 列 | 含义 |
|---|---|
| `gene_id`、`peak_id` | gene 和 peak（`chr:start-end`，hg38） |
| `n_variant` | 候选阶段 caQTL 与 eQTL 共享的 variant 数 |
| `ivw_n_snp` | 去重和 LD 过滤后真正进入 GLS 的 SNP 数 |
| `ivw_n_eff` | 特征分解截断后保留的方向数 |
| `ivw_df` | IVW 自由度 = n_snp − 1 |
| `ivw_sigma2` | 残差方差 Q / df（过滤阈值 ≤ 10） |
| `ivw_beta`、`ivw_se`、`ivw_p` | GLS-IVW 因果效应、SE（已乘过度离散因子）、p |
| `egger_beta`、`egger_se`、`egger_p` | GLS-Egger 斜率、SE、p |
| `smr_variant_id` | SMR 的 top SNP |
| `smr_b`、`smr_se`、`smr_chisq`、`smr_p` | SMR 的 Wald ratio、SE、χ²、p |
| `wnn_n`、`wnn_rho`、`wnn_p` | metacell 数、Spearman rho、p |
| `pseudo_n`、`pseudo_rho`、`pseudo_p` | 个体数、Spearman rho、p |
| `cct_n`、`cct_p` | 参与合并的 p 值个数（应都为 3）、CCT p |
| `ivw_fdr`、`egger_fdr`、`smr_fdr`、`wnn_fdr`、`pseudo_fdr`、`cct_fdr` | 各方法 gene 内 BH 的 FDR |

### 2.7 运行
SLURM 配置：36 核、360G、4 h，`--array=0-2`（`CELLS=(Bn_TCL1A CD4_Tn_CCR7 CD8_Tn_CCR7)`），
日志 `Scripts/Model/logs/%x_%A_%a.out`；脚本内完整写 `module purge; module load ...`，再 `Rscript Scripts/Model/MR_AllPairs.R $CT`。

```bash
cd /scratch/user/tfcui2025/CIMA/Claude_project
mkdir -p Scripts/Model/logs Results/MR tmp
# 1) 冒烟测试（只跑 Bn 的前 200 个 pair）
sbatch --array=0 --export=ALL,TEST_N=200 --time=01:00:00 Scripts/Model/MR_AllPairs_slurm.sh
# 2) 正式运行 3 个 CT
sbatch Scripts/Model/MR_AllPairs_slurm.sh
```

步骤：
1. **冒烟**：检查日志无 `error` status；输出 32 列、无 NA、`ivw_sigma2 <= 10`、`cct_n == 3`。检查完删除 TEST 输出。
2. **正式**：3 个 CT 一起提交。
3. **对账**（`tmp/MR_check_vs_old.R`，登录节点）：旧结果 `$WD/MR_results_formal/<CT>_AllPairs_results.tsv.gz` 按 2.5 节同样规则过滤后，
   与新结果比较 pair 集合、5 个基础 p 值和 cct_p 的最大差异（应≈0）；FDR 因 BH 池变化会不同，只报告显著数变化。

旧日志中的参考值（过滤前，36 核）：

| CT | 候选 pair（QTL） | WNN + pseudobulk 均可用 | gene | peak | 个体数 | 耗时 |
|---|---|---|---|---|---|---|
| Bn_TCL1A | 46,190 | 40,389 | 8,166 | 2,970 | 335 | 12 min |
| CD4_Tn_CCR7 | 109,961 | 86,204 | 10,710 | 5,364 | 372 | 24 min |
| CD8_Tn_CCR7 | 62,136 | 51,649 | 9,443 | 3,351 | 346 | 14 min |

---

## 3. Validation（Valid1–Valid6）

> 原则：**沿用旧脚本的逻辑**，只改输入、路径、命名；这一阶段**不画图**，只产出结果表。
> 环境、数据路径见 `CLAUDE.md`；旧脚本在 `Scripts_reference/`。

### 3.0 共同约定

#### 3.0.1 输入
- MR 结果：`Results/MR/<CT>_AllPairs_results.tsv.gz`（已删 NA、`ivw_sigma2 <= 10`、6 个 p 值在 gene 内 BH）。
  各模块直接用其中的 `*_fdr`，**不再**做 `na.omit`、sigma2 过滤或重新 BH。
- 候选 pair / 背景 peak：`All_GenePeak_<CT>.rds`（peak 写作 `chr-start-end`，用时转成 `chr:start-end`）。
- 6 个方法及 FDR 列：IVW `ivw_fdr`、Egger `egger_fdr`、SMR `smr_fdr`、CCT `cct_fdr`、WNN `wnn_fdr`、Pseudobulk `pseudo_fdr`。
- 细胞类型对应的金标准列：Bn → HiChIP `Bcell` / PCHiC `nB`；CD4 → `CD4` / `nCD4`；CD8 → `CD8` / `nCD8`。

#### 3.0.2 命名与位置
模块名统一为 `Valid<编号>_<方法简写>`：

| 模块 | 内容 | 旧目录 | 旧脚本 |
|---|---|---|---|
| `Valid1_Hyper` | 金标准准备 + gene 内超几何富集 | `Valid1_2_Hyper` | `Hyper1and2_v1.sh` |
| `Valid2_SLDSC` | S-LDSC 遗传力富集 | `Valid3_SLDSC_version2/planD` | `D1_prep.R`、`D2_sldsc_slurm.sh`、`D3_merge_plot.R` |
| `Valid3_Motif` | Motif 富集 + TF 表达检查 | `Valid4_Motif` | `run_motif_enrichment.R`(+`.slurm`)、`Motif1_step2_RunAll.sh` |
| `Valid4_AUC` | 逐 gene ROC AUC | `Valid5_AUC` | `AUC2_all.sh`、`AUC2_all_slurm.sh` |
| `Valid5_GRN` | 与 CIMA SCENIC+ GRN 的重叠 | `Valid6_CIMA_GRN` | `GRN1_peak_gene.sh` |
| `Valid6_GWAS` | 只被 MR 支持的强 link 的 peak 内是否有 GWAS 显著 SNP（Fig8 用） | `Valid8_matchGWAS` | `20260525_QC1/Scripts/source/G1_find_gwas.R` |

| 内容 | 位置 |
|---|---|
| 代码 | `Scripts/Validation/Valid<k>_<名>*.R / *_slurm.sh` |
| 日志 | `Scripts/Validation/logs/` |
| 结果表 | `Results/Validation/Valid<k>_<名>/` |
| 大的中间文件（S-LDSC 的 annot、ldscore） | `tmp/Valid2_SLDSC/` |

#### 3.0.3 运行
- 一律 `sbatch`，每个模块（或模块的每一步）有自己的 `*_slurm.sh`；不同模块之间可以同时投。
- **只有计算量大的步骤才拆成多个 job**；轻的步骤一个 job、脚本内循环 3 个 CT（`CELLS=(Bn_TCL1A CD4_Tn_CCR7 CD8_Tn_CCR7)`）：

| 步骤 | 计算量 | 提交方式 |
|---|---|---|
| Valid1_Hyper | 轻 | 1 个 job（循环 3 个 CT） |
| Valid2_SLDSC `1_prep` | 轻 | 1 个 job |
| Valid2_SLDSC `2_ldsc` | 重（每个 annotation 22 个 LD score + 9 个 h2） | array 0–7，每个 (CT, 集合) 一个 task |
| Valid2_SLDSC `3_merge` | 轻 | 1 个 job |
| Valid3_Motif `1_enrich` | 重（AddMotifs，150G） | array 0–2，每个 CT 一个 task |
| Valid3_Motif `2_summary` | 轻 | 1 个 job |
| Valid4_AUC | 重（6 个金标准 × 几千个 gene 的 ROC） | array 0–2，每个 CT 一个 task |
| Valid5_GRN | 轻 | 1 个 job |
| Valid6_GWAS | 轻（读 9 个 sumstats + 22 个 bim，约 3 min） | 1 个 job |

- **唯一的跨模块依赖**：`Valid4_AUC` 读取 `Valid1_Hyper` 保存的金标准 → 用 `--dependency=afterok:<Valid1 的 jobid>` 提交。
- **模块内部的先后**（同样用 `--dependency=afterok` 串起来）：Valid2 `1_prep` → `2_ldsc` → `3_merge`；Valid3 `1_enrich` → `2_summary`。
- 提交示例：
  ```bash
  cd /scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation && mkdir -p logs
  j1=$(sbatch --parsable Valid1_Hyper_slurm.sh)
  sbatch --dependency=afterok:$j1 Valid4_AUC_slurm.sh
  j2a=$(sbatch --parsable Valid2_SLDSC_1_prep_slurm.sh)
  j2b=$(sbatch --parsable --dependency=afterok:$j2a Valid2_SLDSC_2_ldsc_slurm.sh)
  sbatch --dependency=afterok:$j2b Valid2_SLDSC_3_merge_slurm.sh
  j3=$(sbatch --parsable Valid3_Motif_1_enrich_slurm.sh)
  sbatch --dependency=afterok:$j3 Valid3_Motif_2_summary_slurm.sh
  sbatch Valid5_GRN_slurm.sh
  sbatch Valid6_GWAS_slurm.sh
  ```

---

### 3.1 Valid1_Hyper：金标准准备 + gene 内超几何富集

**目的**：看每个方法判为显著的 peak，是否比该 gene 的全部候选 peak 更多地落在已知的调控区域里。

**输入**：MR 结果表；All_GenePeak；`$BASE/Data/PBMC_validation/` 下的 6 个金标准；EnsDb.Hsapiens.v86。

**提交**：1 个 job，脚本内依次处理 3 个 CT（金标准文件只读一次，再按 CT 取对应的 HiChIP / PCHiC 列）。

**步骤 1：金标准准备**（每个 CT 一份，结果保存给 Valid4_AUC 用）
- `gene_universe` = MR 结果表中的 gene ∩ All_GenePeak 中的 gene。
- 6 个库都整理成"按 gene 分组的 GRanges（GRangesList，names = gene），hg38"，只保留 `gene_universe` 中的 gene：

| 库 | 来源 | 处理 |
|---|---|---|
| EnhancerAtlas | `enhanceratlas.combined.grl.rds` | 直接用 |
| FANTOM5 | `fantom.grl.rds` | 直接用 |
| 4DGenome | `four.d.genome.grl.rds` | 直接用 |
| GTEx | `gtex.rda`（对象 `gtex`） | `meta.pval < 0.05`；去掉与 EnsDb v86 转录本（`seqlevelsStyle = "UCSC"`）重叠的 SNP（即 gene body 内的）；去掉过滤后为空的 gene |
| HiChIP | `41588_2020_745_MOESM3_ESM.xlsx`，sheet "Supplementary Table 3"，`skip = 3`，指定 27 个列名，删掉第 1 行残留表头 | 取该细胞 `_Q < 0.1` 的 loop；区域用 `Interact_Chr/Start/End`，gene 用 `GeneName`；去重 |
| PCHiC | `PCHiC/PCHiC_peak_matrix_cutoff5.txt` | 该细胞 score ≥ 5；区域用 other end（`oeChr` 加 `chr` 前缀、`oeStart`、`oeEnd`）；`baitName` 按 `;` 拆成多个 gene，去掉 `.`；去重 |

- HiChIP 和 PCHiC 是 hg19，用 `hg19ToHg38.over.chain` liftOver：没映射上的丢掉；映射到多条染色体的丢掉；映射成多段的取首尾 `range`。日志里打印保留比例。
- 保存：`Results/Validation/Valid1_Hyper/<CT>_gold_standard.rds`（list，6 个元素，名字为 `enhancerAtlas`、`FANTOM5`、`4DGenome`、`GTEx`、`HiChIP`、`PCHiC`）。
- 日志打印每个库：覆盖的 gene 数（占 universe 比例）、每个 gene 中位区间数、中位宽度。

**步骤 2：超几何检验**
- 每个方法的显著 peak：`*_fdr < 0.1`，按 gene 分组。
- 只检验"至少被一个方法判为显著"的 gene。
- 对每个 gene、每个库（该 gene 在库中有记录才测）、每个方法（该 gene 在此方法下有显著 peak 才测）：
  - 罐子 N = 该 gene 在 All_GenePeak 中的全部候选 peak；K = 其中与该 gene 的金标准区域重叠的 peak 数；
  - n_sig = 该方法在该 gene 中的显著 peak 数；k = 其中重叠的数；
  - `pval = phyper(k − 1, K, N − K, n_sig, lower.tail = FALSE)`；K = 0、n_sig = 0 或 N = K 时 pval = 1。
  - 同时算 MR-tested 版本：N_qtl = 该 gene 的候选 peak 中进入了 MR 表的数量，K_qtl = 其中重叠的数，`pval_qtl` 同样公式。
- peak 坐标：`chr:start-end` 直接按数字建 GRanges（沿用旧代码）。
- 日志检查：MR 表的 peak 落在自己 gene 的 All_GenePeak 背景中的比例（应 ≈ 100%）。

**汇总口径（画图/做表时使用）**：`hyper_long.tsv` 是 gene 级结果。Fig3 对每个 CT × 方法 × 金标准，把各 gene 的 N / K / n_sig / k 分别求和，用合并计数做一次超几何检验得 **global P**，并给 **pooled fold enrichment** = (Σk / Σn_sig) / (ΣK / ΣN)（显著 peak 的金标准命中率 ÷ 候选 peak 的金标准背景率）。只用 N > 0、K > 0、n_sig > 0、N > K 的行。完整定义和读表注意见 `Scripts/Visualization.md` 3.3 节。

**输出**（`Results/Validation/Valid1_Hyper/`）
- `<CT>_hyper_long.tsv`：`gene_id, method, database, pval, pval_qtl, N, K, n_sig, k, N_qtl, K_qtl`（本模块的结果）
- `<CT>_gold_standard.rds`：只为 Valid4_AUC 保存的金标准

**与旧脚本的差异**：读新表，不再 `na.omit`；金标准另存一份给 Valid4；不画直方图。

---

### 3.2 Valid2_SLDSC：S-LDSC 遗传力富集

**目的**：看 MR-Only peak（只有 MR 认为显著的 link）是否富集免疫疾病的遗传力，与全部候选 peak 比较。

**分 3 步**（旧 D1 → D2 → D3）：

**步骤 1 `Valid2_SLDSC_1_prep.R`（准备 annotation；1 个 job，3 个 CT 一起，22 条染色体用 mclapply 并行）**
- 每个 CT 2 个集合（2026-10-14 起）：
  - `MROnly`：IVW / Egger / SMR 任一 `*_fdr < 0.1`，且 WNN、Pseudo 的 `*_fdr` 都 ≥ 0.1 的 pair → 收集 peak（link 层定义，跨 gene 去重）；
  - `pool`：该 CT 在 All_GenePeak 中的全部 peak（Fig5 的 All peaks 对照；与 FDR 无关）。
- 另加第 4 个“细胞类型” `Combined`（2026-10-14）：MROnly = 3 个 CT 各自 MROnly peak 的并集，pool = 3 个 CT pool 的并集（hg19 id 取并集，扩展后 `reduce` 合并跨 CT 的重叠区间）。
- 历史：2026-10-13 曾跑 6 个方法各自 `*_fdr < 0.1` 的集合 + pool，结果表暂存为 `SLDSC_6methods_FDR0.1.tsv`、`set_stats_6methods_FDR0.1.tsv`（jobs 20138809/10/11；对应的 h2 文件和 annot 已删）。2026-10-14 测试过 MR-All（MR 任一显著）和 peak 层 MR-Only，未采用。
- hg38 → hg19 liftOver（`hg38ToHg19.over.chain`）：peak_id 起点按 0-based 处理（+1 后建 GRanges）；只保留一对一映射且在 chr1–22 上的。
- **只做一种扩展**：每个 peak 两侧各扩 500 bp（501 bp 的 peak → 约 1501 bp），然后 `reduce` 合并重叠区间。不跑不扩展的版本。
- 用 1000G EUR plink 的 `.bim`，为每条染色体写单列 thin-annot（SNP 落在区间内 = 1）。
- 共 4（3 CT + Combined）× 2 集合 = 8 个 annotation，每个 22 个文件 → `tmp/Valid2_SLDSC/annot_ld/<CT>/<CT>_<集合>.<chr>.annot.gz`。
- 输出 `Results/Validation/Valid2_SLDSC/set_stats.tsv`（3 个 CT 合并）：每个集合的 peak 数、bp、覆盖 SNP 数和比例、在 pool 中的比例。

**步骤 2 `Valid2_SLDSC_2_ldsc_slurm.sh`（array 0–7，每个 task 一个 (CT, 集合)，task id = 2 × CT + 集合；CT 含 Combined）**
- LD score：22 条染色体并行，`ldsc.py --l2 --bfile 1000G.EUR.QC.<chr> --ld-wind-cm 1 --annot ... --thin-annot --print-snps listHM3.txt`；22 个都在就跳过。
- h2：9 个 GWAS 并行，`ldsc.py --h2 <gwas>.sumstats.gz --ref-ld-chr <annot>.,baselineLD. --w-ld-chr weights.hm3_noMHC. --overlap-annot --frqfile-chr 1000G.EUR.QC. --print-coefficients`；`.results` 已存在就跳过。
- 即：MROnly + baselineLD 跑一次；pool + baselineLD 单独跑一次（pool 的 ldscore 和 `.results` 已存在时自动跳过）。
- 9 个 GWAS：

| 标签 | 文件前缀 |
|---|---|
| AS | `CUSTOM.AS` |
| CD | `CUSTOM.CD` |
| IBD | `PASS.Inflammatory_Bowel_Disease.deLange2017` |
| RA | `PASS.Rheumatoid_Arthritis.Saevarsdottir2022` |
| SLE | `PASS.Lupus.Bentham2015` |
| MS | `PASS_Multiple_sclerosis` |
| T1D | `PASS.Type_1_Diabetes.Chiou2021` |
| Asthma | `GBMI.Asthma` |
| PBC | `PASS.Primary_Biliary_Cirrhosis.Cordell2015` |

- 资源（沿用）：22 核、150G、12 h。`.results` 和 `.log` 写到 `Results/Validation/Valid2_SLDSC/h2/<CT>/`。

**步骤 3 `Valid2_SLDSC_3_merge.R`（合并）**
- 读每个 `.results` 的第 1 行（`L2_0`，即我们的 annotation）：Prop._SNPs、Prop._h2(±se)、Enrichment(±se, p)、Coefficient(±se, z)。
- 从 `.log` 解析总 h2g(±se) 和 intercept(±se)。
- 派生：`tau_p`（双侧）、`tau_p1`（单侧，τ > 0）、`sd_annot = sqrt(p(1−p))`、`tau_star = τ · sd_annot · M / h2g`，M = 5,961,159（`tau_star_se` 同理，忽略 h2g 的不确定性）。
- 检查缺失组合（应为 4 CT × 2 集合 × 9 GWAS = 72 行）。
- 输出 `Results/Validation/Valid2_SLDSC/SLDSC_all.tsv`。

**与旧脚本的差异**：读新表的 `*_fdr`，不再自己过滤和重新 BH；只跑两侧各扩 500 bp 的版本（旧脚本还跑了不扩展的版本）；GWAS 从 13 个减为 9 个，标签直接用短名（IBD、RA、SLE）；旧文件名里的 `0911` 和扩展后缀去掉；不画网格图。

---

### 3.3 Valid3_Motif：Motif 富集 + TF 表达检查

**目的**：看各方法显著 peak 中富集的 TF motif，以及这些 TF 在该细胞类型中是否表达。

**步骤 1 `Valid3_Motif_1_enrich.R`（array 0–2；沿用 8 核、150G、12 h）**
- 前景：每个方法 `*_fdr < 0.1` 的 peak（转成 `chr-start-end`）。
- universe：All_GenePeak 的 peak ∩ ATAC 对象（`atac_<CT>_qc.rds`）中的 peak；前景也取与 universe 的交集。
- 可及性：ATAC 对象 meta.features 中的 `count`（没有就用全部细胞的 counts 行和），取 `log_count = log1p(count)`；在 `AddMotifs` 之前取。
- Motif：JASPAR2020 CORE vertebrates（`all_versions = FALSE`），`AddMotifs`（BSgenome.Hsapiens.UCSC.hg38）。
- 背景：每个方法单独用 `MatchRegionStats` 从 universe 中按 `GC.percent` + `log_count` 匹配 5 万个 peak，种子 = 2026 + 方法序号（IVW 1、Egger 2、SMR 3、CCT 4、WNN 5、Pseudo 6）。
- `FindMotifs(features = 前景, background = 背景)`。
- 输出 `<CT>_motif_enrichment.tsv.gz`：Signac 的全部列 + method、rank、n_fg、n_bg、前景/背景的 GC 和 log_count 中位数。

**步骤 2 `Valid3_Motif_2_summary.R`（3 个 CT 一起）**
- 显著 motif：`p.adjust < 0.05` 且 `fold.enrichment > 1`。
- 每个 CT × motif 按 6 个方法是否显著分组：`Shared_all6`（4 个 MR + 2 个 marginal 全部显著）、`MR_only`、`Marginal_only`、`Partial`。
- motif → TF：去掉 `(var.N)`，按 `::` 拆异二聚体，转大写。
- TF 表达：RNA 对象（`rna_<CT>_qc.rds`）的 counts layer 实为 log1p(CP10k)；算每个 TF 的表达细胞比例 `pct_expr`、表达细胞中的平均 CP10k、全部细胞平均 CP10k。
- 异二聚体取表达比例最低的亚基（limiting TF）；`expressed = 全部亚基都在 RNA 中 且 pct_expr ≥ 0.05`。
- 输出（`Results/Validation/Valid3_Motif/`）：
  - `Motif_sig_3ct.tsv`：显著 motif
  - `Motif_expression_motif.tsv`：每个 CT × motif 的分组和表达

**与旧脚本的差异**：读新表，删掉旧代码里的 sigma2 过滤（旧代码保留了 sigma2 为 NA 的 pair，新表里已经没有这种 pair）；不画 MotifPlot；汇总部分从交互式代码改为脚本。

---

### 3.4 Valid4_AUC：逐 gene ROC AUC

**目的**：在同一个 gene 内，看各方法的 −log10(p) 能否把"与金标准重叠的 peak"排在前面。

**依赖**：先跑 `Valid1_Hyper`，读取 `Results/Validation/Valid1_Hyper/<CT>_gold_standard.rds`。

**步骤**（array 0–2；沿用 16 核、64G）
- 金标准：读 Valid1 保存的 6 个库，展开成带 `gene` 列的 GRanges。
  **GTEx 额外加一步过滤**（沿用旧 AUC 脚本）：只保留与 validation multiome 的 peak（`1_pbmc_multiomics.reprocessed.rda` 中 `pbmc@assays$ATAC@ranges`）重叠的 SNP。
- MR 结果的 peak 用 `StringToGRanges(sep = c(":", "-"))` 转成 GRanges。
- 打标签：在同一个 gene 内，与该 gene 金标准区域重叠的 peak 为正例（TRUE），其余为负例。
- 对每个 gene（正负例都有才算）、每个方法：score = −log10(p)（只用有限且 > 0 的 p），`pROC::roc(response, predictor)`，默认 `direction = "auto"`；记录 AUC、peak 数、正例数、负例数。
- gene 层用 `mclapply` 并行。
- 输出 `Results/Validation/Valid4_AUC/<CT>_AUC_long.tsv.gz`：`cell_type, gold_standard, gene, method, auc, n_peaks, n_pos, n_neg`。

**与旧脚本的差异**：金标准不再自己构建，改读 Valid1 的结果（因此 GTEx 的 p 阈值统一为 Valid1 的 `< 0.05`，旧 AUC 脚本是 `<= 0.05`；HiChIP/PCHiC 的 gene 范围统一为 Valid1 的 `gene_universe`）；每个 CT 一个长表，不再每个金标准一个 rds。

---

### 3.5 Valid5_GRN：与 CIMA SCENIC+ GRN 的重叠

**目的**：看各方法的显著 pair 是否更多地出现在 CIMA 自己的 GRN 中。

**输入**：3 个 CT 的 MR 结果表；`Valid6_CIMA_GRN/Data/` 下的 3 个 csv（路径见 `CLAUDE.md`）。GRN 的 `Region` 格式为 `chr:start-end`、`Gene` 为基因名，与 MR 表一致（已核实）。

**步骤**（1 个 job，3 个 CT 一起）
- `CIMA_eRegulons_Metadata.csv`：去掉第一列行号；`eRegulon = Region_signature_name` 去掉 `_(Nr)` 后缀；保留 `TF, eRegulon, Region, Gene, Consensus_name, R2G_importance, R2G_rho`，去重。
- `CIMA_eRegulons_AUC_Large_than_1_in_cell_type_L4.csv`：取 `cell_type_l4` 为 3 个 CT 的行，`eRegulon_list` 按逗号拆开（Activator 和 Repressor 都要）。
- `CIMA_Selected_High_Quality_eRegulons.csv`：宽表转长表，去掉空值。
- 两套 GRN：
  - `CT` = 该 CT 中 active（AUC > 1）的 eRegulon；
  - `CT_HQ` = `CT` ∩ 对应 lineage 的高质量 eRegulon（Bn → `selected_eRegulons_B`，CD4 → `selected_eRegulons_CD4T`，CD8 → `selected_eRegulons_CD8T`）。
- GRN 的 (Region, Gene) 与 MR pair 对齐，pair 层去重（同一个 pair 可能对应多个 TF）。
- 每个 GRN 集合 × CT × 方法：显著 = `*_fdr < 0.05`；背景 = 该 CT 全部 MR pair；2×2 表（在 GRN 中 / 不在 × 显著 / 不显著），`fisher.test` 得 OR 和 p；另记 `prop_GRN`、`prop_All`、`fold_vs_All`。
- 日志检查：AUC 表和 HQ 表中的 eRegulon 名能在 metadata 中找到的比例；每套 GRN 的 eRegulon 数、TF 数、pair 数、与 MR 重叠的 pair 数。

**输出**（`Results/Validation/Valid5_GRN/`）
3 个 CT 合并在同一张表里，用 `cell_type` 列区分：
- `GRN_sig_counts_by_method.tsv`
- `GRN_summary.tsv`
- `GRN_pairs_inMR.tsv.gz`：与 MR 重叠的 GRN pair（含 TF 列表和 MR 结果）

**与旧脚本的差异**：读新表；去掉旧代码里读 Motif 结果但未使用的部分；不画柱状图；交互式代码改为脚本（`sbatch`）。

---

### 3.6 Valid6_GWAS：MR-only 强 link 与 GWAS 信号的重合

**目的**：找出只被 MR 支持（marginal 方法不显著）的强 peak-gene link 中，peak 本身含有免疫疾病 GWAS 全基因组显著 SNP 的那些，作为 Fig8 的位点图对象。

**输入**：3 个 CT 的 MR 结果表；`$BASE/R_package/S_LDSC/GWAS_data/selected_sumstats/` 中的 9 个 GWAS（与 Valid2 相同的 9 个，标签 AS、CD、IBD、RA、SLE、MS、T1D、Asthma、PBC）；
1000G EUR bim（hg19，`S_LDSC/reference_genome/1000G_EUR_Phase3_plink/`）；`Data/PBMC_validation/hg19ToHg38.over.chain`；`R_package/my_plink/plink_hg38/EAS_hg38.bim`。

**步骤**（`Valid6_GWAS.R`，1 个 job，3 个 CT 一起）
1. **MR 强的 pair**：IVW / Egger / SMR / CCT 任一 `*_fdr < 0.001`；其中 WNN 和 Pseudo 的 FDR 都 ≥ 0.05 记为 `MR_only`，否则 `Both`。
2. **GWAS SNP**：sumstats 只有 rsID 和 Z，P = 2Φ(−|Z|)；先取 P < 1e-5 建缓存：
   rsID → hg19（1000G EUR bim）→ liftOver 到 hg38（只保留一对一、chr1–22）→ 按位置 + 等位基因（两种方向）匹配 `EAS_hg38.bim` 得到 `CHR:BP:REF:ALT`（匹配不上的用 `CHR:BP`，并记 `in_EAS_ref = FALSE`）。
   缓存 `tmp/Valid6_GWAS/GWAS_SNPs_P1e-5_hg38.tsv.gz`，存在即跳过（改 trait 或阈值要删缓存）。抽查 rs2476601（PTPN22）hg38 应为 chr1:113834946。
3. **重合**：只留 P < 5e-8 的 SNP；SNP 的 hg38 位置落在 peak `[start+1, end]` 内（不扩展）。
4. **输出**：`MR_only` pair × 重合 SNP × trait。
5. **对照**：与旧 `Valid8_matchGWAS/G1_Focus_MRonly_GWAS.tsv` 在 pair 层比较（只打印）。

**输出**（`Results/Validation/Valid6_GWAS/`）
- `GWAS_MRonly_overlap.tsv`：每行 = MR_only pair × 重合 SNP × trait（列与旧 `G1_Focus_MRonly_GWAS.tsv` 相同；Fig8 读这张表）
- `GWAS_pair_counts.tsv`：每个 CT × category（Both / MR_only）：MR 强的 pair 数、其中 peak 含 GWAS SNP 的 pair 数

**与旧脚本的差异**：读新表，不再做 sigma2 过滤；GWAS 从 13 个减为项目的 9 个（去掉 CUSTOM.IBD、CUSTOM.SLE、Celiac、RA Ishigaki；IBD / SLE / RA 改用 deLange2017 / Bentham2015 / Saevarsdottir2022，标签用短名）；输出改名并加计数表；R 的临时目录放项目 `tmp/`。

**2026-10-09 结果**（job 20080635，2.5 min）：MR 强的 pair Bn 1,541 / CD4 4,459 / CD8 2,473，其中 MR_only 359 / 1,130 / 609；
MR_only 且 peak 含 GWAS SNP 的 pair 共 27 个（Bn 7 / CD4 10 / CD8 10；21 个 gene、14 个 peak；37 行）。
与旧结果相比：旧 29 个中 27 个保留、没有新增；少掉的 CD4 CERS2、CD4 TNNT3 在旧版中只由 CUSTOM.IBD 支持。

---

## 4. QC

### 4.1 QC_OneSampleMR：null SNP 上 caQTL 与 eQTL z-score 的相关性

**目的**：CIMA 的 caQTL 和 eQTL 来自同一批个体（one-sample MR）。选出两边都不显著（p > 0.05）的近独立 SNP，看每个 peak–gene pair 内 caQTL z 与 eQTL z 在 null 下是否仍然相关。
如果整体相关明显偏离 0，说明样本重叠或共享的技术因素会给 MR 带来额外相关，需要在解释结果时考虑。
复现旧脚本 `Scripts_reference/Server2_ApplyMR0_Prune.sh`（旧 Bn 值：全基因组整体 ρ = 0.0029，p ≈ 9.6e-31）。

**输入**（细胞类型：`Bn_TCL1A`、`CD4_Tn_CCR7`、`CD8_Tn_CCR7`）
- caQTL 全量：`$BASE/Data/OMIX009786/cis_caQTL_Summary_Statics/cis_caQTL/<CT>.txt.gz`（列：`phenotype_id`（`chr1:629700-630201`）、`variant_id`（`chr1_266093`）、`pval_nominal`、`slope`、`slope_se`、`A1(ALT/effect allele)`、`A2(REF)` 等）。
- eQTL 全量：`$BASE/Data/OMIX009786/cis_eQTL_Summary_Statics/cis_eQTL/<CT>.txt.gz`（列同上，`phenotype_id` 为 gene）。
- 候选 pair：`$BASE/Data/OMIX009786/gene_peak/All_GenePeak_<CT>.rds`（`gene_id`、`peak_id`，peak 写作 `chr-start-end`，用时转成 `chr:start-end`）。
- LD 参考：`data/ref_chr/hg38/chr*.{bed,bim,fam}`（ID 为 `CHR:BP:REF:ALT`）；PLINK：`$BASE/R_package/my_plink/plink`。

**步骤**
1. **LD pruning**（`QC_OneSampleMR_prune.R` + PLINK，每个 CT 一个 task）
   - 取该 CT caQTL 的全部 unique `variant_id`，用 `A1/A2` 转成面板 ID `CHR:BP:REF:ALT`（两种等位顺序都尝试匹配，匹配不上的丢掉）。
   - 逐染色体 `plink --bfile <chr> --extract <ID 列表> --maf 0.01 --indep-pairwise 500kb 50 0.1`，合并 `.prune.in`，再换回 `variant_id`。PLINK 临时文件放 `tmp/QC_OneSampleMR/`。
   - 三个 CT 一律用 hg38 面板，不再用旧脚本的 rsID + hg19 EAS 面板。
2. **选 null SNP 并配对**（`QC_OneSampleMR.R`）
   - caQTL、eQTL 各自取 pruned 且 `pval_nominal > 0.05` 的 SNP。
   - 按 `All_GenePeak_<CT>.rds` 的 (gene, peak) 配对，再按 `variant_id` 求 caQTL 与 eQTL 的共享 SNP。
3. **去重**：按 (peak, gene, caQTL slope, eQTL slope) 去重；按 slope 排序后，相邻 SNP 的 caQTL slope 与 eQTL slope 都相差 ≤ 1e-3 的视为重复，只留一个（旧脚本为 0.01，与 MR 步骤的 1e-3 统一）。
4. **相关性**
   - z = slope / slope_se。每个 pair 在 SNP 数 > 2 时算 ρ̂ = cor(caQTL z, eQTL z)，t = ρ̂·√((n−2)/(1−ρ̂²))，双侧 p，df = n − 2。
   - 全基因组整体：把所有 pair 的 SNP 合在一起算一个 ρ 和 p。
5. **图**：按 pair 内 SNP 数分组（1–10、11–20、21–30、31–40、41–50、>50）画 ρ̂ 和 P 值的箱线图，外加 ρ̂ 分布的概览页；每个 CT 一个 3 页的 PDF。

**输出**
- `Results/QC/<CT>_pair_corr.tsv.gz`：每行一个 pair（`peak`、`gene`、`snp_count`、`rho_hat`、`t_stat`、`p_value`）。
- `Results/QC/<CT>_global_rho.tsv`：整体 ρ、p、参与的 SNP 数和 pair 数（另附去重阈值 0.01 的对照）。
- `Results/Figures/Sup1_OneSampleMR_<CT>.pdf`（3 页）。

**位置**：代码 `Scripts/QC/`（`QC_OneSampleMR_prune.R`、`QC_OneSampleMR.R`、对应 `*_slurm.sh`，日志 `Scripts/QC/logs/`）；中间文件 `tmp/QC_OneSampleMR/`。

**运行**：prune 和相关性都是按 CT 的重计算，用 `--array=0-2`；相关性 job 用 `--dependency=afterok` 接在 prune 之后；画图是轻量步骤。
```bash
cd /scratch/user/tfcui2025/CIMA/Claude_project/Scripts/QC && mkdir -p logs
j1=$(sbatch --parsable QC_OneSampleMR_prune_slurm.sh)
sbatch --dependency=afterok:$j1 QC_OneSampleMR_slurm.sh
```

---

## 5. Update History
每次改动流程或代码，在这里加一行：日期 + 改了什么。
- 2026-10-09：建立 `Pipeline.md`，整合 `00_README_Pipeline_Summary.md` 与 MR 实施计划；确定 MR 代码 / 日志 / 结果 / 中间文件的位置；合并步骤为删 NA → `ivw_sigma2 <= 10` → CCT → gene 内 BH；输出只保留 32 列的 tsv.gz。
- 2026-10-09：写好 `Scripts/Model/MR_AllPairs.R`、`MR_AllPairs_slurm.sh` 和对账脚本 `tmp/MR_check_vs_old.R`；输出改为真正的 tab 分隔。
- 2026-10-09：MR 正式运行完成（job 20074030，3 个 CT，13–26 min）；最终 pair 数 Bn 40,087 / CD4 85,527 / CD8 50,931；与旧结果对账：pair 集合相同，所有 p 值和效应量差异为 0。
- 2026-10-09：Validation 部分逐个审阅定稿并写入第 3 节：模块改名为 Valid1_Hyper / Valid2_SLDSC / Valid3_Motif / Valid4_AUC / Valid5_GRN；沿用旧逻辑和旧阈值；不画图；Valid4 读 Valid1 保存的金标准；S-LDSC 只跑两侧扩 500 bp、9 个 GWAS；只有重计算步骤（SLDSC ldsc、Motif enrich、AUC）拆成 array。
- 2026-10-09：写好并运行 5 个验证模块（`Scripts/Validation/`），全部完成，结果在 `Results/Validation/`；Valid1 首次运行因 `seqlevelsStyle<-` 需联网取 UCSC chromInfo 而失败，改为本地重命名染色体（1→chr1、MT→chrM）后重跑成功。
- 2026-10-09：新增 Valid6_GWAS（移植 `G1_find_gwas.R`，GWAS 用项目的 9 个），写入 3.6 节；运行完成，27 个 MR_only + GWAS pair（旧版 29 个中保留 27 个）。
- 2026-10-10：新增第 4 节 QC，加入 QC_OneSampleMR（移植 `Server2_ApplyMR0_Prune.sh`：null SNP 上 caQTL 与 eQTL z-score 的相关性）；三个 CT 统一用 hg38 面板做 LD pruning，去重阈值 1e-3；Update History 顺延为第 5 节。尚未写代码。
- 2026-10-13：Fig3 改为 global test 表；Valid1 结果的汇总口径（global P、pooled fold enrichment）写入 3.1 节，完整定义见 `Scripts/Visualization.md` 3.3 节。
- 2026-10-13：Valid2_SLDSC 显著 peak 阈值由 FDR < 0.05 改为 < 0.1 并重跑（jobs 20138809/10/11）；pool 集合沿用旧 ldscore 和 h2；0.05 版本备份为 `SLDSC_all_FDR0.05.tsv`、`set_stats_FDR0.05.tsv`。
- 2026-10-14：Valid2_SLDSC 的集合改为 MROnly（IVW/Egger/SMR 任一 FDR < 0.1 且 WNN、Pseudo 都不显著的 link 的 peak）+ pool；6 个方法的 FDR 0.1 结果表暂存为 `SLDSC_6methods_FDR0.1.tsv`、`set_stats_6methods_FDR0.1.tsv`，FDR 0.05 备份删除。
- 2026-10-14：Valid2_SLDSC 新增 `Combined`（3 个 CT 的 MROnly 并集 / pool 并集）一组，原 3 个 CT 的结果不变；array 0–7，`SLDSC_all.tsv` 72 行。
