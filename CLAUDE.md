# CIMA peak–gene 因果连接项目（Claude_project）

## 1. 项目目标
用 CIMA（OMIX009786）细胞类型特异的 caQTL 和 eQTL 做 MR，推断 peak → gene 的因果调控连接，
并和 marginal 相关方法比较，用外部证据验证。
- 细胞类型只做：`Bn_TCL1A`、`CD4_Tn_CCR7`、`CD8_Tn_CCR7`
- 6 个方法：IVW、Egger、SMR、CCT（Causal），WNN、Pseudobulk（Marginal）

## 2. 文件安全规则（必须遵守）
- **只允许在 `/scratch/user/tfcui2025/CIMA/Claude_project/` 内创建、修改、删除文件。**
- `Claude_project` 以外的一切数据和文件**只能读取**，禁止任何写操作：不能修改、删除、移动、重命名、覆盖，也不能在其中新建文件或目录。
  这包括 `/scratch/user/tfcui2025/CIMA/` 下的 `Data/`、`Results/`、`R_package/` 以及旧结果目录。
- `data/ref_chr/hg38` 是指向外部目录的软链接，同样**只读**，不能通过它写入。
- 所有输出写到 `Claude_project` 内：结果放 `Results/`，中间文件和临时文件（包括 PLINK 临时文件）放 `tmp/`，SLURM 日志也放在项目内。
  不使用 `$TMPDIR` 或其他外部路径存放输出。
- 如果某一步看起来必须写到项目外，先停下来问用户。

## 3. 协作方式
- 和用户一步一步推进。每一步先说明要做什么，做完汇报结果；不要一次把整个流程跑完。
- 重计算用 `sbatch` 提交到 SLURM，不使用交互节点。登录节点只做轻量的检查和汇总。
- 旧代码和详细说明：`Scripts_reference/`，先看 `Scripts_reference/00_README_Pipeline_Summary.md`。
  旧代码不一定完美；和用户讨论后确定的新做法以总结文档里标【更新】的部分为准。

## 4. 目录
| 目录 | 用途 |
|---|---|
| `Scripts/Model/` | MR 流程代码 |
| `Scripts/Validation/` | 验证代码 |
| `Scripts_reference/` | 旧代码 + 总结 |
| `Results/` | 新结果 |
| `tmp/` | 中间文件、缓存、临时文件 |
| `data/ref_chr/hg38` | LD 参考面板（软链接，只读） |

外部数据直接从 `/scratch/user/tfcui2025/CIMA/` 读取，不复制。

## 5. 运行环境
- R：
  ```bash
  module purge
  module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2
  ```
  - 在 Claude Code 的非交互 shell 中要用 `bash -lc '...'`，并写 `module load`。`ml` 简写在这里不生效，R 会不在 PATH 中。
  - R 代码写成 `.R` 文件，再用 `Rscript file.R` 运行。复杂的 `Rscript -e "..."` 容易因为引号转义问题而静默失败，没有任何输出。
- PLINK 1.9：`/scratch/user/tfcui2025/CIMA/R_package/my_plink/plink`
- LDSC：
  ```bash
  module load Anaconda3/2025.12-2
  eval "$(conda shell.bash hook)"
  conda activate /scratch/user/tfcui2025/CIMA/R_package/S_LDSC/ldsc/ldsc_env
  ```
  然后用 `python .../S_LDSC/ldsc/ldsc.py`。
- SLURM 习惯：
  - R 中核数用 `SLURM_CPUS_PER_TASK`；多个细胞类型用 `--array=0-2`。
  - 提交前先 `mkdir -p` 日志目录（日志目录建在项目内）。

## 6. 关键数据（`BASE=/scratch/user/tfcui2025/CIMA`，全部只读）
| 数据 | 路径 |
|---|---|
| eQTL | `$BASE/Data/OMIX009786/cis_eQTL_Summary_Statics/cis_eQTL/<CT>.txt.gz` |
| caQTL（pair 级 clumping，r²=0.99，20 Mb） | `$BASE/Data/OMIX009786/cis_caQTL_Summary_Statics/cis_caQTL/LD_clump/<CT>_Clumped_r2_0.99_kb_20000.txt.gz` |
| 候选 pair | `$BASE/Data/OMIX009786/gene_peak/All_GenePeak_<CT>.rds`（data.table：`gene_id`、`peak_id`（`chr1-9928-10429`）、`chr`、`TSS`、`peak_center`、`distance`） |
| Pseudobulk | `$BASE/Data/OMIX009786/cell_rna/QC_pseudobulk/{pseudobulk_rna,pseudobulk_atac}_<CT>.tsv.gz`、`libsize_<CT>.csv` |
| 单细胞对象 | `$BASE/Data/OMIX009786/cell_rna/QC_filtered/{rna,atac}_<CT>_qc.rds` |
| WNN metacell（外部 10x PBMC，不分细胞类型） | `$BASE/Data/Valid_Multiome/1_pbmc_multiomics.reprocessed.rda` |
| LD 参考（hg38，1000G EAS，ID 为 `CHR:BP:REF:ALT`） | `data/ref_chr/hg38/chr*.{bed,bim,fam}` |
| 金标准 | `$BASE/Data/PBMC_validation/`（EnhancerAtlas、FANTOM5、4DGenome、GTEx、HiChIP xlsx、PCHiC、liftOver chain） |
| S-LDSC 资源 | `$BASE/R_package/S_LDSC/reference_genome/`（EUR plink、HM3、baselineLD v2.2、weights、frq） |
| GWAS | `$BASE/R_package/S_LDSC/GWAS_data/selected_sumstats/` |
| CIMA GRN | `$BASE/Results/RunMR/Materials_runMR/20260430_LinkPeakGene_more/MRresult_20260911/Valid6_CIMA_GRN/Data/` |
| 旧结果（对账用） | `$BASE/Results/RunMR/Materials_runMR/20260430_LinkPeakGene_more/MRresult_20260911/` |

GWAS 只用 9 个：AS、CD、IBD(deLange2017)、RA(Saevarsdottir2022)、SLE(Bentham2015)、MS、T1D(Chiou2021)、Asthma(GBMI)、PBC(Cordell2015)。

## 7. Model 要点
1. **候选 pair**：All_GenePeak ∩ caQTL ∩ eQTL，共享 variant ≥ 3，且 gene 和 peak 在 pseudobulk 与 WNN 中都存在。
2. **每个 pair**：
   - 两边等位必须一致（不一致直接丢）；两边 slope 都相差 ≤ 1e-3 的 SNP 视为重复，只留一个。
   - 用 PLINK 算带符号的 LD（`--keep-allele-order`），行列顺序以 PLINK 输出为准。
   - Flip，让 beta_x > 0。
3. **Causal 方法**：
   - GLS-IVW / Egger：特征分解截断，保留 99.9% 方差、至少 5 个方向；df = n − p；SE × max(1, σ)。
   - SMR：取 caQTL p 最小的 SNP，由 p 值反推 z。
4. **Marginal 方法**：Spearman，不用任何 covariate。
   - pseudobulk：log1p(CPM)，个体层；
   - WNN：log1p，metacell 层。
5. **合并**：
   - 5 个基础方法中任一 p 为 NA 的 pair 删除；只保留 `ivw_sigma2 <= 10`。
   - 再算 CCT（SMR + IVW + Egger），然后 6 个 p 值在每个 gene 内部做 BH。
6. **输出**：`Results/.../<CT>_AllPairs_results.tsv.gz`（不含截距和状态列）。

## 8. Validation 要点（细节逐个审阅后再更新）
| 模块 | 内容 |
|---|---|
| V1+V2 | 每个 gene 内做超几何检验：显著 peak 与 6 个金标准（EnhancerAtlas、FANTOM5、4DGenome、GTEx、HiChIP、PCHiC）的重叠 |
| V3 | S-LDSC：各方法显著 peak（hg19，±0/500 bp）+ baselineLD，9 个 GWAS，pool 作背景；看 enrichment 和 τ |
| V4 | Motif 富集：Signac + JASPAR2020，背景按 GC 和可及性匹配；检查对应 TF 是否表达 |
| V5 | 逐 gene ROC AUC，用同样的 6 个金标准 |
| V6 | 与 CIMA SCENIC+ GRN 的重叠，做 Fisher 检验 |
| Valid6_GWAS | MR-only 强 link（MR 任一 FDR < 0.001，WNN/Pseudo 不显著）的 peak 内是否有 GWAS P < 5e-8 的 SNP（9 个 GWAS）；Fig8 用 |

## 9. Others（杂项）
- 不好归到前面章节的小知识点放这里。
- 登录节点上用 `bash -lc` 加 `module load` 加载 R 后可以运行 `Rscript`，适合轻量检查（例如查看某个 rds 的列名）。

## 10. Update History（更新备忘）
每次对脚本或整体代码做了更新，在这里用一句话记录：日期 + 改了什么。
- 2026-10-09：建立 `CLAUDE.md` 和 `Scripts_reference/00_README_Pipeline_Summary.md`；确定只做 3 个细胞类型、S-LDSC 只用 9 个 GWAS；MR 合并步骤改为先删 NA、再过滤 `ivw_sigma2 <= 10`，然后算 CCT 和按 gene 的 BH；只写 `Claude_project/` 内的文件。
- 2026-10-09：开始可视化：按 `20261006_paper` 的思路复现论文图，代码在 `Scripts/Visualization/`，说明在 `Scripts/Visualization.md`，图在 `Results/Figures/`；已完成 Fig1、Fig2（v1 + v2）、Fig3、Fig4、Fig5（第 1 页）、Fig6（Motif_all.tsv + 复制 xlsx）、Fig7、Fig8（GWAS 模式，基于 Valid6_GWAS）。
- 2026-10-09：新增验证模块 `Valid6_GWAS`（`Scripts/Validation/Valid6_GWAS.R`，移植旧 `G1_find_gwas.R`，只用 9 个 GWAS），输出 `Results/Validation/Valid6_GWAS/`；Fig8 将基于它复现。
- 2026-10-12：建立 git 仓库（分支 main，远程 `git@github.com:TFCui/MR_PeakGene.git`），白名单 `.gitignore`；新增第 11 节 Git 规则。

## 11. Git 规则
- 只提交 `Scripts/`、`Results/`、`CLAUDE.md`（外加白名单 `.gitignore` 本身）；`data/`、`tmp/`、`Scripts_reference/` 不提交。
- 提交前先检查 `git status`（以及是否有超过 50MB 的大文件）。
- 不 force push。
