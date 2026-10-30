# CIMA peak–gene 因果连接项目（Claude_project）

## 1. 项目目标
用 CIMA（OMIX009786）细胞类型特异的 caQTL 和 eQTL 做 MR，推断 peak → gene 的因果调控连接，
并和 marginal 相关方法比较，用外部证据验证。
- 细胞类型只做：`Bn_TCL1A`、`CD4_Tn_CCR7`、`CD8_Tn_CCR7`
- **细胞类型名称（画图和论文统一）**：图中标签写 `Bn`、`CD4 T cell`、`CD8 T cell`（代码用 `Scripts/Visualization/00_common.R` 的 `CELL_LAB`）。论文在首次出现处定义：naive B cells（Bn_TCL1A，简称 Bn）、naive CD4⁺ T cells（CD4_Tn_CCR7，简称 CD4 T cells）、naive CD8⁺ T cells（CD8_Tn_CCR7，简称 CD8 T cells）；之后正文写 Bn / CD4 T cells / CD8 T cells，作定语或表头用单数（“CD4 T cell links”）。不再单独用 CD4、CD8、CD4 Tn 或 CIMA 原始标签指代细胞类型。
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
| `Manuscript/` | 论文 LaTeX 源；笔记见 `Notes_Claude_Code.md`（Claude Code）和 `CODEX_METHODS_NOTES.md`（CodeX），规则见 `AGENTS.md`（见第 11 节） |
| `paper_reference/` | 核心参考论文 PDF（MR-Egger：Burgess & Thompson 2017 + 附录、GRIP、SMR：Zhu 2016、CIMA Science）；**只读**，Claude Code 和 CodeX 写论文 / 核对方法时可查阅 |

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
| S-LDSC 资源 | `$BASE/R_package/S_LDSC/reference_genome/`（EUR plink、HM3、baselineLD v2.2、weights、frq；hg19）。来源：Gazal 的 Zenodo “S-LDSC reference files” v4（doi:10.5281/zenodo.10515792，https://zenodo.org/records/10515792），5 个 tgz 的大小与 Zenodo 一致 |
| GWAS | `$BASE/R_package/S_LDSC/GWAS_data/selected_sumstats/`（项目用的 9 个）。107 套独立 GWAS 在 `GWAS_data/sumstats_107/`，下载自同一 Zenodo 记录：https://zenodo.org/records/10515792/files/sumstats_indep107.tgz（md5 已核对一致；LDSC 格式，`SNP A1 A2 N Z`）。`selected_sumstats/` 中的 PASS/GBMI 文件与它完全相同；`CUSTOM.*`（AS、CD、IBD、SLE）和 `PASS_Multiple_sclerosis` 是另外处理的 |
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
| V1+V2 | 每个 gene 内做超几何检验：显著 peak 与 6 个金标准（EnhancerAtlas、FANTOM5、4DGenome、GTEx、HiChIP、PCHiC）的重叠；Fig3 汇总为 global test 表（各 gene 的 N/K/n_sig/k 求和后给 global P 和 pooled fold enrichment，定义见 `Scripts/Visualization.md` 3.3 节） |
| V3 | S-LDSC：MR-Only peak（IVW/Egger/SMR 任一 FDR < 0.1 且 WNN、Pseudo 都不显著的 link；hg19，±500 bp）+ baselineLD，9 个 GWAS，pool（全部候选 peak）作对照；另加 3 个 CT 合并（Combined，MR-Only 和 pool 各取并集）一组；看 enrichment |
| V4 | Motif 富集：Signac + JASPAR2020，背景按 GC 和可及性匹配；提取 MR 特异 motif（IVW/Egger/SMR 任一 BH < 0.05 且 fold > 1，WNN/Pseudo 都不显著），记录 % obs / % bg / fold；文献注释（PMID 已核对）+ logo 汇总为 `Results/Figures/Figures_20261026/Table_2.xlsx`；不再检查 TF 表达 |
| V5 | 逐 gene ROC AUC，用同样的 6 个金标准 |
| V6 | 与 CIMA SCENIC+ GRN 的重叠，做 Fisher 检验 |
| Valid6_GWAS | MR-only 强 link（MR 任一 FDR < 0.001，WNN/Pseudo 不显著）的 peak 内是否有 GWAS P < 5e-8 的 SNP（9 个 GWAS）；Fig8 用 |

## 9. Others（杂项）
- 不好归到前面章节的小知识点放这里。
- 登录节点上用 `bash -lc` 加 `module load` 加载 R 后可以运行 `Rscript`，适合轻量检查（例如查看某个 rds 的列名）。

## 10. Git 规则
- 只提交 `Scripts/`、`Results/`、`Manuscript/`（排除规则见 `.gitignore`）、`CLAUDE.md`（外加白名单 `.gitignore` 本身）；`data/`、`tmp/`、`Scripts_reference/` 不提交。
- 提交前先检查 `git status`（以及是否有超过 50MB 的大文件）。
- 不 force push。
- `paper_reference/` 不在提交白名单内（PDF 共约 14 MB，不提交）。


## 11. Manuscript 论文稿件
- 目录 `Manuscript/`：论文 LaTeX 源（主入口 `main.tex`，与 Overleaf 同步）。本地编译：`bash -l tmp/build_ms.sh`（同时编译 `main.tex` 正文和独立的 `supplementary.tex` 补充文件）。
- **Claude Code 与 CodeX 协同工作**，各有分工和各自的笔记文件（共同规则见 `Manuscript/AGENTS.md`，CodeX 的入口是它）：
  - **Claude Code**：整理 LaTeX、修改格式、检查交叉引用、编译论文；维护 `Manuscript/Notes_Claude_Code.md`。它记录论文里每个 figure / table 的来源、进入论文前的修改、含义（C 节），每个 section 的 LaTeX 现状和 TODO（D 节），Claude Code 的更新记录（E 节）。
  - **CodeX**：阅读 Methods、检查公式、讨论科学逻辑、拟定段落结构；记录在 `Manuscript/CODEX_METHODS_NOTES.md`（按论文章节组织）。
  - 两边互相读对方的文件，各自只写自己的；开始论文工作前先读这两个文件。
- 规则：
  1. 论文里任何图表新增 / 替换 / 修改，必须同步更新 `Notes_Claude_Code.md` 的 C 节台账（来源、改动、含义、hash）。
  2. LaTeX 的某个 section 有改动（包括 CodeX 改的段落），核对实际 `.tex` 后更新 D 节；每次 Claude Code 的工作在 E 节记一行：日期 | 工具 | 内容。
  3. 科学判断、公式核对、段落结构的讨论记在 `CODEX_METHODS_NOTES.md`，不记在 Notes_Claude_Code.md。
  4. 论文工作只改 `Manuscript/`（生成表格的辅助脚本放 `tmp/Manuscript_tables/`）；`Scripts/`、`Results/`、`paper_reference/` 只读。
  5. 不 commit / push，除非用户明确要求；编译中间文件已在 `.gitignore` 中排除。

## 12. Update History（更新备忘）
每次对脚本或整体代码做了更新，在这里用一句话记录：日期 + 改了什么。
- 2026-10-09：建立 `CLAUDE.md` 和 `Scripts_reference/00_README_Pipeline_Summary.md`；确定只做 3 个细胞类型、S-LDSC 只用 9 个 GWAS；MR 合并步骤改为先删 NA、再过滤 `ivw_sigma2 <= 10`，然后算 CCT 和按 gene 的 BH；只写 `Claude_project/` 内的文件。
- 2026-10-09：开始可视化：按 `20261006_paper` 的思路复现论文图，代码在 `Scripts/Visualization/`，说明在 `Scripts/Visualization.md`，图在 `Results/Figures/`；已完成 Fig1、Fig2（v1 + v2）、Fig3、Fig4、Fig5（第 1 页）、Fig6（Motif_all.tsv + 复制 xlsx）、Fig7、Fig8（GWAS 模式，基于 Valid6_GWAS）。
- 2026-10-09：新增验证模块 `Valid6_GWAS`（`Scripts/Validation/Valid6_GWAS.R`，移植旧 `G1_find_gwas.R`，只用 9 个 GWAS），输出 `Results/Validation/Valid6_GWAS/`；Fig8 将基于它复现。
- 2026-10-12：建立 git 仓库（分支 main，远程 `git@github.com:TFCui/MR_PeakGene.git`），白名单 `.gitignore`；新增第 11 节 Git 规则。
- 2026-10-13：可视化更新（Fig1、Fig2、Fig8），新图在 `Results/Figures/Figures_20261019/`，详见 `Scripts/Visualization.md`。
- 2026-10-13：Fig3 改为 global test 表（Global P + pooled fold enrichment）；Hyper 相关文档（`Scripts/Visualization.md` 3.3、`Scripts/Pipeline.md` 3.1）已更新。
- 2026-10-13：Valid2_SLDSC 的显著 peak 阈值改为 FDR < 0.1 并重跑（0.05 版本备份为 `*_FDR0.05.tsv`）；Fig5 去掉 s.e.，raw enrichment P < 0.1 标 `*`。
- 2026-10-14：Valid2_SLDSC 改为只跑 MR-Only + pool（6 方法 FDR 0.1 结果表暂存为 `SLDSC_6methods_FDR0.1.tsv`，FDR 0.05 备份已删）；Fig5 改为 MR-Only vs All peaks。
- 2026-10-14：Valid2_SLDSC 新增 Combined（3 个 CT 合并）一组；Fig5 加第 2 页。
- 2026-10-15：论文稿件目录 `Manuscript/`（CodeX 生成）纳入 git 备份；`.gitignore` 放行该目录，排除编译中间文件（aux/bbl/blg/log/out）、`plos_draft.*` 重复文件、`overleaf_project.zip`、`.plos-review/` 和 `*.before` 备份。
- 2026-10-20：Valid3_Motif 去掉 TF 表达检查；新增 MR 特异 motif 表（47 个）、文献注释和带 logo 的论文 Table 2（`Valid3_Motif_3_table.R` → `Figures_20261026/Table_2.xlsx`），详见 `Scripts/Pipeline.md` 3.3。
- 2026-10-20：Manuscript 记录文档由 6 个合并为 `AGENTS.md` + 笔记文件（含逐图表的来源 / 修改台账）；论文加入 Figure_S2（S2–S4 Fig）和 S3 Table（MR 特异 motif）。
- 2026-10-20：新增 `paper_reference/`（核心参考论文）；`CLAUDE.md` 增加第 11 节 Manuscript，Git 规则改为第 10 节，Update History 移到最后（第 12 节）。
- 2026-10-20：`Manuscript/MANUSCRIPT_CONTEXT.md` 改名为 `Notes_Claude_Code.md`；明确 Claude Code 与 CodeX 协同分工，CodeX 的笔记在 `Manuscript/CODEX_METHODS_NOTES.md`。
- 2026-10-20：论文补充图表拆成独立文件 `Manuscript/supplementary.tex`（`supplementary.pdf`）；正文末尾只保留 Supporting information 标题列表，用 `\nameref` 引用；补充表格统一编号为 S1 / S2 / S3 Table。
- 2026-10-29：第 6 节记录 S-LDSC reference 和 107 套 GWAS 的来源（Zenodo 10515792，已核对文件大小 / md5）；`Scripts/Pipeline.md` 3.2 同步。
- 2026-10-29：细胞类型名称统一为 Bn / CD4 T cell / CD8 T cell（第 1 节，图和论文全文已改）；新一轮论文图放 `Results/Figures/Figures_20261102/`，除 `Figure_1B_C` 只出 SVG 外所有图只出 PDF（见 `Scripts/Visualization.md`）；Fig 1B–C / S3 重画（环形图恢复图例和中心文字）。
