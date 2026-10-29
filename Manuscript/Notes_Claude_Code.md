# MR_Peak_Gene 论文笔记（Claude Code）

最后更新：2026-10-23。本文件是 **Claude Code 的笔记**：论文的概述、期刊要求、图表来源台账（C）、各章节的 LaTeX 现状（D）和 Claude Code 的更新记录（E）。工作规则见同目录 `AGENTS.md`。
**分工**：Claude Code 负责整理 LaTeX、修改格式、检查交叉引用、编译论文，维护本文件；CodeX 负责阅读 Methods、检查公式、讨论科学逻辑、拟定段落结构，记录在 `CODEX_METHODS_NOTES.md`（按论文章节组织）。两边互相读对方的文件，各自只写自己的。
Overleaf 导出不含这些 md。**两个根文档**：`main.tex`（正文，pdfLaTeX + BibTeX）和 `supplementary.tex`（补充图表，独立编译成 `supplementary.pdf`，pdfLaTeX）；本地一次编译两个：`../tmp/build_ms.sh`（`bash -l`）。Overleaf 里一次只能编译一个根文档，要编译补充文件时在 Menu → Main document 切到 `supplementary.tex`。
路径约定：不带前缀的路径相对 `Manuscript/`；`../` 开头的是项目根 `/scratch/user/tfcui2025/CIMA/Claude_project/` 下的路径。

**更新规则**（详见 `AGENTS.md`）：图表有任何变动 → 改 C 节；LaTeX 的某个 section 有改动 → 核对实际 `.tex` 后改 D 节（包括 CodeX 改过的段落）；每次 Claude Code 的工作 → 在 E 节记一行（日期 | 工具 | 内容）。科学判断、公式核对和段落结构的讨论不记在这里，去 `CODEX_METHODS_NOTES.md`。

---

## A. 研究概述

**目的**：检验遗传证据能否帮助把可及染色质 peak 连到靶基因，并把 MR 方法与 marginal 相关方法比较，不预设 MR 一定更好。
**细胞类型**：`Bn_TCL1A`（Bn）、`CD4_Tn_CCR7`（CD4）、`CD8_Tn_CCR7`（CD8）。
**数据**：CIMA（OMIX009786）细胞类型特异 caQTL / eQTL 汇总统计、候选 peak–gene 表、个体层 pseudobulk RNA/ATAC、外部 PBMC multiome 的 WNN metacell。数据路径见 `../CLAUDE.md` 第 6 节。
**6 个方法**：GLS-IVW、GLS-MR Egger、SMR、CCT（SMR + IVW + Egger 的等权 Cauchy 合并），以及 WNN metacell Spearman、个体层 pseudobulk Spearman（marginal）。过滤在 CCT 之前，BH 在每个 gene 内、每个方法单独做；发现阈值 q < 0.1。
**论文使用的验证 / 刻画**：CD4 / Bn / CD8 方法比较和方向一致性（Fig 1、S3）；位点示例（Fig 2）；调控参考的超几何富集（S1 / S2 Table）；GO BP / MF（Fig 3）；S-LDSC 与 GWAS 重叠示例（Fig 4、S4）；null-SNP 诊断（S1 Fig）；MR-Egger flip 敏感性（S2 Fig）；MR 特异 motif（S3 Table）。
**v2 起不在论文里的内容**（仍在 `../Results/` 和 git 历史）：ROC-AUC、旧 BP GSEA、旧的两张主文计数表；motif 以 S3 Table 形式回来，SCENIC+ GRN 重叠以 S5 Fig 形式回来（2026-10-23，仅图和 caption，Results / Methods 文字仍是 TODO）。
**证据地图**：`../Scripts/Pipeline.md`（分析流程）、`../Scripts/Visualization.md`（图）、`../Scripts_reference/00_README_Pipeline_Summary.md`（旧代码）；代码在 `../Scripts/{Model,QC,Validation,Visualization}/`；结果在 `../Results/{MR,Validation,QC,Figures}/`。写 Methods 时要回到代码核对细节。

**核心参考论文（`../paper_reference/`，5 个 PDF，只读）**

| 文件 | 内容 | 用在论文哪里 | bib key |
|---|---|---|---|
| `Burgess and Thompson - 2017 - ...MR-Egger method.pdf` + `MR_Egger ... Appdex.pdf` | MR-Egger 的解释与附录（InSIDE、方向性多效性） | Methods：MR-Egger 回归；Introduction | `burgess2017interpreting` |
| `GRIP_MR_Egger.pdf` | MR-Egger 的扩展（GRIP） | Methods：Egger 的 InSIDE / 多效性讨论 | `dudbridge2025grip` |
| `Zhu et al. - 2016 - ... SMR ...pdf` | SMR | Methods：SMR；Introduction | `zhu2016smr` |
| `CIMA_Science.pdf` | CIMA 队列（数据来源） | Methods：Data sources / Data Preprocessing；Introduction | `yin2026cima` |

---

## B. 期刊与投稿

**目标**：PLOS Genetics，Research Article。官方 URL（未联网核对）：https://journals.plos.org/plosgenetics/s/submission-guidelines 、https://journals.plos.org/plosgenetics/s/latex 。
**模板**：`template_reference/` 中的 PLOS LaTeX 模板 v3.8（2026-04）；`plos2025.bst` 与模板中的文件逐字节一致，不改。`main.tex` 保留模板的宏包、行号（右侧连续）、caption 样式、页脚；双倍行距；单栏；数字 Vancouver 引用（`\citep` 别名到 `\cite`）。
**标题**：全标题 ≤ 200 字符，短标题 ≤ 70 字符（通用模板写 250，本稿按 200，短标题暂定）。
**章节顺序**：Title Page → Abstract（< 300 词，不带引用）→ Author Summary（150–200 词，第一人称）→ Introduction → Results → Discussion → Materials and Methods → Acknowledgments → References → Supporting Information。标题层级最多三级（现在只用 section / subsection）。
**草稿期约定（投稿前要改）**：`\input` 分文件、PDF 图可见（`\draftfigurestrue`）、TODO 可见（`\showtodostrue`）。投稿版要求单一 `.tex`、不嵌入图片：另做一版，内联各 section、`\draftfiguresfalse`（保留 caption 和 label）、删除 TODO 和占位符（只隐藏 TODO 命令不能去掉全部占位符），保持连续行号；supporting 条目只留 caption，文件另交；图表 caption 跟在第一次文中引用之后；正文引用图用 “Fig”，supporting 条目用 PLOS 的名称。出版用图（TIFF/EPS）以后单独准备，目前保留提供的 PDF。
**版面约定**：全文竖版 Letter（8.5×11 in），图表用普通 float；不用横版页，不强制分页。

**STROBE-MR**：需提交已完成的清单作为 supporting information（含稿件页码和对应文字）；目前**没有**清单，也没有空占位。

**投稿系统字段（以下均为 TODO，不能推断或编造；funding 不放 Acknowledgments）**

| 字段 | 仍需要的信息 |
|---|---|
| Data Availability Statement | 作者认可的声明；accession / DOI / 仓库链接；已核实的访问限制和申请流程；与 Methods 的数据来源一致 |
| Author Contributions | 作者认可的 CRediT 角色，与最终作者名单一致 |
| Funding / Financial Disclosure | 核实过的资助方、grant 编号、受助人、资助方作用；不能默认“无资助” |
| Competing Interests | 所有作者认可的声明；不能默认“无” |
| Author metadata | 作者名、单位、通讯信息、需要时的 ORCID；与 Title Page 一致 |

投稿前还要核对：期刊当前指南与上述字段；参考文献的完整作者信息（现有条目有以“三位作者 + others”结尾的）和期刊缩写。

---

## C. 图表台账

每个条目：**位置 / 含义 / 来源 / 进入论文前的修改 / caption 与引用 / 校验**。**FIG** = `../Results/Figures/Figures_20261027/`（用户的最新图和编号，按图逐张更新到论文）。2026-10-23 起 Fig 1–3 已按 `Figures_20261027` 更新；下面其余条目（Fig 4、S 图、表）仍写旧的 `Figures_20261026` 版本，各自更新时再改，在那之前它们对应的 FIG 文件可能已改名或拆分（见 `../Scripts/Visualization.md`）。
旧的 7 个 PDF 图的论文内副本与当时的 FIG 同名文件已核对**逐字节一致**（2026-10-20）；hash 为 SHA-256 前 12 位，刷新用 `sha256sum figures/*.pdf`。
FIG 里每个图另有 `.svg`（文字可编辑）；论文只用 PDF。
编号总览：主文 Fig 1–4（在 `main.tex` / `sections/results.tex`）；Supporting 在独立文件 `supplementary.tex` 中：S1–S5 Fig，S1 / S2 / S3 / S4 Table（2026-10-20 起统一编号，不再有 “S1.2 Tables”；2026-10-23 S7 / S8、S2–S4、S5 / S6 各自合并后为 S1–S4，再新增 S5 Fig（GRN）后为 S1–S5，图编号由 `supplementary.tex` 的计数器自动生成，手写的列表和 `\nameref` 要同步）。
**正文与补充的引用方式**：正文用 PLOS 的 `\nameref{S1_Fig}` 等引用补充条目（label 在 `sections/supporting_info.tex` 的 Supporting information 列表里，只含编号和标题；补充 PDF 里的 caption label 另用 `fig:supp-*` / `tab:supp-*`，两套互不依赖）。补充 caption 里指向正文的 6 处是**固定文字**，正文图号变化时要手动同步：S3 Fig 的 “Fig 1B–C”（1 处）、S4 Fig 的 “Fig 4A”（1 处）。新增或重排补充条目时，要同时改 `supplementary.tex`、`sections/supporting_info.tex`（列表和 label）以及正文里的 `\nameref`。

### Fig 1 — 研究框架、方法重叠与方向比较（CD4）
- **位置**：主文 Fig 1（`sections/results.tex`，label `fig:overview`，第 4 页），`figures/Figure_1.pdf`，`width=\linewidth, height=0.7\textheight, keepaspectratio`（1706×1427 pt，宽 6.5 in 时约 5.4 in 高）。panel 字母 A / B / C 已画在图里，LaTeX 不再另写。
- **含义**：A = 研究框架：MR 用 caQTL（exposure）/ eQTL（outcome）汇总统计；左下 “Correlation across individuals”（CIMA pseudobulk），右下 “Correlation across cell types”（外部 PBMC multiome 的 WNN metacell）。B = CD4 的 UpSet（前 15 个显著组合，6 个方法，p.adj < 0.1，右侧条形图 “Significant pairs”，标题只写细胞类型）。C = 六个方法都显著的 pair 的环形图：扇区是 IVW 方向与两个 marginal 相关方向的关系（都一致 / 只与 Multiomics 不一致 / 只与 Pseudobulk 不一致 / 都不一致），图内只标 73% 和 24% 两个扇区并用短线注明含义，n、SMR / IVW / Egger 内部一致比例、两个小扇区占比都不在图里，写在 caption。图里 “Causal Joint” = CCT。
- **来源**（2026-10-25 更新）：`Figure_1B_C` 已缩小为 9×3.6 in 并简化后（`Plot1_MRresults.R`），用户重新在 Figma 里拼接；原写：用户在 Figma 里把 `Figures_20261027/Figure_1A.pdf`（手工示意图，非脚本生成）和 `Figure_1B_C.pdf`（`../Scripts/Visualization/Plot1_MRresults.R`，CD4）拼接排版，得到 `../Results/Figures/Figures_20261027/Figure_1.pdf`（手工，非脚本生成）；输入 `../Results/MR/<CT>_AllPairs_results.tsv.gz`；`Visualization.md` 3.1。
- **修改**：论文内直接复制该文件，不改画面（2026-10-25 换成新拼接的版本；caption 的 B / C 重写：写明 50,704 个 ≥1 方法显著的 pair、前 15 个组合、四种扇区颜色含义及数量 918 / 304 / 17 / 20（73% / 24% / 1.4% / 1.6%）、n = 1,259 和 SMR/IVW/Egger 一致 1,258/1,259；数量由 `Results/MR` 表直接计数，与图一致）。（2026-10-23）原先 LaTeX 里 1A 和 1B–C 上下叠放并另写 A / B / C 标签，已改为单个 `\includegraphics`。旧的 “示意图 pseudobulk 分支写 across metacells” 的问题在新图中已修正（现在是 across individuals / across cell types）。
- **正文数字**：六方法都显著的 pair 数 Bn 469 / CD4 1,259 / CD8 759；三个 MR 方向一致 469/469、1,258/1,259、758/759；IVW 与两个 marginal 方向都一致的比例约 75% / 73% / 74%（取自图上标注）。
- **引用**：`results.tex` 的 A common analysis set 小节（Fig 1A）；Shared links 小节已于 2026-10-29 被用户删除。caption 里的 Bn / CD8 指向 S3 Fig（原 S5 / S6）。
- **caption（2026-10-29 Overleaf 改写）**：用户重写 A / B / C 三句，去掉了 50,704、前 15 个组合、四个扇区的数量 / 百分比和 SMR/IVW/Egger 一致 1,258/1,259；上面 “修改” 里记录的 2026-10-25 caption 内容已不在论文中。
- **注意**：Bn / CD8 的方法重叠已放进 S3 Fig（`figures/Figure_S3.pdf`）；论文内旧的 `figures/Figure_1A.pdf`、`Figure_1B_C.pdf` 已无引用。
- **hash**：`c23ab1bc59ab`（`figures/Figure_1.pdf` = `Figures_20261027/Figure_1_flat.pdf`；Figma 原文件 `Figure_1.pdf` 的 hash 为 `485ffbcd9a73`）。

### Fig 2 — 四个示例 peak–gene link
- **位置**：主文 Fig 2（`results.tex`，label `fig:examples`），`figures/Figure_2.pdf`，`width=\linewidth, height=0.65\textheight, keepaspectratio`（22×20 in 的图缩放进竖版页）。
- **含义**：4 行，A = SIDT2（Bn），B = RPL11（CD4），C = ERAP2（CD4），D = ARPC1B（CD4）。每行 4 列：caQTL / eQTL −log10 P 轨道（红三角 = MR 保留的工具变量，金三角 = SMR lead SNP）、MR 散点和 IVW / Egger 线、WNN metacell Spearman、个体层 pseudobulk Spearman。
- **来源**：`../Scripts/Visualization/Plot2_examples_v2.R` + `Plot2_examples_slurm.sh`；pair 由 `../Scripts/Visualization/Example_pairs.tsv`（`pair_set = Figure_2`）指定；每个 pair 重跑 QTL 对齐 → 去重 → PLINK LD → GLS IVW / Egger；`Visualization.md` 3.2。
- **修改**：项目内由上一版 6 页 `Fig2_selected` 改为 1 页 4 行，行标签 A–D；论文内不改画面。
- **引用**：**2026-10-29 起正文无引用**（原 Shared links 小节及四个例子的描述被用户删除；float 仍在 A common analysis set 小节）。
- **caption（2026-10-29 Overleaf 改写）**：按列说明 4 列内容（caQTL/eQTL 轨道、MR 拟合 + SMR/IVW/Egger/Causal Joint p 值、80 个 WNN metacell、个体层 pseudobulk），不再写 A–D 对应的基因名。
- **版面**（2026-10-23）：`Figures_20261027/Figure_2.pdf` 与论文内副本逐字节一致，未改；受宽度限制（6.5 in），图内文字在竖版页上偏小。
- **hash**：`f360de84b125`。

### Fig 3 — GO 富集（BP / MF）
- **位置**：主文 Fig 3（`results.tex`，label `fig:go`），`figures/Figure_3.pdf`，`width=\linewidth`。
- **含义**：A = GO Biological Process，B = GO Molecular Function，dotplot；每个 panel 按 DAVID FDR 取前 10 个 term；x = GeneRatio，点大小 = Count，颜色 = p.adj（< 1e-10 截为最深色）。图上标注 A 4,031 个基因、B 4,185 个基因。
- **来源**：`../Scripts/Visualization/Plot9_GO.R`；输入是 3 个 CT 合并的 MR-Only gene（`../Results/Validation/Valid7_GO/DAVID_input/MROnly_genes_all.txt`，4,509 个）在 DAVID 网站的结果 `DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv`（只含 GO BP / MF）；`Visualization.md` 3.10。
- **修改**：由旧 `Fig9_GO_DAVID.pdf` 重画为 13×5.5 in 的 A / B 并排；图中 FDR 一律写 p.adj。论文内不改画面。
- **待确认**：图上 4,031 / 4,185 与输入 4,509 的关系（推测是 DAVID 识别的基因数，未核实）；背景集、DAVID 版本、term 筛选规则要写进 Methods。2026-10-29：用户在 Methods “Gene ontology analysis” 写的是 clusterProfiler `gseGO`（GSEA），与本图的 DAVID 来源不一致，待用户核对；caption 去掉了色标截断（≤1e-10）和 “remain to be documented” 句。
- **版面**（2026-10-23）：`Figures_20261027/Figure_3.pdf` 与论文内副本逐字节一致，未改。
- **hash**：`0d10d08a12c3`。

### Fig 4 — S-LDSC 与 GWAS 重叠的位点示例
- **位置**：主文 Fig 4（`results.tex`，label `fig:disease`，第 7 页），`figures/Figure_4.pdf`（2026-10-23 起只有 1 页）。
- **含义**：A = S-LDSC enrichment，3 个 CT 合并（Combined）的 MR-Only peak vs 全部候选 peak，9 个 GWAS（AS、CD、IBD、RA、SLE、MS、T1D、Asthma、PBC）；MR-Only = IVW / Egger / SMR 任一 FDR < 0.1 且 WNN、Pseudobulk 都不显著的 link 的 peak（hg19，两侧各扩 500 bp，baselineLD）；柱上 `*` 标 enrichment raw P < 0.1，不画 s.e.。B = 4 个位点：RGL2（Bn）、ERBB3、POLR2E、UBAC2（CD4；peak 内 2 个 GWAS SNP，均为 asthma），每个位点顶部写 Gene / Peak（hg38 坐标），底部写 细胞类型 | GWAS 性状，中间是 Gene / SNPs / Peak / 5 个方法的 link / 分段 x 轴；紫色菱形 = GWAS SNP，红 / 蓝 link = 正 / 负效应，虚线 = 不显著，星号 = p.adj < 0.05 / 0.01 / 0.001。
- **来源**：`../Scripts/Visualization/Plot10_Combined.R` + `Plot10_Combined_slurm.sh`；panel A 数据 `../Results/Validation/Valid2_SLDSC/SLDSC_all.tsv`；panel B 的 pair 来自 `../Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv`；`Visualization.md` 3.9，`Pipeline.md` 3.2 / 3.6。
- **修改**：项目内由 Fig5 的 Combined 页和 Fig8 的单 pair 位点图合并而成（原 Fig10，改名 Figure_4）；panel B 行名只在最左列写一次，弧线旁标星号；pair 每页 4 个，共 3 页。**2026-10-23 重画 panel B**：去掉页底说明（星号、红 / 蓝、线宽 ∝ −log10 p.adj、灰虚线）和黑色边框，改成相邻 pair 之间浅灰细竖线；顶部 `Gene: … Peak: …`，底部 `CT | GWAS trait`；左侧行名放大；正文的 GNAS (CD4) 与 S7 的 UBAC2 (CD4) 互换。去掉的解释文字已写进论文 Fig 4 caption（`results.tex`）和 S7 / S8 caption。项目内输出拆成 `Figure_4.pdf`（第 1 页，正文）和 `Figure_S4.pdf`（2 页 = S4.1 / S4.2，补充）；论文内 Figure_4.pdf 放主文，Figure_S4.pdf 的第 1 / 2 页放 S7 / S8 Fig（补充图编号本步未改）。
- **2026-10-23 重画 panel A**：由 9 个独立 y 轴的小面板改为单个面板，9 个 GWAS 共用同一个 y 轴，GWAS 名写在 x 轴，每个 GWAS 的 MR-Only / All peaks 两根柱紧贴（数据、星号标准不变，仍是 Combined 3 个 CT）；`results.tex` Fig 4 caption 的 A 句加入 “for each GWAS trait, shown on a common y axis”。panel B 不变。
- **caption（2026-10-29 Overleaf 改写）**：A 句写明星号 = p < 0.1；B 只保留轨道 / 颜色 / 线宽 / 星号说明，去掉疾病缩写表、四个位点列表、位点标题说明和 “//” 轴断说明。正文只剩 S-LDSC 小节的两句（引用 Fig 4A、4B），原 GWAS 小节（`res:gwas`）已删除。
- **待确认**：MR-only 选择规则现写在 Methods “Heritability enrichment analysis from S-LDSC”（IVW/Egger/SMR 任一 p.adj < 0.1 且两个 marginal 都不显著）；项目文档的说法如上，写定前要用代码核对。
- **hash**：`58fa3582037e`（`Figure_4.pdf`，2026-10-23 重画 panel A，与 FIG 逐字节一致）。

### S1 Fig — null-SNP 诊断
- **位置**：Supporting（`supplementary.tex`，label `fig:supp-qc`；正文用 `\nameref{S1_Fig}`），`figures/Figure_S1.pdf`，`width=\linewidth`（9.5×4.8 in）。
- **含义**：2 个 panel：A = caQTL 与 eQTL z 的 pair 内相关 ρ̂（按 pair 内共享 null SNP 数分组），B = 相关性检验的 p 值；每组 3 根并排箱子 = Bn / CD4 / CD8（紫 / 绿 / 橙）；虚线 = ρ = 0 和 p = 0.5。各 CT 的全局 ρ / P 写在图例里；pair 数（547,452 / 1,051,953 / 686,829）和 SNP–pair 观测数（17,074,938 / 33,263,695 / 21,537,822）写在 caption。
- **来源**：`../Scripts/Visualization/Plot_Sup1_OneSampleMR.R`；输入 `../Results/QC/1_OneSampleMR/<CT>_pair_corr.tsv.gz` 和 `<CT>_global_rho.tsv`；上游 `../Scripts/QC/QC1_1_OneSampleMR_prune.R`、`QC1_2_OneSampleMR.R`；说明 `Pipeline.md` 4.1，`Visualization.md` 3.11。
- **修改**：项目内由旧的每 CT 3 页 PDF 合并为 1 页（原 Sup1_OneSampleMR）；**2026-10-23 重画**：原 3 行 × 2 列（panel A–F，每 CT 一行）合并为 2 个 panel，三个细胞类型在同一个 panel 内并排，换 Dark2 配色，去掉每个箱子的 n 和行标题，说明移到 caption（旧版备份在 `../tmp/Plot_S1_test/old_20261027/`）；数据和分组不变（各组 pair 数、中位数与旧图一致）。方法：hg38 面板 LD pruning（500 kb / 50 / r² 0.1，MAF 0.01），两边 nominal P > 0.05，共享 SNP，slope 差 ≤ 1e-3 去重，ρ̂ 仅对 SNP 数 > 2 的 pair 计算。论文内不改画面；caption 改为 A / B 两个 panel。
- **待补**：2026-10-29 用户删除了 Methods “Null-SNP diagnostic” 小节，S1 Fig 目前正文无引用、Methods 无对应文字。
- **hash**：`92e1b359ae5a`（与 FIG 逐字节一致）。

### S2 Fig — MR-Egger 的 flip 前后对比（Bn / CD4 / CD8，原 S2 / S3 / S4 Fig 合并）
- **位置**：Supporting，label `fig:supp-flip`（原 `fig:supp-flip-bn/-cd4/-cd8` 已删除），`figures/Figure_S2.pdf`（单页，13×9.6 in），`width=\linewidth`。
- **含义**：3 行 = Bn / CD4 / CD8，每行 4 个 panel：左半 = flip 前后差别小的 pair 的 [原始编码 | flip]，右半 = 差别大的 pair 的 [原始编码 | flip]；flip = 让 caQTL 效应 > 0（与正式 MR 一致）；同一 pair 两图坐标范围相同；红 / 蓝线 = IVW / Egger，panel 上方写 IVW P、Egger P；金三角 = SMR lead SNP，灰线 = 95% CI。图内不再写 SNP 数、β、intercept 和说明文字（说明在 caption）。
  - 差别小：Bn SIDT2（194 SNP，Egger p 2.9e-171 → 7.0e-171）、CD4 RPL11（25，1.7e-11 → 2.5e-7）、CD8 LCLAT1（195，6.9e-142 → 1.3e-140）。
  - 差别大（2026-10-23 换成 SNP 数多的例子，原 GSTM4 10 / RGS12 32 / YBEY 40 SNP 弃用）：Bn **PIGX**（chr3:195747820-195748321，136 SNP，Egger p 1.3e-11 → 0.175，gene 内 FDR 6.6e-11 → 0.44）、CD4 **PDE4DIP**（chr1:149608426-149608927，74 SNP，1.5e-10 → 0.675，FDR 6e-10 → 0.90）、CD8 **KLRC3**（chr12:10553007-10553508，134 SNP，8.1e-9 → 0.267，FDR 2.3e-8 → 0.53）。IVW 对 flip 不变。
- **来源**：`../Scripts/QC/QC2_2_FlipPoints_plot.R` + `_slurm.sh`；上游全基因组扫描 `QC2_1_FlipScan.R`（176,545 个 pair，不 flip 的 Egger 在 gene 内 BH，flip 版与正式表一致）；候选挑选 `QC2_1_FlipScan_select.R`（不 flip p.adj < 0.05 且 flip 后 ≥ 0.1：Bn 700 / CD4 1,763 / CD8 1,056；其中 SNP 数 ≥ 50 的 187 / 507 / 305 个另存 `FlipScan_candidates_nsnp50.tsv`，大差别例子从这里挑）；数值 `../Results/QC/2_FlipPoints/FigureS2_estimates.tsv`；说明 `Pipeline.md` 4.2。
- **修改**：项目内由 3 页（每 CT 一页，上下两个 pair、左右 flip 前后）重排为单页 3 × 4，只出 `Figure_S2.pdf` + `Figure_S2.svg`（旧 3 页版和旧估计表备份在 `../tmp/QC2_FlipPoints_test/old_20261027/`）。论文内一个 float，caption 写明行 / 列含义、SNP 数、线和符号。2026-10-23 补充图重新编号：S2–S4 合并为 S2 后，原 S5 / S6 / S7 顺延为 S3 / S4 / S5。
- **引用**：目前只有 `results.tex` 末尾 TODO（`\nameref{S2_Fig}`）；Methods 的 “Allele-orientation sensitivity of MR-Egger” 小节 2026-10-29 已被用户删除，正文还没写。
- **panel 字母**（2026-10-23）：每个 peak–gene pair 一个字母，按行从左到右：A SIDT2、B PIGX（Bn），C RPL11、D PDE4DIP（CD4），E LCLAT1、F KLRC3（CD8）；字母画在每个 pair 标题行的左侧。caption 同步写 A–F。
- **hash**：`723913309d81`（`Figure_S2.pdf`，与 FIG 逐字节一致）。

### S3 Fig — Bn / CD8 的方法重叠与方向比较（原 S5 / S6 Fig 合并）
- **位置**：Supporting，label `fig:supp-bn-cd8`（原 `fig:supp-bn` / `fig:supp-cd8` 已删除；正文用 `\nameref{S3_Fig}`），`figures/Figure_S3.pdf`（单页，15×11 in），`width=\linewidth`。
- **含义**：A = Bn（上），B = CD8（下）（2026-10-23 加 panel 字母，每个细胞类型一个，画在每行左上角；`Plot1_MRresults.R` 默认不再重写 `Figure_1B_C`）；每个 CT 一行：左 UpSet（前 25 个显著组合，6 个方法，p.adj < 0.1，右侧条形图 “Significant pairs”），右环形图（六方法都显著的 pair：IVW 方向与两个 marginal 相关方向的关系，中心 = n 和 SMR / IVW / Egger 内部方向一致比例）。图内没有 A / B 字母，上下由 caption 说明。数字：六方法都显著 Bn 469 / CD8 759；一致比例 75% / 74%。
- **来源**：`../Scripts/Visualization/Plot1_MRresults.R`（同 Fig 1B–C，一个脚本同时出 `Figure_1B_C`（CD4）和 `Figure_S3`（Bn + CD8）；`Visualization.md` 3.1）；`../Results/Figures/Figures_20261027/Figure_S3.pdf`。
- **修改**：2026-10-23 图本身未改，直接复制为 `figures/Figure_S3.pdf`；论文内由原来的两个 float（用旧 `Figure_1B_C.pdf` 第 1 / 3 页）合并为一个 float。S3 / S4 Fig 合并后，原 S7 → S5（再顺延）：最终补充图编号为 S1–S4。
- **引用**：Results 的 A common analysis set、Shared links 小节和 Fig 1 caption（三处 `\nameref{S3_Fig} and \nameref{S4_Fig}` 改为 `\nameref{S3_Fig}`）。
- **注意**：论文内 `figures/Figure_1A.pdf`、`figures/Figure_1B_C.pdf`（旧版）已没有任何引用（Fig 1 用 `Figure_1.pdf`，S3 用 `Figure_S3.pdf`）；暂留未删。
- **hash**：`e93ab3b469f0`（`Figure_S3.pdf`，加字母后重生成，与 FIG 逐字节一致）。

### S4 Fig — 更多 GWAS 重叠位点（原 S7 / S8 Fig 合并；因 S2–S4、S5–S6 合并而顺延为 S4）
- **位置**：Supporting，label `fig:supp-disease-2`（原 `fig:supp-disease-3` 已删除），`figures/Figure_S4.pdf`，单页，`width=\linewidth`（24×18.7 in，宽 6.5 in 时约 5.1 in 高）。2026-10-23 起；之前是 S7 / S8 两个 float，用 `Figure_4.pdf` 第 2 / 3 页。
- **含义**：panel A 与 Fig 4A 相同；S7 的 panel B = GLS2（CD8，T1D / asthma）、GNAS（CD4，IBD / Crohn's；2026-10-23 与正文 UBAC2 互换）、UBA7（Bn，IBD）、TMBIM1（CD8，Crohn's）；第二行 = GNAS（CD8，IBD / Crohn's）、ERAP1（Bn，Crohn's）、TMPRSS3（CD4，T1D）、ZSCAN16（Bn，RA）。
- **来源 / 修改**：同 Fig 4（版式同 2026-10-23 的新 panel B）；项目内 `Figure_S4.{pdf,svg}` 为单页，panel A 同 Fig 4A，panel B 两行 × 4 个 pair（共 8 个），页高 18.7 in；caption 写明两行的 pair 和符号说明。`sections/supporting_info.tex` 删除 S8 Fig 条目，`results.tex` 的 `\nameref` 只剩 S7；`supplementary.tex` 的 “S1–S8 Figs” 改为 “S1–S7 Figs”。
- **2026-10-23**：panel A 随 Fig 4A 重画（共用 y 轴，GWAS 名在 x 轴，两根柱紧贴），panel B 不变；caption 不变（“repeated from Fig 4A”）。
- **hash**：`7464e354d049`（`Figure_S4.pdf`，与 FIG 逐字节一致）；`Figure_4.pdf` 为 `58fa3582037e`。

### S5 Fig — 与 CIMA SCENIC+ 调控网络（GRN）的重叠（2026-10-23 新增）
- **位置**：Supporting，label `fig:supp-grn`（正文用 `\nameref{S5_Fig}`，目前只在 `results.tex` 末尾 TODO 里；`methods.tex` 的 GRN 小节 2026-10-29 已被用户删除），`figures/Figure_S5.pdf`，`width=\linewidth`（10×9.5 in）。
- **含义**：A = 各方法显著（gene 内 BH p.adj < 0.05，比主图的 0.1 严）且在 GRN 中的 pair 数（行 = 细胞类型 eRegulon / 高质量 eRegulon，列 = Bn / CD4 / CD8，柱 = 6 个方法）；B = 显著 pair 在 GRN 中的 odds ratio（2×2 Fisher，背景 = 该 CT 全部 MR pair），CT 与 CT_HQ 并排，虚线 OR = 1。所有比较 Fisher P < 1e-30。与 MR 重叠的 GRN pair 数：CT 864 / 2,520 / 1,364，CT_HQ 817 / 458 / 382（eRegulon 数 47 / 40 / 33 和 29 / 2 / 8）。OR 范围：CT 集 2.3–9.6（Bn）、2.8–7.2（CD4）、3.4–7.5（CD8），CT_HQ 集更高（CD4 / CD8 的 SMR 达 10.8 / 11.6）；WNN 的 OR 最低，SMR 最高。
- **来源**：`../Results/Validation/Valid5_GRN/`（`../Scripts/Validation/Valid5_GRN.R`，2026-10-09 跑，MR 表此后未变）→ `../Scripts/Visualization/Plot_Sup2_GRN.R`（旧版 `Plot7_GRN.R` 的 `Fig7_GRN.pdf` 为 14×11 in 的 facet 图，已重画）；`Visualization.md` 3.7，`Pipeline.md` 3.5。
- **修改**：新图，不是旧图复用；p.adj 写法、方法配色与全文一致；GRN 集合用浅 / 深灰蓝（避免与 6 个方法的配色冲突）；两个 panel 的三列对齐。
- **待补**：Results / Methods 文字（`results.tex` 末尾 TODO；Methods 目前没有 GRN 小节，写时要对照 `Valid5_GRN.R` 并补 CIMA GRN 的来源）；eRegulon 集合、背景、阈值要写进 Methods。
- **hash**：`09621f1760ed`（`Figure_S5.pdf`，与 FIG 逐字节一致）。

### S1 Table — 超几何富集 p 值
- **位置**：Supporting，`tables/supplementary_hyper_p.tex`，label `tab:supp-hyper-p`。2026-10-29 起正文无引用（Results 的相关句和 Methods “Regulatory-reference enrichment” 被用户删除；新的 Methods “Hypergeometric test” 未引用）。
- **含义**：18 行（3 CT × 6 方法）× 6 个参考（EnhancerAtlas、FANTOM5、4DGenome、GTEx、HiChIP、PCHiC）的 pooled 超几何 p 值。每个 gene 内：N = 候选 peak 数，K = 与参考重叠的候选 peak 数，n_sig = 该方法显著 peak 数（FDR < 0.1），k = 其中与参考重叠的数；只用 N > 0、K > 0、n_sig > 0、N > K 的 gene；Global P = `phyper(Σk−1, ΣK, ΣN−ΣK, Σn_sig, lower.tail = FALSE)`，不额外多重校正。
- **来源**：`../Results/Validation/Valid1_Hyper/<CT>_hyper_long.tsv`（`../Scripts/Validation/Valid1_Hyper.R`）→ `../Scripts/Visualization/Plot3_Hyper.R` → `../Results/Figures/Figures_20261027/Table_1_2.xlsx` + `Table_1_2_global.tsv`；说明 `Visualization.md` 3.3，`Pipeline.md` 3.1。**不是**把各 gene 的 P 值合并：Valid1 对每个 gene 做超几何检验，Table_1_2 把各 gene 的 N / K / n_sig / k 求和后做一次检验（Global P）并算 pooled fold，共 3 CT × 6 方法 × 6 参考 = 108 个格子。
- **修改**：项目内由旧的 gene 级 p 值直方图改为 global test 表；论文内由 xlsx 的数值**逐格转写成 LaTeX**（`booktabs`，p 值字符串保留低于浮点下限的指数，如 1.0×10^-416）；工作簿不随论文提交。2026-10-23：`Table_1_2.xlsx` 由每个 CT 一个 sheet 合并为**单个 sheet**（上块 Global P = S1 Table，下块 pooled fold = S2 Table，各 18 行，CT 单元格纵向合并，每 CT 内每列最高 fold 高亮）；S1 / S2 Table 的 216 个数值与 `Table_1_2_global.tsv` 逐格核对**一致，数值没有变化**；唯一改动是方法列的名称与图统一（SMR / GLS-IVW / GLS-MR Egger / CCT / WNN / Pseudobulk → Causal SMR / Causal IVW / Causal Egger / Causal Joint / Marginal Multiomics / Marginal Pseudobulk，脚注说明 Causal Joint = CCT 等），并把 `\tabcolsep` 由 3 pt 调为 2.5 pt 以免长方法名溢出。旧 xlsx 备份在 `../tmp/Tables_test/old_20261027/`。
- **hash**：`Table_1_2.xlsx` = `721d5095aeb8`。

### S2 Table（原 “S1.2 Tables”）— 超几何富集 fold enrichment
- **位置**：Supporting，`tables/supplementary_hyper_fold.tex`，label `tab:supp-hyper-fold`；2026-10-29 起正文无引用。2026-10-20 去掉了 `\thetable = S1.2` 的特殊写法，改成统一的 S2 Table。
- **含义**：同 S1 Table 的 18 × 6 布局，值为 pooled fold = (Σk/Σn_sig)/(ΣK/ΣN)（显著 peak 落在参考中的比例 / 候选 peak 落在参考中的比例），保留两位小数；每个细胞类型内每个参考的最大值加粗。注意它没有置信区间，且 fold 要和显著 peak 数一起看（显著 peak 越多 P 越小但 fold 常更低）。
- **来源 / 修改**：同 S1 Table（`Table_1_2.xlsx` 下半块逐格转写；方法名已与图统一）。
- **数字（2026-10-29 起不在正文中）**：SMR 在 18 个 CT / 参考组合中的 17 个 fold 最大，CD4 PCHiC 是 IVW 最大；Bn 中 SMR 与 WNN 的 fold：FANTOM5 6.18 vs 3.56，HiChIP 4.46 vs 3.08。
- **hash**：来源同上。

### S3 Table — MR 特异 TF motif
- **位置**：Supporting 末尾，`tables/supplementary_motif.tex`，label `tab:supp-motif`，`longtable`（竖版，`\scriptsize`，跨页重复表头，按细胞类型分组的整行标题）。目前正文只有 TODO 引用。
- **含义**：13 个 MR 特异 motif（Bn 10 / CD4 2 / CD8 1）。2026-10-23 起“MR 特异”= IVW / MR-Egger / SMR 任一 Signac p.adj < **0.01** 且 fold > 1，同时 WNN 和 Pseudobulk 都不显著（p.adj ≥ 0.05），并且有明确的文献解释（CCT 不参与筛选，只列出）。Valid3 结果里的完整定义（p.adj < 0.05）仍得到 47 个（Bn 23 / CD4 14 / CD8 10）：其中 4 个没有明确解释（Bn CREB3L4、Bn Pax2、CD4 ZNF460、CD8 HINFP），其余 30 个 MR p.adj ≥ 0.01，共 34 个被去掉（见 `Table_3.xlsx` 的 `Excluded` sheet）。阈值对比（去掉无解释行后）：0.05 → 43、0.02 → 22、0.01 → 13、0.005 → 5、0.001 → 1。Motif 富集：Signac + JASPAR2020 CORE vertebrates，前景 = 各方法 FDR < 0.1 的 peak，背景从全部候选 peak 中按 GC 和可及性（log1p count）匹配 5 万个。
- **来源**：`../Results/Validation/Valid3_Motif/Motif_MR_specific.tsv`（统计量，`Valid3_Motif_2_summary.R`）+ `Motif_MR_specific_literature.tsv`（文献注释；`../tmp/Valid3_Motif_lit/` 起草、`curate.py` 去掉限定语，86 个 PMID 均用 NCBI eutils 核对过题目、年份和摘要）+ `logos/<ID>.png`（`Valid3_Motif_3_table.R`，PNG 里标题带 TF 名和 Motif ID）；说明 `Pipeline.md` 3.3。完整版（含 % obs / % bg / fold / FDR 高亮 / immune relevance / 完整解释 / Excluded）在 `../Results/Figures/Figures_20261027/Table_3.xlsx`（hash `aeab6d63fc36`）。
- **进入论文前的修改**：由 `../tmp/Manuscript_tables/make_motif_table.py` 从两个 TSV 生成 `.tex`（筛选规则与 `Table_3.xlsx` 相同，两边的 motif 集合已核对一致；PMID 取自 `reference_note`，不手抄）。相对 Table_3.xlsx：**删除** TF motif 名、Motif ID、Immune relevance 三列和 % obs / % bg / fold，只留各方法的 p.adj；MR 一侧 p.adj < 0.01 且 fold > 1 的数值加粗，PMID 只作文本（不进 `references.bib`）。**版式（2026-10-23）**：去掉 Cell type 列，改为按 CT 分组的整行标题（`Bn (10 motifs)` 等）；logo 宽由 0.78 in 放大到 1.74 in（约 2.2 倍，高约 0.93 in），p.adj 六列各 0.4 in，interpretation 列左对齐、宽 1.4 in；表只保留用到的 13 个 logo（`figures/motif_logos/`，其余副本已删，可从 `../Results/Validation/Valid3_Motif/logos/` 重新复制）。表下说明写明筛选定义和 **Reference (PMID) 的来源**：按 TF 在 PubMed 中检索同一或相关免疫细胞类型、免疫介导疾病中的功能文献，PMID 逐个核对 PubMed 记录（题目、年份、摘要）——这段只写了文档里有记录的流程，具体检索词没有保存，**待用户确认措辞**。
- **注意**：`Table_3.xlsx` 同样是 13 行（`Motif_table`、`Details`）+ `Excluded`（34 行）三个 sheet；logo 放大到 1.8 × 0.96 in。剩下的 13 个里 Creb5、CENPB、POU3F1 的 immune relevance 是 Weak（解释句仍明确）。页脚写的是“Signac 报告的 adjusted p-value”；p.adj 的校正方式（项目文档写 BH）还没用 Signac 源码核对。
- **待补**：Results 的 motif 小节还没写。2026-10-29 用户新写的 Methods “Motif enrichment analysis” 用的是 p.adj < 0.05 且 fold > 1，而本表用 p.adj < 0.01（并按文献解释筛选），两者待统一。

### S4 Table — 三个细胞类型的数据概览
- **位置**：Supporting 末尾（S3 Table 之后），`tables/supplementary_data_overview.tex`，label `tab:supp-data`，竖版 `table`（booktabs，4 列：项目 / Bn / CD4 / CD8，按 B–E 四块分组）；2026-10-29 起正文无引用（用户改写 Methods “Data sources” 后去掉了 `\nameref{S4_Table}`）。
- **含义**：给出分析用到的 Bn / CD4 / CD8 的基本数字：个体数（335 / 372 / 346，同时通过 RNA 和 ATAC 过滤）、scRNA / scATAC 细胞数（及每个体中位数）、pseudobulk 矩阵的基因 / peak 数、eQTL 样本量（405 / 406 / 405）和 caQTL 样本量（385 / 385 / 384）、“至少一个 caQTL 的 peak 数”（P < 1e-8）、候选 pair 漏斗（All_GenePeak → 共享变异 ≥ 3 → 两个边际数据均有 → 去 NA → ivw_sigma2 ≤ 10；最终 40,087 / 85,527 / 50,931 pair，基因 8,149 / 10,701 / 9,427，peak 2,951 / 5,338 / 3,327）、每个 pair 的 IVW 变异数（中位数 22）、WNN 的 80 个 metacell。表注说明三套样本数的区别（QTL 样本量 / pseudobulk 个体数 / WNN metacell）、caQTL 与 eQTL 样本重叠、CIMA 全队列数字（428 人、73 个 L4 类型、±1 Mb、5,051,397 个 MAF > 0.1 变异，引自 CIMA 论文）。
- **来源**：`../Scripts/Visualization/Table4_DataOverview.R`（slurm：`Table4_DataOverview_slurm.sh`；输出 `../Results/Figures/Figures_20261027/Table_4.xlsx`（单 sheet）+ `Table_4.tsv`（长表））；漏斗数字解析自 `../Scripts/Model/logs/MR_AllPairs_20074030_{0,1,2}.out`，最后一行与 `../Results/MR/<CT>_AllPairs_results.tsv.gz` 的 pair / 基因 / peak 数逐项核对（脚本内 `stopifnot`）。
- **进入论文前的修改**：`../tmp/Manuscript_tables/make_data_table.py` 从 `Table_4.tsv` 生成 `.tex`（不手抄）。
- **注意**：(1) “peaks with ≥ 1 caQTL at P < 1e-8” 是名义 P 阈值（即选工具变量的阈值），不是 CIMA 论文 FDR 定义的 caPeak，表注已写明；(2) 没有 eGene 数（本地只有名义 P 的汇总统计），如需要要从 CIMA 补充表 S8 取数并标注出处；(3) 表注里的 CIMA 全队列数字来自 `paper_reference/` 的 CIMA 论文，请作者最终核对；(4) Methods “Data sources” 的 TODO 保留（外部 PBMC 引用、协变量、数据获取声明仍缺）。
- **hash**：`Table_4.xlsx` 见 E 节 / `sha256sum`。

---

## D. 各章节写作记录

（这里记录 LaTeX 里**现在有什么**、引用了哪些图表、还有哪些 TODO；CodeX 对 Methods / 公式的科学审阅和段落设计见 `CODEX_METHODS_NOTES.md`。若 CodeX 改过某一节，以实际 `.tex` 为准更新本节。）

状态词：Not started / Outline / Initial draft / Under revision。TODO 数是 `\TODO{}` 的个数（`sections/*.tex` 中 2026-10-29 统计，Overleaf 导入后共 2 个）。
2026-10-29 起，各节文字以用户在 Overleaf 上的版本为准（`tmp/Manuscript_overleaf.zip` 原样导入，未改任何文字）；导入前的本地版本备份在 `../tmp/Manuscript_backup_20261029/sections/`。

### Title Page（`sections/title_page.tex`）— Initial draft
已写：全标题 “Mendelian randomization and multi-omic evidence for regulatory peak–gene links in human immune cells”，作者 Tengfei Cui¹、Yang Ni²\*、Yuchao Jiang¹\*，两个单位（Texas A&M 统计系；UT Austin 统计与数据科学系），通讯邮箱 yang.ni@austin.utexas.edu、yuchaojiang@tamu.edu。2026-10-29：用户删去了标题页上的 “Short title” 行（`main.tex` 里 `\ShortTitle` 的定义仍在，目前未使用）。缺：作者顺序和单位需作者确认；ORCID（如需要）。

### Abstract（`abstract.tex`）— Not started
只有一句占位（“will be written after the Results narrative…”），无 TODO。≤ 300 词。

### Author Summary（`author_summary.tex`）— Not started
空（只有标题）。150–200 词，非技术性，第一人称。

### Introduction（`introduction.tex`）— Not started
空（2026-10-29 用户删去了五段提纲和 TODO，只剩标题和 label）。

### Results（`results.tex`）— Under revision（1 个 TODO）
- **A common analysis set enables comparison of six linking methods**：只有 6 个方法定义一段，引用 Fig 1A。Fig 1、Fig 2 两个 float 放在这一小节。（原先的 q 定义、发现数比较段落和 “Shared links can have different MR and marginal effect directions” 小节已被用户删除。）
- **Fig 1 caption**（用户改写）：A = 框架（MR 用 caQTL/eQTL；Marginal multiomics 用外部 PBMC metacell；Marginal pseudobulk 用 CIMA 个体；均为 Spearman），B = CD4 UpSet（BH p.adj < 0.1），C = 1,259 个六方法都显著 link 的方向扇区；Bn / CD8 见 S3 Fig。不再写 50,704、扇区数量和 SMR/IVW/Egger 一致比例。
- **Fig 2 caption**（用户改写）：按列说明（第 1 列 caQTL/eQTL 轨道，第 2 列 MR 拟合与 4 个 p 值，第 3 列 80 个 WNN metacell，第 4 列个体层 pseudobulk）；不再写 A–D 的基因名。**正文目前没有引用 Fig 2**（`fig:examples` 无 `\ref`）。
- **Selected links are enriched for regulatory-reference support**：只剩一句（列出 6 个参考），不再引用 S1 / S2 Table，也不再有 fold 数字。
- **GO summaries…**：一句，引用 Fig 3。Fig 3 caption 去掉了色标截断和 “remain to be documented” 句。
- **MR-only regions show exploratory disease-heritability enrichment**：一句引用 Fig 4A，一句列出 Fig 4B 的 RGL2（Bn）、ERBB3、POLR2E、UBAC2（CD4）。原 “GWAS-overlapping peaks…” 小节（`res:gwas`）已删除，S4 Fig 正文不再引用。
- **Fig 4 caption**（用户改写）：A 星号 = p < 0.1；B 只写轨道 / 颜色 / 线宽 / 星号说明，去掉了疾病缩写表、位点列表、标题说明和 “//” 轴断说明。
- **末尾 TODO**：补 S2 Fig（flip 对比）、S5 Fig（GRN）和 S3 Table（MR 特异 motif）的结果文字。

### Discussion（`discussion.tex`）— Outline
三条空提纲标题：Principal findings / Complementary biological information / Limitations and future directions（2026-10-29 用户清空内容，删去 “Scope of inference” 条和 TODO）。

### Materials and Methods（`methods.tex` + `data_preprocessing.tex`）— Under revision（1 个 TODO）
2026-10-29 Overleaf 版本的 15 个小节：
1. **Data sources**：用户新写：CIMA（OMIX009786）提供的 caQTL / eQTL 汇总统计与单细胞计数；三个 naive 细胞类型（Bn_TCL1A、CD4_Tn_CCR7、CD8_Tn_CCR7）及简称；“External multiome data for WNN” 只有标题句。占位 `(Citation)`。不再引用 S4 Table。
2. **Data Preprocessing**（`data_preprocessing.tex`）：用户重写为提纲（6 条空 bullet）+ 五段草稿：候选 pair 与过滤（caQTL p < 1e-8、去重复 IV（差 < 1e-3）、≥ 3 个 IV、两边数据可用；pair 数 3 个 `XXX`，段末自注要移到最后）；LD clumping（clump_r2 = 0.99，clump_kb = 200000，hg38 EAS）；LD 矩阵（PLINK `--keep-allele-order`、`--r`，按 exposure 方向翻转）；Pseudobulk（≥ 100 细胞、median + 3 scaled MAD，335 / 372 / 346 个体，CPM + log1p）；WNN metacell（Signac 重新计数 CIMA peak，RNA PC 1–50 + ATAC LSI 2–50，80 个 metacell；含 “Add how to get … later” 自注）。原先的 LD 参数冲突 TODO 已删除（用户写的是 0.99；项目文档为 20 Mb = 20000 kb，正文写的 200000 待核对）。
3. **MR Framework**：用户重写：三条 IV 假设（enumerate）、水平多效性与 IVW / MR-Egger（InSIDE）、个体层模型（带截距 μ）、汇总统计、Γ_j = βγ_j + α_j（带编号公式）。
4. **Inverse-variance weighted method**、5. **MR-Egger regression**、6. **GLS-IVW with correlated variants**、7. **GLS MR-Egger with correlated variants**：公式未改；第 7 节末尾新增过滤规则段（先删任一 p 为 NA，再删 IVW 残差方差 > 10；gene 内按方法 BH；q < 0.1），段首有用户自注 “( We can explain the filteration rule here… )”。
8. **Summary-data-based Mendelian randomization**：lead SNP、p 值反推 z、SE 与 χ² 统计量、无 HEIDI；无 TODO。
9. **Marginal pseudobulk correlations**（label `sec:pseudobulk-correlations`）：CPM + log(1+CPM)，Spearman，`cor.test(exact=FALSE)`，无协变量，gene 内 BH，FDR ≤ 0.1。
10. **WNN metacell correlations**（label `sec:wnn-correlations`）：外部 PBMC 80 个 metacell，Signac `FeatureMatrix` 重新计数，log(1+x)，Spearman，gene 内 BH，FDR ≤ 0.1。**唯一 TODO**：外部数据描述与引用、WNN 方法引用、metacell 构建细节。
11. **Causal joint**：CCT 等权合并 IVW / Egger / SMR，给出 T 统计量公式（未引用 `liu2020cauchy`）。
12. **Motif enrichment analysis**：Signac AddMotifs + JASPAR2020（746 motif），FindMotifs，motif 间 BH，p.adj < 0.05 且 fold > 1；MR 特异定义；`(table number)` 占位。注意：S3 Table 现用 p.adj < 0.01 的定义。
13. **Hypergeometric test**：6 个参考资源的描述、liftOver、`phyper(q-1, m, n, k, lower.tail=F)`；多个 `(citation)`。不再写 pooled fold 定义，也不引用 S1 / S2 Table。
14. **Gene ontology analysis**：写的是 clusterProfiler `gseGO`（GSEA，pos 打分，10–500 基因）；注意 Fig 3 来自 DAVID（`Plot9_GO.R`，见 C 节 Fig 3），两者待用户核对。
15. **Heritability enrichment analysis from S-LDSC**：MR-only 定义、背景 peak、baselineLD v2.2、1000G EUR、HapMap3、9 个 GWAS 与下载来源；`(Version later)`、`(Later)`；T1D 写成 “T1B”。
原先的 Candidate peak–gene pairs、Instrument selection and LD preparation、Filtering and multiple testing、Regulatory-reference enrichment、Gene Ontology summaries、Disease heritability and GWAS overlap、Null-SNP diagnostic、Allele-orientation sensitivity of MR-Egger、Overlap with the CIMA regulatory network 小节已被用户删除（S1 / S2 / S5 Fig 在 Methods 中不再有对应小节）。`methods.tex` 末尾约 30 行空行。

### Acknowledgments（`acknowledgments.tex`）— Not started
空（只有标题；占位 TODO 已删除）。

### References（`references.bib` + `plos2025.bst`）— Initial draft
9 个条目：`benjamini1995fdr`、`bowden2015egger`、`zhu2016smr`、`liu2020cauchy`、`finucane2015sldsc`、`burgess2017interpreting`、`dudbridge2025grip`、`yin2026cima`、`burgess2016combining`。待补：软件 / 资源 / GWAS 引用；部分条目作者列表以 “others” 结尾，需补全作者信息和期刊缩写。

### Supporting Information（`supplementary.tex` + `tables/`，正文列表 `sections/supporting_info.tex`）— Initial draft
现有：独立补充文件 `supplementary.tex`（首页有 Supporting figures / tables 目录，单倍行距、无行号，页眉 “Supporting Information”）：S1–S5 Fig、S1 / S2 / S3 / S4 Table，caption 已按当前图表更新（见 C 节）。正文末尾只有 “Supporting information” 列表（编号 + 标题，8 条）；正文用 `\nameref` 引用。待做：STROBE-MR 清单；投稿时把条目拆成单独文件上传的安排；GO 设置和 enrichment 阈值确认后更新相关 caption；S3 Table 的 Methods / Results 文字。

---

## E. 更新记录

格式：日期 | 工具 | 内容。

- 2026-10-15 | CodeX | 生成初稿 `Manuscript/`（PLOS 模板适配、全部章节、旧图 Fig1–Fig9 + Sup1、旧主文计数表、AUC / GRN / motif / BP GSEA 的结果段落）；纳入 git 备份。该稿可从 git 提交 `7cf0b1a` 恢复。
- 2026-10-19 | Overleaf（手改）→ 2026-10-20 导入 | `Manuscript_MR_v2`（Overleaf 导出，2026-10-19）整体替换 CodeX 稿：新图 Figure_1A / 1B_C / 2 / 3 / 4 / S1；删除旧 AUC、GRN、motif、BP GSEA 展示和两张主文计数表；Results、caption 和 Methods 下游描述按新图重写；S1 / S1.2 Table 由 `Table_1.xlsx` 逐格转写；新增 Data Preprocessing 小节（从 Word 稿导入，LD 参数冲突标 TODO，不编造内容）；未重跑任何分析。
- 2026-10-20 | Claude Code | 加入 Figure_S2 → S2 / S3 / S4 Fig（原 S2–S5 顺延为 S5–S8）和 S3 Table（native LaTeX，47 个 MR 特异 motif，logo 复制到 `figures/motif_logos/`，`main.tex` 加 `longtable` 包和 logo 的 graphicspath）；Methods 加两个仅含 TODO 的小节，Results 末尾加 TODO；编译通过（32 页，无未定义引用；`supplementary.tex` 第 77–111 行附近有 1 处 14 pt 的 overfull，未处理）。
- 2026-10-20 | Claude Code | 记录文档由 6 个合并为 `AGENTS.md` + 本文件（删除 `README.md`、`UPDATE_NOTES.md`、`submission_requirements.md`、`figure_manifest.tsv`，内容并入 B / C / D 节）；新增 `../paper_reference/`（核心参考论文）；`../CLAUDE.md` 增加 Manuscript 一节并把 Update History 移到最后。
- 2026-10-20 | CodeX | 只读梳理项目与 SMR 写作计划（未改 LaTeX、未重跑分析）：核对 `MR_AllPairs.R` 的 lead SNP 选择、p-derived z、ratio/SE/chi-square 和无 HEIDI 实现，并阅读本地 Zhu 2016 Online Methods Eq. 4–5；后续需明确忽略 exposure/outcome 估计协方差的近似，不能把 CIMA 两类 QTL 默认写成独立样本。直接核对结果表：Bn/CD4/CD8 为 40,087/85,527/50,931 pairs，pseudobulk n=335/372/346，WNN n=80（不等同 QTL 样本量）。另发现现有 GLS 文字的截断集合尚未体现代码中优先保留前 5 个合格特征方向的规则，留待对应小节修订。
- 2026-10-20 | Claude Code | `MANUSCRIPT_CONTEXT.md` 改名为 `Notes_Claude_Code.md`；明确 Claude Code（LaTeX 整理 / 格式 / 交叉引用 / 编译）与 CodeX（Methods 阅读 / 公式检查 / 科学逻辑 / 段落结构，记录在 `CODEX_METHODS_NOTES.md`）的分工；同步修改 `AGENTS.md`、`../CLAUDE.md`、`CODEX_METHODS_NOTES.md` 和 `acknowledgments.tex` 中的文件名。
- 2026-10-20 | Claude Code | 补充图表从正文拆出：新增独立文档 `supplementary.tex`（自带序言，单倍行距、无行号，目录列表；`sections/supplementary.tex` 删除），`main.tex` 去掉补充相关的序言设置（`plossupport*` caption 格式、`longtable`、motif logo 的 graphicspath）并在参考文献后 `\input{sections/supporting_info}`（PLOS 格式的 11 条编号 + 标题列表）；正文 11 处引用改成 `\nameref{S?_Fig/Table}`；表格编号统一为 S1 / S2 / S3 Table（去掉 S1.2）；补充 caption 指向正文的 6 处改成固定文字（Fig 1B–C、Fig 4A）；`tmp/build_ms.sh` 改成同时编译两个文档。编译通过：正文 23 页、补充 13 页，均无未定义引用，原先的 14 pt overfull 随之消失。
- 2026-10-23 | Claude Code | 主文 Fig 1–3 按 `Figures_20261027` 更新：新增 `figures/Figure_1.pdf`（用户在 Figma 拼接的 1A + 1B–C），`results.tex` 的 Fig 1 float 改为单个 `\includegraphics`（去掉 LaTeX 写的 A / B / C 标签）；Fig 2、Fig 3 与 `Figures_20261027` 逐字节一致，未改；旧 `Figure_1B_C.pdf` 保留给 S5 / S6 Fig，待补充图步骤更新。C 节 Fig 1A 与 1B–C 合并为 Fig 1，D 节同步。编译通过（正文 23 页、补充文件无未定义引用，无 overfull）。
- 2026-10-23 | Claude Code | Fig 4 / Figure_S4 重画 panel B 并更新到论文：`Plot10_Combined.R` 去掉页底说明和黑框，改浅灰竖线分隔，顶部 Gene / Peak、底部 CT | trait，行名放大，UBAC2 (CD4) 与 GNAS (CD4) 互换；新增 `figures/Figure_S4.pdf`，`supplementary.tex` 的 S7 / S8 改用它的第 1 / 2 页；Fig 4、S7、S8 caption 和 Results 的 GWAS 示例段（GNAS → UBAC2）同步；编译通过（正文 23 页、补充 13 页，无未定义引用、无 overfull）。
- 2026-10-23 | Claude Code | Figure_S4 改为单页（panel A 不变，panel B 两行 × 4 = 8 个 pair，24×18.7 in）：`Plot10_Combined.R` 的 `PAGES` 改成 2 页（正文 4 个 / 补充 8 个），页高按行数计算，只输出 `Figure_S4.{pdf,svg}`（删除 `Figure_S4_1/2.svg`）；论文里 S7 / S8 Fig 合并为 S7 Fig（`supplementary.tex` 单个 float、`supporting_info.tex` 去掉 S8、`results.tex` 的引用只留 S7、列表计数 S1–S7）；编译通过（正文 23 页、补充 13 页，无未定义引用、无 overfull）。
- 2026-10-23 | Claude Code | Figure_S1 重画并更新到论文：`Plot_Sup1_OneSampleMR.R` 把 6 个 panel 合并为 2 个（A = ρ̂，B = P 值；三个 CT 并排箱，Dark2 配色，图例带各 CT 全局 ρ / P，去掉箱上 n 和行标题），9.5×4.8 in；`supplementary.tex` 的 S1 Fig 改 `width=\linewidth`，caption 改为 A / B 并写入各 CT 的 pair 数和 SNP–pair 观测数；编译通过（正文 23 页、补充 13 页，无未定义引用、无 overfull）。
- 2026-10-23 | Claude Code | Figure_S2 重画并更新到论文：`QC2_2_FlipPoints_plot.R` 把 3 页合并为单页 3 行 × 4 panel（每行 = 一个 CT：差别小的 pair [原始 | flip] + 差别大的 pair [原始 | flip]），去掉 SNP 数 / β / intercept / 说明文字，panel 上方写 IVW P、Egger P；“差别大”的例子从全基因组 flip 扫描（`<CT>_Egger_FlipScan.tsv.gz`）中换成 SNP 数多的 PIGX (136) / PDE4DIP (74) / KLRC3 (134)；`QC2_1_FlipScan_select.R` 新增 `FlipScan_candidates_nsnp50.tsv`。论文：`supplementary.tex` 的 3 个 flip float 合并为 S2 Fig，S2–S4 合并后补充图重新编号为 S1–S5（原 S5 / S6 / S7 → S3 / S4 / S5），同步 `supporting_info.tex`、`results.tex`、`methods.tex` 的 `\nameref` 和列表；编译通过（正文 23 页、补充 10 页，无未定义引用、无 overfull）。
- 2026-10-23 | Claude Code | Figure_S3 更新到论文：`Figures_20261027/Figure_S3.pdf`（Bn 上、CD8 下，图本身未改）复制为 `figures/Figure_S3.pdf`；`supplementary.tex` 中 Bn、CD8 两个 float（旧 `Figure_1B_C.pdf` 第 1 / 3 页）合并为一个 S3 Fig；补充图编号变为 S1–S4（原 S5 → S4），同步 `supporting_info.tex`、`results.tex` 的 `\nameref`（4 处）和列表计数；编译通过（正文 23 页、补充 11 页，无未定义引用、无 overfull）。
- 2026-10-23 | Claude Code | Figure_S2 / S3 重画加 panel 字母（S2：每个 peak–gene pair 一个，A–F；S3：每个细胞类型一个，A = Bn、B = CD8）并更新到论文，caption 同步；新增 S5 Fig（Valid5_GRN：`Plot_Sup2_GRN.R` 画 2 个 panel，A = 显著且在 GRN 中的 pair 数，B = OR），`supplementary.tex` 加 S5 Fig float，`supporting_info.tex` 加 S5 条目，`results.tex` 末尾 TODO 和 `methods.tex` 新小节（仅 TODO）引用 `\nameref{S5_Fig}`；`Plot1_MRresults.R` 默认不再重写已用于 Figma 的 `Figure_1B_C`；编译通过（正文 23 页、补充 12 页，无未定义引用、无 overfull）。
- 2026-10-23 | Claude Code | 表格：`Table_1_2.xlsx` 合并为单 sheet（`Plot3_Hyper.R`），S1 / S2 Table 方法名与图统一（216 个数值逐格核对一致）；MR 特异 motif 表收紧到 p.adj < 0.01、去掉无明确解释的行（47 → 13，`Valid3_Motif_3_table.R` 在建表阶段筛选，`Table_3.xlsx` 加 `Excluded` sheet），S3 Table 重排（按 CT 分组、logo 放大、加 Reference 来源说明；`make_motif_table.py` 重写）；编译通过（正文 23 页、补充 10 页，无未定义引用，仅有 <1 pt overfull）。
- 2026-10-23 | Claude Code | Fig 4A / Figure_S4A（S-LDSC）重画并更新到论文：`Plot10_Combined.R` 把 9 个独立 y 轴的小面板改为单个面板（共用 y 轴，GWAS 名在 x 轴，MR-Only / All peaks 两根柱紧贴），panel B 不变；复制 `figures/Figure_4.pdf`、`figures/Figure_S4.pdf`，Fig 4 caption 的 A 句加入 “common y axis”；C 节 Fig 4 / S4 的修改和 hash 同步。
- 2026-10-23 | Claude Code | 新增 S4 Table（三个细胞类型的数据概览）：`Table4_DataOverview.R`（+ slurm）算出个体 / 细胞 / QTL 样本量 / caQTL peak / 漏斗 / 每 pair 变异数，输出 `Figures_20261027/Table_4.xlsx` + `Table_4.tsv`（漏斗末行与 `Results/MR` 表核对一致）；`tmp/Manuscript_tables/make_data_table.py` 生成 `tables/supplementary_data_overview.tex`，`supplementary.tex`（S3 Table 后）、`supporting_info.tex`（S4 Table 条目）、Methods “Data sources”（`\nameref{S4_Table}`）同步；编译无未定义引用。Table_4.xlsx hash a8253961edc3。
- 2026-10-25 | Claude Code | Fig 1 更新到论文：复制用户在 Figma 重新拼好的 `Figure_1.pdf`（B 缩小为前 15 个组合、C 只标两个大扇区）到 `figures/Figure_1.pdf`；`results.tex` 的 Fig 1 caption 的 B / C 重写（见 C 节 Fig 1 条目）；编译通过，第 4 页检查排版正常。S3 Fig 尚未换成新缩小版。
- 2026-10-25 | Claude Code | Fig 1 透明度拍平：Figma 导出的 `Figure_1.pdf` 含约 158 张带 SMask 的位图，部分 PDF 阅读器显示为黑底；用 Ghostscript（PDF 1.3 兼容、整页 720 ppi 无损）重写为 `Figures_20261027/Figure_1_flat.pdf` 并放入 `figures/Figure_1.pdf`（无 SMask、无透明度；代价：图内文字变为位图，不再可选 / 可搜）。Figma 原文件保留不动。
- 2026-10-29 | Overleaf（手改）→ Claude Code 导入 | 导入 `tmp/Manuscript_overleaf.zip`（用户在 Overleaf 上只改了文字）：42 个文件中 33 个（`main.tex`、`supplementary.tex`、`references.bib`、全部 `figures/`、`tables/`、`supporting_info.tex`）与本地逐字节一致，未改；只有 `sections/` 下 9 个文件（abstract、acknowledgments、author_summary、data_preprocessing、discussion、introduction、methods、results、title_page）不同，按用户要求逐字节覆盖（文字、拼写、占位一律未改）。10-25 旧导出包与导入前本地 sections 一致，无需合并；导入前版本备份在 `../tmp/Manuscript_backup_20261029/sections/`。TODO 21 → 2。C 节 Fig 1–4、S1/S2/S5 Fig、S1–S4 Table 的引用 / caption 说明和 D 节同步。编译通过：正文 21 页（原 23 页）、补充 11 页，无未定义引用；正文无 overfull，补充仍是原有 13 处 <1 pt 表格 overfull。
