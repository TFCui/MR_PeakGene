# 2_BOARD：论文修改看板

记录要做什么、已经定了什么。流程见 `1_WORKFLOW.md`，数字见 `3_LEDGER.md` 第 1 节。路径相对 `Manuscript/`。

## 0. 怎么用

- 每次会话读第 1、2 节，再读当前小节（第 4 节）。
- 工作卡只出现在正在做的小节下面。每个问题都有 `建议：` 和 `作者：` 两行；作者在 `作者：` 后面写一句就行（写“同意”也可以），结构上的想法写在卡末的 `补充：`。
- 收卡后删掉工作卡。该小节的“已定”里只留以后还会约束写作的决定，并更新 2.2 的状态和 2.3。
- 第 1 节定下的内容不再讨论；有新的决定就追加进去。全文层面的问题放到第 3 节，由作者以后统一处理。

## 1. 统一表达

### 1.1 术语与写法

| 项             | 写法                                                                                                                                                                                                                                                                                                                        |
| -------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 细胞类型       | 首次出现写 naive B cells (Bn cells)、naive CD4⁺ T cells (CD4 T cells)、naive CD8⁺ T cells (CD8 T cells)，之后写 Bn cells / CD4 T cells / CD8 T cells；作定语时用单数（"CD4 T cell links"）。不单独写 "Bn"、"CD4"；CIMA 原始标签只在 Methods 的 Data sources 小节出现                                                      |
| 方法名（正文） | causal IVW、causal Egger（不写 causal MR-Egger）、causal SMR、causal joint、marginal multiomics、marginal pseudobulk。Methods 小节标题保留 GLS-IVW 等全称；图表标签不改（Causal IVW … Marginal Pseudobulk）                                                                                                                |
| 类别名         | Methods 和 Results 的开头各定义一次："We refer to the four MR-based tests (…) collectively as the causal methods, and to the two correlation-based methods as the marginal methods." 不写 "our causal methods"（这些方法不是本文提出的）；本文的贡献落在框架上（"an MR framework for peak–gene linking"），框架不起短名字 |
| 显著性         | $P_{\mathrm{adj}}$ = 每个 gene 内部、每个方法分别做 BH 后的 P 值，用宏 `\q`；全文写 `$\q<0.1$`；名义 P 写斜体 *P*。不写 p.adj、FDR、q；DAVID 自己的 FDR 写 "DAVID FDR"                                                                                                                                              |
| peak / enhancer | 生物学概念和小节标题写 "enhancer–gene"（与论文标题、Abstract、Intro、Discussion 一致，2026-11-20 作者决定）；被检验的单位写 "peak–gene pair"，显著的写 "link"；首次出现写 "chromatin accessibility peak" |
| pair / link    | 被检验的单位写 "peak–gene pair"，显著的 pair 写 "link"；marginal 方法估计的量写 "peak–gene association"；eRegulon 沿用 "region–gene pair"。不写 linkage（易与 genetic linkage / LD 混淆）、edge、causal link（2026-11-22）                                                                                                                                                                                                                                                                     |
| GWAS           | 只报告 7 个疾病：AS、CD、IBD、RA、SLE、MS、T1D，统称 "immune-mediated diseases"（不写 "blood traits"）；S-LDSC 实际跑了 9 个，asthma 和 PBC 不报告                                                                                                                                                                          |
| GO             | 只写 DAVID 过表示分析                                                                                                                                                                                                                                                                                                       |
| 时态           | 默认一般现在时；我们做过的操作、得到的结果、前人研究用过去时（细则见 1_WORKFLOW 4.4；2026-11-22 全文核对） |
| 格式           | 基因符号斜体；科学计数写 `$6.1\times10^{-4}$`；美式拼写（center）；图写 `Fig~\ref{fig:overview}B`，表写 `Table~\ref{tab:motif}`，补充条目写 `\nameref{S1_Fig}`                                                                                                                                                      |
| 参考文献       | PLOS 数字编号（`\cite{}`）；正文不写作者名（不写 "X et al." 和年份），引用跟在论断 / 工具名 / 研究名后（如 `Signac \cite{…}`、`the TRIPOD study \cite{…}`）；只有 S4 Table 来源列保留 "Author et al.~\cite{…}"；缺的地方写 `[CITE: 论断]` |

### 1.2 定义与阈值（主定义在 Methods 的 Study design and terminology；Results 开头用一两句话再说一次）

| 名称                | 定义                                                                                                                                                        | 用在                                          |
| ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| MR-only link        | causal IVW / Egger / SMR 中至少一个 $\q<0.1$，且两个 marginal 方法都 $\ge 0.1$；causal joint 不参与                                                     | GO、S-LDSC、S5 Fig / Fig 3D；3 个细胞类型合并 |
| strong MR-only link | 4 个 causal 方法（含 joint）中任一 $\q<0.001$，且两个 marginal 方法都 $\ge 0.05$                                                                        | GWAS 重叠（Fig 4B）                           |
| MR-specific motif   | 前景是每个方法各自 $\q<0.1$ 的 peak；motif 层面要求 $P_{\mathrm{adj}}<0.05$ 且 fold > 1；IVW / Egger / SMR 中至少一个富集，且两个 marginal 方法都不富集 | Table 1                                       |
| GRN overlap         | 每个方法各自的显著 pair，阈值 $\q<0.05$（比别处严格）                                                                                                     | Fig 3E、S8 Fig                                |

### 1.3 说法边界

| 主题               | 可以写                                                                                                                                                                            | 不能写                                                                            |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------- |
| causal             | "MR estimates causal effects under the instrumental-variable assumptions, whereas correlation-based methods estimate associations"；"MR-supported link"；"putative causal effect" | identify / confirm causal links；true causal effect；prove；marginal 方法“错了” |
| 多个方法一致       | "consistent evidence across methods with different assumptions"                                                                                                                   | independent evidence（工具变量和数据都相同）                                      |
| 新颖性             | "Unlike existing peak–gene linking methods, which rely on correlations …, we estimate …"                                                                                       | "the first"（文献核对之前不写）；novel                                            |
| 比较结果           | MR 的 OR / fold 更高；causal methods 找到的 link 更少                                                                                                                             | MR 优于 marginal 方法；outperform；validate / confirm                             |
| 外部参考和 GRN     | "overlap"、"were enriched for"；fold 是描述性的                                                                                                                                   | 金标准；验证了因果关系；fold 的差异显著（没有做方法之间的检验）                   |
| S-LDSC             | 点估计更高；nominal evidence（P < 0.1）                                                                                                                                           | 显著富集；MR-only 比全部 peak “更富集”（没有检验两组之间的差异）                |
| GWAS 重叠          | 提示候选靶基因                                                                                                                                                                    | 变异通过这个 peak 调控该基因；该基因是致病基因（没有做共定位）                    |
| motif              | 只在 MR 方法中富集；如实写出 6 个都由 Egger 选出、其中 4 个只由 Egger 选出                                                                                                        | TF 结合这些 peak；marginal 方法“看不到”这些机制                                 |
| GO                 | 富集项（背景是全基因组，fold 1.16–1.41）                                                                                                                                         | MR-only 基因有“MR 特有的功能”                                                   |
| Egger 截距、多效性 | —                                                                                                                                                                                | "no pleiotropy"（结果表里没有 Egger 截距）                                        |
| 样本重叠诊断       | 相关的数值很小                                                                                                                                                                    | 没有相关；证明了独立；证明了 MR 的校准                                            |

### 1.4 内容取舍（已定）

- 不放进论文的分析：ROC-AUC、旧的 BP GSEA、旧的两张主文计数表。
- GRN 只报 OR（Fig 3E；Bn 和 CD8 T cells 放在 S8 Fig）。
- Fig 4A 和全文只报告 7 个 GWAS。
- S2 Fig 中“flip 前后差别大”的例子选 SNP 数多的 pair（PIGX、PDE4DIP、KLRC3）。
- Table 1 放在正文：每个细胞类型 2 个 motif，每个 motif 引一篇文献；注释不加 "exploratory" 这类限定语。
- S1 / S2 Table 中的方法名与图一致。版面：竖版 Letter，普通 float。

## 2. 全文地图（架构由作者决定，Claude 负责记录和维护）

### 2.1 主线: 

作者的一句话论点：（待作者写）

### 2.2 各节一览

状态：未动 / 已调结构 / 已润色（日期）/ 有新改动。

**Abstract**（`sections/abstract.tex`，1 段 10 句，292 词）：背景 → 现有方法靠相关（correlation does not imply causation）→ MR 框架和四个检验 → CIMA 三个细胞类型、两种对照方法 → link 数量和方向 → eRegulon OR、6 个外部参考 → MR-only peak 的 HiChIP 3.28 倍和 7 个疾病的遗传力点估计 → GWAS 位点例子（*ERAP1*、*ERBB3*、*THBS3*）→ 结论句。状态：已重写（11-18），待作者 polish。

**Introduction**（`sections/introduction.tex`，7 段，约 1,230 词）：P1 背景（调控元件、GWAS 非编码变异、enhancer–gene map 难以实验获得）→ P2 现有 peak–gene linking 方法（跨细胞 / metacell、跨个体的相关）→ P3 marginal association 的问题：跨细胞类型（Xie 2023 的例子、TRIPOD）、跨个体（scMultiMap、本文关于 pseudobulk 的推论、GRaNIE Box 1），段末 “correlation does not imply causation” → P4 MR、分子性状、CIMA 提供同一细胞类型的 caQTL 和 eQTL → P5 三个难点（工具变量少且相关、多效性、单一方法的假设）→ P6 本文框架（causal IVW / Egger / SMR / joint）→ P7 应用与结果（定性，不放数字）。状态：已改写（11-18），P3 新增（11-18），待作者 polish。

**Author Summary**：空。

**Results**（`sections/results.tex`）

| 小节 | 标题 | 段 | 每段一句话 | 状态 |
| ---- | ---- | -- | ---------- | ---- |
| R1 | A Mendelian randomization framework for causal enhancer–gene inference | P1–P7 | IV 理由 → causal IVW / Egger（模型式、GLS、一致性、InSIDE、φ）→ NOME 与 sample overlap（S1 Fig）→ SMR / joint → CIMA 数据 → 两种 marginal 方法的思路 → 术语 → Fig 1A | 已扩写（11-22） |
| R2 | Concordance and discordance between causal effects and marginal associations | P1 | 候选 pair 数；CD4 各方法 link 数；multiomics 多（回扣 S1 Sec 1，可能有细胞类型混杂）；MR-only 定义与 4,749；870 组合；S4 Fig 共享比例 → Fig 1B、S3 Fig、S4 Fig | 已润色（11-22），待作者看 |
| | | P2 | 1,259 个 link 的方向：MR 之间一致；不一致集中在 metacell 相关（25.7% vs 2.9%，Bn / CD8 比例）→ Fig 1C、S3 Fig | |
| | | P3 | S6 Fig：全部 pair 上 IVW 与 pseudobulk r 0.24–0.30、与 multiomics 0.03–0.04；MR 之间 0.55–0.87 | 11-22 新增 |
| | | P4 | *SIDT2* / *RPL11* / *ERAP2* / *ARPC1B* 四个例子（只描述）→ Fig 2 | |
| | | P5 | 机制：跨细胞类型 / donor 的共同变异；混杂掩盖（*RPL11*、*ARPC1B*，S1 Sec 1）；marginal 相关不能确立因果 | 11-22 从 P4 拆出 |
| R3 | Orthogonal regulatory evidence for MR-inferred enhancer–gene links | P1 | Fig 3A 描述（信号强度、头部基因差异，不写含义）→ S7 Fig | 已润色（11-22），按 panel 顺序重排 |
| | | P2 | GO（重写）：转录调控、RNA 加工（*SRSF7* / *DDX5*）、泛素化（*UBE2N*）+ 全基因组背景限定 → Fig 3B | |
| | | P3 | 6 个外部参考 × 各方法（S1 Table 都富集；S2 Table SMR fold 17/18 最大、multiomics 18/18 最小；P 值受 link 数影响） | |
| | | P4 | MR-only 与参考：HiChIP 8 组都在对角线上方，3.28 倍；其余 5 个参考 1.89–2.87 → Fig 3C、S5 Fig | |
| | | P5 | CIMA eRegulon OR（$\q<0.05$，Fisher）→ Fig 3D、S8 Fig | |
| | | P6–P7 | Motif：47 个 MR-specific，Table 1 全部由 Egger 选出（等位方向句已删）；MNT、ATF3、FOSL2::JUND；R3 收尾句（Claude 草拟，待作者确认） | |
| R4 | Heritability enrichment and candidate target genes for immune-mediated diseases | P1 | S-LDSC：全部 peak 本身富集；MR-only 点估计更高，4 个 P < 0.1；覆盖 0.25% vs 15.1%；收尾句（Claude 草拟）→ Fig 4A | 已润色（11-22） |
| | | P2–P4 | 选例依据（24 个 strong MR-only + GWAS link）；*ERAP1* / *ERBB3* / *THBS3* 各一段，统一为 位点 → 已有证据 → MR β、$\q$、marginal 不显著 → 含义 → Fig 4B | |

**Discussion**（`sections/discussion.tex`）：3 段：D1 研究概述 → D2 验证结果总结 → 局限 First：样本重叠与 NOME 合并论述（同一个协方差、偏倚 ∝ 1/F、强工具变量为主、S1 Fig 为辅；末尾 LIML 两句和 independent cohorts 一句）。原 “allele orientation” 段 11-20 移到 Methods 的 Allele orientation 小节。状态：3 段已润色（11-22，待作者看）：D1 加核心发现（MR 与 metacell 相关基本无关、与 donor 相关较一致）；D2 按 Results 改正 GO 主题、去掉 “strong biological support” / “therapeutic intervention”，补 S-LDSC；局限段只改两句。作者决定暂不补写新段落。

**Methods**（`sections/methods.tex` + `data_preprocessing.tex`）

| 小节                                 | 段 | 内容                                                                  | 状态               |
| ------------------------------------ | -- | --------------------------------------------------------------------- | ------------------ |
| M1 Data sources                      | 2  | CIMA 数据；外部 PBMC multiome；三个细胞类型及标签                     | 已润色（11-02）    |
| M0b Study design and terminology     | 2  | 框架；类别名、$P_{\mathrm{adj}}$、pair / link、MR-only 的定义       | 已润色（11-02）    |
| M2 Data Preprocessing                | 4  | 候选 pair → 工具变量 → 等位与去重 → LD 矩阵                        | 已润色（11-02）    |
| M3 Mendelian randomization           | 5  | 三条假设 → 为什么能避开混杂（S1 Text）→ 记号、NOME 与 overlap 偏倚 → MR-Egger 式 → IVW 偏倚条件、InSIDE | 已补充（11-22）    |
| M4a Causal IVW                       | 4  | LD → GLS；记号、Σ 只含 Γ̂ 的误差；截断（参数、Eckart–Young）；推断、φ̂ 只放大 SE、φ̂ > 10 | 已补充（11-22）    |
| M4b Causal Egger                     | 3  | 加截距的 GLS → 等位对齐、相关变异的 InSIDE、对测量误差更敏感 → 推断 | 已补充（11-22）    |
| M5 Causal SMR                        | 2  | 基本 SMR（比值、T_SMR、不依赖 NOME、假设两样本独立、不能区分多效性）→ 本研究的算法（lead SNP、P → z、χ²₁）、与 IVW 的关系 | 已重写（11-22）    |
| M6 Causal joint                      | 1  | CCT 的公式；相关 P 值下仍有效；由最小 P 值主导 | 已补充（11-22）    |
| M7 Marginal pseudobulk               | 3  | count 矩阵 → 个体过滤 → ln(1 + CPM) → Spearman → 边际关联、无协变量、含混杂 | 已补充（11-22）    |
| M8 Marginal multiomics               | 3  | 外部数据（TRIPOD：7,790 细胞、14 类型、80 个 WNN metacell）；重新计数；Spearman | 已补充（11-22）    |
| M9 Filtering and multiple testing    | 1  | 删 NA → $\hat\phi \le 10$ → joint → gene 内 BH                   | 已调结构           |
| M10 Hypergeometric enrichment        | 7  | 6 个参考；计数定义；pooled 检验；三种分析（各方法 × CT → S1/S2 Table；MR-only 3 CT 合并 → Fig 3C / S5；按 K/N 分 8 组 + bootstrap CI） | 已重写（11-09），待润色 |
| M11 Gene ontology analysis           | 1  | DAVID 输入、类别、作图规则                                            | 已调结构           |
| M12 GRN overlap                      | 1  | 两套 eRegulon；2×2 Fisher 检验                                       | 已调结构           |
| M13 S-LDSC                           | 2  | 注释的构造、参数；7 个疾病和来源                                      | 已调结构           |
| M14 GWAS overlap                     | 1  | strong MR-only link + GWAS P < 5e-8                                   | 已调结构           |
| M15 Motif enrichment                 | 1  | 前景、背景、motif 集、检验和阈值                                      | 已调结构           |
| M16 Allele orientation               | 1  | flip 操作；IVW 不受影响、Egger 受影响                                 | 已调结构           |
| M17 Error correlation due to sample overlap | 3 | 动机（同一 cohort、偏倚 ∝ 1/F̄）→ null SNP 的 z 相关估计 ρ_j → pooled 与逐 pair 两种估计 | 已重写（11-21）    |

**Supplementary：S1 Text**（`sections/s1_text.tex`，承接 M3–M4、M17，2026-11-02 重组，2026-11-21 改写）

| 小节                                                     | 承接正文       | 状态            |
| -------------------------------------------------------- | -------------- | --------------- |
| S1.1 Marginal association and causal effect（Fig A）       | M3、M7、M8、R1 | 已写（11-20），符号改为 G / ξ（11-21） |
| S1.2 MR model（样本集合、残差相关）                        | M3             | 已改写（11-21） |
| S1.3 Inverse-variance weighted method                      | M3 / M4a       | 已润色（11-02） |
| S1.4 Consistency of the IVW estimator                      | M3             | 已润色（11-02） |
| S1.5 Sample overlap and measurement error（含 NOME）       | M3、M17、R1-P3、Discussion | 新增（11-21） |
| S1.6 MR-Egger regression and the InSIDE assumption         | M3 / M4b / M16 | 已润色（11-02），恢复 hats（11-21） |
| S1.7 Correlated variants（含 LD 交叉项）                   | M4             | 已润色（11-02），补一句（11-21） |
| S1.8 Multiplicative random effects                         | M4a–b / M9    | 已润色（11-02） |
| S1.9 Stabilizing GLS estimation via eigenvalue truncation  | M4a            | 已润色（11-02） |

### 2.3 信息在哪里

**定义和关键数字**（首次给出的位置 → 其他出现的位置）：

- 类别名、$P_{\mathrm{adj}}$、MR-only：Methods Terminology → R2-P1（MR-only 的短版定义和 4,749）。strong MR-only：M14 → R4-P2。MR-specific motif：M15 → R3-P6。GRN 阈值 0.05：M12 → R3-P5、Fig 3D caption。
- 候选 pair 数 40,087 / 85,527 / 50,931：R2-P1（Methods 不再写）、S3 Table。
- pseudobulk 335 / 372 / 346 人：M7、S3 Table。80 个 metacell：R1-P6、M8。
- caQTL $P<10^{-8}$、clumping 参数：M2；$P<10^{-8}$ 也在 R1-P1 / P3、M17、S1.5（χ² > 32.8、1/F̄ ≲ 3%）、Discussion。$\hat\phi \le 10$：M4a、M9、S1.8。
- sample overlap：不写人数（R1、M17、Discussion、S3 Table 脚注）；全局 ρ̂ 0.0035 / 0.0060 / 0.0053 只在 S1 Fig 图注。

**图表引用**（2026-11-22 更新；`\TODO` 里的提及不算引用）

| 条目             | 被引用的位置                                 |
| ---------------- | -------------------------------------------- |
| Fig 1A / 1B / 1C | R1-P1、P5 / R2-P1、R3-P5 / R2-P2、R2-P5 |
| Fig 2A–D        | R2-P4、R2-P5（D） |
| Fig 3A / 3B / 3C / 3D | R3-P1 / R3-P2 / R3-P4 / R3-P5；M11、M12 |
| Fig 4A / 4B      | R4-P1 / R4-P2；M14 |
| Table 1          | R3-P6、R3-P7；M15 |
| S1 Fig           | R1；Discussion 局限段；M17 |
| S2 Fig           | M16（Results 中的引用 11-22 随 Egger 等位方向句删除） |
| S3 Fig           | R2-P1、R2-P2、R3-P5；Fig 1 图注 |
| S4 Fig           | R2-P1 |
| S5 Fig           | R3-P4；Fig 3C 图注 |
| S6 Fig           | R2-P3 |
| S7 Fig           | R3-P1；Fig 3A 图注 |
| S8 Fig           | R3-P5；M12；Fig 3D 图注 |
| S1 / S2 Table    | R3-P3；M10 |
| S3 Table         | R2-P1 |
| S1 Text          | R1（多处）、R2-P1、R2-P5；M3–M8、M17；Discussion 局限段 |

### 2.4 每一部分该写什么、不该写什么

- **Abstract**：只概括 Results 和 Discussion 中已定的内容；不引文献；< 300 词；说法强度不能超过 Results。
- **Introduction**：背景 → 现有 peak–gene linking 方法（跨个体或跨细胞的相关，估计的是 association）→ MR 能补充什么 → 本文做了什么。不写结果的细节。
- **Results**：
  - 每段写一个发现：先说做了什么（一句话，并写 "see Methods"），再写发现和数字，然后引图。
  - 定义只引用，不重写（R1-P1 的短版除外）。
  - 解释和意义留到 Discussion；参数和软件留在 Methods。
- **Discussion**：主要发现 → 与已有工作的关系 → 局限（causal 的含义、样本重叠、工具变量数和 LD 参考、比较与验证的局限）→ 今后工作。不重复 Results 的数字。
- **Methods**：每个小节按 目的 → 输入 → 步骤 → 理由 → 参数 的顺序写；推导放 S1 Text；阈值和定义与 1.2 一致。
- **S1 Text**：放正文 Methods 没有写的推导和假设说明；不写 Results。

## 3. 全文层面的问题（作者以后统一处理）

- **标题**：`main.tex` 里的标题是 "Demystifying enhancer-gene linkage …"。正文已不再用 linkage（2026-11-22），标题暂按作者意见不动，以后由作者定（候选：enhancer–gene links / regulation；连字符应改为 en dash）；`\ShortTitle` 还没用上。
- ~~**未引用的图**~~：2026-11-22 已全部在 Results 中引用（见 2.3）。
- **文献**（2026-11-18 两轮补引已完成）：占位符全部换成 `\cite`；bib 新增 19 条（PMID / DOI 已核对），删 `fang2001cblb`；Methods 新加短小节 “Terminology”（`sec:terminology`，$\q$ 和 MR-only 的定义）。剩余：部分 bib 作者以 "others" 结尾，投稿前确认 PLOS 格式。
- **软件版本**（3_LEDGER 1.8）：放进 S1 Text 的软件表，还是写在正文各处。
- **暂缓的数据来源**：
  - WNN 外部 PBMC 数据：只剩重新计数时的 Signac 版本（M8）；来源、QC、WNN 聚合 11-09 已按 TRIPOD（jiang2022tripod）写入 M1、M8。10x 数据集页面的引用、sctransform、Signac 文献仍是 `[CITE]`（10x 官网拒绝自动访问，需作者从浏览器取 URL）
  - ~~GTEx、EnhancerAtlas、FANTOM5、4DGenome 对象的来源和基因组版本（M10）~~：2026-11-09 已写入，来源为 TRIPOD 整理的版本（FANTOM5 / 4DGenome 经 HACER；GTEx 引 2017 Nature），正文写 “without coordinate conversion”，不写 hg38
  - ~~AS、CD、MS 三个 GWAS 的来源（M13）~~：2026-11-05 已由作者提供并写入 Methods 和 S4 Table
- **S3 Table 表注**：CIMA 全队列数字（428 人、73 个 L4 类型、±1 Mb、5,051,397 个变异）待作者核对。
- **投稿相关**：选哪个 section（Integrative Omics / Gene Expression and Regulation / Methods，选 Methods 需要模拟验证）；STROBE-MR 核对表；AI 使用披露的范围；Table 1 的 motif logo（PLOS 不允许表格里放图）；代码存档要不要 DOI；投稿系统的各个字段（见 3_LEDGER 第 3 节）。
- **稿件以外**：`Scripts/Pipeline.md:82` 写的“pseudobulk 未按个体细胞数过滤”容易误读，需要改写（需要单独开会话）。
- 从原始 h5ad 重新验证 pseudobulk（不影响正文）。

## 4. 小节工作区

### Results R1

- 已定：
  - 2026-11-22 扩写为 7 段（约 1,000 词）：P1 IV 理由（保留 "alleles are randomly allocated at meiosis"）→ P2 causal IVW / Egger（模型式、GLS、截断、一致性条件、InSIDE、等位对齐、φ）→ P3 NOME 与 sample overlap（偏倚 ∝ 1/F；Results 里可以报告 S1 Fig 的结论 "was small"）→ P4 causal SMR / joint → P5 CIMA 数据 → P6 两种 marginal 方法的思路 → P7 术语（类别名、$\q$、link、MR-only，"see Methods"）。
  - 公式只留 IVW / Egger 的模型式和 SMR 比值，其余性质用文字。
  - R1、Methods、S1 Text 之间可以重复，当作强调，但要换说法。
  - sample overlap 不写具体人数（只说同一 cohort、大部分重叠），R1、Methods、S3 Table 都一样。
- 遗留问题：TensorQTL 待补 [CITE]（Methods M1）。

### Results R2–R4（2026-11-22 审读 + 润色，待作者看；diff：`git diff --no-index ../tmp/ms_rounds/20261122_Results/sections/results.tex sections/results.tex`）

- 已定（作者 11-22 批注）：
  - R2-P1 放总体结果：各方法 link 数、MR-only 数（4,749）、870 组合、S4 Fig 的共享比例（只报数字，不解释 multiomics 的 91%）；multiomics link 多，只回扣 Introduction / S1 Sec 1 的细胞类型混杂（加限定语）。
  - MR 方法之间一致写 “consistent across different MR methods”。
  - S6 Fig 单独成段（R2-P3），支撑 “MR 估计与 marginal 关联基本无关”；不写 “differ substantially”（Abstract 11-18 的决定）。
  - Fig 2：*RPL11* 的表达范围、*ERAP2* 的双峰分布都不解释；*RPL11* / *ARPC1B* 的 pseudobulk 不显著只回扣 S1 Sec 1 的“掩盖”，不做生物学推测。
  - Fig 3A 只描述，不写含义。
  - GO：*SRSF7* / *DDX5* / *UBE2N* 的 term 归属已核对，属实，保留；段末加全基因组背景的限定。
  - S1 / S2 Table：S1（P 值）中 MR **不**更显著（multiomics link 最多，P 最小），只说 S2 Table 的 fold 更大，并注明 P 值受 link 数影响。
  - Fig 3D：阈值 0.05（图注已改）；Egger 的两处例外（S8 Fig 中 Bn HQ、CD8 HQ）不单独点出，正文用 “generally”。
  - Table 1：删去 “Egger 依赖等位方向” 一句（S2 Fig 在 Results 中不再引用）；仍如实写出 6 个都由 Egger 选出、4 个只由 Egger 选出。
  - Fig 4B：β 来自结果表（LEDGER 1.5）；三个例子各一段、结构平行。
- 收尾句（Claude 草拟）：R3 末句、R4-P1 的 “likely reflects their small size” 和末句，作者 11-22 已确认。
- 遗留问题：
  - *ERAP1* 的 peak 同时与 *ERAP2* 相连（β +2.28，强得多且方向相反）；*ERBB3* 的 peak 与 *RPS26* 相连（β −1.23）。审稿人可能质疑 LD / 多效性。作者 11-22 决定暂不处理；写 Discussion 局限段（D4，相邻基因共享 caQTL 变异）时再定。
  - **待办 M14**：Methods GWAS overlap 只写了筛选标准，需补外部注释来源（HiChIP / PCHi-C / EnhancerAtlas / GTEx 同 M10；FAVOR、RegulomeDB、GeneHancer、ENCODE cCRE 用 2026-11 查询；Open Targets 共定位 H4；文献中的 LD 用 1000G EUR，MR SNP 的 LD 用 1000G EAS）。
  - Fig 1A 图内 “Median randomization” 由作者在 Figma 里改（图注已改）。
  - MNT（Sézary 综合征）和 TEF（B-ALL 融合）的人类证据偏弱，作者已知。

### Discussion

- 遗留问题（原 D0–D8 的候选要点，每段一行；写的时候由作者勾选）：
  - 11-22 润色：D2 的 “consistently demonstrated strong biological support”、“immune-related functions”（与 GO 结果不符）、“therapeutic intervention” 已改；作者如果想保留治疗靶点的意思，可以写 “may inform future functional studies”。
  - 局限段以 “First” 开头，但还没有后续；最后如果不补写其他局限（D4 / D6 / D7），就删掉 “First”。
  - 结构：D3 单独成段还是并入 D1？局限分成 D4–D7 四段，还是压缩成两段？
  - D1 主要发现：框架；causal 方法找到的 link 更少；方向一致性（与 pseudobulk 不一致 1.0–3.0%，与 multiomics 不一致 24–26%）；方法之间一致（措辞见 1.3）；SMR 的 OR 和 fold 最高（解释由作者给）。
  - D2 与已有工作的关系：peak–gene linking 方法估计的是 association；分子 QTL 的 SMR / MR 已有先例（对比式写法，不写 first）；与共定位的区别；CIMA 的 GRN；方向不一致主要来自跨细胞类型的相关（解释由作者给）。
  - D3 MR-only link 的生物学信息：GO、motif（依赖 Egger）、S-LDSC（nominal）、GWAS 重叠（没有共定位）。
  - D4 局限：IV 假设；相邻 peak 共享 caQTL 变异（exclusion restriction）；反向因果（没有做 Steiger）；IVW / Egger 的假设。
  - D5 样本重叠 + NOME（2026-11-21 已合并为一段，已定）：两者作用于同一个协方差，只有违反 NOME 时 overlap 才产生偏倚，约 ∝ 1/F；主要论据是强工具变量，null SNP 的 z 相关很小作辅助；对真正有因果的 pair，残差本身相关，所以不写 “negligible”；不写人数；LIML 压缩成两句；今后工作（independent cohorts）一句暂留段末，写 D8 时再移。
  - D6 工具变量数、EAS LD 参考、SMR 没有做 HEIDI、等位方向与 Egger："risk factor" 改为 exposure；补 S2 Fig 的规模（LEDGER 1.7）；方向惯例待补 [CITE]。
  - D7 比较与验证的局限：WNN 是外部、跨细胞类型的数据；pseudobulk 没有协变量；S-LDSC 证据弱且祖先不匹配；GRN 不是金标准；GWAS 重叠没有共定位；只看了 3 个细胞类型和 cis；BH 只在 gene 内做。
  - D8 今后工作：两样本 MR；多变量 MR 或 fine-mapping；HEIDI 或共定位；与等位编码无关的 Egger（GRIP）；细胞类型匹配的 multiome；更多细胞类型；模拟研究。结尾句由作者给。

### Abstract

- 已定（2026-11-18）：
  - 按全文重写，旧稿强于证据的 4 句已替换：不写 identify causal / pleiotropic、differs substantially、validate、significant heritability enrichment；也不再写 428 individuals（那是 CIMA 全队列人数）。
  - 细胞类型写全称（naive B cells、naive CD4⁺ T cells、naive CD8⁺ T cells），不用 Bn cells 等简称；方法写通用名（inverse-variance weighted、MR-Egger、summary-data-based MR、joint test），不写 causal IVW 等本文的方法名。
  - 数字只留 3.28（MR-only peak 与 HiChIP）、six references、seven diseases；Results 或这些数字改动时，Abstract 要同步。
- 遗留问题：
  - 292 词，离 300 词上限很近；作者要加内容时先删方向一致性那一句或 eRegulon 那一句。
  - 备选句（方向，三个细胞类型都成立，数字见 3_LEDGER 1.2）：“…, whereas 24–26% disagreed in direction with the correlation across metacells, compared with 1–3% for the correlation across individuals.”
  - 结论句 “Our results show that MR with population-scale single-cell data complements correlation-based mapping of enhancers to target genes” 由作者定稿；motif 结果没有写进 Abstract。
  - Results 的 “each MR-based test identified fewer links” 只和 marginal multiomics 比；Abstract 写的是 “than either correlation-based method”（3_LEDGER 1.2 的 link 数表在三个细胞类型都成立）。

### Introduction / Author Summary

- 已定（2026-11-18）：
  - 主旨句是作者给的：correlation does not imply causation；用 MR 和 population-scale single-cell QTL 推断因果（保留为文件顶部的注释，P3 段末点题）。
  - 用词：生物学概念写 “enhancer–gene”（P1），被检验的单位写 “peak–gene pair”，显著的写 “link”；三个细胞类型在 P7 首次定义。
  - 不与 scMultiMap / TRIPOD / FusioMR 的 Introduction 雷同：改写后共有的 5 词片段只剩专有名词（scRNA-seq 全称等），以后改动这一节时保持。
  - 最后一段只写定性结论，不放数字；说法按 1.3（不写 outperform / validate / first）。
- 遗留问题：
  - P2 跨个体相关的引文作者已换成 `kamal2023granie`、`zhang2026ra`：GRaNIE 是 bulk RNA-seq + ATAC-seq 跨个体（31–45 人的 iPSC-derived macrophages）的 Pearson 相关，可以用；`zhang2026ra`（Front Immunol 2026）在 PubMed 按标题和关键词都检索不到，bib 无 PMID，待作者核对 DOI 以及它是否做了跨个体的 peak–gene 相关。
  - ~~P3 Xie 2023 的基因名~~：2026-11-22 已写 *Stim2*（PMC10363915 全文确认）。
  - P3 第 2 轮（11-22，已改，待作者看）：首句保留作者原中心句（改为 “biased or masked”）→ 理论 1 句指向 S1 Sec 1 → 后果句（arise / masked / biased / reversed）→ 三个例子，每个点明 U：Xie 2023 *Stim2*（U = shared upstream regulation，原文用语；不写 co-expression，作者不希望突出）、TRIPOD Fig 2C（U = TF，表达为基因所必需、与 peak 可及性负相关 → 掩盖；不写 trio / conditional 细节）、scMultiMap（U = donor；按原文写 pooled cells，旧稿 “donor-level pseudobulk profiles” 的归属错误已改，作者认可）。例子细节（TRIPOD 两类细胞、ψ_Y > 0、Cov(X,U) < 0）放到 S1 Sec 1 新增段；符号对应是我们的解读。Methods 两个 marginal 小节已指向 S1 Sec 1，不改。
  - P4 “few studies have used MR to estimate the effect of chromatin accessibility at a peak on …”：作者的论断，没有引文；强度由作者定（Discussion D2 要交代分子 QTL 之间的 SMR / MR 先例）。
  - P5 多效性的例子 “through another regulatory element in the same region” 取自 Discussion D4 的要点，作者确认是否保留。
  - P7 末句 “Together, these analyses show how MR … can complement correlation-based peak–gene linking” 是收尾句，由作者定稿。
  - 篇幅约 1,230 词（scMultiMap 约 900、FusioMR 约 1,100）；要压缩时先压 P3 和 P7。
- Author Summary 留到最后写：150–200 词，非技术性，用第一人称。

### Methods M1 Data sources

- 遗留问题：
  - "three major naive lymphocyte populations" 里的 "major" 要不要保留？
  - 作者以后会补充更多数据细节（见第 3 节“暂缓的数据来源”）。

### Methods M2 Data Preprocessing

- 已定：删掉 Final analysis set 段，候选 pair 数只放在 R1。
- 遗留问题：Methods 里不再引用 S3 Table（Results 仍引用），这样可以吗？

### Methods M3 Mendelian randomization

- 已定（2026-11-22）：M3–M8 各补一两句 S1 Text 的性质（边际斜率含混杂、NOME 与 overlap 偏倚、IVW 偏倚条件、InSIDE、Σ 不含 γ̂ 的误差、截断参数、φ̂ 只放大 SE、Egger 对测量误差更敏感、CCT 在相关下有效、marginal 方法估计的是含混杂的边际关联）。
- 遗留问题：MR-Egger 式里的误差项写成 $\varepsilon_j \sim N(0, \phi s_{Yj}^2)$（与 S1 Text 一致）。正文要不要只留回归式，去掉误差分布？

### Methods M4 Causal IVW / Causal Egger

- 已定：数学下标已统一为 IVW / Egger；Eckart–Young 已引 `eckart1936approximation`；表示估计值的粗体向量一律带 hat（2026-11-21，与 Results、S1 Text 一致）。

### Methods M5 Causal SMR

- 已定（2026-11-22 按 Zhu 2016 Online Methods 重写）：第 1 段基本 SMR（单工具变量比值 → $\hat\beta_{\mathrm{SMR}}$ 与 $T_{\mathrm{SMR}}$ 同一行 → $1/T$ 形式 → 计入 $\hat\gamma$ 的误差、不依赖 NOME，但假设 cov($\hat\gamma$, $\hat\Gamma$) = 0，指向 overlap 一节 → 单个 SNP 不能区分因果和多效性）；第 2 段本研究：lead SNP、由 P 反推 z、χ²₁（与 `Scripts/Model/MR_AllPairs.R` 一致）、与 IVW 的关系。不写 Var、不写 P 值下限。

### Methods M6–M8

- M6：下标改为 causal 名，补 CCT 在相关 P 值下有效的性质（2026-11-22）。
- M7、M8：外部数据的细节见第 3 节，润色时不补；相关系数下标为 individual / metacell（2026-11-22）。

### Methods M9 Filtering and multiple testing

- 遗留问题：gene 内做 BH 的理由由作者给出；$\hat\phi \le 10$ 的理由按作者的原意写（$\hat\phi$ 约为 10 时 GLS 拟合失败，通常是因为离群的工具变量或 LD 矩阵接近奇异）。

### Methods M10–M15

- M11：DAVID 知识库的版本待确认。
- M14：已定 7 个 GWAS（2026-11-06，去掉 *UBAC2* / Asthma，Methods 写 “seven GWAS”）。

### Methods M16 Allele orientation

- 遗留问题：单独成节，还是并入 M4？

### Methods M17 Error correlation due to sample overlap

- 已定（2026-11-21）：节名改为 “Error correlation due to sample overlap”（`\label{sec:overlap}` 不变）；三段：动机（同一 cohort、接近 one-sample，不写人数；偏倚 ∝ 1/F̄，指向 S1 Text）→ null SNP 的 z 相关估计 ρ_j → 两种估计（pooled 观测不独立；逐 pair t 检验）。只写方法，结果数字放 S1 Fig 图注。

### S1 Text

- **工作卡（2026-11-22，S1 Text 重排；W2 已完成，等作者看）**
  - 作者的决定：
    - 理论按同一样本 N 推导，开头说明实际样本因各细胞类型、各模态单独 QC 略有不同，同一样本是 overlap 偏倚最大的情况；
    - 新增第 3 节 "Assumptions: NOME and InSIDE"，放在 MR model 之后、overlap 之前；
    - S1 Text 标题改为 "Statistical framework of the MR analysis"，各节编号 1–9，正文引用写成 (S1 Text, Section k)；
    - 目录合成一个列表，不分 Text / Figures / Tables 三组；按 PLOS 惯例保留 S1 Text、S# Fig、S# Table 的名字。
  - 作者第二轮意见：
    - 残差相关只写 ρ_ε 一个式子，删去 σ_XY ⇒ β̃ 的推导；
    - NOME 一节加上强工具变量的论证（\|z\| > 5.73、χ² > 32.8，SE < 18%，并说明这不是严格检验）；
    - 偏倚部分改为定性解释，定量结论引 Burgess 2016（已核对原文：Eq 3 和摘要），删去自推的 (κβ̃ − β)/F̄。
  - 节号对应：
    - 1 Marginal
    - 2 MR model
    - 3 Assumptions
    - 4 Overlap
    - 5 IVW（含一致性）
    - 6 MR-Egger
    - 7 Correlated variants
    - 8 Random effects
    - 9 GLS truncation
  - S1 Text 内部用 `\label{s1sec:*}` 和 `\ref`；正文的节号是手写的。
  - 第三轮意见（11-22，已改）：
    - 第 2、5、6 节开头写明这里介绍的是 standard IVW / MR-Egger，causal IVW / causal Egger 是 GLS 推广（Section 7）；
    - Assumptions 开头改为：NOME 关于数据，作用是让 γ̂ 可当作已知；IVW（工具变量全部有效或平衡多效性）和 MR-Egger（允许方向性多效性）都依赖 InSIDE；
    - InSIDE 的两种表述并列写出，各给一式：population（α_j ⫫ γ_j，J → ∞）、fixed-instrument-set（cov_v = 0）；
    - IVW 的 Consistency 段补充分解式 Σv αγ/Σv = cov_v + ᾱ_v γ̄_v，对应这两种表述；
    - 第 4 节改名 "Sample overlap"，删去 Residual variance 段，Egger 段只说推导类似（多了截距），不再提 SMR；
    - IVW 和随机效应两节的权重统一为 v_j（删去 w_j）；
    - 目录中的图表改用短标题（`\caption[...]`），图注正文不变。
  - 第四轮意见（11-22，已改）：
    - S1 Section 2 新增 Fig B（TikZ 绘制的经典 MR 图：U 在上方，α_j 为下方虚线弧，G_j–U 之间用打 × 的点线表示独立性）；Fig A 同步改为 U 在上方。
    - 正文 Methods 的 "Mendelian randomization" 小节（sec:mr-overview）删除，三条 IV 假设的编号列表（lawlor2008mendelian）移到 S1 Section 2，Section 3 开头改为指回它。
    - 正文 Causal IVW 第一段补定义：exposure / outcome、IV 假设、γ̂_j / s_Xj / Γ̂_j / s_Yj / J、NOME 和偏倚，以及“standard 估计量见 S1 Text 5–6”。
    - Causal IVW 模型式后补 β、φ 的定义和一致性条件；Causal Egger 补 α_0（平均多效性）、InSIDE 和方向性多效性的说明。
    - R1 的占位 "(ref method sec)" 改为 (S1 Text, Section 2)。
  - 第五轮意见（11-22，已改）：S1 Section 4 开头重写 Section 2 的两个边际模型（X_i、Y_i）和残差相关 ρ_ε，承接下文；"Correlation between the estimation errors" 一段补上最小二乘误差的显式写法 δ = Σ_i c_ij ε_ij（c_ij 只取决于基因型，caQTL 和 eQTL 共用）。
  - 第六轮意见（11-22，已改，正文）：
    - Methods 的 Causal SMR 删去 "does not rely on NOME" 一句（作者：不提）。
    - Data preprocessing、Filtering、Terminology 合并为 **Analysis settings**（`data_preprocessing.tex`，在 Data sources 之后，三个 label 都保留），BH 只写一遍。
    - 阈值按代码统一为 q < 0.1：原 Filtering 写的 "≤ 0.1" 是错的，作者已确认。
    - "causal / marginal methods" 的定义移进 Analysis settings。
    - Results R1：
      - 公式后紧接 α_0 和 IVW / Egger 的区别；LD、InSIDE、φ、NOME 简写，细节指向 Methods / S1；
      - 删去 SMR 多效性句、"rest on different assumptions"、"We refer to these four tests"、细胞类型句和最后一段；
      - causal joint 的名字并入 CCT 句；q 在 R2 第一次用阈值处定义。
    - 细胞类型统一为 Bn / CD4 T / CD8 T cells（Results、图注、Methods 第 128 行、supplementary 图注）；外部数据和文献中的 naive 不改。
  - 第七轮（11-22，已改）：
    - S6 Fig 图注语法改正；S3 Fig 标题改为 "Overview of results from the six methods in Bn cells and CD8 T cells."；supplementary 图注里单独写的 "Bn" 都补全为 "Bn cells"。
    - R1 给 "balanced pleiotropy" 加了一句解释。
    - Supporting information 列表的标题改为与 supplementary 图表的完整 caption 第一句一致；S1 Text 的描述按现在的 9 节重写；S1 Table 的 "p-values" 改为 "$P$ values"。
- 已定（2026-11-21）：
  - 小节顺序：Marginal association and causal effect → MR model → IVW → Consistency of the IVW estimator → **Sample overlap and measurement error**（开头讲 NOME；独立的 NOME 小节已删）→ MR-Egger / InSIDE → Correlated variants → Multiplicative random effects → GLS 截断。
  - 符号：工具变量 G；Fig A 结构误差 ξ_X、ξ_Y；边际回归残差 ε^X_ij、ε^Y_ij；MR 回归残差 ε_j；"MR model" 结构方程为 ψ_Y U_i + ξ_Yi（U 与 Fig A 相同）；带 hat 为估计值、不带为总体值（全篇，含 Egger 节）；样本集合 𝒮_X、𝒮_Y，N_X、N_Y、N_s，κ = N_s/N_Y，F_j、F̄。
  - 新小节的结论：overlap 只通过违反 NOME 产生偏倚，严格 NOME 下无偏；bias ≈ (κβ̃ − β)/F̄，本数据 ≲ 3%。Egger 只写受 I²_GX 支配（`bowden2016nome`），不给数字；SMR 只写 lead SNP 也满足 χ² > 32.8（单工具变量比值没有有限矩，不套 1/F̄ 公式）；残差方差额外项的大小写成 eQTL / caQTL 期望 z² 之比（不是 O(1/F)）。
  - 不写 CIMA QTL 的协变量（PEER / PC），作者自己核对后再定。
  - S1 Fig 图注写全局 ρ̂（0.0035 / 0.0060 / 0.0053）。
- 遗留问题：
  - “Both cases agree with previous analyses of bias due to participant overlap in two-sample MR \cite{burgess2016overlap}” 一句来自对该文结论的记忆，作者核对原文后再保留。
  - Overdispersion 小节 "In our setting, many instruments may affect gene expression through pathways other than …"，这个措辞合适吗？
  - Correlated variants 小节中 causal IVW 的极限式（β + γᵀΨα / γᵀΨγ）是 Claude 推出来的，A.4 原文只给了 Egger 的式子，请作者核对。
