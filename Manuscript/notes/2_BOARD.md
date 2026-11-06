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
| peak           | 全文写 "peak"（首次出现写 "chromatin accessibility peak"），不写 "enhancer"（标题也要改，见第 3 节）                                                                                                                                                                                                                        |
| pair / link    | 被检验的单位写 "peak–gene pair"，显著的 pair 写 "link"                                                                                                                                                                                                                                                                     |
| GWAS           | 只报告 7 个疾病：AS、CD、IBD、RA、SLE、MS、T1D，统称 "immune-mediated diseases"（不写 "blood traits"）；S-LDSC 实际跑了 9 个，asthma 和 PBC 不报告                                                                                                                                                                          |
| GO             | 只写 DAVID 过表示分析                                                                                                                                                                                                                                                                                                       |
| 格式           | 基因符号斜体；科学计数写 `$6.1\times10^{-4}$`；美式拼写（center）；图写 `Fig~\ref{fig:overview}B`，表写 `Table~\ref{tab:motif}`，补充条目写 `\nameref{S1_Fig}`                                                                                                                                                      |
| 参考文献       | 暂时不处理，最后统一处理（缺的地方写 `[CITE: 论断]`）                                                                                                                                                                                                                                                                     |

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

**Abstract**（1 段，6 句，旧稿）：背景 → 现有方法靠相关 → MR 框架 → CIMA 应用 → Hi-C / HiChIP 验证、motif 和遗传力。状态：未动（有 4 句强于证据，见第 4 节）。

**Introduction / Author Summary**：空。

**Results**（`sections/results.tex`）

| 小节 | 标题                                                                              | 段 | 每段一句话                                                                                                                | 状态     |
| ---- | --------------------------------------------------------------------------------- | -- | ------------------------------------------------------------------------------------------------------------------------- | -------- |
| R1   | Study overview and comparison of peak–gene linking methods                       | P1 | CIMA 数据、三个细胞类型、六个方法分两类、类别名、$P_{\mathrm{adj}}$ 和 MR-only 的定义 → Fig 1A                         | 已调结构 |
|      |                                                                                   | P2 | 同一组候选 pair；marginal multiomics 最多（CD4 37,181），causal 方法更少；重叠 → Fig 1B、S3 Fig                          | 已调结构 |
|      |                                                                                   | P3 | 六个方法都显著的 1,259 个 link：MR 之间几乎全部同号，73% 与两个 marginal 方法同号，24.1% 只与 multiomics 不同号 → Fig 1C | 已调结构 |
| R2   | Concordance and discordance between MR estimates and marginal correlations        | P1 | MR 与两种相关分别估计什么；多数一致，也有不一致，下面看 4 个例子                                                          | 已调结构 |
|      |                                                                                   | P2 | *SIDT2*（Bn）六个方法都一致；*RPL11*（CD4）个体层不显著 → Fig 2A–B                                                  | 已调结构 |
|      |                                                                                   | P3 | *ERAP2*（CD4）MR 与个体层同为负，跨 metacell 为正 → Fig 2C                                                             | 已调结构 |
|      |                                                                                   | P4 | *ARPC1B*（CD4）只有 MR 显著 → Fig 2D                                                                                   | 已调结构 |
|      |                                                                                   | P5 | 小结：关联可能与因果效应不同（结论偏强）                                                                                  | 已调结构 |
| R3   | Functional characterization and regulatory network support for peak–gene links   | P0 | `\TODO`：Fig 3A、3D、S5、S7 还没有正文                                                                                  | 有新改动 |
|      |                                                                                   | P1 | MR-only 基因的 DAVID GO（4,509 个基因）→ Fig 3B–C                                                                       | 已调结构 |
|      |                                                                                   | P2 | 与 CIMA SCENIC+ eRegulon 的重叠，多数比较中 causal 方法的 OR 更高 → Fig 3E、S8 Fig                                       | 有新改动 |
|      |                                                                                   | P3 | 6 个外部参考的超几何富集：SMR 在 17/18 个组合中 fold 最大 → S1 / S2 Table                                                | 已调结构 |
| R4   | Disease heritability enrichment and GWAS overlap of MR-supported peak–gene links | P1 | S-LDSC：MR-only peak 的点估计更高，4 个疾病 P < 0.1 → Fig 4A                                                             | 已调结构 |
|      |                                                                                   | P2 | strong MR-only link 的 GWAS 重叠例子（*ERBB3* 等）→ Fig 4B                                                             | 已调结构 |
| R5   | Motif enrichment highlights candidate regulators of MR-supported peak–gene links | P1 | 前景、背景、47 个 MR-specific motif；Table 1 依赖 Egger                                                                   | 已调结构 |
|      |                                                                                   | P2 | 例子与文献：MNT、ATF3、FOSL2::JUND 等                                                                                     | 已调结构 |

**Discussion**（`sections/discussion.tex`）：只有 D1–D8 的注释骨架，加上两段旧文字：D5（样本重叠，S1 Fig）和 D6（等位方向，S2 Fig）。状态：未动。

**Methods**（`sections/methods.tex` + `data_preprocessing.tex`）

| 小节                                 | 段 | 内容                                                                  | 状态               |
| ------------------------------------ | -- | --------------------------------------------------------------------- | ------------------ |
| M1 Data sources                      | 2  | CIMA 数据；外部 PBMC multiome；三个细胞类型及标签                     | 已润色（11-02）    |
| M0b Study design and terminology     | 2  | 框架；类别名、$P_{\mathrm{adj}}$、pair / link、MR-only 的定义       | 已润色（11-02）    |
| M2 Data Preprocessing                | 4  | 候选 pair → 工具变量 → 等位与去重 → LD 矩阵                        | 已润色（11-02）    |
| M3 Mendelian randomization           | 4  | 一般定义和三条假设 → 记号 → MR-Egger 式（α₀ = 0 即 IVW）→ 本研究 | 已润色（11-02）    |
| M4a Causal IVW                       | 4  | LD → GLS；记号；截断（Eckart–Young）；推断                          | 已润色（11-02）    |
| M4b Causal Egger                     | 1  | 加截距的 GLS、推断                                                    | 已润色（11-02）    |
| M5 Causal SMR                        | 2  | lead SNP、比值、由 P 反推 z、T_SMR、不做 HEIDI                        | 未动（有作者意见） |
| M6 Causal joint                      | 1  | CCT 的公式与尾部近似                                                  | 已调结构           |
| M7 Marginal pseudobulk               | 2  | count 矩阵 → 个体过滤 → ln(1 + CPM) → Spearman                     | 已调结构           |
| M8 Marginal multiomics               | 3  | 80 个 WNN metacell；重新计数；Spearman                                | 已调结构           |
| M9 Filtering and multiple testing    | 1  | 删 NA → $\hat\phi \le 10$ → joint → gene 内 BH                   | 已调结构           |
| M10 Hypergeometric enrichment        | 2  | 6 个参考；N / K / n_sig / k；合并后做一次检验                         | 已调结构           |
| M11 Gene ontology analysis           | 1  | DAVID 输入、类别、作图规则                                            | 已调结构           |
| M12 GRN overlap                      | 1  | 两套 eRegulon；2×2 Fisher 检验                                       | 已调结构           |
| M13 S-LDSC                           | 2  | 注释的构造、参数；7 个疾病和来源                                      | 已调结构           |
| M14 GWAS overlap                     | 1  | strong MR-only link + GWAS P < 5e-8                                   | 已调结构           |
| M15 Motif enrichment                 | 1  | 前景、背景、motif 集、检验和阈值                                      | 已调结构           |
| M16 Allele orientation               | 1  | flip 操作；IVW 不受影响、Egger 受影响                                 | 已调结构           |
| M17 Sensitivity: overlapping samples | 4  | 两个边际模型 → pruning 和 null SNP → 全局和逐 pair 相关             | 已调结构           |

**Supplementary：S1 Text**（`sections/s1_text.tex`，承接 M3–M4，2026-11-02 重组）

| 小节                                                     | 承接正文       | 状态            |
| -------------------------------------------------------- | -------------- | --------------- |
| S1.1 MR model and summary statistics                     | M3             | 已润色（11-02） |
| S1.2 The NOME assumption                                 | M3             | 已润色（11-02） |
| S1.3 Inverse-variance weighted method                    | M3 / M4a       | 已润色（11-02） |
| S1.4 Consistency of the IVW estimator                    | M3             | 已润色（11-02） |
| S1.5 MR-Egger regression and the InSIDE assumption       | M3 / M4b / M16 | 已润色（11-02） |
| S1.6 Correlated variants                                 | M4             | 已润色（11-02） |
| S1.7 Multiplicative random effects and overdispersion    | M4a–b / M9    | 已润色（11-02） |
| S1.8 Truncated eigendecomposition for the GLS estimators | M4a            | 已润色（11-02） |

### 2.3 信息在哪里

**定义和关键数字**（首次给出的位置 → 其他出现的位置）：

- 类别名、$P_{\mathrm{adj}}$、MR-only：M0b → R1-P1（短版）。strong MR-only：M14 → R4-P2。MR-specific motif：M15 → R5-P1。GRN 阈值 0.05：M12 → R3-P2、Fig 3E caption。
- 候选 pair 数 40,087 / 85,527 / 50,931：R1-P2（Methods 不再写）、S3 Table。
- pseudobulk 335 / 372 / 346 人：R1-P1、M7。80 个 metacell：R1-P1、M8。
- caQTL $P<10^{-8}$、clumping 参数：M2。$\hat\phi \le 10$：M9、S1.7。

**图表引用**（2026-11-04；`\TODO` 里的提及不算引用）

| 条目             | 被引用的位置                                 |
| ---------------- | -------------------------------------------- |
| Fig 1A / 1B / 1C | R1-P1 / R1-P2、R3-P2 / R1-P3                 |
| Fig 2A–D        | R2-P2（A、B）、R2-P3（C）、R2-P4（D）        |
| Fig 3A           | **未引用**                             |
| Fig 3B / 3C      | R3-P1；M11                                   |
| Fig 3D           | **未引用**                             |
| Fig 3E           | R3-P2；M12                                   |
| Fig 4A / 4B      | R4-P1 / R4-P2；M14                           |
| Table 1          | R5-P1、R5-P2；M15                            |
| S1 Fig           | D5；M17                                      |
| S2 Fig           | R5-P1；D6；M16                               |
| S3 Fig           | R1-P2、R1-P3、R3-P2；Fig 1 caption           |
| S4 Fig           | **未引用**                             |
| S5 Fig           | **未引用**（只有 Fig 3D caption 提到） |
| S6.1–S6.3 Fig   | **未引用**                             |
| S7.1–S7.3 Fig   | **未引用**                             |
| S8 Fig           | R3-P2；M12；Fig 3 caption                    |
| S1 / S2 Table    | R3-P3；M10                                   |
| S3 Table         | R1-P2                                        |
| S1 Text          | M3、M4a、M4b                                 |

### 2.4 每一部分该写什么、不该写什么

- **Abstract**：只概括 Results 和 Discussion 中已定的内容；不引文献；< 300 词；说法强度不能超过 Results。
- **Introduction**（以后写）：背景 → 现有 peak–gene linking 方法（跨个体或跨细胞的相关，估计的是 association）→ MR 能补充什么 → 本文做了什么。不写结果的细节。
- **Results**：
  - 每段写一个发现：先说做了什么（一句话，并写 "see Methods"），再写发现和数字，然后引图。
  - 定义只引用，不重写（R1-P1 的短版除外）。
  - 解释和意义留到 Discussion；参数和软件留在 Methods。
- **Discussion**：主要发现 → 与已有工作的关系 → 局限（causal 的含义、样本重叠、工具变量数和 LD 参考、比较与验证的局限）→ 今后工作。不重复 Results 的数字。
- **Methods**：每个小节按 目的 → 输入 → 步骤 → 理由 → 参数 的顺序写；推导放 S1 Text；阈值和定义与 1.2 一致。
- **S1 Text**：放正文 Methods 没有写的推导和假设说明；不写 Results。

## 3. 全文层面的问题（作者以后统一处理）

- **标题**：`main.tex` 里的标题是 "Demystifying enhancer-gene linkage …"，要按 1.1 改成 peak 的说法；`\ShortTitle` 还没用上。
- **未引用的图**：Fig 3A、3D，以及 S4、S5、S6.1–S6.3、S7.1–S7.3 Fig。放在哪里、怎么写由作者决定（见 R1、R3 的遗留问题）。
- **文献**：正文有 27 处 `[CITE: …]`；bib 中有些条目以 "others" 结尾，要补全；DAVID、SCENIC+、TensorQTL 等软件也要引用。
- **软件版本**（3_LEDGER 1.8）：放进 S1 Text 的软件表，还是写在正文各处。
- **暂缓的数据来源**：
  - WNN 外部 PBMC 数据：来源（作者说来自 UPenn）、细胞 QC、原始的 WNN 聚类、RNA 聚合方式；笔记本上运行时的 Signac 版本（M1、M8）
  - GTEx、EnhancerAtlas、FANTOM5、4DGenome 对象的来源和基因组版本（M10）
  - ~~AS、CD、MS 三个 GWAS 的来源（M13）~~：2026-11-05 已由作者提供并写入 Methods 和 S4 Table
- **S3 Table 表注**：CIMA 全队列数字（428 人、73 个 L4 类型、±1 Mb、5,051,397 个变异）待作者核对。
- **投稿相关**：选哪个 section（Integrative Omics / Gene Expression and Regulation / Methods，选 Methods 需要模拟验证）；STROBE-MR 核对表；AI 使用披露的范围；Table 1 的 motif logo（PLOS 不允许表格里放图）；代码存档要不要 DOI；投稿系统的各个字段（见 3_LEDGER 第 3 节）。
- **稿件以外**：`Scripts/Pipeline.md:82` 写的“pseudobulk 未按个体细胞数过滤”容易误读，需要改写（需要单独开会话）。
- 从原始 h5ad 重新验证 pseudobulk（不影响正文）。

## 4. 小节工作区

### Results R1

- 已定：
  - R1-P1 用一两句话定义类别名、$P_{\mathrm{adj}}$ 和 MR-only，并写 "see Methods"（作者要求：读者在 Results 开头就要明白这些术语）。
  - R1-P3 用 24.1%（单个扇区），写清是“只与 marginal multiomics 不一致”。
- 遗留问题：
  - P1 有 7 句，是把定义放在段末单独成一组句子，还是把数据介绍缩短？TensorQTL 待补 [CITE]。
  - P2 要不要给出各方法的发现数（现在只给了 multiomics 的 37,181，数字见 LEDGER 1.2）？"for a fair comparison" 要不要保留？S4 Fig（细胞类型之间的共享）和 S6.1–S6.3 Fig（方法两两的相关）要不要在 R1 引用？
  - P3：Bn 和 CD8 T cells 的对应比例（LEDGER 1.2）要不要写进正文？

### Results R2

- 已定：*ERAP2* 的 MR 方向是负（已改）；*ARPC1B* 写成 "effect of peak … on *ARPC1B*"（已改）。
- 遗留问题：
  - P1 里 "which is reasonable" 想表达什么？为什么一致是合理的？
  - P2 的 "no clear pleiotropy"（有 `\TODO`）：删掉，还是补算 Egger 截距（属于补分析）？结论句里的 "within the same cell type" 指的是 pseudobulk 吗？
  - P3 的 "surprisingly" 要不要保留？
  - P5（有 `\TODO`）：保留这一段，还是并入 P4 的结尾？最后一句强于证据（见 1.3）；“未观测混杂”这个解释是不是作者的观点？

### Results R3

- 已定：R3-P3（超几何）2026-11-01 新增，目前只有事实。
- 遗留问题：
  - 小节开头的 `\TODO`：Fig 3A（Miami 图）、Fig 3D（MR-only 与 HiChIP 的校准）、S5 Fig、S7.1–S7.3 Fig（Manhattan 图）需要作者给出要点后再写。小节标题要不要随之调整？
  - P1：举 *SRSF7*、*DDX5*、*ZNF143*、*IKZF3* 想说明什么？
  - P2：要不要点出 causal SMR 在 6 个比较中 OR 都最高？要不要提 causal Egger 的例外？
  - P3：想从这段得出什么结论？要不要和 P2 合并？

### Results R4

- 遗留问题：
  - P1 的 OneK1K 一句（有 `\TODO`）：删掉，还是由作者提供出处？目前没有结论句，作者想说什么？
  - P2：~~*UBAC2*（asthma）~~ 2026-11-06 已从 Fig 4B 和正文去掉，换成 *ERAP1*（Bn，CD）；Methods 改为 “seven GWAS”。最后一句（“Although we cannot directly say…”）语法仍待改写。

### Results R5

- 已定：收紧说法，如实写出“6 个都在 causal Egger 中富集、4 个只由 Egger 选出”，并指向 M16 和 S2 Fig。
- 遗留问题：
  - P1：要不要加作者的限定句（例如“应谨慎解读”）？
  - P2：*ATF3* 只点了名，要不要补一句文献（Table 1 已引用）？CD8 T cells 的两个例子和 CENPB 要不要换（文献证据偏弱）？（11-06：Table 1 已全部改引人类研究，motif 不变；MNT（Sézary 综合征）和 TEF（B-ALL 融合）的人类证据偏弱，作者已知。）

### Discussion

- 遗留问题（原 D0–D8 的候选要点，每段一行；写的时候由作者勾选）：
  - 结构：D3 单独成段还是并入 D1？局限分成 D4–D7 四段，还是压缩成两段？
  - D1 主要发现：框架；causal 方法找到的 link 更少；方向一致性（与 pseudobulk 不一致 1.0–3.0%，与 multiomics 不一致 24–26%）；方法之间一致（措辞见 1.3）；SMR 的 OR 和 fold 最高（解释由作者给）。
  - D2 与已有工作的关系：peak–gene linking 方法估计的是 association；分子 QTL 的 SMR / MR 已有先例（对比式写法，不写 first）；与共定位的区别；CIMA 的 GRN；方向不一致主要来自跨细胞类型的相关（解释由作者给）。
  - D3 MR-only link 的生物学信息：GO、motif（依赖 Egger）、S-LDSC（nominal）、GWAS 重叠（没有共定位）。
  - D4 局限：IV 假设；相邻 peak 共享 caQTL 变异（exclusion restriction）；反向因果（没有做 Steiger）；IVW / Egger 的假设。
  - D5（现有第 1 段）样本重叠：只能说相关的数值很小；“SNP 越多越独立”一句不准确，要删掉或改写；今后工作那句移到 D8。
  - D6（现有第 2 段）工具变量数、EAS LD 参考、SMR 没有做 HEIDI、等位方向与 Egger："risk factor" 改为 exposure；补 S2 Fig 的规模（LEDGER 1.7）；方向惯例待补 [CITE]。
  - D7 比较与验证的局限：WNN 是外部、跨细胞类型的数据；pseudobulk 没有协变量；S-LDSC 证据弱且祖先不匹配；GRN 不是金标准；GWAS 重叠没有共定位；只看了 3 个细胞类型和 cis；BH 只在 gene 内做。
  - D8 今后工作：两样本 MR；多变量 MR 或 fine-mapping；HEIDI 或共定位；与等位编码无关的 Egger（GRIP）；细胞类型匹配的 multiome；更多细胞类型；模拟研究。结尾句由作者给。

### Abstract

- 遗留问题：
  - 旧稿中强于证据的 4 句：
    - "directly identify causal and pleiotropic …"（T6 禁用；pleiotropic 没有对应的结果）
    - "our method … differs substantially"（只有一部分成立）
    - "orthogonally validate …"（参考数据集共有 6 个；validate 太强）
    - "significant enrichment of … heritability"（只有 nominal 证据）
  - 428 人是 CIMA 全队列的人数，不是本文 QTL 的样本量。
  - 逐句要点（等 Results 和 Discussion 定稿后再写）：背景 → 现有方法估计 association → MR 框架（承认方法不是本文首创）→ CIMA 的 3 个细胞类型和候选 pair → 主要比较结果 → 外部数据的支持（GRN / 超几何）→ MR-only link 的特征 → 结论句（由作者给）。

### Introduction / Author Summary

- 留到最后写。Author Summary 150–200 词，非技术性，用第一人称。

### Methods M1 Data sources

- 遗留问题：
  - "three major naive lymphocyte populations" 里的 "major" 要不要保留？
  - 作者以后会补充更多数据细节（见第 3 节“暂缓的数据来源”）。

### Methods M2 Data Preprocessing

- 已定：删掉 Final analysis set 段，候选 pair 数只放在 R1。
- 遗留问题：Methods 里不再引用 S3 Table（Results 仍引用），这样可以吗？

### Methods M3 Mendelian randomization

- 遗留问题：MR-Egger 式里的误差项写成 $\varepsilon_j \sim N(0, \phi s_{Yj}^2)$（与 S1 Text 一致）。正文要不要只留回归式，去掉误差分布？

### Methods M4 Causal IVW / Causal Egger

- 遗留问题：
  - 数学下标 `\text{GLS-IVW}` / `\text{GLS-Egger}` 要不要改成 IVW / Egger？
  - Eckart–Young 定理要引用 Eckart & Young 1936（并入文献统一处理）。

### Methods M5 Causal SMR

- 遗留问题（作者意见，还没有落实）：参照 `paper_reference/` 中 Zhu et al. 2016 的 Online Methods（"Summary data–based Mendelian randomization analysis"）重写：先讲基本的 SMR，再讲本研究的 causal SMR 怎么算。$\hat\beta_{\mathrm{SMR}}$ 的公式之后，在同一行给出 $T_{\mathrm{SMR}}$；不写 $p=\max(\cdot)$ 这类细节，也不写 Var($\hat\beta$)；直接写 z 统计量和 $T_{\mathrm{SMR}}$ 怎么构造。

### Methods M6–M8

- M6：作者意见是“暂时这样”。
- M7、M8：外部数据的细节见第 3 节，润色时不补。

### Methods M9 Filtering and multiple testing

- 遗留问题：gene 内做 BH 的理由由作者给出；$\hat\phi \le 10$ 的理由按作者的原意写（$\hat\phi$ 约为 10 时 GLS 拟合失败，通常是因为离群的工具变量或 LD 矩阵接近奇异）。

### Methods M10–M15

- M11：DAVID 知识库的版本待确认。
- M14：已定 7 个 GWAS（2026-11-06，去掉 *UBAC2* / Asthma，Methods 写 “seven GWAS”）。

### Methods M16 Allele orientation

- 遗留问题：单独成节，还是并入 M4？

### Methods M17 Sensitivity: overlapping samples

- 遗留问题："MR study is usually base on Two-Sample data …"、"Let's recall" 这类口语，润色时改写。

### S1 Text

- 遗留问题（2026-11-02 润色后留下）：
  - NOME 小节最后一句 "All instruments in this study had caQTL nominal P < 10⁻⁸" 只陈述事实，要不要保留？
  - Correlated variants 小节最后转述了原文的提醒：变异之间都相关时，相关变异的 InSIDE 条件不太可能满足。这对本文不利，保留还是删掉？
  - Overdispersion 小节 "In our setting, many instruments may affect gene expression through pathways other than …"，这个措辞合适吗？
  - Correlated variants 小节中 causal IVW 的极限式（β + γᵀΨα / γᵀΨγ）是 Claude 推出来的，A.4 原文只给了 Egger 的式子，请作者核对。
