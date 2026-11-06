# 1_WORKFLOW：论文怎么改、怎么润色

路径相对 `Manuscript/`。每次论文会话都读本文件和 `2_BOARD.md`（第 1、2 节和当前小节）；需要数字、图表、投稿或同步信息时，再查 `3_LEDGER.md`。只有这三个文件加 `AGENTS.md`，不新建零散的笔记。

## 1. 分工

- **作者**：决定文章架构、段落顺序、解释和结论；回答工作卡里的问题；看过改动、合格就留下。
- **Claude**：核心任务是**润色**，把非母语的英文改得流畅、地道、准确、学术。同时联系全文其他部分的信息，提醒“这一段该写什么、不该写什么”。架构由作者决定，Claude 照着实施，再更新 2_BOARD 第 2 节。
- **CodeX**：只在投稿前做一次独立审阅（第 6 节 P-CodeX）。
- 科学事实以 `../Scripts/Pipeline.md` 和代码为准，数字以 `../Results/` 和 3_LEDGER 第 1 节为准。

## 2. 一轮流程（一个小节一轮，一个会话）

Results 和 Methods 一个小节一轮；Discussion、Abstract 整节一轮；S1 Text 可以几个小节合成一轮。

1. **W1 审读（不改 .tex）**：在 2_BOARD 第 4 节对应的小节下写工作卡（模板见下）。
2. **作者回答**：在每个 `作者：` 后写一句；结构上的想法写在 `补充：`。
3. **W2 润色（一次完成）**：
   - 先把要改的 .tex 按原路径复制到快照目录，作为本轮的快照（例：`mkdir -p ../tmp/ms_rounds/20261105_R2/sections && cp sections/results.tex ../tmp/ms_rounds/20261105_R2/sections/`）。
   - 实施作者定下的结构改动，并按第 4 节逐句润色。
   - 跑 `python tools/ms_check.py invariants --snapshot ../tmp/ms_rounds/<日期>_<小节>` 和 `bash -l ../tmp/build_ms.sh`。
   - 在聊天里给出：
     - change log：只列实质改动；
     - **表达要点**：挑 3–5 处“原句 → 新句 → 为什么这样更地道”；
     - 自检清单的结果；
     - 本轮 diff 命令：`git diff --no-index ../tmp/ms_rounds/<…>/sections/<文件> sections/<文件>`。
   - 改写中发现的新问题，追加到工作卡里，记为 Q(n+1)【新】。
4. **作者查看**：在 VS Code 里看改动。合格的就留在正文；有意见就提出，Claude 再改一轮（同一小节最多两轮，第二轮只处理指出的问题）。
5. **W3 收卡**：
   - 删掉工作卡。
   - 该小节的“已定”里只留以后仍会约束写作的决定。
   - 更新 2_BOARD 2.2 的状态和 2.3。
   - 在 3_LEDGER 第 5 节记一行。

**commit 和同步**：不按轮次 commit。每隔一段时间，或者论文有较大进展时，等作者说一声，再 commit 到 GitHub，并用 `bash tmp/overleaf_sync.sh push`（在项目根目录运行）同步到 Overleaf。这一步也记进 3_LEDGER 第 5 节。

**架构调整**：作者随时提出，Claude 实施后更新 2_BOARD 第 2 节（2.2、2.3）。

### 工作卡模板

```
#### 【工作卡】Results R2（W1 2026-11-05）

**A. 变动提示**
- 上次收卡后的变化：……（Overleaf pull 进来的改动 / 新图、新数字 / 无）
- 涉及本节的全文问题（第 3 节）：……
- 作者更新提示：

**B. 段落地图**
| 段 | 现在说什么（一句中文） | 提示（保留 / 合并 / 移动 / 拆分；由作者决定） |
|---|---|---|

**C. 全文联系提醒**（逐段：该写什么、不该写什么，依据 2.3、2.4、1.3 和其他章节）
- P1：……

**D. scMultiMap 的启发**（第 4.4 节对照表中对应的那一行；附页码）
- ……

**E. 问题**
Q1 ……（背景一行：原文位置或证据）
   建议：……
   作者：

**F. 润色预告**：这一节语言上最主要的 2–4 个问题。

补充：
```

## 3. 每轮的固定动作

1. 在项目根目录运行 `bash tmp/overleaf_sync.sh status`。如果 Overleaf 上有新改动，先告诉作者，`pull` 之后再开始。
2. 定下这一轮要改的 .tex。在这一轮完成之前，作者不要在 Overleaf 上改同一个文件（同步时整个文件覆盖）。
3. 检查：`python tools/ms_check.py invariants --snapshot DIR`（或 `--base REF`）、`placeholders`，然后 `bash -l ../tmp/build_ms.sh`。
4. 在 3_LEDGER 第 5 节记一行：日期 | 工具和模型 | 改了什么 | 怎么核对。

## 4. 润色规范

### 4.1 目标

流畅、地道、准确、学术。读者是 PLOS Genetics 的读者：懂遗传学，但不一定懂 MR 或单细胞 multiome。意思不变，作者的限定语保留，清楚的句子不动；作者自己重写过的句子，只查拼写和语法。

### 4.2 句子层面

- 一句话只讲一个主要意思；超过约 35 个词就拆开。
- 已知信息放在前面，新信息放在后面；主语和谓语尽量靠近。
- 少用名词化（"performed an analysis of" → "analyzed"）。
- 同一个概念全文只用一个词（按 2_BOARD 1.1），不为了避免重复而换同义词。
- 先说发现，再给数字。限定语要和证据强度对应（nominal、tended to）；没有检验支撑时不加 significantly、robust、novel。不把 "suggest" 改成 "demonstrate"。
- **确实是错误或不够正式，要改**：
  - 冠词和单复数
  - 错误的搭配（"related with" → "related to"）
  - 口语化的评论（"which is reasonable"、"surprisingly"、"Let's recall"）
  - 指代不清（"the results of X were …"）
  - 推论没有依据就写 "Therefore / This demonstrates"（这里要改的是推论本身，不是这个词）
- **可以用，但不要过度**：
  - Notably、Interestingly、Moreover、Furthermore，每节一两次就够（scMultiMap 也用）。
  - 破折号和嵌套括号能拆就拆。
  - delve、leverage、underscore、pivotal 本身没有错，但被写滥了，能换成平实的词就换。
- **格式**：基因符号斜体；`Fig~\ref{}`、`\nameref{}`；美式拼写；按 1.1 写科学计数。

### 4.3 段落层面

- 首句给出本段论点，接着是证据，最后一句是作者的解释或过渡（如果有）。
- 段与段之间有一句承接。
- 解释性、结论性的句子只来自作者（见第 5 节红线）。

### 4.4 参照 scMultiMap（Su et al. 2025, Nat Commun）

`../paper_reference/Su et al. - 2025 - scMultiMap ….pdf`（只读）。用 `/sw/eb/sw/poppler/24.04.0-GCC-13.2.0/bin/pdftotext <pdf> -` 读，输出到 stdout，不写文件。页码是期刊页码。

**只学写法**：句子怎么写、段落里的内容按什么顺序展开、版面怎么安排、什么位置写什么。不照搬它的架构和观点，也不学它的结论强度；措辞强度一律按 2_BOARD 1.3（例如不写 outperform）。

**对照表**：本文很多分析和 scMultiMap 相似，写到这些位置时，直接看它对应的段落。

| 本文                                              | scMultiMap 对应段落                                                                                                                                                                                                                  | 可借鉴的写法                                                                                                                                                                              | 不能照搬的地方                                                 |
| ------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------- |
| R1-P1 研究概述                                    | Results "Overview of scMultiMap"（p2–3）                                                                                                                                                                                            | 先交代输入和模型，每个设计选择都说明它相对现有方法的理由，再写 "More details … can be found in the Methods section"；用 "In summary, scMultiMap takes … as input and generates …" 收尾 | 它用公式介绍新模型；我们的方法不是新提出的，只说明框架和类别名 |
| R1-P2–P3 发现数、重叠和方向                      | benchmark 段（p3–4）、reproducibility 段（p4–5）                                                                                                                                                                                   | 每块按 "To evaluate X, we …" → 设计为什么能检验 X →（Methods）→ "Fig. 1a shows that …" → "In comparison, …" 展开；数字放在句末                                                     | 我们没有模拟，也没有 type I error 的结果；不写“更准”         |
| R2 例子                                           | PICALM 段（p7）                                                                                                                                                                                                                      | 按“发现 → 已知生物学 → 变异和数字 → 方法补充了什么 → 加限定的含义（"This implies that … may …"）”的顺序写                                                                         | 我们的例子是说明方法之间的差异，不是疾病机制                   |
| R3 外部参考的重叠（超几何、Fig 3D、S5 Fig）和 GRN | 一致性验证段（p4："Fig. 2b shows that peak-gene pairs identified by scMultiMap are the most consistent with promoter capture Hi-C, H3K27ac HiChIP, and cell-type-specific eQTL data …"）；与 SCENIC+ 的比较（p4）                   | 先说用了哪几类外部数据（Methods），再写结果；重叠是否显著放在括号里；百分比放到补充图                                                                                                     | 不写 "most consistent" 这种排名式结论，只描述 fold 和 OR       |
| R3 GO                                             | GO 段（p4–6）；Methods "GO enrichment analysis"（p10）                                                                                                                                                                              | Results 只写最主要的几个 term 和它们的功能主题；Methods 一句话写清工具、检验方法和阈值                                                                                                    | 我们的背景是全基因组，fold 小，不说“MR 特有的功能”           |
| R4 S-LDSC                                         | p5："we hypothesize that a more powerful and accurate method … should yield higher heritability enrichment"；"It is important to note that …"；"Fig. 4a shows that …"；Methods "LDSC analysis for heritability enrichment"（p10） | 先写预期，再交代适用范围，最后写结果；Methods 按 目的 → S-LDSC 算什么（一句话）→ 输入 → 协变量（baselineLD）及理由 的顺序写                                                            | 我们只有 nominal 证据，不写 "higher and more significant"      |
| R4 GWAS 位点例子                                  | PICALM 段（p7）                                                                                                                                                                                                                      | 先写已知位点，再写变异和数字，然后写本文补充了什么，最后加限定（"may contribute to … by …"）                                                                                            | 我们没有做共定位，也没有做差异分析                             |
| Methods：对照方法（M7、M8）                       | "Other methods for statistical inference of peak-gene associations"（p9）                                                                                                                                                            | 写明对照方法的输入、版本和参数，以及为了公平比较做了什么（"to ensure a fair comparison"）                                                                                                 | —                                                             |
| Methods：外部参考（M10）                          | 一致性分析的 Methods（p9："for promoter capture Hi-C data, we used interactions with a CHiCAGO score ≥5; for H3K27ac HiChIP data, we used the significant interactions based on FitHiChIP (FDR <0.1)"）                             | 每个数据集一句话：用了哪部分、阈值多少、阈值的出处                                                                                                                                        | —                                                             |
| Methods：数据预处理（M2）                         | "Preprocessing of fragment counts"（p10）                                                                                                                                                                                            | 先写选择及其理由（附引用），再写步骤，工具带版本号（"MACS2 (v2.2.9.1)"）                                                                                                                  | —                                                             |

**通用写法**（作为参考，不是硬性要求）：

*Results*

- 每个分析的开头写目的和做法（"To evaluate X, we …"），后面跟 "(see Methods)"。
- 发现句可以让图作主语（"Fig 2A shows that …"）；补充图放在括号里。
- 与对照方法比较时用 "In comparison" / "In contrast"；原因加限定（"possibly because"）；打平或不显著的结果也如实写出，并给出阈值。
- 正文只留关键的效应量和 P 值，完整的统计量放到图和表里。
- 衔接用 "First / We then / We next / Finally"、"We further …"。
- 小节标题：scMultiMap 第一个小节用名词短语，其余用结论句。本文用哪种由作者决定。

*Methods*

- 每个小节按 目的 → 输入 → 步骤 → 理由（"This is because …"）→ 参数 → "Details are given in S1 Text" 的顺序写。
- 公式前先说依据（"Based on …, we have"）；记号在第一次使用前用 "Let … denote" 定义；推导和算法放进 S1 Text。
- 软件写成 "R package X (v1.14.0)"；每个阈值都写明出处。

*时态*

- 用 "we"，主动语态。
- Results：做了什么用过去时；图在说明什么、一般事实用现在时。
- Methods：操作用过去时；数学用现在时。
- PLOS 写法："Fig 1A"、"S1 Fig"。

### 4.5 自检清单（W2 结束前逐条检查，在聊天里报告）

- [ ] 每段的首句给出论点，段与段之间有承接
- [ ] 术语和定义与 2_BOARD 1.1、1.2 一致；数字与 3_LEDGER 第 1 节一致，与 2.3 记录的首次出现处一致
- [ ] 说法不超过 1.3 的边界；没有新增作者没给过的解释或结论
- [ ] invariants 没有计划外的差异；编译通过，log 里没有 undefined
- [ ] 没有拼写错误；美式拼写；基因符号斜体；图表引用格式正确

## 5. 红线

1. **不改事实**：数字、P 值、阈值、样本量、`\cite` key、`\label`、公式和图表内容都不改。发现事实错误只报告，由作者决定。
2. **不凭空添加**：结果、方法细节、引用、结论都不编。缺的地方写 `\TODO{}` 或 `[CITE: 要支持的论断]`。不凭记忆添加参考文献。
3. **解释和结论只来自作者**。
4. **只写 `Manuscript/`**：`../Scripts/`、`../Results/`、`../paper_reference/` 只读；辅助脚本放在 `../tmp/Manuscript_tables/`，每轮快照放在 `../tmp/ms_rounds/`。
5. **不 commit、不 push**（GitHub 和 Overleaf 都一样），除非作者在本次会话里明确要求。

## 6. Prompt 模板

作者一句话就能触发，例如 `W1 Results R2`、`W2 Results R2`、`W3 Results R2`。Claude 按下面的内容执行。

**W1（审读）**

```
Read notes/1_WORKFLOW.md, notes/2_BOARD.md (sections 1-2 and the <unit> entry) and the <unit> text.
Check notes/3_LEDGER.md section 1 for numbers and the scMultiMap row in 1_WORKFLOW 4.4 (read those pages).
Write the work card (template in 1_WORKFLOW section 2) under the <unit> entry of 2_BOARD; fold in its
遗留问题. Do not edit any .tex. Questions only where meaning or a decision is needed.
```

**W2（润色）**

```
Snapshot the target .tex into ../tmp/ms_rounds/<date>_<unit>/ (mirrored paths, e.g. sections/results.tex). Apply the author's answers and 补充,
then polish <unit> following 1_WORKFLOW section 4 (native, fluent, precise, academic English; keep
numbers, citations, labels, equations and hedging). Run invariants --snapshot and build_ms.sh.
Report: change log, 3-5 expression notes (before -> after -> why), checklist 4.5, diff command.
Append new questions to the card as Q(n+1)【新】.
```

**W3（收卡）**

```
Delete the <unit> work card; keep only lasting decisions under 已定; update 2_BOARD 2.2 status and 2.3;
append one line to 3_LEDGER section 5.
```

**局部检查**

```
Check only spelling and grammar in lines <a>-<b> of sections/<file>.tex. Do not change wording,
word order, terms, numbers or citations.
```

**P-CodeX（投稿前一次）**

```
(CodeX) You get only: git diff <base> -- Manuscript/sections. List every change that alters meaning:
numbers, citations, claim strength, terminology, "causal" wording, added or removed qualifiers,
negations, method details. Problems only, with file:line and before/after.
```
