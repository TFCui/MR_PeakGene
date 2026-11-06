# Figure 3 (2026-11-04 起 4 个 panel): A = CD4 T cell 的 Miami 图 (Causal Joint 向上 / Marginal Multiomics 向下; 原 S7 第 1 页);
#           B = MR-Only gene 的 GO 富集 (DAVID; BP、MF 两个 facet 合为一个 panel);
#           C = MR-Only link 与 HiChIP 的分箱校准图 (S5 panel E 的同一口径, 由 Plot13_MROnly_ObsExp.R 存 rds; S5 仍保留 6 个库);
#           D = CD4 T cell 各方法显著 pair 与 CIMA SCENIC+ GRN 的重叠 (odds ratio); Bn / CD8 T cell 的移到 Figure S8 (本脚本一并输出)
# 版式 (2026-11-04): 各 panel 不写 title / subtitle, 统一用浅薄荷色小标签条 (theme_strip_mint, 00_common.R): A "CD4 T cell" / B 两个 GO 类别 / C "HiChIP" / D "CD4 T cell"
# --- A: 由 Plot15_Manhattan.R 构建并存为 tmp/Plot15_Manhattan/miami_CD4.rds, 这里只读入拼图 (须先跑 Plot15) ---
# 配色 (2026-11-04, 薄荷色板 PAL_MINT, 00_common.R): A causal 钢蓝 / marginal 薄荷绿 / 例子珊瑚 (见 Plot15); B GO p.adj 浅薄荷 → 青绿 → 深蓝;
#   C HiChIP 深蓝分箱 + 薄荷 gene 点 (见 Plot13); D / S8 GRN 浅薄荷 (cell-type) / 钢蓝 (high-quality)
# 由 Plot9_GO.R (原 Figure 3) 与 Plot_Sup2_GRN.R (原 Figure S5) 合并 (2026-10-29); S5 的 panel A (重叠 pair 数) 不再画, S5 Fig 已从论文删除
#
# --- B / C: GO (DAVID Functional Annotation Chart; 3 个细胞类型的 MR-Only gene 合并为一个 list) ---
# gene list: Results/Validation/Valid7_GO/DAVID_input/MROnly_genes_all.txt (3 个 CT 并集, 4,509 个 gene), 在 DAVID 网站上运行
#   (OFFICIAL_GENE_SYMBOL, Homo sapiens, 默认 background; 只导出 GOTERM_BP_DIRECT / GOTERM_MF_DIRECT)
# DAVID 只导出 P < 0.1 的 term; FDR 用 DAVID 的 FDR 列 (BH, 每个类别内部)
# B = GO BP, C = GO MF; 每个 panel 按 FDR 取前 10 个 term (并列按 P-Value)
#   dotplot (仿 enrichplot::dotplot): x = GeneRatio (Count / List Total), 点大小 = Count, 颜色 = FDR (图例写作 p.adj, 图中不再出现 "FDR")
# 两个 panel 共用颜色标尺 (FDR, log10 刻度; FDR < 1e-10 的统一显示为最深色, 否则 GO MF "protein binding" FDR ~1e-85 会把其余颜色压淡)
#   和点大小标尺 (Count, log10 刻度), 图例合并放右侧
# 输入: Results/Validation/Valid7_GO/DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv (DAVID 以 list 第一个 gene 命名)
#
# --- E (+ S8): GRN odds ratio ---
# 显著 (p.adj < 0.05, 与 Valid5_GRN 一致, 比正文的 0.1 严格) pair 在 GRN 中的 odds ratio (GRN vs 非 GRN, 2x2 Fisher; 每个细胞类型一个小面板,
#   CT 与 CT_HQ 两种 eRegulon 集合并排); CT = 该细胞类型中 active (AUC > 1) 的 eRegulon; CT_HQ = CT 与对应 lineage 高质量 eRegulon 的交集
#   (流程见 Scripts/Pipeline.md 3.5); 统计量已在 Valid5_GRN 中算好, 这里只画图; 说明 (阈值 / GRN 集合 / pair 数) 写在论文 caption
# 输入: Results/Validation/Valid5_GRN/GRN_sig_counts_by_method.tsv, GRN_summary.tsv; A 的输入: tmp/Plot15_Manhattan/miami_CD4.rds; D 的输入: tmp/Plot13_MROnly/hichip_fig3.rds
#
# 输出: Results/Figures/Figures_20261102/Figure_3.pdf + Figure_S8.pdf (只出 PDF)
# 运行: Rscript Scripts/Visualization/Plot9_GO_GRN.R   (小表, 登录节点即可)

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages({ library(patchwork); library(scales) })

OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261102")
HEIGHTS <- c(1.05, 2)           # A (Miami) : 下半 2×2 网格 的高度比 (2026-11-04 第 2 轮)
GRID_W  <- c(0.85, 1)            # 下半网格: 左列 (GO BP / MF) : 右列 (HiChIP / GRN) 宽度比
LEG_H   <- 2.75                  # C 图例: p.adj 色条高度 (cm), 与 Count 图例同高
BASE3   <- 11                    # Figure 3 基础字号 (第 2 轮加大, 减少空白感)

theme_pub <- function(base = BASE3) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      plot.title         = element_text(face = "bold", size = base + 1.5, colour = "#111827"),
      plot.subtitle      = element_text(size = base - 1, colour = "#6B7280", margin = margin(b = 4)),
      axis.line          = element_line(colour = COL_AX3, linewidth = 0.5),
      axis.ticks         = element_line(colour = COL_AX3, linewidth = 0.45),
      axis.text.x        = element_text(size = base - 1, colour = COL_AX3),
      axis.text.y        = element_text(size = base - 1.5, colour = "#1F2937", lineheight = 0.82),
      axis.title.x       = element_text(face = "bold", size = base, colour = "#1F2937", margin = margin(t = 3)),
      axis.title.y       = element_blank(),
      panel.grid.major.x = element_line(colour = "#EEF0F3", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "#F5F6F8", linewidth = 0.3),
      legend.title       = element_text(face = "bold", size = base - 0.5),
      legend.text        = element_text(size = base - 1),
      plot.title.position = "plot",
      plot.tag           = element_text(face = "bold", size = 18),
      plot.margin        = margin(3, 6, 3, 3)
    ) + theme_strip_mint(base)
}

# 带完整边框的 panel (GO / GRN 与 HiChIP 一致): 边框加粗, 不用 L 形轴线
TH_BOX <- theme(axis.line = element_blank(), panel.border = element_rect(colour = COL_AX3, fill = NA, linewidth = 0.9),
                axis.ticks = element_line(colour = COL_AX3, linewidth = 0.6), axis.ticks.length = unit(1.6, "mm"))

###================================================================================================
### B / C: GO
DAVID_FILE <- file.path(VAL_DIR, "Valid7_GO/DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv")
CATS       <- c("GOTERM_BP_DIRECT", "GOTERM_MF_DIRECT")
TOP_N      <- 10
CAT_LAB    <- c(GOTERM_BP_DIRECT = "GO biological process", GOTERM_MF_DIRECT = "GO molecular function")   # 标签条文字
COL_GO     <- unname(PAL_MINT[c("m7", "m6", "m5", "m4", "m2")])   # FDR 越小颜色越深: 深蓝 → 钢蓝 → 青绿 → 薄荷 → 浅薄荷 (2026-11-04 薄荷色板)

res <- fread(DAVID_FILE)
setnames(res, c("P-Value", "List Total", "Fold Enrichment", "Full Term"), c("pvalue", "list_total", "fold", "full_term"))
msg("DAVID 结果: 每个类别的 term 数 (只含 P < 0.1) / FDR < 0.05 / FDR < 0.1")
print(res[, .(n = .N, fdr_lt_005 = sum(FDR < 0.05), fdr_lt_01 = sum(FDR < 0.1)), by = Category])
res <- res[Category %in% CATS]

## 每个 panel 前 10 个 term
top <- res[order(match(Category, CATS), FDR, pvalue)][, head(.SD, TOP_N), by = Category]
stopifnot(uniqueN(top$Category) == length(CATS))
top[, gene_ratio := Count / list_total]
top[, term := stringr::str_wrap(Term, 24)]   # 长 term 折成两行 (第 2 轮: 40 → 24), 让 GO panel 变窄
print(top[, .(Category, Term, Count, list_total, gene_ratio = round(gene_ratio, 3), fold, FDR)])

# 共用标尺的范围; FDR 颜色下限截在 1e-10
FDR_FLOOR  <- 1e-10
lim_fdr    <- c(FDR_FLOOR, max(top$FDR))
lim_count  <- range(top$Count)
fdr_breaks <- 10^c(-10, -7, -4, -1)
fdr_breaks <- fdr_breaks[fdr_breaks <= lim_fdr[2]]
fdr_labels <- ifelse(fdr_breaks == FDR_FLOOR, "≤ 1e-10", format(fdr_breaks, scientific = TRUE))

# BP / MF 各一个 panel (2026-11-04 第 2 轮: 竖排在左列), 各带标签条; 不写 title / "n = … genes" (写在 caption)
# 两个 panel 共用颜色 / 大小标尺; 图例只放在 MF panel 内右下角 (MF 只有 protein binding 在右侧, 右下空白)
go_panel <- function(cat, legend = FALSE) {
  d <- copy(top[Category == cat])
  d[, ykey := factor(full_term, levels = rev(full_term))]   # FDR 最小的在最上面
  d[, cat_lab := CAT_LAB[[cat]]]
  ggplot(d, aes(y = ykey)) +
    geom_point(aes(x = gene_ratio, size = Count, fill = FDR), shape = 21, colour = PAL_MINT[["m7"]], stroke = 0.35) +
    facet_wrap(~ cat_lab) +
    scale_size_continuous(name = "Count", limits = lim_count, range = c(1.8, 6.5), trans = "log10", breaks = c(50, 100, 300, 1000, 3000)) +
    scale_fill_gradientn(name = "p.adj", colours = COL_GO, trans = "log10", limits = lim_fdr,
                         oob = scales::squish, breaks = fdr_breaks, labels = fdr_labels) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    scale_y_discrete(labels = setNames(d$term, d$full_term)) +
    guides(size = guide_legend(override.aes = list(fill = PAL_MINT[["m4"]]), order = 2),
           fill = guide_colourbar(order = 1, barheight = unit(LEG_H, "cm"), barwidth = unit(0.45, "cm"))) +   # 色条与 Count 图例同高
    labs(x = "GeneRatio") +
    theme_pub() + TH_BOX +
    theme(legend.position = if (legend) "inside" else "none", legend.position.inside = c(0.985, 0.03), legend.justification = c(1, 0),
          legend.box = "horizontal", legend.box.just = "top", legend.background = element_rect(fill = alpha("white", 0.9), colour = NA),
          legend.key.height = unit(0.45, "cm"), legend.spacing.x = unit(4, "mm"))
}

###================================================================================================
### E: GRN odds ratio
GRN_FDR <- 0.05   # 与 Valid5_GRN 一致 (比正文主图的 0.1 严格)
sig_tab <- fread(file.path(VAL_DIR, "Valid5_GRN/GRN_sig_counts_by_method.tsv"))
grn_sum <- fread(file.path(VAL_DIR, "Valid5_GRN/GRN_summary.tsv"))
stopifnot(nrow(sig_tab) == 2 * length(CELLS) * length(METHOD_KEY), all(sig_tab$method %in% METHOD_KEY))

# 日志: 写 caption 用的数字
cat("\n===== 每个 GRN 集合 x 细胞类型: eRegulon 数 / 与 MR 重叠的 GRN pair 数\n")
print(grn_sum[, .(grn_set, cell_type, n_eRegulon, n_pair_inMR, n_MR, pct_MR_in_GRN = round(pct_MR_in_GRN, 2))])
cat("\n===== 显著且在 GRN 中的 pair 数\n")
print(dcast(sig_tab, grn_set + cell_type + n_pairs_GRN ~ method, value.var = "n_sig_GRN"))
cat("\n===== OR\n")
print(dcast(sig_tab[, .(grn_set, cell_type, method, OR = round(OR, 2))], grn_set + cell_type ~ method, value.var = "OR"))
cat(sprintf("\nFisher P: 最大 = %.1e (所有比较 P < 1e-30: %s)\n", max(sig_tab$fisher_p), all(sig_tab$fisher_p < 1e-30)))

GRN_LAB <- c(CT = "Cell-type eRegulons", CT_HQ = "High-quality eRegulons")
COL_GRN <- c(CT = PAL_MINT[["m3"]], CT_HQ = "#3492B2")   # 浅薄荷 / 亮青蓝 (第 2 轮: 钢蓝 → 更鲜亮的 #3492B2, 取自同一参考色板) (CT_HQ 是 CT 的子集; 2026-11-04 薄荷色板)
sig_tab[, `:=`(method  = factor(method, levels = METHOD_KEY),
               grn_set = factor(grn_set, levels = names(GRN_LAB), labels = GRN_LAB),
               ct      = factor(CELL_LAB[cell_type], levels = CELL_LAB[CELLS]))]

# strip_fmt: 标签条文字格式 (Figure 3E: "CIMA GRN overlap (CD4 T cell)"; S8 只写细胞类型)
grn_plot <- function(cts, strip_fmt = "%s") ggplot(sig_tab[cell_type %in% cts], aes(method, OR, fill = grn_set)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_col(position = position_dodge(0.78), width = 0.72) +
  facet_wrap(~ ct, nrow = 1, labeller = as_labeller(function(x) sprintf(strip_fmt, x))) +
  scale_fill_manual(values = setNames(COL_GRN, GRN_LAB), name = NULL) +
  scale_x_discrete(labels = sub(" ", "\n", lab_paper(lab_key))) +   # 方法名两行水平 (Causal / SMR); 论文写法 Causal joint 等
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "Odds ratio") +   # 不写 title / subtitle (阈值、GRN 集合写在 caption)
  theme_pub() + TH_BOX +
  theme(axis.title.y = element_text(face = "bold", size = BASE3, colour = "#1F2937", angle = 90, margin = margin(r = 3)),
        axis.text.x = element_text(size = BASE3 - 1, colour = "#1F2937", lineheight = 0.9),
        panel.grid.major.x = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.985, 0.985), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA), legend.margin = margin(2, 4, 2, 4),
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = BASE3 - 1))
p_or <- grn_plot("CD4_Tn_CCR7", "CIMA GRN overlap (%s)")                       # Figure 3E: 只画 CD4 T cell (2026-11-04)
p_or_s8 <- grn_plot(c("Bn_TCL1A", "CD8_Tn_CCR7"))     # Figure S8: Bn 和 CD8 T cell

###================================================================================================
### A: Miami 图 (Plot15 构建)
MIAMI_RDS <- file.path(TMP_DIR, "Plot15_Manhattan/miami_CD4.rds")
if (!file.exists(MIAMI_RDS)) stop("缺少 ", MIAMI_RDS, ": 先运行 Scripts/Visualization/Plot15_Manhattan.R")
miami <- readRDS(MIAMI_RDS)
# 字号在这里加大 (& 作用于 Miami 的所有子图; 不影响 S7); y 轴标题已在 Plot15 中放在左侧单独一列 (两半对齐)
miami <- miami & theme(axis.text = element_text(size = BASE3 - 1.5, colour = "black"), axis.text.x = element_text(size = BASE3 - 1.5),
                       strip.text = element_text(face = "bold", size = BASE3, colour = "#1F2937", margin = margin(3, 4, 3, 4)))
p_miami <- wrap_elements(full = miami) +   # 标签条 "CD4 T cell" 在 Plot15 中加在上半轨道上 (2026-11-04 起不再用标题)
  theme(plot.tag = element_text(face = "bold", size = 18), plot.margin = margin(0, 0, 0, 0))   # 作为一个整体 panel, 上下两半不分别打 tag

###================================================================================================
### D: HiChIP 校准图 (Plot13 构建; 薄荷色板, 标签条 "HiChIP", 与 S5 panel E 同一口径)
HICHIP_RDS <- file.path(TMP_DIR, "Plot13_MROnly/hichip_fig3.rds")
if (!file.exists(HICHIP_RDS)) stop("缺少 ", HICHIP_RDS, ": 先运行 Scripts/Visualization/Plot13_MROnly_ObsExp.R")
p_hichip <- readRDS(HICHIP_RDS) + theme(plot.tag = element_text(face = "bold", size = 18), plot.margin = margin(3, 8, 3, 3))

###================================================================================================
### 拼图 (2026-11-04 第 2 轮): A (Miami, 通栏) / 下半 2×2 网格: 左列 B = GO BP、C = GO MF, 右列 D = HiChIP、E = GRN (CD4 T cell)
# 扁平 2×2 网格: 同一行两图的绘图区上下对齐、同一列左右对齐; 网格按行排 (BP, HiChIP, MF, GRN), 字母按列走 → 自定义 tag 顺序
grid4 <- wrap_plots(go_panel(CATS[1]), p_hichip, go_panel(CATS[2], legend = TRUE), p_or, ncol = 2) + plot_layout(widths = GRID_W)
fig <- (p_miami / grid4) + plot_layout(heights = HEIGHTS) + plot_annotation(tag_levels = list(c("A", "B", "D", "C", "E"))) &
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(face = "bold", size = 18, hjust = 0, vjust = 1))   # tag 画在图的左上角, 不额外占一行 (减少 panel 间空白)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_W <- 13; FIG_H <- 14
out <- file.path(OUT_DIR, "Figure_3.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H)
print(fig)
invisible(dev.off())
msg("Saved: ", out)

### Figure S8: GRN odds ratio, Bn 和 CD8 T cell (从 Figure 3 移出)
out_s8 <- file.path(OUT_DIR, "Figure_S8.pdf")
cairo_pdf(out_s8, width = 10, height = 5.2)
print(p_or_s8)
invisible(dev.off())
msg("Saved: ", out_s8)
