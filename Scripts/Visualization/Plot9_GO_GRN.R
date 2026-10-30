# Figure 3: A / B = MR-Only gene 的 GO 富集 (DAVID); C = 各方法显著 pair 与 CIMA SCENIC+ GRN 的重叠 (odds ratio)
# 由 Plot9_GO.R (原 Figure 3) 与 Plot_Sup2_GRN.R (原 Figure S5) 合并 (2026-10-29); S5 的 panel A (重叠 pair 数) 不再画, S5 Fig 已从论文删除
#
# --- A / B: GO (DAVID Functional Annotation Chart; 3 个细胞类型的 MR-Only gene 合并为一个 list) ---
# gene list: Results/Validation/Valid7_GO/DAVID_input/MROnly_genes_all.txt (3 个 CT 并集, 4,509 个 gene), 在 DAVID 网站上运行
#   (OFFICIAL_GENE_SYMBOL, Homo sapiens, 默认 background; 只导出 GOTERM_BP_DIRECT / GOTERM_MF_DIRECT)
# DAVID 只导出 P < 0.1 的 term; FDR 用 DAVID 的 FDR 列 (BH, 每个类别内部)
# A = GO BP, B = GO MF; 每个 panel 按 FDR 取前 10 个 term (并列按 P-Value)
#   dotplot (仿 enrichplot::dotplot): x = GeneRatio (Count / List Total), 点大小 = Count, 颜色 = FDR (图例写作 p.adj, 图中不再出现 "FDR")
# 两个 panel 共用颜色标尺 (FDR, log10 刻度; FDR < 1e-10 的统一显示为最深色, 否则 GO MF "protein binding" FDR ~1e-85 会把其余颜色压淡)
#   和点大小标尺 (Count, log10 刻度), 图例合并放右侧
# 输入: Results/Validation/Valid7_GO/DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv (DAVID 以 list 第一个 gene 命名)
#
# --- C: GRN odds ratio ---
# 显著 (p.adj < 0.05, 与 Valid5_GRN 一致, 比正文的 0.1 严格) pair 在 GRN 中的 odds ratio (GRN vs 非 GRN, 2x2 Fisher; 每个细胞类型一个小面板,
#   CT 与 CT_HQ 两种 eRegulon 集合并排); CT = 该细胞类型中 active (AUC > 1) 的 eRegulon; CT_HQ = CT 与对应 lineage 高质量 eRegulon 的交集
#   (流程见 Scripts/Pipeline.md 3.5); 统计量已在 Valid5_GRN 中算好, 这里只画图; 说明 (阈值 / GRN 集合 / pair 数) 写在论文 caption
# 输入: Results/Validation/Valid5_GRN/GRN_sig_counts_by_method.tsv, GRN_summary.tsv
#
# 输出: Results/Figures/Figures_20261102/Figure_3.pdf (只出 PDF)
# 运行: Rscript Scripts/Visualization/Plot9_GO_GRN.R   (小表, 登录节点即可)

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages({ library(patchwork); library(scales) })

OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261102")

theme_pub <- function(base = 10) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      plot.title         = element_text(face = "bold", size = base + 1.5, colour = "#111827"),
      plot.subtitle      = element_text(size = base - 1, colour = "#6B7280", margin = margin(b = 4)),
      axis.line          = element_line(colour = "#9CA3AF", linewidth = 0.3),
      axis.ticks         = element_line(colour = "#9CA3AF", linewidth = 0.3),
      axis.text.x        = element_text(size = base - 1.5, colour = "#4B5563"),
      axis.text.y        = element_text(size = base - 1, colour = "#1F2937", lineheight = 0.85),
      axis.title.x       = element_text(face = "bold", size = base - 0.5, colour = "#1F2937", margin = margin(t = 4)),
      axis.title.y       = element_blank(),
      panel.grid.major.x = element_line(colour = "#EEF0F3", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "#F5F6F8", linewidth = 0.3),
      legend.title       = element_text(face = "bold", size = base - 0.5),
      legend.text        = element_text(size = base - 1),
      plot.title.position = "plot",
      plot.tag           = element_text(face = "bold", size = base + 6),
      plot.margin        = margin(8, 10, 8, 8)
    )
}

###================================================================================================
### A / B: GO
DAVID_FILE <- file.path(VAL_DIR, "Valid7_GO/DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv")
CATS       <- c("GOTERM_BP_DIRECT", "GOTERM_MF_DIRECT")
CAT_LAB    <- c(GOTERM_BP_DIRECT = "GO Biological Process", GOTERM_MF_DIRECT = "GO Molecular Function")
TOP_N      <- 10
COL_DARK   <- "#1F5A96"; COL_LIGHT <- "#CFE3F3"   # FDR 越小颜色越深

res <- fread(DAVID_FILE)
setnames(res, c("P-Value", "List Total", "Fold Enrichment", "Full Term"), c("pvalue", "list_total", "fold", "full_term"))
msg("DAVID 结果: 每个类别的 term 数 (只含 P < 0.1) / FDR < 0.05 / FDR < 0.1")
print(res[, .(n = .N, fdr_lt_005 = sum(FDR < 0.05), fdr_lt_01 = sum(FDR < 0.1)), by = Category])
res <- res[Category %in% CATS]

## 每个 panel 前 10 个 term
top <- res[order(match(Category, CATS), FDR, pvalue)][, head(.SD, TOP_N), by = Category]
stopifnot(uniqueN(top$Category) == length(CATS))
top[, gene_ratio := Count / list_total]
top[, term := stringr::str_wrap(Term, 40)]
print(top[, .(Category, Term, Count, list_total, gene_ratio = round(gene_ratio, 3), fold, FDR)])

# 共用标尺的范围; FDR 颜色下限截在 1e-10
FDR_FLOOR  <- 1e-10
lim_fdr    <- c(FDR_FLOOR, max(top$FDR))
lim_count  <- range(top$Count)
fdr_breaks <- 10^c(-10, -7, -4, -1)
fdr_breaks <- fdr_breaks[fdr_breaks <= lim_fdr[2]]
fdr_labels <- ifelse(fdr_breaks == FDR_FLOOR, "≤ 1e-10", format(fdr_breaks, scientific = TRUE))

go_panel <- function(cat) {
  d <- copy(top[Category == cat])
  d[, ykey := factor(full_term, levels = rev(full_term))]   # FDR 最小的在最上面
  ggplot(d, aes(y = ykey)) +
    geom_point(aes(x = gene_ratio, size = Count, colour = FDR)) +
    scale_size_continuous(name = "Count", limits = lim_count, range = c(2, 8), trans = "log10", breaks = c(50, 100, 300, 1000, 3000)) +
    scale_colour_gradient(name = "p.adj", low = COL_DARK, high = COL_LIGHT, trans = "log10", limits = lim_fdr,
                          oob = scales::squish, breaks = fdr_breaks, labels = fdr_labels) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    scale_y_discrete(labels = setNames(d$term, d$full_term)) +
    labs(x = "GeneRatio", title = CAT_LAB[[cat]], subtitle = sprintf("n = %s genes", format(d$list_total[1], big.mark = ","))) +
    theme_pub()
}

###================================================================================================
### C: GRN odds ratio
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
COL_GRN <- c(CT = "#A6CEE3", CT_HQ = "#1F78B4")      # ColorBrewer Paired 浅蓝 / 亮蓝 (同色系浅 / 深: CT_HQ 是 CT 的子集; 不与 6 个方法的配色冲突)
sig_tab[, `:=`(method  = factor(method, levels = METHOD_KEY),
               grn_set = factor(grn_set, levels = names(GRN_LAB), labels = GRN_LAB),
               ct      = factor(CELL_LAB[cell_type], levels = CELL_LAB[CELLS]))]

p_or <- ggplot(sig_tab, aes(method, OR, fill = grn_set)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_col(position = position_dodge(0.78), width = 0.72) +
  facet_wrap(~ ct, nrow = 1) +
  scale_fill_manual(values = setNames(COL_GRN, GRN_LAB), name = NULL) +
  scale_x_discrete(labels = lab_key) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "Odds ratio", title = "Overlap with the CIMA regulatory network",
       subtitle = sprintf("Significant pairs (p.adj < %g) in the GRN vs. non-GRN pairs", GRN_FDR)) +
  theme_pub() +
  theme(axis.title.y = element_text(face = "bold", size = 9.5, colour = "#1F2937", angle = 90, margin = margin(r = 4)),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 9, colour = "#1F2937"),
        panel.grid.major.x = element_blank(), panel.border = element_rect(colour = "#9CA3AF", fill = NA, linewidth = 0.3),
        strip.background = element_rect(fill = "#F5F6F8", colour = "#9CA3AF", linewidth = 0.3), strip.text = element_text(face = "bold", size = 10),
        legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
        legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 9))

###================================================================================================
### 拼图: 上行 A / B (GO, 图例收在右侧), 下行 C (GRN, 通栏)
top_row <- wrap_plots(lapply(CATS, go_panel), nrow = 1) + plot_layout(guides = "collect")
fig <- (top_row / free(p_or)) + plot_layout(heights = c(1, 1.1)) + plot_annotation(tag_levels = "A")   # free(): C 不与上行的轴对齐, 否则左侧留出 A 的 term 标签宽度的空白

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_W <- 13; FIG_H <- 11
out <- file.path(OUT_DIR, "Figure_3.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H)
print(fig)
invisible(dev.off())
msg("Saved: ", out)
