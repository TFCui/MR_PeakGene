# Figure 3: MR-Only gene 的 GO 富集 (DAVID Functional Annotation Chart; 3 个细胞类型的 MR-Only gene 合并为一个 list)
# gene list: Results/Validation/Valid7_GO/DAVID_input/MROnly_genes_all.txt (3 个 CT 并集, 4,509 个 gene), 在 DAVID 网站上运行
#   (OFFICIAL_GENE_SYMBOL, Homo sapiens, 默认 background; 只导出 GOTERM_BP_DIRECT / GOTERM_MF_DIRECT)
# DAVID 只导出 P < 0.1 的 term; FDR 用 DAVID 的 FDR 列 (BH, 每个类别内部)
# 2 个 panel 并排: A = GO BP, B = GO MF; 每个 panel 按 FDR 取前 10 个 term (并列按 P-Value)
#   dotplot (仿 enrichplot::dotplot): x = GeneRatio (Count / List Total), 点大小 = Count, 颜色 = FDR (图例写作 p.adj, 图中不再出现 "FDR")
# 两个 panel 共用颜色标尺 (FDR, log10 刻度; FDR < 1e-10 的统一显示为最深色, 否则 GO MF "protein binding" FDR ~1e-85 会把其余颜色压淡)
#   和点大小标尺 (Count, log10 刻度), 图例合并放右侧
# 旧版 (每个 CT 一份 DAVID 结果, 4 类别 x 3 CT, dotplot + barplot) 见 git 历史和 Figures_20261019/Fig9_GO_DAVID.pdf
# 输入: Results/Validation/Valid7_GO/DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv (DAVID 以 list 第一个 gene 命名)
# 输出: Results/Figures/Figures_20261027/Figure_3.pdf + Figure_3.svg
# 运行: Rscript Scripts/Visualization/Plot9_GO.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages(library(patchwork))

DAVID_FILE <- file.path(VAL_DIR, "Valid7_GO/DAVID_output/DAVIDChartReport_A1BG_2026-10-19.csv")
OUT_DIR    <- file.path(PROJ, "Results/Figures/Figures_20261027")
CATS       <- c("GOTERM_BP_DIRECT", "GOTERM_MF_DIRECT")
CAT_LAB    <- c(GOTERM_BP_DIRECT = "GO Biological Process", GOTERM_MF_DIRECT = "GO Molecular Function")
TOP_N      <- 10
COL_DARK   <- "#1F5A96"; COL_LIGHT <- "#CFE3F3"   # FDR 越小颜色越深

## ---- 读 DAVID 结果 ----
res <- fread(DAVID_FILE)
setnames(res, c("P-Value", "List Total", "Fold Enrichment", "Full Term"), c("pvalue", "list_total", "fold", "full_term"))
msg("DAVID 结果: 每个类别的 term 数 (只含 P < 0.1) / FDR < 0.05 / FDR < 0.1")
print(res[, .(n = .N, fdr_lt_005 = sum(FDR < 0.05), fdr_lt_01 = sum(FDR < 0.1)), by = Category])
res <- res[Category %in% CATS]

## ---- 每个 panel 前 10 个 term ----
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

## ---- 单个 panel (dotplot) ----
panel_plot <- function(cat) {
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

fig <- wrap_plots(lapply(CATS, panel_plot), nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

## ---- 输出 ----
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_W <- 13; FIG_H <- 5.5
out <- file.path(OUT_DIR, "Figure_3.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H)
print(fig)
invisible(dev.off())
msg("写出: ", out)
out <- file.path(OUT_DIR, "Figure_3.svg")
svglite::svglite(out, width = FIG_W, height = FIG_H)   # svglite: 文字保留为可编辑文本, 便于拼图
print(fig)
invisible(dev.off())
msg("写出: ", out)
