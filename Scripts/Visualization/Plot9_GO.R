# Fig9: MR-Only gene 的 GO BP GSEA (Valid7_GO; exponent = 1)
# 3 页, 每组 gene set (Bn / CD4 / CD8 / Combined) 取前 5 个 term:
#   第 1 页: 按 GSEA pvalue 取前 5, 4 组的 -log10(P), 一列 4 行
#   第 2 页: 第 1 页中 Across Cell Types (Combined) 的 5 个 term 的 -log10(FDR), 虚线 = FDR 0.1
#   第 3 页: 按经验 P 取前 5 (并列时按 GSEA pvalue; 前 5 基本都并列在下限 1/5001), -log10(emp_p), 点线 = 5000 次随机能达到的上限
# 输入: Results/Validation/Valid7_GO/GO_GSEA_MROnly.tsv.gz
# 输出: Results/Figures/Figures_20261019/Fig9_GO_GSEA.pdf
# 运行: Rscript Scripts/Visualization/Plot9_GO.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages(library(patchwork))

res <- fread(file.path(VAL_DIR, "Valid7_GO/GO_GSEA_MROnly.tsv.gz"))

GSETS    <- c(CELLS, "Combined")
SET_LAB  <- c(CELL_DESC, Combined = "Across Cell Types")
TOP_N    <- 5
N_PERM   <- unique(res$n_perm); stopifnot(length(N_PERM) == 1)
COL_FILL <- "#6FA8DC"; COL_LINE <- "#4F8CC4"   # 与 Fig5 的 MR-Only 同色

theme_pub <- function(base = 9) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      plot.title         = element_text(face = "bold", size = base + 3, colour = "#111827"),
      plot.subtitle      = element_text(size = base - 0.5, colour = "#6B7280", margin = margin(b = 6)),
      strip.background   = element_rect(fill = "#F5F6F8", colour = NA),
      strip.text         = element_text(face = "bold", size = base, colour = "#1F2937", hjust = 0, margin = margin(t = 3, b = 3, l = 4)),
      axis.line          = element_line(colour = "#9CA3AF", linewidth = 0.3),
      axis.ticks         = element_line(colour = "#9CA3AF", linewidth = 0.3),
      axis.text.x        = element_text(size = base - 1.5, colour = "#4B5563"),
      axis.text.y        = element_text(size = base - 1, colour = "#1F2937", lineheight = 0.85),
      axis.title.x       = element_text(face = "bold", size = base, colour = "#1F2937", margin = margin(t = 5)),
      axis.title.y       = element_blank(),
      panel.grid.major.x = element_line(colour = "#EEF0F3", linewidth = 0.3),
      panel.spacing.y    = unit(0.9, "lines"),
      plot.title.position = "plot",
      plot.margin        = margin(10, 14, 8, 10)
    )
}

## ---- 每组前 5 个 term (ord = 排序列, 依次作为并列时的次序) ----
pick_top <- function(ord) {
  d <- res[gene_set %in% GSETS][order(match(gene_set, GSETS), get(ord[1]), get(ord[length(ord)]))][, head(.SD, TOP_N), by = gene_set]
  stopifnot(nrow(d) == TOP_N * length(GSETS))
  lev <- sprintf("%s  (n = %s MR-Only genes)", SET_LAB[GSETS], format(d[match(GSETS, gene_set), n_genes], big.mark = ","))
  d[, panel := factor(sprintf("%s  (n = %s MR-Only genes)", SET_LAB[gene_set], format(n_genes, big.mark = ",")), levels = lev)]
  d[, term := stringr::str_wrap(Description, 45)]
  # y 轴: 面板内按排序从上到下; 同一 term 可能出现在多个面板, 用 gene_set 前缀区分, 标签只显示 term
  d[, ykey := factor(paste(gene_set, ID, sep = "___"), levels = rev(paste(gene_set, ID, sep = "___")))]
  d
}
show_top <- function(d, title) {
  msg(title)
  print(d[, .(gene_set, Description = substr(Description, 1, 60), setSize, NES = round(NES, 2), pvalue = signif(pvalue, 2),
              p.adjust = signif(p.adjust, 2), emp_p = signif(emp_p, 2))])
}
top     <- pick_top("pvalue")             # 第 1、2 页
top_emp <- pick_top(c("emp_p", "pvalue")) # 第 3 页
show_top(top,     "第 1 页: 每组 GSEA P 最小的 5 个 term")
show_top(top_emp, "第 3 页: 每组经验 P 最小的 5 个 term (并列时按 GSEA P)")

## ---- 画图函数: 水平条形图 ----
bar_plot <- function(dt, xcol, xlab, vline, vline_lty = "dashed", ceiling = NULL) {
  dt <- copy(dt)[, x := -log10(get(xcol))]
  p <- ggplot(dt, aes(x = x, y = ykey)) +
    geom_col(fill = COL_FILL, colour = COL_LINE, linewidth = 0.3, width = 0.68) +
    geom_vline(xintercept = vline, linetype = vline_lty, colour = "#9CA3AF", linewidth = 0.4) +
    scale_y_discrete(labels = setNames(dt$term, as.character(dt$ykey))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = xlab) +
    theme_pub()
  if (!is.null(ceiling)) p <- p + geom_vline(xintercept = ceiling, linetype = "dotted", colour = "#6B7280", linewidth = 0.5)
  if (uniqueN(dt$panel) > 1) p <- p + facet_wrap(~ panel, ncol = 1, scales = "free")
  p
}

p1 <- bar_plot(top, "pvalue", expression(bold(-log[10]~italic(P))), -log10(0.05)) +
  labs(title = "GO BP enrichment of MR-Only genes (GSEA)",
       subtitle = "Top 5 terms per gene set; dashed line: P = 0.05")

cmb <- top[gene_set == "Combined"]
p2 <- bar_plot(cmb, "p.adjust", expression(bold(-log[10]~FDR)), -log10(0.1)) +
  labs(title = "GO BP enrichment of MR-Only genes: Across Cell Types",
       subtitle = sprintf("n = %s MR-Only genes; BH FDR\ndashed line: FDR = 0.1", format(cmb$n_genes[1], big.mark = ",")))

p3 <- bar_plot(top_emp, "emp_p", expression(bold(-log[10]~italic(P)[empirical])), -log10(0.05), ceiling = -log10(1 / (N_PERM + 1))) +
  labs(title = "GO BP enrichment of MR-Only genes: empirical P",
       subtitle = sprintf("Top 5 terms per gene set by empirical P (ties broken by GSEA P)\n%s random gene sets of the same size drawn from tested genes\ndashed line: P = 0.05; dotted line: minimum attainable P = 1/%s",
                          format(N_PERM, big.mark = ","), format(N_PERM + 1, big.mark = ",")))

## ---- 输出 ----
out <- file.path(FIG_DIR, "Fig9_GO_GSEA.pdf")
cairo_pdf(out, width = 7.5, height = 10, onefile = TRUE)
print(p1)
print(p2 / plot_spacer() + plot_layout(heights = c(1, 1.7)))
print(p3)
invisible(dev.off())
msg("写出: ", out)
