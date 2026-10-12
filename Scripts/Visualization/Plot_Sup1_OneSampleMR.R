# Sup1: QC_OneSampleMR — null SNP (caQTL、eQTL 都不显著) 上 caQTL z 与 eQTL z 的相关性
# 一个 PDF, 3 页 (每页一个细胞类型); 每页 2 个箱线图: 左 = rho_hat, 右 = 相关性检验的 P 值; 横轴为 pair 内共享 SNP 数分组
# subtitle: 全基因组整体 rho 和 P
# 输入: Results/QC/<CT>_pair_corr.tsv.gz, <CT>_global_rho.tsv
# 输出: Results/Figures/Sup1_OneSampleMR.pdf
# 运行: Rscript Scripts/Visualization/Plot_Sup1_OneSampleMR.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages(library(gridExtra))
QC_DIR <- file.path(PROJ, "Results/QC")

BRK <- c(2, 10, 20, 30, 40, 50, Inf)
LAB <- c("3–10", "11–20", "21–30", "31–40", "41–50", ">50")
COL_GRP <- colorRampPalette(c("#BFD7EA", "#3C5488"))(length(LAB))
BASEF <- 9
FONT  <- "sans"   # cairo_pdf 嵌入字体, 支持 ρ 等符号
fmt_p <- function(p) sprintf("%.1e", p)
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

theme_qc <- function() theme_classic(base_size = BASEF, base_family = FONT) + theme(
  text = element_text(color = "black"), axis.text = element_text(color = "black", size = BASEF - 1),
  axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
  plot.caption = element_text(size = BASEF - 2, color = "grey40", hjust = 0, margin = margin(t = 6)),
  legend.position = "none", plot.margin = margin(6, 12, 6, 8))

box_plot <- function(d, yvar, ylab, ref, cap, ylim) {
  nn <- d[, .(n = .N), by = grp]
  ggplot(d, aes(grp, .data[[yvar]], fill = grp)) +
    geom_hline(yintercept = ref, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.4, color = COL_INK, fatten = 2.2) +
    geom_text(data = nn, aes(grp, ylim[2], label = paste0("n = ", fmt_n(n))), inherit.aes = FALSE,
              size = 2.3, vjust = 1, color = "grey25", family = FONT) +
    scale_fill_manual(values = COL_GRP) +
    coord_cartesian(ylim = ylim) +   # 只缩放视野, 不删数据 (scale limits 会在算箱线图前丢点)
    labs(caption = cap, x = "Number of shared null SNPs per peak–gene pair", y = ylab) +
    theme_qc()
}

make_page <- function(ct) {
  d <- fread(file.path(QC_DIR, paste0(ct, "_pair_corr.tsv.gz")))[!is.na(rho_hat)]
  g <- fread(file.path(QC_DIR, paste0(ct, "_global_rho.tsv")))
  d[, grp := cut(snp_count, BRK, labels = LAB)]
  ttl <- sprintf("%s (%s)", ct, CELL_DESC[[ct]])
  sub <- sprintf("Global correlation ρ = %.4f, P = %s  (%s null SNPs, %s peak–gene pairs)", g$rho, fmt_p(g$p_value), fmt_n(g$n_snp), fmt_n(nrow(d)))
  rl <- max(abs(quantile(d$rho_hat, c(0.01, 0.99)))) * 1.15   # y 范围取 1%-99% 分位数
  p1 <- box_plot(d, "rho_hat", expression("Correlation of z-scores (" * hat(rho) * ")"), 0, "Dashed line: expected under the null (ρ̂ = 0)", c(-rl, rl))
  p2 <- box_plot(d, "p_value", expression(italic(P) * " value of correlation test"), 0.5, "Dashed line: expected median under the null (P = 0.5)", c(0, 1))
  hdr <- grid::gTree(children = grid::gList(
    grid::textGrob(ttl, x = 0.015, y = 0.72, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = BASEF + 3, fontfamily = FONT)),
    grid::textGrob(sub, x = 0.015, y = 0.25, hjust = 0, gp = grid::gpar(fontsize = BASEF, col = "grey30", fontfamily = FONT))))
  grid.arrange(hdr, arrangeGrob(p1, p2, ncol = 2), nrow = 2, heights = c(0.13, 0.87))
}

out <- file.path(FIG_DIR, "Sup1_OneSampleMR.pdf")
cairo_pdf(out, width = 11, height = 4.8, onefile = TRUE, family = FONT)
for (ct in CELLS) make_page(ct)
invisible(dev.off())
msg("Saved: ", out)
