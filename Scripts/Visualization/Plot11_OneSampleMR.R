# Plot11 (论文当前 S1 Fig): QC_OneSampleMR — null SNP (caQTL、eQTL 都不显著) 上 caQTL z 与 eQTL z 的相关性
# Figure S1: 1 页, 2 个 panel: A = pair 内相关 rho_hat, B = 相关性检验的 P 值; 横轴为 pair 内共享 null SNP 数分组,
#            每组 3 根并排箱子 = 3 个细胞类型 (图例写各细胞类型的全基因组整体 rho 和 P); 说明 (虚线含义 / pair 数) 写在论文 caption
# 输入: Results/QC/1_OneSampleMR/<CT>_pair_corr.tsv.gz, <CT>_global_rho.tsv
# 输出: Results/Figures/Figures_20261102/Figure_S1.pdf (只出 PDF)
# 运行: Rscript Scripts/Visualization/Plot11_OneSampleMR.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages(library(patchwork))
QC_DIR <- file.path(PROJ, "Results/QC/1_OneSampleMR")

BRK <- c(2, 10, 20, 30, 40, 50, Inf)
LAB <- c("3–10", "11–20", "21–30", "31–40", "41–50", ">50")
CT_SHORT <- CELL_LAB[CELLS]   # Bn / CD4 T cell / CD8 T cell (00_common.R)
COL_CT   <- setNames(c("#7570B3", "#1B9E77", "#D95F02"), CT_SHORT)   # ColorBrewer Dark2 (色盲友好); 与 6 个方法的配色不同
COL_CT_D <- setNames(c("#4B4785", "#116B50", "#963F01"), CT_SHORT)   # 同色加深: 箱子描边 / 中位线
BASEF <- 11
FONT  <- "sans"   # cairo_pdf 嵌入字体, 支持 ρ 等符号
fmt_p <- function(p) sprintf("%.1e", p)
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

# 三个细胞类型合并读入
D <- rbindlist(lapply(CELLS, function(ct) {
  d <- fread(file.path(QC_DIR, paste0(ct, "_pair_corr.tsv.gz")))[!is.na(rho_hat)]
  d[, `:=`(ct = factor(CT_SHORT[[ct]], levels = CT_SHORT), grp = cut(snp_count, BRK, labels = LAB))]
}))
G <- rbindlist(lapply(CELLS, function(ct) fread(file.path(QC_DIR, paste0(ct, "_global_rho.tsv")))[, ct := CT_SHORT[[ct]]]))
G[, lab := sprintf("%s: ρ = %.4f, P = %s", ct, rho, fmt_p(p_value))]
LEG <- setNames(G$lab, G$ct)

# 日志: 每个 CT 的 pair 数 / SNP 数 / 全局 rho, P (写论文 caption 用); 每组中位数 (与旧图对账)
print(G[, .(ct, n_pair = n_pair_cor, n_snp, rho = signif(rho, 4), p_value = signif(p_value, 3))])
print(D[, .(n_pair = .N, median_rho = round(median(rho_hat), 4), median_p = round(median(p_value), 3)), by = .(ct, grp)][order(ct, grp)])

theme_qc <- function() theme_classic(base_size = BASEF, base_family = FONT) + theme(
  text = element_text(color = "black"), axis.text = element_text(color = "black", size = BASEF - 1),
  axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
  legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = BASEF),
  legend.key.size = unit(0.4, "cm"), plot.margin = margin(6, 12, 6, 8),
  plot.tag = element_text(face = "bold", size = BASEF + 8, family = FONT))

box_plot <- function(yvar, ylab, ref, ylim) {
  ggplot(D, aes(grp, .data[[yvar]], fill = ct, color = ct)) +
    geom_hline(yintercept = ref, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_boxplot(width = 0.72, position = position_dodge(0.8), outlier.shape = NA, linewidth = 0.45, fatten = 2) +
    scale_fill_manual(values = alpha(COL_CT, 0.55), labels = LEG) +
    scale_color_manual(values = COL_CT_D, labels = LEG) +
    coord_cartesian(ylim = ylim) +   # 只缩放视野, 不删数据 (scale limits 会在算箱线图前丢点)
    labs(x = "Number of shared null SNPs per peak–gene pair", y = ylab) +
    theme_qc()
}

rl <- max(abs(quantile(D$rho_hat, c(0.01, 0.99)))) * 1.15   # y 范围取 1%-99% 分位数 (三个 CT 合并)
p1 <- box_plot("rho_hat", expression("Correlation of z-scores (" * hat(rho) * ")"), 0, c(-rl, rl))
p2 <- box_plot("p_value", expression(italic(P) * " value of correlation test"), 0.5, c(0, 1))
fig <- (p1 | p2) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "top")

FIG_W <- 9.5; FIG_H <- 4.8
OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261102")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- file.path(OUT_DIR, "Figure_S1.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H, family = FONT)
print(fig)
invisible(dev.off())
msg("Saved: ", out)
