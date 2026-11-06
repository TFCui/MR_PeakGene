# Plot12 (论文当前 S4 Fig): 每个方法在 3 个细胞类型之间显著 peak–gene link 的 Venn 图
# Figure S4: 1 页, 2×3 个 panel (按 METHODS 顺序: SMR / IVW / Egger; Joint / Multiomics / Pseudobulk), 每个 panel 一个三圆 Venn
#   2026-11-04 起只用 3 个 CT 都检验过的候选 pair (14,890 个; 避免“另一个 CT 没检验”造成的假 unique)
#   call set = 某 CT 中该方法 p.adj (gene 内 BH, 各 CT 自己的 MR 表) < 0.1 且属于共有 pair 的 link; 用 peak_id + gene_id 匹配
#   区域内: count + (占三者并集的 %; 并集数只写入 tsv, 2026-11-04 起图中不显示); 圆外标签: 该 CT 在共有 pair 中的显著数 / 该 CT 全部显著数 (占比)
#   细胞类型用页底共用图例; 圆大小固定 (不按比例); 低饱和配色 (雾蓝 / 鼠尾草绿 / 杏色)
# 输入: Results/MR/<CT>_AllPairs_results.tsv.gz
# 输出: Results/Figures/Figures_20261102/Figure_S4.pdf (只出 PDF) + Figure_S4_counts.tsv (长表)
# 运行: Rscript Scripts/Visualization/Plot12_Venn.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages({ library(patchwork); library(ggforce) })

FDR_CUT  <- 0.1
CT_SHORT <- CELL_LAB[CELLS]   # Bn / CD4 T cell / CD8 T cell
COL_CT   <- setNames(c("#9DBEDF", "#A9D3A2", "#F3C29B"), CT_SHORT)   # 低饱和、清新的雾蓝 / 鼠尾草绿 / 杏色 (2026-11-04; 与 Figure S1 的 Dark2 不同), 保证圆内黑字清楚
COL_CT_D <- setNames(c("#6F95BD", "#74A86B", "#D49564"), CT_SHORT)   # 同色加深: 圆边
BASEF <- 11
FONT  <- "sans"
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

res_list <- setNames(lapply(CELLS, read_mr), CELLS)
keys <- lapply(res_list, function(d) paste(d$peak_id, d$gene_id, sep = "__"))
common <- Reduce(intersect, keys); n_common <- length(common)
stopifnot(n_common == 14890)
cat("\n===== 各 CT 检验的 pair 数; 3 个 CT 都检验过的 pair 数 =", n_common, "\n")
print(data.table(cell_type = CELLS, n_tested = lengths(keys)))

###------------------------------------------------------------------------------------------------
### 1. 每个方法 (只看共有 pair): 7 个区域的 count / 占并集 %; 每个 CT: 全部显著数, 共有 pair 中的显著数和占比, % unique
REGIONS <- c("100", "010", "001", "110", "101", "011", "111")   # 位 = Bn, CD4, CD8
venn_counts <- function(mm) {
  col <- paste0(METHOD_KEY[[mm]], "_fdr")
  S_all <- lapply(CELLS, function(ct) { d <- res_list[[ct]]; x <- d[[col]]; unique(keys[[ct]][!is.na(x) & x < FDR_CUT]) })
  S <- lapply(S_all, function(z) z[z %chin% common])   # 只保留共有 pair
  U <- Reduce(union, S)
  code <- paste0(as.integer(U %in% S[[1]]), as.integer(U %in% S[[2]]), as.integer(U %in% S[[3]]))
  reg <- data.table(method = mm, region = REGIONS, n = as.integer(table(factor(code, levels = REGIONS))))
  reg[, `:=`(n_union = length(U), pct = 100 * n / length(U))]
  stopifnot(sum(reg$n) == length(U))
  ct <- data.table(method = mm, cell_type = CELLS, n_sig_all = lengths(S_all), n_sig = lengths(S),
                   n_unique = reg$n[match(c("100", "010", "001"), reg$region)])
  ct[, `:=`(pct_common = 100 * n_sig / n_sig_all, pct_unique = 100 * n_unique / n_sig)]
  list(reg = reg, ct = ct)
}
VC  <- setNames(lapply(METHODS, venn_counts), METHODS)
REG <- rbindlist(lapply(VC, `[[`, "reg"))
CTN <- rbindlist(lapply(VC, `[[`, "ct"))

cat("\n===== 7 个区域 (位 = Bn, CD4 T cell, CD8 T cell)\n")
print(dcast(REG, method ~ region, value.var = "n")[match(METHODS, method), c("method", REGIONS), with = FALSE])
cat("\n===== 每个 CT: 全部显著数 / 共有 pair 中的显著数和占比 / 只在本 CT 显著的数和 % (共有 pair 口径)\n")
print(CTN[, .(method, cell_type, n_sig_all, n_sig_common = n_sig, pct_common = round(pct_common, 1), n_unique, pct_unique = round(pct_unique, 1))])

###------------------------------------------------------------------------------------------------
### 2. 画图 (固定几何: 半径 1 的 3 个圆)
CIRC <- data.table(ct = factor(CT_SHORT, levels = CT_SHORT), x0 = c(-0.55, 0.55, 0), y0 = c(0.32, 0.32, -0.63), r = 1)
POS  <- data.table(region = REGIONS,
                   x = c(-0.98, 0.98, 0, 0, -0.60, 0.60, 0),
                   y = c(0.55, 0.55, -1.12, 0.78, -0.30, -0.30, 0.04))
LABP <- data.table(ct = CT_SHORT, x = c(-1.05, 1.05, 0), y = c(1.40, 1.40, -1.70), vjust = c(0, 0, 1))   # 圆外标签: Bn 左上 / CD4 右上 / CD8 正下

venn_panel <- function(mm) {
  reg <- merge(REG[method == mm], POS, by = "region")
  reg[, lab := sprintf("%s\n(%.1f%%)", fmt_n(n), pct)]
  ctl <- merge(LABP, CTN[method == mm][, ct := CT_SHORT[cell_type]], by = "ct")
  ctl[, lab := sprintf("%s / %s\n(%.1f%%)", fmt_n(n_sig), fmt_n(n_sig_all), pct_common)]   # 两行, 避免上方两个标签相撞
  ggplot() +
    geom_circle(data = CIRC, aes(x0 = x0, y0 = y0, r = r, fill = ct, color = ct), alpha = 0.45, linewidth = 0.7) +
    geom_text(data = reg, aes(x, y, label = lab), size = 3.5, lineheight = 0.88, fontface = "bold", family = FONT, color = "grey10") +
    geom_text(data = ctl, aes(x, y, label = lab, color = ct, vjust = vjust), size = 3.2, lineheight = 0.9, fontface = "bold", family = FONT) +
    scale_fill_manual(values = COL_CT, name = NULL, guide = guide_legend(override.aes = list(alpha = 0.6))) +
    scale_color_manual(values = COL_CT_D, guide = "none") +
    coord_fixed(xlim = c(-1.75, 1.75), ylim = c(-2.15, 1.85), expand = FALSE, clip = "off") +
    labs(title = mm) +
    theme_void(base_size = BASEF, base_family = FONT) +
    theme(plot.title = element_text(face = "bold", size = BASEF + 1, hjust = 0.5, color = pal_method[[mm]], margin = margin(b = 4)),
          plot.margin = margin(4, 6, 4, 6),
          plot.tag = element_text(face = "bold", size = BASEF + 8, family = FONT))
}

fig <- wrap_plots(lapply(METHODS, venn_panel), nrow = 2, byrow = TRUE) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.text = element_text(size = BASEF, family = FONT),
        legend.key.size = unit(0.45, "cm"), legend.spacing.x = unit(0.3, "cm"))

FIG_W <- 10; FIG_H <- 8
OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261102")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- file.path(OUT_DIR, "Figure_S4.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H, family = FONT)
print(fig)
invisible(dev.off())
msg("Saved: ", out)

### 3. 长表 (写 caption / 正文用)
tab <- rbind(
  REG[, .(method, item = paste0("region_", region), n, denom = n_union, pct = round(pct, 2))],
  CTN[, .(method, item = paste0("common_", cell_type), n = n_sig, denom = n_sig_all, pct = round(pct_common, 2))],
  CTN[, .(method, item = paste0("unique_", cell_type), n = n_unique, denom = n_sig, pct = round(pct_unique, 2))])
tab[, method := factor(method, levels = METHODS)]
setorder(tab, method)
fwrite(tab, file.path(OUT_DIR, "Figure_S4_counts.tsv"), sep = "\t")
msg("Saved: ", file.path(OUT_DIR, "Figure_S4_counts.tsv"))
