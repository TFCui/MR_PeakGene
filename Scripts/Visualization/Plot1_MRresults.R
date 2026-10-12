# Fig1: MR 总体结果 (FDR < 0.1)
#   page 1   每个方法的显著 pair 数
#   page 2-4 每个细胞类型一张 UpSet (前 25 个组合)
#   page 5   六方法全显著 pairs 的效应方向一致性
# 另外: 在新 FDR 下检查 24 x 3 个 example pairs 是否仍满足各自类别 (结果只打印到日志, 并写 tmp/example_pairs_check.tsv)
# 移植自 20261006_paper/Scripts/Plot1_MRresults.sh
# 输入: Results/MR/<CT>_AllPairs_results.tsv.gz; (对照) 旧 MR_results_formal/<CT>_AllPairs_results.rds
# 输出: Results/Figures/Fig1_MR_overview.pdf
# 运行: Rscript Scripts/Visualization/Plot1_MRresults.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/01_example_pairs.R")
suppressPackageStartupMessages({ library(patchwork); library(ComplexHeatmap); library(grid); library(scales) })

FDR_CUT <- 0.1
res_list <- setNames(lapply(CELLS, read_mr), CELLS)
print(data.table(cell_type = CELLS, n_pairs = sapply(res_list, nrow)))

###------------------------------------------------------------------------------------------------
### 0. Example pairs: 新 FDR 下是否仍满足类别 (旧 FDR 作对照)
OLD_MR_DIR <- file.path(BASE, "Results/RunMR/Materials_runMR/20260430_LinkPeakGene_more/MRresult_20260911/MR_results_formal")
ex_chk <- rbindlist(lapply(CELLS, function(ct) {
  new <- check_examples(EXAMPLE_PAIRS[[ct]], res_list[[ct]])
  old_mr <- readRDS(file.path(OLD_MR_DIR, paste0(ct, "_AllPairs_results.rds")))[!is.na(ivw_sigma2) & ivw_sigma2 <= 10]
  old <- check_examples(EXAMPLE_PAIRS[[ct]], as.data.table(old_mr))
  new[, `:=`(cell_type = ct, pass_old = old$pass)]
}))
setcolorder(ex_chk, c("cell_type", "ord", "cat_id", "gene_id", "peak_id", "in_table", "pass_old", "pass"))
cat("\n===== Example pairs: pass (new FDR) vs pass_old (paper FDR)\n")
print(ex_chk[, .(n = .N, in_table = sum(in_table), pass_old = sum(pass_old, na.rm = TRUE), pass_new = sum(pass, na.rm = TRUE)),
             by = .(cell_type, cat_id)])
cat("\n----- pairs that changed or fail:\n")
print(ex_chk[!(pass %in% TRUE) | !(pass_old %in% TRUE) | pass != pass_old], digits = 3)
fwrite(ex_chk, file.path(TMP_DIR, "example_pairs_check.tsv"), sep = "\t")

###------------------------------------------------------------------------------------------------
### 1. 显著矩阵
sig_mat <- function(dt) {
  m <- sapply(METHODS, function(mm) { x <- dt[[paste0(METHOD_KEY[[mm]], "_fdr")]]; !is.na(x) & x < FDR_CUT })
  rownames(m) <- paste(dt$peak_id, dt$gene_id, sep = "__"); m
}

### Page 1: 每个方法的显著 pairs 数
cnt <- rbindlist(lapply(CELLS, function(ct) { m <- sig_mat(res_list[[ct]]); data.table(cell_type = ct, n_total = nrow(m), method = METHODS, n_sig = colSums(m)) }))
cnt[, method := factor(method, levels = METHODS)][, pct := n_sig / n_total]
print(dcast(cnt, cell_type + n_total ~ method, value.var = "n_sig"))

plot_count <- function(ct) {
  d <- cnt[cell_type == ct]
  ggplot(d, aes(method, n_sig, fill = method)) +
    geom_col(width = 0.68, color = NA) +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)", comma(n_sig), 100 * pct)), vjust = -0.25, size = 2.7, lineheight = 0.9) +
    scale_fill_manual(values = pal_method, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = comma) +
    labs(title = sprintf("%s  (n = %s pairs)", ct, comma(d$n_total[1])), x = NULL, y = sprintf("Significant pairs (FDR < %g)", FDR_CUT)) +
    theme_paper() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
}
p1 <- wrap_plots(lapply(CELLS, plot_count), nrow = 1) + plot_annotation(tag_levels = "a")

### Page 2-4: UpSet (单方法特有的组合用该方法颜色, 其余深灰)
plot_upset <- function(ct, n_top = 25) {
  m  <- sig_mat(res_list[[ct]]); m <- m[rowSums(m) > 0, , drop = FALSE]
  cm <- make_comb_mat(m, mode = "distinct")
  cm <- cm[rank(-comb_size(cm), ties.method = "first") <= n_top]
  cs <- comb_size(cm); cd <- comb_degree(cm)
  comb_cols <- unname(ifelse(cd == 1, pal_method[set_name(cm)][regexpr("1", comb_name(cm))], COL_INK))
  ht <- UpSet(cm, set_order = intersect(METHODS, set_name(cm)), comb_order = order(cs, decreasing = TRUE),
    pt_size = unit(2.6, "mm"), lwd = 1.6, comb_col = comb_cols, bg_col = c("#F3F5F8", "#FFFFFF"), bg_pt_col = "#D9DEE3",
    row_names_gp = gpar(fontsize = 8, fontface = "bold", col = pal_method[set_name(cm)]),
    top_annotation = upset_top_annotation(cm, add_numbers = TRUE, numbers_gp = gpar(fontsize = 7), gp = gpar(fill = comb_cols, col = NA),
      height = unit(5, "cm"), annotation_name_gp = gpar(fontsize = 9), axis_param = list(gp = gpar(fontsize = 8))),
    right_annotation = upset_right_annotation(cm, add_numbers = TRUE, numbers_gp = gpar(fontsize = 7), gp = gpar(fill = pal_method[set_name(cm)], col = NA),
      width = unit(3.2, "cm"), annotation_name_gp = gpar(fontsize = 9), axis_param = list(gp = gpar(fontsize = 8))))
  draw(ht, column_title = sprintf("%s: significant pairs across methods (FDR < %g; %s pairs significant in >= 1 method)", ct, FDR_CUT, comma(nrow(m))),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"), padding = unit(c(2, 12, 2, 2), "mm"))   # 左侧留白, 防止长方法名被截
  invisible(NULL)
}

### Page 5: 六方法全显著 pairs 的方向一致性
dir_dt <- rbindlist(lapply(CELLS, function(ct) {
  dt <- res_list[[ct]]; d <- dt[rowSums(sig_mat(dt)) == length(METHODS)]
  if (!nrow(d)) return(NULL)
  S5 <- sign(as.matrix(d[, .(ivw_beta, egger_beta, smr_b, wnn_rho, pseudo_rho)]))
  d[, .(cell_type = ct, gene_id, peak_id,
        five = fifelse(abs(rowSums(S5)) == 5, "Concordant", "Discordant"),
        mr3  = fifelse(abs(rowSums(S5[, 1:3, drop = FALSE])) == 3, "Concordant", "Discordant"),
        obs  = fcase(sign(ivw_beta) == sign(wnn_rho) & sign(ivw_beta) == sign(pseudo_rho), "Concordant with both",
                     sign(ivw_beta) != sign(wnn_rho) & sign(ivw_beta) != sign(pseudo_rho), "Discordant with both",
                     sign(ivw_beta) != sign(wnn_rho), "Discordant with Marginal Multiomics", default = "Discordant with Marginal Pseudobulk"))]
}))
print(dir_dt[, .(n_all6 = .N, five_disc = sum(five == "Discordant"), mr3_disc = sum(mr3 == "Discordant"), obs_disc = sum(obs != "Concordant with both")), by = cell_type])

pal_dir2 <- c(Concordant = COL_CONC, Discordant = COL_DISC)
pal_dir4 <- c(`Concordant with both` = COL_CONC, `Discordant with Marginal Multiomics` = pal_method[["Marginal Multiomics"]],
              `Discordant with Marginal Pseudobulk` = pal_method[["Marginal Pseudobulk"]], `Discordant with both` = COL_DISC)

plot_dir <- function(col, title, pal, nrow_leg = 1) {
  d <- dir_dt[, .N, by = c("cell_type", col)]; setnames(d, col, "status")
  d[, status := factor(status, levels = names(pal))][, pct := N / sum(N), by = cell_type]
  d[, txt_col := fifelse(status %in% names(pal)[1], "#333333", "white")]
  tot <- d[, .(N = sum(N)), by = cell_type]
  ggplot(d, aes(cell_type, N, fill = status)) +
    geom_col(width = 0.62, color = "white", linewidth = 0.3, position = position_stack(reverse = TRUE)) +
    geom_text(data = d[pct >= 0.05], aes(label = sprintf("%.0f%%", 100 * pct), color = txt_col), position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.8, fontface = "bold") +
    geom_text(data = tot, aes(cell_type, N, label = comma(N)), inherit.aes = FALSE, vjust = -0.4, size = 3) +
    scale_fill_manual(values = pal, drop = FALSE) + scale_color_identity() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = comma) +
    guides(fill = guide_legend(nrow = nrow_leg)) +
    labs(title = title, x = NULL, y = sprintf("Pairs significant in all 6 methods (FDR < %g)", FDR_CUT)) +
    theme_paper() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1))
}
p5 <- (plot_dir("five", "Causal SMR, IVW, Egger\n+ Marginal Multiomics, Pseudobulk", pal_dir2) |
       plot_dir("mr3",  "Causal SMR, IVW, Egger", pal_dir2) |
       plot_dir("obs",  "Causal IVW vs.\nMarginal Multiomics / Pseudobulk", pal_dir4, nrow_leg = 4)) +
      plot_annotation(tag_levels = "a", title = "Effect direction among pairs significant in all six methods",
                      theme = theme(plot.title = element_text(face = "bold", size = 12)))

###------------------------------------------------------------------------------------------------
### 输出 (5 页)
out <- file.path(FIG_DIR, "Fig1_MR_overview.pdf")
cairo_pdf(out, width = 11, height = 5.2, onefile = TRUE)
print(p1)
for (ct in CELLS) plot_upset(ct)
print(p5)
invisible(dev.off())
msg("Saved: ", out)
