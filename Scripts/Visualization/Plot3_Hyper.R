# Fig3: gene 内超几何富集 p 值直方图 (6 方法 x 6 金标准), 每个细胞类型一页
# 移植自 20261006_paper/Scripts/Plot3_val_Hyper.sh, 只改输入/输出路径
# 输入: Results/Validation/Valid1_Hyper/<CT>_hyper_long.tsv
# 输出: Results/Figures/Fig3_Hyper.pdf
# 运行: Rscript Scripts/Visualization/Plot3_Hyper.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")

in_dir <- file.path(VAL_DIR, "Valid1_Hyper")
hyper_list <- setNames(lapply(CELLS, function(ct) fread(file.path(in_dir, paste0(ct, "_hyper_long.tsv")))), CELLS)

m_lbls   <- c(SMR = "Causal SMR", IVW = "Causal IVW", Egger = "Causal Egger",
              CCT = "Causal Joint", WNN = "Marginal Multiomics", Pseudobulk = "Marginal Pseudobulk")   # 数据中的 method -> 显示名
db_names <- c("enhancerAtlas", "FANTOM5", "4DGenome", "GTEx", "HiChIP", "PCHiC")

base_th <- theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95", color = "grey80"),
        strip.text = element_text(face = "bold", size = 9),
        axis.text = element_text(color = "black", size = 8.5),
        panel.spacing = unit(3, "mm"))

# ---------- 单个细胞类型 -> 6x6 图; 各格子 N/K/n/k 求和后算一个 global P ----------
global_p <- function(hyper_long) {
  p_dt <- hyper_long[!is.na(pval) & N > 0 & K > 0 & n_sig > 0 & N > K][
    , `:=`(method_lbl = factor(m_lbls[as.character(method)], levels = METHODS),
           database   = factor(database, levels = db_names))]
  stopifnot(!anyNA(p_dt$method_lbl), !anyNA(p_dt$database))
  g_dt <- p_dt[, .(n_gene = .N, sum_N = sum(N), sum_K = sum(K), sum_nsig = sum(n_sig), sum_k = sum(k)),
               by = .(method_lbl, database)][
    , gp := phyper(sum_k - 1L, sum_K, sum_N - sum_K, sum_nsig, lower.tail = FALSE)][
    , p_txt := ifelse(gp < 1e-300, "Global P < 1e-300", sprintf("Global P = %.2e", gp))]
  list(p_dt = p_dt, g_dt = g_dt)
}

plot_hyper <- function(x, TAG) {
  ggplot(x$p_dt, aes(x = pval, fill = method_lbl)) +
    geom_histogram(breaks = seq(0, 1, 0.05), color = "grey25", linewidth = 0.15, alpha = 0.85) +
    geom_text(data = x$g_dt, aes(x = 0.98, y = Inf, label = p_txt),
              hjust = 1, vjust = 1.4, size = 2.7, fontface = "italic", inherit.aes = FALSE) +
    facet_grid(method_lbl ~ database, scales = "free_y") +
    scale_fill_manual(values = pal_method) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    labs(title = sprintf("Hypergeometric Enrichment (%s)", TAG),
         x = expression(italic(P) * "-value"), y = "Gene Count") +
    base_th
}

res <- lapply(hyper_list, global_p)

# 日志: 每个格子的 gene 数和 global P
for (ct in CELLS) {
  cat("\n=====", ct, "\n")
  print(dcast(res[[ct]]$g_dt, method_lbl ~ database, value.var = "n_gene"))
  print(dcast(res[[ct]]$g_dt[, .(method_lbl, database, gp = signif(gp, 3))], method_lbl ~ database, value.var = "gp"))
}

out <- file.path(FIG_DIR, "Fig3_Hyper.pdf")
pdf(out, width = 13, height = 10.5, onefile = TRUE)
for (ct in CELLS) print(plot_hyper(res[[ct]], ct))
invisible(dev.off())
msg("Saved: ", out)
