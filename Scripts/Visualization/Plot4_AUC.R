# Fig4: 逐 gene ROC AUC (小提琴 + 箱线), 每个细胞类型一页, 每页 6 个 panel = 6 个金标准
# 移植自 20261006_paper/Scripts/Plot4_val_AUC.sh
# 输入: Results/Validation/Valid4_AUC/<CT>_AUC_long.tsv.gz
# 输出: Results/Figures/Fig4_AUC.pdf
# 运行: Rscript Scripts/Visualization/Plot4_AUC.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")

gold_lbls <- c(enhancerAtlas = "EnhancerAtlas", FANTOM5 = "FANTOM", `4DGenome` = "4DGenome",
               GTEx = "GTEx cisQTL", HiChIP = "HiChIP", PCHiC = "PCHiC")   # 数据中的金标准名 -> 显示名 (沿用 paper)

auc <- rbindlist(lapply(CELLS, function(ct) fread(file.path(VAL_DIR, "Valid4_AUC", paste0(ct, "_AUC_long.tsv.gz")))))
auc[, `:=`(gold_lbl     = factor(gold_lbls[gold_standard], levels = gold_lbls),
           method_label = factor(lab_key[method], levels = METHODS))]
stopifnot(!anyNA(auc$gold_lbl), !anyNA(auc$method_label))

auc_sum <- auc[, .(n_genes = .N, mean_auc = mean(auc, na.rm = TRUE), median_auc = median(auc, na.rm = TRUE),
                   sd_auc = sd(auc, na.rm = TRUE)), by = .(cell_type, gold_lbl, method_label)][, se_auc := sd_auc / sqrt(n_genes)]
setorder(auc_sum, cell_type, gold_lbl, method_label)
cat("\n----- mean AUC (rows = cell_type x gold standard)\n")
print(dcast(auc_sum[, .(cell_type, gold_lbl, method_label, m = round(mean_auc, 3))], cell_type + gold_lbl ~ method_label, value.var = "m"))

plot_cell <- function(TAG) {
  d  <- auc[cell_type == TAG]
  nn <- auc_sum[cell_type == TAG]
  ggplot(d, aes(x = method_label, y = auc, fill = method_label)) +
    geom_violin(trim = TRUE, alpha = 0.8, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black") +
    geom_hline(yintercept = 0.5, linetype = 2, color = "grey50") +
    geom_text(data = nn, aes(x = method_label, y = 0.27, label = paste0("n=", n_genes)),
              inherit.aes = FALSE, size = 2.5) +
    facet_wrap(~ gold_lbl, nrow = 2, drop = FALSE) +
    scale_fill_manual(values = pal_method) +
    scale_y_continuous(limits = c(0.25, 1), breaks = seq(0.25, 1, by = 0.25)) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1),
          strip.background = element_rect(fill = "grey95", color = "grey80"),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(3, "mm")) +
    labs(title = sprintf("%s (%s)", TAG, CELL_DESC[[TAG]]), y = "Gene-wise ROC AUC")
}

out <- file.path(FIG_DIR, "Fig4_AUC.pdf")
pdf(out, width = 13, height = 8.5, onefile = TRUE)
for (TAG in CELLS) print(plot_cell(TAG))
invisible(dev.off())
msg("Saved: ", out)
