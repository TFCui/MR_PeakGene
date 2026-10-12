# Fig7: 与 CIMA SCENIC+ GRN 的重叠 (上: 显著且在 GRN 中的 pair 数; 下: OR, GRN vs 非 GRN)
# 移植自 20261006_paper/Scripts/Plot7_val_GRN.sh 第 6 节; 统计量已在 Valid5_GRN 中算好, 这里只画图
# 输入: Results/Validation/Valid5_GRN/GRN_sig_counts_by_method.tsv (FDR < 0.05)
# 输出: Results/Figures/Fig7_GRN.pdf
# 运行: Rscript Scripts/Visualization/Plot7_GRN.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
library(patchwork)

FDR_CUT <- 0.05   # 与 Valid5_GRN 一致
sig_tab <- fread(file.path(VAL_DIR, "Valid5_GRN/GRN_sig_counts_by_method.tsv"))
stopifnot(nrow(sig_tab) == 2 * length(CELLS) * length(METHOD_KEY), all(sig_tab$method %in% METHOD_KEY))

print(dcast(sig_tab, grn_set + cell_type + n_pairs_GRN ~ method, value.var = "n_sig_GRN"))
print(dcast(sig_tab[, .(grn_set, cell_type, method, OR = round(OR, 2))], grn_set + cell_type ~ method, value.var = "OR"))

sig_tab[, `:=`(method  = factor(method, levels = METHOD_KEY),
               grn_set = factor(grn_set, levels = c("CT", "CT_HQ")))]

p_n <- ggplot(sig_tab, aes(method, n_sig_GRN, fill = method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_sig_GRN), vjust = -0.4, size = 3.5) +
  facet_grid(grn_set ~ cell_type) +
  scale_fill_manual(values = pal_key, guide = "none") +
  scale_x_discrete(labels = lab_key) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "No. of significant pairs\noverlapping CIMA GRN",
       subtitle = sprintf("FDR < %.2f; rows = GRN set (CT / CT_HQ), columns = cell type", FDR_CUT)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey92"), axis.text.x = element_text(angle = 35, hjust = 1))

p_or <- ggplot(sig_tab, aes(method, OR, fill = grn_set)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ cell_type) +
  scale_fill_manual(values = c(CT = "#0055AA", CT_HQ = "#E69F00"), name = "GRN") +
  scale_x_discrete(labels = lab_key) +
  labs(x = NULL, y = "Odds ratio (GRN vs non-GRN)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"), axis.text.x = element_text(angle = 35, hjust = 1))

out <- file.path(FIG_DIR, "Fig7_GRN.pdf")
pdf(out, width = 14, height = 11)
print(p_n / p_or + plot_layout(heights = c(2, 1)))
invisible(dev.off())
msg("Saved: ", out)
