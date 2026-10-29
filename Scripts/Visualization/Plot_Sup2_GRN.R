# Sup2 (Figure S5): 各方法的显著 pair 与 CIMA SCENIC+ GRN 的重叠
# 2 个 panel: A = 显著 (p.adj < 0.05) 且在 GRN 中的 pair 数 (行 = GRN 集合, 列 = 细胞类型, 柱 = 6 个方法);
#             B = 显著 pair 在 GRN 中的 odds ratio (GRN vs 非 GRN, 2x2 Fisher; 每个细胞类型内 CT 与 CT_HQ 并排)
# GRN 集合: CT = 该细胞类型中 active (AUC > 1) 的 eRegulon; CT_HQ = CT 与对应 lineage 的高质量 eRegulon 的交集 (流程见 Scripts/Pipeline.md 3.5)
# 统计量已在 Valid5_GRN 中算好, 这里只画图 (沿用 Plot7_GRN.R 的数据读取和检查); 说明 (阈值 / GRN 集合 / pair 数) 写在论文 caption
# 输入: Results/Validation/Valid5_GRN/GRN_sig_counts_by_method.tsv, GRN_summary.tsv
# 输出: Results/Figures/Figures_20261027/Figure_S5.pdf + Figure_S5.svg
# 运行: Rscript Scripts/Visualization/Plot_Sup2_GRN.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages({ library(patchwork); library(scales) })

FDR_CUT <- 0.05   # 与 Valid5_GRN 一致 (本图的阈值与正文主图的 0.1 不同)
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

CT_SHORT <- c(Bn_TCL1A = "Bn", CD4_Tn_CCR7 = "CD4", CD8_Tn_CCR7 = "CD8")
GRN_LAB  <- c(CT = "Cell-type\neRegulons", CT_HQ = "High-quality\neRegulons")   # 换行: 右侧行标签放得下; 图例里换回空格
COL_GRN  <- c(CT = "#A9B8C8", CT_HQ = "#44546A")      # 同色系浅 / 深: CT_HQ 是 CT 的子集
sig_tab[, `:=`(method  = factor(method, levels = METHOD_KEY),
               grn_set = factor(grn_set, levels = names(GRN_LAB), labels = GRN_LAB),
               ct      = factor(CT_SHORT[cell_type], levels = CT_SHORT))]
BASEF <- 11

theme_s5 <- function() theme_bw(base_size = BASEF) + theme(
  text = element_text(color = "black"), axis.text = element_text(color = "black", size = BASEF - 1),
  panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
  panel.border = element_rect(colour = "grey35", linewidth = 0.5),
  strip.background = element_rect(fill = "#F5F6F8", colour = "grey35", linewidth = 0.5), strip.text = element_text(face = "bold", size = BASEF),
  axis.text.x = element_text(angle = 40, hjust = 1), plot.tag = element_text(face = "bold", size = BASEF + 9),
  plot.margin = margin(4, 8, 4, 4))

p_n <- ggplot(sig_tab, aes(method, n_sig_GRN, fill = method)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = comma(n_sig_GRN)), vjust = -0.4, size = (BASEF - 3) / .pt) +
  facet_grid(grn_set ~ ct) +
  scale_fill_manual(values = pal_key, guide = "none") +
  scale_x_discrete(labels = lab_key) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)), labels = comma) +
  labs(x = NULL, y = "Significant pairs overlapping the CIMA GRN") +
  theme_s5()

p_or <- ggplot(sig_tab, aes(method, OR, fill = grn_set)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_col(position = position_dodge(0.78), width = 0.72) +
  facet_wrap(~ ct, nrow = 1) +
  scale_fill_manual(values = setNames(COL_GRN, GRN_LAB), labels = function(x) gsub("\n", " ", x), name = NULL) +
  scale_x_discrete(labels = lab_key) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = NULL, y = "Odds ratio (GRN vs non-GRN pairs)") +
  theme_s5() + theme(legend.position = "top", legend.text = element_text(size = BASEF), legend.key.size = unit(0.4, "cm"),
                     plot.margin = margin(4, 8 + 37, 4, 4))   # 右侧留出 panel A 行标签的宽度, 使 A、B 的三列对齐

fig <- (p_n / p_or) + plot_layout(heights = c(2, 1.05)) + plot_annotation(tag_levels = "A")

FIG_W <- 10; FIG_H <- 9.5
OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261027")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- file.path(OUT_DIR, "Figure_S5.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H)
print(fig)
invisible(dev.off())
msg("Saved: ", out)
out <- file.path(OUT_DIR, "Figure_S5.svg")
svglite::svglite(out, width = FIG_W, height = FIG_H)   # svglite: 文字保留为可编辑文本, 便于拼图
print(fig)
invisible(dev.off())
msg("Saved: ", out)
