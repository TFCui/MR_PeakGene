# Fig5: S-LDSC 遗传力富集 (ext500; 9 traits x 3 CT; MR-Only vs All peaks)
# 2026-10-13: 去掉 SE 误差棒; enrichment raw P < 0.1 的柱子上加粗体黑色 *
# 2026-10-14: 第 2 页 = Combined (3 个 CT 并集)
# 2026-10-14: 只画 MR-Only (IVW/Egger/SMR 任一 FDR < 0.1 且 WNN、Pseudo 都不显著的 link 的 peak) + All peaks; 不再画 6 个方法
# 移植自 20261006_paper/Scripts/Plot5_val_SLDSC.sh Part 2 (版式); 2 页: 第 1 页 3 个 CT, 第 2 页 Combined
# 输入: Results/Validation/Valid2_SLDSC/SLDSC_all.tsv
# 输出: Results/Figures/Figures_20261019/Fig5_SLDSC.pdf
# 运行: Rscript Scripts/Visualization/Plot5_SLDSC.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")

res <- fread(file.path(VAL_DIR, "Valid2_SLDSC/SLDSC_all.tsv"))

KEYS      <- c("MROnly", "pool")
TRAITS    <- c("AS", "CD", "IBD", "RA", "SLE", "MS", "T1D", "Asthma", "PBC")
tag_lab   <- c(CD4_Tn_CCR7 = "CD4", Bn_TCL1A = "Bn", CD8_Tn_CCR7 = "CD8", Combined = "Bn + CD4\n+ CD8")   # 行顺序沿用 paper: CD4 / Bn / CD8
lab_all   <- c(MROnly = "MR-Only", pool = "All peaks")
pal_fill  <- c(`MR-Only` = "#6FA8DC", `All peaks` = COL_CONC)
pal_line  <- c(`MR-Only` = "#4F8CC4", `All peaks` = "#A3AFBC")   # 柱子加一圈深一点的边

theme_pub <- function(base = 9) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      plot.title         = element_text(face = "bold", size = base + 3, colour = "#111827"),
      plot.subtitle      = element_text(size = base - 0.5, colour = "#6B7280", margin = margin(b = 6)),
      strip.background   = element_rect(fill = "#F5F6F8", colour = NA),
      strip.text.x       = element_text(face = "bold", size = base, colour = "#1F2937", margin = margin(t = 3, b = 3)),
      strip.text.y       = element_text(face = "bold", size = base, colour = "#1F2937", angle = 0, margin = margin(l = 4, r = 4)),
      axis.line          = element_line(colour = "#9CA3AF", linewidth = 0.3),
      axis.ticks         = element_line(colour = "#9CA3AF", linewidth = 0.3),
      axis.ticks.length.x = unit(0.08, "cm"),
      axis.text.x        = element_text(face = "bold", size = base - 2, colour = "#374151", angle = 45, hjust = 1, vjust = 1),
      axis.text.y        = element_text(size = base - 1.5, colour = "#4B5563"),
      axis.title.y       = element_text(face = "bold", size = base, colour = "#1F2937", margin = margin(r = 5)),
      panel.grid.major.y = element_line(colour = "#EEF0F3", linewidth = 0.3),
      panel.spacing.x    = unit(0.7, "lines"),
      panel.spacing.y    = unit(0.8, "lines"),
      legend.position    = "top",
      legend.justification = "left",
      legend.box         = "horizontal",
      legend.title       = element_text(face = "bold", size = base - 0.5),
      legend.text        = element_text(size = base - 1),
      legend.key.size    = unit(0.32, "cm"),
      legend.margin      = margin(b = -2),
      plot.margin        = margin(10, 12, 8, 10)
    )
}

facet_layer <- if (requireNamespace("ggh4x", quietly = TRUE)) {
  ggh4x::facet_grid2(tag ~ trait, scales = "free_y", independent = "y")
} else {
  message("ggh4x 未安装, 使用 facet_wrap 代替")
  facet_wrap(tag ~ trait, nrow = 3, scales = "free_y", labeller = labeller(.multi_line = FALSE))
}

# 一页的数据 (cts = 要画的 CT, 行顺序按 cts)
page_dt <- function(cts) {
  dt <- res[trait %in% TRAITS & set %in% KEYS & cell_type %in% cts]
  stopifnot(nrow(dt) == length(TRAITS) * length(cts) * length(KEYS))
  dt[, `:=`(method = factor(lab_all[set], levels = lab_all[KEYS]),
            tag    = factor(tag_lab[cell_type], levels = tag_lab[cts]),
            trait  = factor(trait, levels = TRAITS))]
  dt[, `:=`(val = enrich, star_y = pmax(enrich, 0))]
  cat("----- enrichment (enrich_p < 0.1 marked *):", paste(cts, collapse = ", "), "\n")
  print(dcast(dt[, .(tag, trait, method, e = sprintf("%.1f%s (p=%.3f)", enrich, ifelse(enrich_p < 0.1, "*", ""), enrich_p))], tag + trait ~ method, value.var = "e"))
  dt
}

plot_page <- function(dt, title, subtitle) {
  ggplot(dt, aes(x = method, y = val)) +
    geom_hline(yintercept = 1, linetype = "22", colour = "#9CA3AF", linewidth = 0.35) +
    geom_col(aes(fill = method, colour = method), width = 0.7, linewidth = 0.3, alpha = 0.95) +
    geom_text(data = dt[enrich_p < 0.1], aes(y = star_y), label = "*", fontface = "bold", colour = "black", size = 5, vjust = 0.3) +
    facet_layer +
    scale_fill_manual(values = pal_fill, name = "Peaks") +
    scale_colour_manual(values = pal_line, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = scales::pretty_breaks(n = 3)) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(colour = pal_line[lab_all[KEYS]]))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Heritability enrichment") +
    theme_pub()
}

mr_def <- "MR-Only = peaks of links significant (FDR < 0.1) in Causal SMR, IVW or Egger but in neither Multiomics nor Pseudobulk"
foot   <- "Each annotation conditioned on 97 baseline-LD annotations. Dashed line: enrichment = 1. * enrichment P < 0.1."

# 第 1 页: 3 个 CT 各自
p1 <- plot_page(page_dt(c("CD4_Tn_CCR7", "Bn_TCL1A", "CD8_Tn_CCR7")),
                "S-LDSC heritability enrichment: MR-Only peaks vs all candidate peaks",
                paste0("Peak extension = 500 bp. ", mr_def, ".\nAll peaks = all candidate peaks of the cell type. ", foot))

# 第 2 页: Combined (3 个 CT 的并集); 只有一行, 用 aspect.ratio 让面板和第 1 页差不多高
p2 <- plot_page(page_dt("Combined"),
                "S-LDSC heritability enrichment: Bn + CD4 + CD8 combined",
                paste0("Peak extension = 500 bp. MR-Only = union of the MR-Only peaks of Bn, CD4 and CD8 (defined per cell type).\n",
                       mr_def, ".\nAll peaks = union of all candidate peaks of the three cell types. ", foot)) +
  theme(aspect.ratio = 2)

out <- file.path(FIG_DIR, "Fig5_SLDSC.pdf")
if (capabilities("cairo")) cairo_pdf(out, width = 12, height = 7.5, onefile = TRUE) else pdf(out, width = 12, height = 7.5, useDingbats = FALSE, onefile = TRUE)
print(p1); print(p2)
invisible(dev.off())
msg("Saved: ", out)
