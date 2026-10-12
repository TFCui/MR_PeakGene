# Fig5: S-LDSC 遗传力富集 (ext500; 9 traits x 3 CT; 6 个方法的显著 peak + All peaks)
# 移植自 20261006_paper/Scripts/Plot5_val_SLDSC.sh Part 2; 只做第 1 页 (paper 第 2/3 页需要 MR-All 集合, 本项目没有跑)
# 输入: Results/Validation/Valid2_SLDSC/SLDSC_all.tsv
# 输出: Results/Figures/Fig5_SLDSC.pdf
# 运行: Rscript Scripts/Visualization/Plot5_SLDSC.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")

res <- fread(file.path(VAL_DIR, "Valid2_SLDSC/SLDSC_all.tsv"))

KEYS      <- c(METHOD_KEY, "pool")
TRAITS    <- c("AS", "CD", "IBD", "RA", "SLE", "MS", "T1D", "Asthma", "PBC")
tag_lab   <- c(CD4_Tn_CCR7 = "CD4", Bn_TCL1A = "Bn", CD8_Tn_CCR7 = "CD8")   # 行顺序沿用 paper: CD4 / Bn / CD8
lab_all   <- c(lab_key, pool = "All peaks")
pal_fill  <- c(pal_method, `All peaks` = COL_CONC)
pal_line  <- c(pal_method, `All peaks` = "#A3AFBC")   # 浅色背景柱加一圈深一点的边

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

dt <- res[trait %in% TRAITS & set %in% KEYS]
stopifnot(nrow(dt) == length(TRAITS) * length(CELLS) * length(KEYS))
dt[, `:=`(method = factor(lab_all[set], levels = lab_all[KEYS]),
          tag    = factor(tag_lab[cell_type], levels = tag_lab),
          trait  = factor(trait, levels = TRAITS))]
dt[, `:=`(val = enrich, lo = pmax(enrich - enrich_se, 0), hi = enrich + enrich_se)]

cat("----- enrichment (enrich_p < 0.05 marked *)\n")
print(dcast(dt[, .(tag, trait, method, e = sprintf("%.1f%s", enrich, ifelse(enrich_p < 0.05, "*", "")))], tag + trait ~ method, value.var = "e"))

p1 <- ggplot(dt, aes(x = method, y = val)) +
  geom_hline(yintercept = 1, linetype = "22", colour = "#9CA3AF", linewidth = 0.35) +
  geom_col(aes(fill = method, colour = method), width = 0.76, linewidth = 0.3, alpha = 0.95) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, colour = COL_INK, linewidth = 0.3) +
  facet_layer +
  scale_fill_manual(values = pal_fill, name = "Peaks") +
  scale_colour_manual(values = pal_line, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = scales::pretty_breaks(n = 3)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(colour = pal_line[lab_all[KEYS]]))) +
  labs(title    = "S-LDSC heritability enrichment of MR-linked peaks",
       subtitle = paste0("Peak extension = 500 bp; significant peak-gene links at FDR < 0.05 per method; ",
                         "All peaks = all candidate peaks of the cell type. Each annotation conditioned on 97 baseline-LD annotations. Dashed line: enrichment = 1."),
       x = NULL, y = "Heritability enrichment (± s.e.)") +
  theme_pub()

out <- file.path(FIG_DIR, "Fig5_SLDSC.pdf")
if (capabilities("cairo")) cairo_pdf(out, width = 18, height = 7.5) else pdf(out, width = 18, height = 7.5, useDingbats = FALSE)
print(p1)
invisible(dev.off())
msg("Saved: ", out)
