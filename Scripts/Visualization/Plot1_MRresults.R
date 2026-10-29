# Fig1: MR 总体结果 (FDR < 0.1)
#   每个细胞类型一行: 左 UpSet (Figure 1B; 前 15 个组合, 已含每个方法的显著 pair 数)
#                    + 右 六方法全显著 pairs 的效应方向一致性环形图 (Figure 1C)
#   CD4 -> Figure_1B_C (正文; 已在 Figma 拼成 Figure_1, 默认不重写); Bn (上, 字母 A) + CD8 (下, 字母 B) -> Figure_S3 (补充, 同一页)
#   图中阈值写作 p.adj (即各方法 gene 内 BH 后的 *_fdr 列)
# 另外: 在新 FDR 下检查 24 x 3 个 example pairs 是否仍满足各自类别 (结果只打印到日志, 并写 tmp/example_pairs_check.tsv)
# 移植自 20261006_paper/Scripts/Plot1_MRresults.sh
# 输入: Results/MR/<CT>_AllPairs_results.tsv.gz; (对照) 旧 MR_results_formal/<CT>_AllPairs_results.rds
# 输出: Results/Figures/Figures_20261027/Figure_1B_C.{pdf,svg} (CD4) + Figure_S3.{pdf,svg} (Bn + CD8)
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

### 每个方法的显著 pairs 数 (只打印)
cnt <- rbindlist(lapply(CELLS, function(ct) { m <- sig_mat(res_list[[ct]]); data.table(cell_type = ct, n_total = nrow(m), method = METHODS, n_sig = colSums(m)) }))
print(dcast(cnt, cell_type + n_total ~ factor(method, levels = METHODS), value.var = "n_sig"))

### 左: UpSet (单方法特有的组合用该方法颜色, 其余深灰); 画在当前 viewport 内
plot_upset <- function(ct, n_top = 15) {
  m  <- sig_mat(res_list[[ct]]); m <- m[rowSums(m) > 0, , drop = FALSE]
  cm <- make_comb_mat(m, mode = "distinct")
  cm <- cm[rank(-comb_size(cm), ties.method = "first") <= n_top]
  cs <- comb_size(cm); cd <- comb_degree(cm)
  comb_cols <- unname(ifelse(cd == 1, pal_method[set_name(cm)][regexpr("1", comb_name(cm))], COL_INK))
  ax_at <- function(x) { b <- pretty(c(0, x), n = 3); b[b <= x] }   # 刻度只取 3-4 个, 避免重叠
  ra <- upset_right_annotation(cm, add_numbers = TRUE, numbers_gp = gpar(fontsize = 6.5), gp = gpar(fill = pal_method[set_name(cm)], col = NA),
    width = unit(2.2, "cm"), annotation_name_gp = gpar(fontsize = 8),
    axis_param = list(at = 0, labels = "", gp = gpar(col = NA)))   # 右侧条形图不显示 X 轴刻度 (每根条已标数字)
  ra@anno_list[[1]]@label <- "Significant pairs"   # upset_right_annotation 把 "Set size" 写死, 这里直接改
  ht <- UpSet(cm, set_order = intersect(METHODS, set_name(cm)), comb_order = order(cs, decreasing = TRUE),
    pt_size = unit(2, "mm"), lwd = 1.4, comb_col = comb_cols, bg_col = c("#F3F5F8", "#FFFFFF"), bg_pt_col = "#D9DEE3",
    row_names_gp = gpar(fontsize = 7.5, fontface = "bold", col = pal_method[set_name(cm)]),
    top_annotation = upset_top_annotation(cm, add_numbers = TRUE, numbers_rot = 0, numbers_gp = gpar(fontsize = 6), gp = gpar(fill = comb_cols, col = NA),
      height = unit(3, "cm"), annotation_name_gp = gpar(fontsize = 8), axis_param = list(at = ax_at(max(cs)), gp = gpar(fontsize = 6.5))),
    right_annotation = ra)
  draw(ht, newpage = FALSE, column_title = ct,
       column_title_gp = gpar(fontsize = 11, fontface = "bold"), padding = unit(c(2, 10, 2, 2), "mm"))   # 左侧留白, 防止长方法名被截
  invisible(NULL)
}

### 右: 六方法全显著 pairs 的方向一致性
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

### 环形图: 扇区 = Causal IVW 与两种 Marginal 方法的方向关系; 环下注明 SMR / IVW / Egger 三者内部的一致性
pal_dir4 <- c(`Concordant with both` = COL_CONC, `Discordant with Marginal Multiomics` = pal_method[["Marginal Multiomics"]],
              `Discordant with Marginal Pseudobulk` = pal_method[["Marginal Pseudobulk"]], `Discordant with both` = COL_DISC)
lab_dir4 <- c("Concordant with both Marginal methods", "Discordant with Marginal Multiomics only",
              "Discordant with Marginal Pseudobulk only", "Discordant with both Marginal methods")

plot_donut <- function(ct) {
  dd <- dir_dt[cell_type == ct]
  d  <- dd[, .N, by = .(status = obs)][, status := factor(status, levels = names(pal_dir4))][order(status)]
  d[, pct := N / sum(N)][, txt_col := fifelse(status == names(pal_dir4)[1], "#333333", "white")]
  # 图内不放图例 / 中心文字 / 小扇区标注; 两个大扇区 (>= 3%) 用短径向线 + 文字说明含义, n / SMR-IVW-Egger 一致比例写进论文 caption
  lab_short <- setNames(c("Concordant with\nboth Marginal", "Discordant with\nMultiomics only", "", ""), names(pal_dir4))
  d[, lab := lab_short[as.character(status)]]
  # position_stack 把第一个水平放在最上 (y 最大), 故 theta 累加顺序与 status 水平相反
  dr <- d[order(-as.integer(status))][, mid := (cumsum(N) - N / 2) / sum(N)]
  dl <- merge(d, dr[, .(status, mid)], by = "status")[lab != ""][, `:=`(hj = fifelse(mid < 0.5, 0, 1), y = mid * sum(dr$N))]
  ggplot(d, aes(x = 2, y = N, fill = status)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    geom_text(data = d[pct >= 0.03], aes(label = sprintf("%.0f%%", 100 * pct), color = txt_col),
              position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
    geom_segment(data = dl, aes(x = 2.5, xend = 2.85, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4, color = "grey35") +
    geom_text(data = dl, aes(x = 3.0, y = y, label = lab, hjust = hj), inherit.aes = FALSE, size = 2.6, color = "grey20", lineheight = 0.9) +
    coord_polar(theta = "y", clip = "off") + xlim(1.1, 3.3) +
    scale_fill_manual(values = pal_dir4, labels = lab_dir4, drop = FALSE) + scale_color_identity() +
    labs(title = "Effect direction (all six significant)") +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9.5), legend.position = "none",
          plot.margin = margin(6, 40, 6, 6))
}

###------------------------------------------------------------------------------------------------
### 输出: 每个细胞类型一行 (左 UpSet + 右环形图)
###   Figure_1B_C : CD4 (正文, 1 页, 9 x 3.6 in), 之后与 Figure_1A 在 Figma 合并为 Fig 1
###   Figure_S3   : Bn (上) + CD8 (下) 同一页 (补充, 9 x 7.2 in)
### PDF + SVG (svglite: 文字保留为可编辑文本, 便于拼图)
OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261027")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ROW_W <- 9; ROW_H <- 3.6
plot_row <- function(ct, row = 1, tag = NULL) {
  pushViewport(viewport(layout.pos.row = row, layout.pos.col = 1))
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.68, 0.32), "npc"))))
  pushViewport(viewport(layout.pos.col = 1)); plot_upset(ct); popViewport()
  print(plot_donut(ct), vp = viewport(layout.pos.col = 2))
  if (!is.null(tag)) grid.text(tag, x = unit(2, "mm"), y = unit(1, "npc") - unit(2, "mm"), just = c("left", "top"), gp = gpar(fontsize = 24, fontface = "bold"))   # panel 字母 (每个细胞类型一个); 必须最后画: ComplexHeatmap 的 draw() 默认铺白色背景, 会盖住先画的字母
  popViewport(2)
}
plot_page <- function(cts, tags = FALSE) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(length(cts), 1)))
  for (i in seq_along(cts)) plot_row(cts[i], i, if (tags) LETTERS[i] else NULL)
  popViewport()
}
FIGS <- list(Figure_1B_C = "CD4_Tn_CCR7", Figure_S3 = c("Bn_TCL1A", "CD8_Tn_CCR7"))
# Figure_1B_C 已在 Figma 中拼成 Figure_1, 默认不重写; 需要重出时设环境变量 PLOT1_WRITE_1BC=1
# Figure_S3 的每个细胞类型行带 panel 字母 (A = Bn, B = CD8); Figure_1B_C 不带 (字母在 Figma 的 Figure_1 里)
for (nm in names(FIGS)) {
  if (nm == "Figure_1B_C" && !nzchar(Sys.getenv("PLOT1_WRITE_1BC"))) { msg("Skip ", nm, " (设 PLOT1_WRITE_1BC=1 才重写)"); next }
  cts <- FIGS[[nm]]; h <- ROW_H * length(cts); tg <- (nm == "Figure_S3")
  out <- file.path(OUT_DIR, paste0(nm, ".pdf"))
  cairo_pdf(out, width = ROW_W, height = h, onefile = TRUE); plot_page(cts, tg); invisible(dev.off())
  msg("Saved: ", out)
  out <- file.path(OUT_DIR, paste0(nm, ".svg"))
  svglite::svglite(out, width = ROW_W, height = h); plot_page(cts, tg); invisible(dev.off())
  msg("Saved: ", out)
}
