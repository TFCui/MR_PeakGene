# Plot14 (论文 S6 Fig, 待定): 5 个方法两两之间的相关性与密度矩阵, 每个细胞类型一页, 每页两个 panel 上下排 (A 上, B 下)
# 5 个方法: Causal SMR / IVW / Egger, Marginal Multiomics / Pseudobulk (CCT 没有 effect, 不画)
#   panel A = 原始 effect (SMR smr_b, IVW ivw_beta, Egger egger_beta; WNN wnn_rho, Pseudobulk pseudo_rho)
#   panel B = signed -log10(P) = sign(effect) * -log10(P) (旧图 Server6_wgen1_plot.sh 的口径)
#   (2026-11-04 起原 S6 / S7 合并为一个 PDF; 页顶只写细胞类型, 不写 pair 数)
# 每个 panel 是 5x5 矩阵:
#   下三角: 二维核密度 (MASS::kde2d, n = 160, 带宽 1.5 x bandwidth.nrd), 填色 density^0.25 (同 smoothScatter), 蓝色系;
#           x = 0 / y = 0 灰点线; y = x 虚线 (A 只在两轴都是 MR beta 时画, B 每格都画)
#   上三角: 该 CT 全部 pair 的 Spearman r (不截断), 背景按 r 低饱和发散色, 字号随 |r|
#   对角线: 方法名 (方法配色) + 该方法自身数值的分布曲线 (density(), 高度缩放到格子下半部, 只看形状)
#   刻度只在最下一行 / 最左一列; 左上角格子不画纵轴刻度 (第一行没有密度图, 刻度易被误读为曲线读数)
# 坐标范围: A 每个方法单独对称范围 = 3 个 CT 中 |x| 的 99% 分位的最大值 (向上取整到 1 / 1.5 / 2 ...) (3 页共用);
#           B 所有方法共用一个范围 = 合并后 |x| 的 99% 分位 (向上取整)
#           超出范围的 pair 不进入密度 (直接舍去, 不压到边上); 显示比例写入 tsv
# 输入: Results/MR/<CT>_AllPairs_results.tsv.gz
# 输出: Results/Figures/Figures_20261102/Figure_S6.pdf (3 页: Bn / CD4 T cell / CD8 T cell; 每页 A effect + B signed -log10 P)
#       + Figure_S6_cor.tsv (type = cor: 每个 刻度 x CT x 方法对的 r 和 n; type = range: 每个方法的轴范围和显示比例)
# 运行: Rscript Scripts/Visualization/Plot14_MethodCorr.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages({ library(patchwork); library(MASS) })

M5     <- setdiff(METHODS, "Causal Joint")
EFFCOL <- c(`Causal SMR` = "smr_b", `Causal IVW` = "ivw_beta", `Causal Egger` = "egger_beta",
            `Marginal Multiomics` = "wnn_rho", `Marginal Pseudobulk` = "pseudo_rho")
PCOL   <- setNames(paste0(METHOD_KEY[M5], "_p"), M5)
IS_MR  <- setNames(grepl("^Causal", M5), M5)
DIAG_LAB <- sub(" ", "\n", M5)   # 两行: Causal\nSMR
names(DIAG_LAB) <- M5

Q_LIM  <- 0.99
KDE_N  <- 160
KDE_HX <- 1.5   # 带宽放大 1.5 倍, 尾部更平滑
FONT   <- "sans"
PAL_DENS <- c("#FFFFFF", "#DCEAF5", "#7FB2D9", "#2F6FA8")
PAL_R    <- colorRampPalette(c("#E39C8A", "#F6D9D0", "#FFFFFF", "#D7E7F4", "#86B4DA"))(201)
COL_YX   <- "#E08A6B"

###------------------------------------------------------------------------------------------------
### 1. 读数据, 两种刻度的值
dat <- setNames(lapply(CELLS, function(ct) {
  d <- read_mr(ct)
  eff <- as.data.table(lapply(EFFCOL, function(cc) d[[cc]]))
  slp <- as.data.table(lapply(M5, function(mm) sign(d[[EFFCOL[[mm]]]]) * -log10(d[[PCOL[[mm]]]])))
  setnames(eff, M5); setnames(slp, M5)
  stopifnot(!anyNA(eff), !anyNA(slp), all(is.finite(as.matrix(slp))))
  list(effect = eff, signedP = slp)
}), CELLS)

nice_up <- function(v) { e <- 10^floor(log10(v)); ceiling(v / e * 2) / 2 * e }   # 向上取到 1 / 1.5 / 2 / 2.5 ... x 10^k
LIMS <- list(
  effect  = sapply(M5, function(mm) nice_up(max(sapply(CELLS, function(ct) quantile(abs(dat[[ct]]$effect[[mm]]), Q_LIM))))),
  signedP = { L <- ceiling(quantile(abs(unlist(lapply(CELLS, function(ct) unlist(dat[[ct]]$signedP)))), Q_LIM)); setNames(rep(L, 5), M5) })
cat("\n===== 轴范围 (对称, ±)\n"); print(sapply(LIMS, identity))

###------------------------------------------------------------------------------------------------
### 2. 每个格子
theme_cell <- function(show_x, show_y) theme_bw(base_size = 8, base_family = FONT) + theme(
  panel.grid = element_blank(), panel.border = element_rect(fill = NA, colour = "grey35", linewidth = 0.4),
  axis.title = element_blank(), axis.text = element_text(size = 7, colour = "grey20"),
  axis.ticks = element_line(linewidth = 0.3, colour = "grey35"), axis.ticks.length = unit(1.5, "pt"),
  axis.text.x = if (show_x) element_text(size = 7, colour = "grey20") else element_blank(),
  axis.ticks.x = if (show_x) element_line(linewidth = 0.3, colour = "grey35") else element_blank(),
  axis.text.y = if (show_y) element_text(size = 7, colour = "grey20") else element_blank(),
  axis.ticks.y = if (show_y) element_line(linewidth = 0.3, colour = "grey35") else element_blank(),
  legend.position = "none", plot.margin = margin(1.5, 1.5, 1.5, 1.5))
brk <- function(L) { b <- pretty(c(-L, L), n = 4); b[abs(b) < L * 0.9] }   # 不放两端刻度, 避免相邻格子的刻度字相撞
sc_xy <- function(Lx, Ly) list(scale_x_continuous(breaks = brk(Lx)), scale_y_continuous(breaks = brk(Ly)),
                               coord_cartesian(xlim = c(-Lx, Lx), ylim = c(-Ly, Ly), expand = FALSE))

lower_cell <- function(x, y, Lx, Ly, yx, show_x, show_y) {
  keep <- abs(x) <= Lx & abs(y) <= Ly
  x <- x[keep]; y <- y[keep]
  k <- kde2d(x, y, n = KDE_N, h = KDE_HX * c(bandwidth.nrd(x), bandwidth.nrd(y)), lims = c(-Lx, Lx, -Ly, Ly))
  g <- CJ(j = seq_along(k$y), i = seq_along(k$x))
  g[, `:=`(x = k$x[i], y = k$y[j], z = k$z[cbind(i, j)])]
  g[, z := (pmax(z, 0) / max(z))^0.25]
  p <- ggplot(g, aes(x, y)) +
    geom_raster(aes(fill = z), interpolate = TRUE) +
    scale_fill_gradientn(colours = PAL_DENS, limits = c(0, 1)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.25, linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.25, linetype = "dotted")
  if (yx) p <- p + geom_abline(slope = 1, intercept = 0, colour = COL_YX, linewidth = 0.35, linetype = "dashed")
  p + sc_xy(Lx, Ly) + theme_cell(show_x, show_y)
}

upper_cell <- function(r) {
  col <- PAL_R[round((r + 1) * 100) + 1]
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = col) +
    annotate("text", x = 0.5, y = 0.5, label = sprintf("r = %.2f", r), family = FONT, fontface = "bold",
             size = 3.4 + 3.2 * abs(r), colour = COL_INK) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_cell(FALSE, FALSE)
}

diag_cell <- function(x, mm, L, show_x, show_y) {
  x <- x[abs(x) <= L]
  dd <- density(x, n = 512, from = -L, to = L)
  dt <- data.table(x = dd$x, y = -L + dd$y / max(dd$y) * 1.05 * L)   # 曲线占下半格, 坐标用本方法的范围 (左列 / 底行刻度才对)
  ggplot(dt, aes(x, y)) +
    geom_ribbon(aes(ymin = -L, ymax = y), fill = pal_method[[mm]], alpha = 0.22) +
    geom_line(colour = pal_method[[mm]], linewidth = 0.45) +
    annotate("text", x = 0, y = 0.52 * L, label = DIAG_LAB[[mm]], family = FONT, fontface = "bold",
             size = 3.6, lineheight = 0.9, colour = pal_method[[mm]]) +
    sc_xy(L, L) + theme_cell(show_x, show_y)
}

###------------------------------------------------------------------------------------------------
### 3. 每页 = 一个 CT; 两个 panel (A effect, B signed -log10 P), 每个 panel 一个 5x5 矩阵
COR <- list(); RNG <- list()
make_matrix <- function(sc, ct) {
  X <- dat[[ct]][[sc]]; L <- LIMS[[sc]]
  cells <- vector("list", 25)
  for (i in 1:5) for (j in 1:5) {   # i = 行 (y), j = 列 (x)
    mi <- M5[i]; mj <- M5[j]
    show_x <- i == 5; show_y <- j == 1 && i > 1   # 左上角格子不画纵轴刻度
    cells[[(i - 1) * 5 + j]] <-
      if (i == j) diag_cell(X[[mi]], mi, L[[mi]], show_x, show_y)
      else if (i > j) lower_cell(X[[mj]], X[[mi]], L[[mj]], L[[mi]], yx = sc == "signedP" || (IS_MR[[mi]] && IS_MR[[mj]]), show_x, show_y)
      else {
        r <- cor(X[[mi]], X[[mj]], method = "spearman")
        COR[[length(COR) + 1]] <<- data.table(scale = sc, cell_type = ct, method1 = mi, method2 = mj, r = r, n = nrow(X))
        upper_cell(r)
      }
  }
  RNG[[length(RNG) + 1]] <<- data.table(scale = sc, cell_type = ct, method1 = M5, axis_limit = unname(L[M5]),
                                        n = nrow(X), n_shown = sapply(M5, function(mm) sum(abs(X[[mm]]) <= L[[mm]])))
  wrap_plots(cells, nrow = 5)
}
PANEL_TITLE <- c(effect = "Effect estimate", signedP = "Signed −log10(P)")
panel_head <- function(tag, txt) ggplot() +   # panel 标题行 (单独一行, 不压在矩阵上)
  annotate("text", x = 0, y = 0.5, hjust = 0, vjust = 0.5, family = FONT, size = 4.8,
           label = sprintf("bold(%s)~~~bold('%s')", tag, txt), parse = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") + theme_void()
make_page <- function(ct) {   # 两个 panel 上下排 (补充材料页面不够宽): 标题行 A / 矩阵 A / 标题行 B / 矩阵 B
  wrap_plots(list(panel_head("A", PANEL_TITLE[["effect"]]),  wrap_elements(full = make_matrix("effect", ct)),
                  panel_head("B", PANEL_TITLE[["signedP"]]), wrap_elements(full = make_matrix("signedP", ct))),
             ncol = 1, heights = c(0.045, 1, 0.045, 1)) +
    plot_annotation(title = CELL_LAB[[ct]],
                    theme = theme(plot.title = element_text(face = "bold", size = 16, family = FONT, margin = margin(b = 2)),
                                  plot.margin = margin(8, 8, 6, 8)))
}

OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261102")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_W <- 8.5; FIG_H <- 17.6
out <- file.path(OUT_DIR, "Figure_S6.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H, family = FONT, onefile = TRUE)
for (ct in CELLS) { print(make_page(ct)); msg(ct, " done") }
invisible(dev.off())
msg("Saved: ", out)

###------------------------------------------------------------------------------------------------
### 4. 日志 + tsv
COR <- rbindlist(COR); RNG <- rbindlist(RNG)
RNG[, pct_shown := 100 * n_shown / n]
for (sc in c("effect", "signedP")) for (ct in CELLS) {
  cat(sprintf("\n===== Spearman r: %s, %s\n", sc, ct))
  m <- matrix(NA_real_, 5, 5, dimnames = list(METHOD_KEY[M5], METHOD_KEY[M5])); diag(m) <- 1
  z <- COR[scale == sc & cell_type == ct]
  m[cbind(match(z$method1, M5), match(z$method2, M5))] <- z$r
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  print(round(m, 3))
}
cat("\n===== 每个方法显示的 pair 比例 (%)\n")
print(dcast(RNG, scale + method1 ~ cell_type, value.var = "pct_shown")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])
tab <- rbind(COR[, .(type = "cor", scale, cell_type, method1, method2, value = round(r, 4), n, pct_shown = NA_real_)],
             RNG[, .(type = "range", scale, cell_type, method1, method2 = NA_character_, value = axis_limit, n, pct_shown = round(pct_shown, 2))])
fwrite(tab, file.path(OUT_DIR, "Figure_S6_cor.tsv"), sep = "\t")
msg("Saved: ", file.path(OUT_DIR, "Figure_S6_cor.tsv"))
