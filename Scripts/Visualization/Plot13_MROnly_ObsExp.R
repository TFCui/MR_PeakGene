# Plot13 (论文当前 S5 Fig): MR-Only link 的金标准命中率, 观测 (k / n_sig) vs 背景 (K / N), 分箱校准图
#   每个 gene 一行 (3 个细胞类型合并; 口径见 Scripts/Validation/Valid1_Hyper_MROnly.R); 只用 K > 0、N > K 的 gene (同 Fig3 / Table S1–S2)
#   MR-Only 在 gene 层很稀疏 (55% 的 gene 只有 1 个 MR-Only peak, K/N 中位 0.01–0.08), 单 gene 的 k/n_sig 几乎只取 0 / 1,
#   所以不画 2D 密度; 改为: 淡灰 gene 点作背景 + 按 K/N 等数量分 8 箱, 每箱画 (ΣK/ΣN, Σk/Σn_sig) 和 95% CI, 对角线 y = x
#   CI (2026-11-04 起): 箱内按 gene 有放回重抽 (bootstrap, B = 2000, 百分位法), 每次重算 Σk/Σn_sig; 保留同一 gene 内 peak 的相关性和 gene 间差异
#   (原来的精确二项 Clopper–Pearson CI 把每个 peak 当独立试验, 偏窄; 仍写在 Figure_S5_bins.tsv 作对照)
#   每个 panel 坐标按数据放大, x、y 同一范围且都用平方根刻度 (同一变换, 对角线仍是 45° 直线; 把挤在低值区的箱拉开); panel 副标题一行: XXX genes | Fold = (pooled fold) (P = 合并超几何 P)
#   淡蓝小点 = 落在 panel 范围内的 gene (y 很高的不画); 坐标两端各留 3% 空白, 使 x = 0 / y = 0 的点不压在边框上 (用带符号的平方根变换, 允许略小于 0)
#   至少命中 1 个的 gene 观测 vs 期望 %、对角线以上 % 只写 Figure_S5_summary.tsv (写 caption 用)
#   6 个 panel (2×3): GTEx / EnhancerAtlas / FANTOM5 / 4DGenome / HiChIP / PCHiC
# 输入: Results/Validation/Valid1_Hyper/MROnly_gene_overlap_combined.tsv
# 输出: Results/Figures/Figures_20261102/Figure_S5.pdf (只出 PDF) + Figure_S5_summary.tsv (每库汇总) + Figure_S5_bins.tsv (分箱)
#       tmp/Plot13_MROnly/hichip_fig3.rds (正文 Figure 3C 的 HiChIP panel: 薄荷色板, 标签条, 无副标题; 2026-11-04)
# 运行: Rscript Scripts/Visualization/Plot13_MROnly_ObsExp.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages(library(patchwork))

IN  <- file.path(VAL_DIR, "Valid1_Hyper/MROnly_gene_overlap_combined.tsv")
OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261102")
DB_ORDER <- c(GTEx = "GTEx", enhancerAtlas = "EnhancerAtlas", FANTOM5 = "FANTOM5", `4DGenome` = "4DGenome", HiChIP = "HiChIP", PCHiC = "PCHiC")
N_BIN <- 8
COL_BIN <- "#1F4E79"; COL_GENE <- "#8FB3D3"; COL_DIAG <- "#C44E52"; COL_ENR <- "#F3F7FB"   # 同一蓝色系: 分箱深蓝 / gene 淡蓝 / 富集侧极浅蓝; 对角线柔和红
N_BOOT <- 2000; set.seed(20261104)
BASEF <- 10; FONT <- "sans"
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")
fmt_lp <- function(lp) {   # log10 P -> "P = 2.3e-73" (不下溢)
  if (lp > -3) return(sprintf("P = %.3f", 10^lp))
  e <- floor(lp); sprintf("P = %.1fe%d", 10^(lp - e), e)
}

D <- fread(IN)[K > 0 & N > K]
D[, `:=`(x = K / N, y = k / n_sig, db = factor(DB_ORDER[database], levels = DB_ORDER),
         p_hit = phyper(0, K, N - K, n_sig, lower.tail = FALSE))]   # 期望: 至少命中 1 个的概率 (超几何, 精确)

### 每库汇总: n / 至少命中 1 个的 gene (观测 vs 期望) / pooled fold / 合并计数的超几何 P (log 计算)
S <- D[, {
  sN <- sum(N); sK <- sum(K); sn <- sum(n_sig); sk <- sum(k)
  lp <- phyper(sk - 1, sK, sN - sK, sn, lower.tail = FALSE, log.p = TRUE) / log(10)
  .(n_gene = .N, pct_nsig1 = 100 * mean(n_sig == 1), pct_hit_obs = 100 * mean(k > 0), pct_hit_exp = 100 * mean(p_hit),
    pct_above = 100 * mean(y > x), sum_N = sN, sum_K = sK, sum_nsig = sn, sum_k = sk,
    fold = (sk / sn) / (sK / sN), log10_p = lp)
}, keyby = db]
cat("\n===== 每个库汇总\n")
print(S[, .(db, n_gene, pct_nsig1 = round(pct_nsig1, 1), pct_hit_obs = round(pct_hit_obs, 1), pct_hit_exp = round(pct_hit_exp, 1),
            pct_above = round(pct_above, 1), fold = round(fold, 3), log10_p = round(log10_p, 1))])

### 分箱: 每个库内按 K/N 排序, 等数量分 N_BIN 箱; 每箱 ΣK/ΣN vs Σk/Σn_sig
### 95% CI = 按 gene 的 bootstrap (百分位法); 对照: 精确二项 CI (lo_binom / hi_binom)
boot_ci <- function(k, n, B = N_BOOT) {
  m <- length(k); idx <- matrix(sample.int(m, m * B, replace = TRUE), nrow = m)
  r <- colSums(matrix(k[idx], nrow = m)) / colSums(matrix(n[idx], nrow = m))
  unname(quantile(r, c(0.025, 0.975)))
}
D[, bin := cut(frank(x, ties.method = "first"), N_BIN, labels = FALSE), by = db]
B <- D[, {
  sk <- sum(k); sn <- sum(n_sig); cb <- binom.test(sk, sn)$conf.int; ci <- boot_ci(k, n_sig)
  .(n_gene = .N, bx = sum(K) / sum(N), by = sk / sn, lo = ci[1], hi = ci[2], lo_binom = cb[1], hi_binom = cb[2], sum_k = sk, sum_nsig = sn)
}, keyby = .(db, bin)]
cat("\n===== 分箱 (每库", N_BIN, "箱)\n")
print(B[, .(db, bin, n_gene, bx = round(bx, 4), by = round(by, 4), lo_boot = round(lo, 4), hi_boot = round(hi, 4),
            lo_binom = round(lo_binom, 4), hi_binom = round(hi_binom, 4), width_ratio = round((hi - lo) / (hi_binom - lo_binom), 2), ratio = round(by / bx, 2))])
cat(sprintf("\nbootstrap CI 下限 <= 背景比 (跨过对角线) 的箱: %d / %d; 精确二项: %d / %d\n",
            B[lo <= bx, .N], nrow(B), B[lo_binom <= bx, .N], nrow(B)))

theme_s5 <- function() theme_bw(base_size = BASEF, base_family = FONT) + theme(
  text = element_text(color = "black"), axis.text = element_text(color = "grey20", size = BASEF - 1),
  panel.border = element_rect(color = "grey35", linewidth = 0.4), panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
  panel.grid.minor = element_blank(), axis.ticks = element_line(color = "grey35", linewidth = 0.3),
  plot.title = element_text(face = "bold", size = BASEF + 1, hjust = 0, margin = margin(b = 2)),
  plot.subtitle = element_text(size = BASEF - 2, color = "grey25", hjust = 0, lineheight = 1.1, margin = margin(b = 4)),
  plot.margin = margin(4, 8, 4, 4), plot.tag = element_text(face = "bold", size = BASEF + 5))

ssqrt <- scales::trans_new("ssqrt", function(x) sign(x) * sqrt(abs(x)), function(x) sign(x) * x^2)   # 带符号平方根: 允许坐标略小于 0
EXP <- 0.03   # 两端留白 (变换后空间的比例)
sq_brk <- function(lim) if (lim > 0.3) c(0, 0.05, 0.1, 0.2, 0.4, 0.6) else c(0, 0.02, 0.05, 0.1, 0.2)   # 平方根刻度下低值区不挤
AX_LAB <- c(x = "Background overlap ratio (K / N)", y = "MR-Only overlap ratio (k / n_sig)")

COLS_S5   <- list(bin = COL_BIN, gene = COL_GENE, enr = COL_ENR, diag = COL_DIAG)
# Figure 3 的 HiChIP panel (2026-11-04 起薄荷色板): 分箱深蓝 / gene 薄荷绿 / 富集侧极浅薄荷 (比 S5 更透明) / 对角线灰
COLS_FIG3 <- list(bin = PAL_MINT[["m7"]], gene = PAL_MINT[["m3"]], enr = "#F5FAF3", diag = "grey45")
# fig3 = TRUE: 不写 title, 改用标签条 (theme_strip_mint); 边框 / 轴 / 分箱线条加粗
# 2026-11-04 起 S5 也用 Figure 3D 的风格 (fig3 = TRUE); cart = TRUE 只给 Figure 3D (coord_cartesian, 便于与 GRN 对齐), S5 保持 coord_equal (对角线 45°)
panel <- function(d, show_x = TRUE, show_y = TRUE, cols = COLS_S5, show_sub = TRUE, fig3 = FALSE, strip = d, cart = FALSE, ann_size = 4.4, ann_sep = " ") {   # ann_sep = "\n": S5 小 panel 中 Fold / P 分两行, 不压到对角线   # strip: Figure 3 的标签条文字
  lw <- if (fig3) list(line = 0.75, err = 0.7, pt = 2.9, stroke = 0.5, gene_a = 0.35) else list(line = 0.45, err = 0.5, pt = 2.3, stroke = 0.4, gene_a = 0.4)
  s <- D[db == d][, strip_lab := strip]; b <- B[db == d][, strip_lab := strip]; st <- S[db == d]   # strip_lab: Figure 3 的标签条
  lim <- max(b$hi, b$bx) * 1.06                                   # x、y 同一范围, 对角线 45°
  lo_l <- -(sqrt(lim) * EXP)^2; hi_l <- (sqrt(lim) * (1 + EXP))^2   # 留白后的坐标范围 (原始单位)
  enr <- data.table(x = c(lo_l, hi_l, lo_l), y = c(lo_l, hi_l, hi_l))   # 对角线以上 (富集一侧) 淡底色, 铺满留白
  sub <- if (show_sub) sprintf("%s genes | Fold = %.2f (%s)", fmt_n(st$n_gene), st$fold, fmt_lp(st$log10_p)) else NULL
  brk <- sq_brk(lim); brk <- brk[brk <= lim]
  if (fig3 && lim > 0.3) brk <- setdiff(brk, 0.05)   # 范围大的 panel (如 HiChIP): 去掉低值区挤在一起的 0.05 刻度
  ggplot() +
    geom_polygon(data = enr, aes(x, y), fill = cols$enr) +
    geom_point(data = s[x <= lim & y <= lim], aes(x, y), color = cols$gene, alpha = lw$gene_a, size = 0.8, shape = 16) +
    geom_abline(slope = 1, intercept = 0, color = cols$diag, linetype = "dashed", linewidth = 0.5) +
    geom_line(data = b, aes(bx, by), color = cols$bin, linewidth = lw$line) +
    geom_errorbar(data = b, aes(x = bx, ymin = lo, ymax = hi), color = cols$bin, width = sqrt(lim) * 0.025, linewidth = lw$err) +
    geom_point(data = b, aes(bx, by), shape = 21, fill = cols$bin, color = "white", size = lw$pt, stroke = lw$stroke) +
    scale_x_continuous(trans = ssqrt, breaks = brk, labels = scales::label_number(drop0trailing = TRUE)) +
    scale_y_continuous(trans = ssqrt, breaks = brk, labels = scales::label_number(drop0trailing = TRUE)) +
    # Figure 3: coord_cartesian (x、y 范围相同), 由 Plot9 的宽高比保持近似正方形; 固定纵横比在 patchwork 网格里会留空白、无法与 GRN 对齐
    (if (cart) coord_cartesian(xlim = c(lo_l, hi_l), ylim = c(lo_l, hi_l), expand = FALSE) else coord_equal(xlim = c(lo_l, hi_l), ylim = c(lo_l, hi_l), expand = FALSE)) +
    # Figure 3: 右下角 (对角线以下的空白处) 写 pooled fold 和合并超几何 P (与 S5 副标题同一口径; 2026-11-04)
    (if (fig3) annotate("label", fill = "white", label.size = 0, label.padding = unit(0.12, "lines"), x = (sqrt(hi_l) * 0.97)^2, y = (sqrt(hi_l) * 0.13)^2, hjust = 1, vjust = 0, size = ann_size, fontface = "bold", colour = "#1F2937",   # 加粗放大 (2026-11-04)
                        label = sprintf("Fold = %.2f%s(%s)", st$fold, ann_sep, fmt_lp(st$log10_p)), lineheight = 0.95) else NULL) +
    labs(title = if (fig3) NULL else d, subtitle = sub, x = if (show_x) AX_LAB[["x"]] else NULL, y = if (show_y) AX_LAB[["y"]] else NULL) +
    theme_s5() +
    if (fig3) list(facet_wrap(~ strip_lab), theme_strip_mint(11),
                   theme(panel.border = element_rect(color = COL_AX3, linewidth = 0.9), axis.ticks = element_line(color = COL_AX3, linewidth = 0.6),
                         axis.ticks.length = unit(1.6, "mm"), axis.text = element_text(color = COL_AX3, size = BASEF),
                         axis.title = element_text(face = "bold", color = "#1F2937", size = BASEF + 1), panel.grid.major = element_line(color = "grey93", linewidth = 0.3)))
}

dbs <- levels(D$db)   # 轴名只在最下一行 (x) 和最左一列 (y) 显示
# S5 (2026-11-04): Figure 3D 风格 (薄荷色板、标签条写数据库名、加粗边框、右下角 Fold / P; 不再写副标题); gene 数只写 summary tsv
fig <- wrap_plots(lapply(seq_along(dbs), function(i) panel(dbs[i], show_x = i > 3, show_y = i %in% c(1, 4), cols = COLS_FIG3, show_sub = FALSE,
                                                           fig3 = TRUE, ann_size = 3.5, ann_sep = "\n")), nrow = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "A") & theme(plot.tag.position = c(0, 1), plot.tag = element_text(face = "bold", size = BASEF + 6, hjust = 0, vjust = 1))

FIG_W <- 10; FIG_H <- 7.2
out <- file.path(OUT_DIR, "Figure_S5.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H, family = FONT)
print(fig)
invisible(dev.off())
msg("Saved: ", out)

### Figure 3 用的 HiChIP panel (2026-11-04): 薄荷色板, 标签条 "HiChIP", 线条加粗, 不写副标题数字; 存 rds 由 Plot9_GO_GRN.R 拼图 (改这里后先跑本脚本, 再跑 Plot9)
TMP13 <- file.path(TMP_DIR, "Plot13_MROnly"); dir.create(TMP13, recursive = TRUE, showWarnings = FALSE)
saveRDS(panel("HiChIP", cols = COLS_FIG3, show_sub = FALSE, fig3 = TRUE, strip = "MR-only link overlap (HiChIP)", cart = TRUE), file.path(TMP13, "hichip_fig3.rds"))
msg("Saved: ", file.path(TMP13, "hichip_fig3.rds"))

fwrite(S[, .(database = db, n_gene, pct_nsig1 = round(pct_nsig1, 2), pct_hit_obs = round(pct_hit_obs, 2), pct_hit_exp = round(pct_hit_exp, 2),
             pct_above = round(pct_above, 2), sum_N, sum_K, sum_nsig, sum_k, fold = round(fold, 4), log10_p = round(log10_p, 2))],
       file.path(OUT_DIR, "Figure_S5_summary.tsv"), sep = "\t")
fwrite(B[, .(database = db, bin, n_gene, bg_ratio = signif(bx, 4), obs_ratio = signif(by, 4), ci_lo = signif(lo, 4), ci_hi = signif(hi, 4),
             ci_lo_binom = signif(lo_binom, 4), ci_hi_binom = signif(hi_binom, 4), sum_k, sum_nsig)],
       file.path(OUT_DIR, "Figure_S5_bins.tsv"), sep = "\t")
msg("Saved: Figure_S5_summary.tsv, Figure_S5_bins.tsv")
