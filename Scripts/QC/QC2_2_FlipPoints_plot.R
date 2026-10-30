#!/usr/bin/env Rscript
# ==============================================================================
# QC2_2_FlipPoints_plot: MR 中 "flip (让 caQTL 效应 beta_x > 0)" 前后的对比图 (流程见 Scripts/Pipeline.md 第 4.2 节)
# pair: 下面 SEL 的 6 个 (每个 CT 各 2 个):
#   flip 前后差别小 = Figure 2 / Fig2 的 example (SIDT2 Bn, RPL11 CD4, LCLAT1 CD8)
#   flip 前后差别大 = 全基因组 flip 扫描 (QC2_1_FlipScan.R) 中 "Egger 不 flip 极显著, flip 后不显著" 且 SNP 数多 (>= 50) 的 pair (PIGX Bn 136, PDE4DIP CD4 74, KLRC3 CD8 134);
#                     候选表见 QC2_1_FlipScan_select.R 输出的 FlipScan_candidates_nsnp50.tsv
# 每个 pair: QTL 对齐 -> 去重 -> PLINK LD (只做一次), 然后
#   flip = TRUE : s_x = sign(beta_x); bx, by, R 同乘符号 (与 Scripts/Model/MR_AllPairs.R 相同)
#   flip = FALSE: 原始等位编码, 不翻转
#   分别拟合 GLS IVW / Egger (函数拷自 Plot2_examples_v2.R, 与 MR_AllPairs.R 一致)
# 输出: Results/QC/2_FlipPoints/FigureS2_estimates.tsv (pair x flip x model 的 beta / se / p / Egger 截距)
#       Results/Figures/Figures_20261102/Figure_S2.pdf (只出 PDF; 单页, 3 行 = Bn / CD4 T cell / CD8 T cell, 每行 4 个 panel;
#         左半 = panel A = 差别小的 pair [No flipping | Flipping], 右半 = panel B = 差别大的 pair [No flipping | Flipping];
#         同一 pair 的两个 panel 坐标范围相同; panel 内只写 IVW P / Egger P; 所有 SNP 同样画红三角, 不单独标 lead SNP)
# 运行: sbatch Scripts/QC/QC2_2_FlipPoints_plot_slurm.sh
# ==============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrepel) })
options(width = 250)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ---------------------------- 路径 ---------------------------------------------
CIMA_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(CIMA_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
OUT_DIR  <- file.path(PROJ, "Results/QC/2_FlipPoints")
FIG_DIR  <- file.path(PROJ, "Results/Figures/Figures_20261102")
PLINK    <- file.path(CIMA_DIR, "R_package/my_plink/plink")
REF_DIR  <- file.path(PROJ, "data/ref_chr/hg38")
LD_DIR   <- file.path(PROJ, "tmp/QC2_FlipPoints", paste0("ld_", Sys.getpid()))   # PLINK 临时文件, 结束自动删除
OMIX     <- file.path(CIMA_DIR, "Data/OMIX009786")
for (d in c(OUT_DIR, FIG_DIR, LD_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

eqtl_file        <- function(ct) file.path(OMIX, "cis_eQTL_Summary_Statics/cis_eQTL", paste0(ct, ".txt.gz"))
caqtl_clump_file <- function(ct) file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL/LD_clump", paste0(ct, "_Clumped_r2_0.99_kb_20000.txt.gz"))

CELLTYPES <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")

# MR 固定参数 (与 MR_AllPairs.R 一致)
MIN_SNP <- 3; DEDUP_TOL <- 1e-3; KEEP_ALLELE <- TRUE
EIG_TRUNCATE <- TRUE; EIG_VAR_KEEP <- 0.999; EIG_MIN_KEEP <- 5; EIG_TOL <- 1e-10; P_FLOOR <- 1e-300; SIGMA_MODE <- "over"

# ---------------------------- 选定的 pairs --------------------------------------
SEL <- data.table(
  group     = rep(c("small", "large"), each = 3),    # flip 前后 Egger 差别小 / 大; 只用于排序 (图里不写)
  cell_type = rep(CELLTYPES, 2),
  gene_id   = c("SIDT2", "RPL11", "LCLAT1", "PIGX", "PDE4DIP", "KLRC3"),
  peak_id   = c("chr11:117185715-117186216", "chr1:23691494-23691995", "chr2:30447755-30448256",
                "chr3:195747820-195748321", "chr1:149608426-149608927", "chr12:10553007-10553508"))
msg("pairs: ", nrow(SEL))

# ---------------------------- MR 函数 (拷自 Plot2_examples_v2.R) ----------------
dedup_keep <- function(x, y, tol = DEDUP_TOL) {
  keep <- logical(length(x)); ref_x <- ref_y <- numeric(0)
  for (j in seq_along(x)) {
    if (length(ref_x) && any(abs(ref_x - x[j]) <= tol & abs(ref_y - y[j]) <= tol)) next
    keep[j] <- TRUE; ref_x <- c(ref_x, x[j]); ref_y <- c(ref_y, y[j])
  }
  keep
}

make_plink_id <- function(variant_id, ref, alt) {
  chr <- sub("^chr", "", sub("_.*$", "", variant_id)); bp <- sub("^[^_]*_", "", variant_id)
  chr <- fcase(chr == "X", "23", chr == "Y", "24", chr %chin% c("M", "MT"), "26", default = chr)
  paste(chr, bp, ref, alt, sep = ":")
}
ref_bfile <- function(peak) file.path(REF_DIR, sub(":.*$", "", peak))

run_plink_ld <- function(ids, bfile, pre) {
  fs <- paste0(pre, c(".snps.txt", ".ld", ".snplist", ".log", ".nosex")); on.exit(unlink(fs), add = TRUE)
  fwrite(data.table(unique(ids[!is.na(ids)])), fs[1], col.names = FALSE)
  kao <- if (KEEP_ALLELE) "--keep-allele-order" else ""
  cmd <- paste(PLINK, "--bfile", bfile, "--extract", fs[1], kao, "--r square --write-snplist --threads 1 --memory 4000 --out", pre, "--silent")
  st  <- suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))
  if (st != 0 || !file.exists(fs[2]) || !file.exists(fs[3])) return(NULL)
  snps <- fread(fs[3], header = FALSE)[[1]]
  R <- as.matrix(fread(fs[2], header = FALSE, na.strings = c("nan", "-nan", "NaN", "NA"))); storage.mode(R) <- "numeric"
  dimnames(R) <- list(snps, snps); R
}

run_gls_mr <- function(bx, by, se_y, R, intercept = FALSE, eg = NULL) {
  n <- length(bx)
  if (is.null(eg)) eg <- eigen(outer(se_y, se_y) * R, symmetric = TRUE)
  lam <- eg$values; V <- eg$vectors
  pos  <- lam > EIG_TOL * lam[1]
  cum  <- cumsum(pmax(lam, 0)) / sum(pmax(lam, 0))
  keep <- if (EIG_TRUNCATE) (c(TRUE, cum[-n] < EIG_VAR_KEEP) | seq_len(n) <= EIG_MIN_KEEP) & pos else pos
  k <- sum(keep); p <- if (intercept) 2L else 1L; df <- n - p
  if (k < p || df < 1) stop(sprintf("k = %d, df = %d", k, df))
  d_half <- numeric(n); d_half[keep] <- lam[keep]^(-0.5)
  W  <- V %*% diag(d_half, n) %*% t(V)
  X  <- if (intercept) cbind(alpha = 1, beta = bx) else cbind(beta = bx)
  ys <- as.numeric(W %*% by); Xs <- W %*% X
  fit <- lm(ys ~ Xs - 1); est <- unname(coef(fit))
  se_fixed <- sqrt(diag(solve(crossprod(Xs)))); sigma <- sqrt(sum(residuals(fit)^2) / df)
  infl <- switch(SIGMA_MODE, over = max(1, sigma), under = 1 / min(1, sigma), none = 1)
  se <- se_fixed * infl; pv <- 2 * pt(abs(est / se), df = df, lower.tail = FALSE)
  data.table(model = if (intercept) "Egger" else "IVW", n_snp = n, beta = est[p], se_beta = se[p], p_beta = pv[p],
             alpha = if (intercept) est[1] else NA_real_, p_alpha = if (intercept) pv[1] else NA_real_)
}

# 对齐 -> 去重 -> LD 只做一次; 然后 flip = FALSE / TRUE 各拟合一次 IVW + Egger
run_pair_both <- function(gene, peak, ca_data, eq_data) {
  pid <- paste(peak, gene, sep = "__")
  ca <- ca_data[.(peak, gene), nomatch = 0L]; eq <- eq_data[.(gene), nomatch = 0L]
  if (!nrow(ca) || !nrow(eq)) return(list(status = "no_qtl"))
  ca[, snp_id := make_plink_id(variant_id, `A2(REF)`, `A1(ALT/effect allele)`)]
  ca <- ca[, .(variant_id, snp_id, ea_x = `A1(ALT/effect allele)`, oa_x = `A2(REF)`, beta_x = slope, se_x = slope_se, p_x = pval_nominal)]
  eq <- eq[, .(variant_id,          ea_y = `A1(ALT/effect allele)`, oa_y = `A2(REF)`, beta_y = slope, se_y = slope_se, p_y = pval_nominal)]
  dat <- merge(ca, eq, by = "variant_id")[ea_x == ea_y & oa_x == oa_y]
  dat <- dat[dedup_keep(beta_x, beta_y)]
  if (nrow(dat) < MIN_SNP) return(list(status = "few_snp_qtl"))
  ld <- run_plink_ld(dat$snp_id, ref_bfile(peak), file.path(LD_DIR, gsub("[^A-Za-z0-9]", "_", pid)))
  if (is.null(ld)) return(list(status = "plink_fail"))
  d <- dat[match(rownames(ld), snp_id)][!is.na(variant_id)]
  R <- ld[d$snp_id, d$snp_id, drop = FALSE]
  if (anyNA(R)) { bad <- which(apply(R, 1, function(x) any(!is.finite(x)))); d <- d[-bad]; R <- R[-bad, -bad, drop = FALSE] }
  if (nrow(d) < MIN_SNP) return(list(status = "few_snp_ld"))

  fits <- lapply(c(FALSE, TRUE), function(flip) {
    s_x <- if (flip) ifelse(d$beta_x < 0, -1, 1) else rep(1, nrow(d))
    bx <- d$beta_x * s_x; by <- d$beta_y * s_x; R_use <- R * outer(s_x, s_x)
    eg <- eigen(outer(d$se_y, d$se_y) * R_use, symmetric = TRUE)
    res <- rbindlist(lapply(c(FALSE, TRUE), function(itc) tryCatch(run_gls_mr(bx, by, d$se_y, R_use, intercept = itc, eg = eg), error = function(e) NULL)))
    list(res = res[, flipped := flip], pts = data.table(variant_id = d$variant_id, x = bx, y = by, se_x = d$se_x, se_y = d$se_y))
  })
  list(status = "ok", fits = fits, n_snp = nrow(d))
}

read_filtered <- function(file, keys) {           # 只读第一列在 keys 里的行
  kf <- tempfile(tmpdir = LD_DIR, fileext = ".txt"); writeLines(unique(keys), kf); on.exit(unlink(kf), add = TRUE)
  fread(cmd = sprintf("zcat %s | awk 'NR==FNR{a[$1]=1; next} FNR==1 || ($1 in a)' %s -", shQuote(file), shQuote(kf)))
}

# ---------------------------- 主循环: 按细胞类型 ---------------------------------
OBJ <- list(); EST <- list()
for (CT in intersect(CELLTYPES, SEL$cell_type)) {
  msg("==================== ", CT, " ====================")
  pc <- SEL[cell_type == CT]
  eq_sub <- read_filtered(eqtl_file(CT), unique(pc$gene_id))[, .(phenotype_id, variant_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se, pval_nominal)]
  setkey(eq_sub, phenotype_id, variant_id)
  ca_sub <- fread(caqtl_clump_file(CT))[unique(pc[, .(phenotype_id = peak_id, gene_id)]), on = .(phenotype_id, gene_id), nomatch = 0L,
                                        .(phenotype_id, variant_id, gene_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se, pval_nominal)]
  setkey(ca_sub, phenotype_id, gene_id)
  mr <- fread(file.path(MR_DIR, paste0(CT, "_AllPairs_results.tsv.gz")))[pc, on = .(gene_id, peak_id), nomatch = NULL]

  for (i in seq_len(nrow(pc))) {
    g <- pc$gene_id[i]; pk <- pc$peak_id[i]; key <- paste(CT, g, pk, sep = "|")
    o <- run_pair_both(g, pk, ca_sub, eq_sub)
    msg(sprintf("  %s__%s: %s", pk, g, o$status))
    if (o$status != "ok") next
    OBJ[[key]] <- c(o, list(cell_type = CT, gene_id = g, peak_id = pk, group = pc$group[i]))
    est <- rbindlist(lapply(o$fits, `[[`, "res"))[, `:=`(cell_type = CT, gene_id = g, peak_id = pk)]
    # 对账: flip 版与正式 MR 表一致
    f <- mr[gene_id == g & peak_id == pk]
    est[flipped == TRUE & model == "IVW",   formal_beta := f$ivw_beta]
    est[flipped == TRUE & model == "Egger", formal_beta := f$egger_beta]
    EST[[key]] <- est
  }
  rm(eq_sub, ca_sub); invisible(gc())
}
EST <- rbindlist(EST)
setcolorder(EST, c("cell_type", "gene_id", "peak_id", "flipped", "model", "n_snp", "beta", "se_beta", "p_beta", "alpha", "p_alpha", "formal_beta"))

cat("\n===== 对账: flip 版 beta 与正式 MR 表的最大差异\n")
print(EST[flipped == TRUE, .(max_abs_diff = max(abs(beta - formal_beta))), by = model])
cat("\n===== flip 前后对比 (不 flip 减 flip)\n")
cmp <- dcast(EST, cell_type + gene_id + peak_id + model ~ flipped, value.var = c("beta", "p_beta", "alpha"))
cmp[, `:=`(d_beta = beta_FALSE - beta_TRUE, ratio_p = p_beta_FALSE / p_beta_TRUE)]
print(cmp[, .(max_abs_d_beta = max(abs(d_beta)), median_p_ratio = median(ratio_p)), by = model])
print(cmp[model == "Egger", .(cell_type, gene_id, beta_noflip = signif(beta_FALSE, 3), beta_flip = signif(beta_TRUE, 3),
                              p_noflip = signif(p_beta_FALSE, 2), p_flip = signif(p_beta_TRUE, 2),
                              int_noflip = signif(alpha_FALSE, 2), int_flip = signif(alpha_TRUE, 2))])
fwrite(EST[, formal_beta := NULL], file.path(OUT_DIR, "FigureS2_estimates.tsv"), sep = "\t")
msg("Saved: ", file.path(OUT_DIR, "FigureS2_estimates.tsv"))

# ---------------------------- 作图 -----------------------------------------------
COL_IVW <- "#E64B35"; COL_EGG <- "#3C5488"                    # 方法配色同 00_common.R
COL_SNP <- "#DC0000"; COL_ERR <- "#C8D1D9"; COL_REF <- "grey60"; SEP_COL <- "#C8CDD3"
BASE <- 12
fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p)

# lim = list(x, y): 同一个 pair 的两张图 (不 flip / flip) 用相同坐标范围; 只画散点和拟合线 (文字放在 panel_header, 避免压住点和线)
mr_core <- function(fit, lim, show_x = TRUE, show_y = TRUE) {
  pts <- fit$pts
  r <- fit$res; iv <- r[model == "IVW"]; eg <- r[model == "Egger"]
  lines <- data.table(series = factor(c("Causal IVW", "Causal Egger"), levels = c("Causal IVW", "Causal Egger")),
                      slope = c(iv$beta, eg$beta), intercept = c(0, eg$alpha))
  ggplot(pts, aes(x, y)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = COL_REF, linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = COL_REF, linewidth = 0.4) +
    geom_errorbar(aes(ymin = y - 1.96 * se_y, ymax = y + 1.96 * se_y), color = COL_ERR, width = 0, linewidth = 0.4) +
    geom_errorbar(aes(xmin = x - 1.96 * se_x, xmax = x + 1.96 * se_x), orientation = "y", color = COL_ERR, width = 0, linewidth = 0.4) +
    geom_abline(data = lines, aes(slope = slope, intercept = intercept, color = series), linewidth = 1.3, alpha = 0.8) +
    geom_point(shape = 24, fill = COL_SNP, color = "white", size = 2.8, stroke = 0.4) +
    scale_color_manual(values = c(`Causal IVW` = COL_IVW, `Causal Egger` = COL_EGG), name = NULL) +
    coord_cartesian(xlim = lim$x, ylim = lim$y) +
    labs(x = if (show_x) "caQTL effect size" else NULL, y = if (show_y) "eQTL effect size" else NULL) +
    theme_classic(base_size = BASE) +
    theme(axis.text = element_text(color = "black", size = BASE - 2),
          legend.position = "none", legend.key.width = unit(10, "mm"), legend.text = element_text(size = BASE),
          plot.margin = margin(2, 8, 2, 4))
}
# panel 上方三行: 标题 / IVW P (红) / Egger P (蓝)
panel_header <- function(ttl, p_ivw, p_egg) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.83, label = ttl, fontface = "bold", size = BASE / .pt, color = "#1F2937") +
    annotate("text", x = 0.5, y = 0.50, label = paste0("IVW P = ", fmt_p(p_ivw)), fontface = "bold", size = (BASE - 1) / .pt, color = COL_IVW) +
    annotate("text", x = 0.5, y = 0.17, label = paste0("Egger P = ", fmt_p(p_egg)), fontface = "bold", size = (BASE - 1) / .pt, color = COL_EGG) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") + theme_void(base_size = BASE) +
    theme(plot.margin = margin(0, 4, 0, 4))
}
mr_panel <- function(fit, ttl, lim, show_x = TRUE, show_y = TRUE) {
  r <- fit$res
  (panel_header(ttl, r[model == "IVW"]$p_beta, r[model == "Egger"]$p_beta) / mr_core(fit, lim, show_x, show_y)) + plot_layout(heights = c(0.3, 1))
}

# 页顶只放一行文字的窄轨: Gene: XXX  Peak: chr:start-end (与 Fig 4 同风格)
text_track <- function(txt, parse = TRUE, size = BASE, tag = NULL) {   # tag = panel 字母 (A / B, 只在第一行), 画在左侧
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = txt, parse = parse, color = "#1F2937", size = size / .pt) +
    { if (!is.null(tag)) annotate("text", x = 0, y = 0.5, label = tag, hjust = 0, fontface = "bold", size = (BASE + 8) / .pt, color = "black") } +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") + theme_void(base_size = BASE) +
    theme(plot.margin = margin(2, 4, 2, 4))
}

# 一个 pair = 标题行 + 2 个 panel (不 flip | flip); 坐标范围取两版的并集 (含 0 和 95% CI)
build_pair <- function(o, show_x, tag = NULL) {
  P <- rbindlist(lapply(o$fits, `[[`, "pts"))
  rng <- function(v, s) { r <- range(c(0, v - 1.96 * s, v + 1.96 * s)); r + c(-1, 1) * 0.05 * diff(r) }
  lim <- list(x = rng(P$x, P$se_x), y = rng(P$y, P$se_y))
  hd  <- text_track(sprintf('"Gene: "~bold("%s")~~"Peak: "~bold("%s")', o$gene_id, o$peak_id), tag = tag)
  pn  <- mr_panel(o$fits[[1]], "No flipping", lim, show_x, TRUE) | mr_panel(o$fits[[2]], "Flipping", lim, show_x, FALSE)
  (hd / pn) + plot_layout(heights = c(0.11, 1))
}

# 一行 = 一个细胞类型: 行标签 | 差别小的 pair | 差别大的 pair (两个 pair 之间一条浅灰竖线); 图例统一放页底
row_label <- function(ct) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = CT_SHORT[[ct]], angle = 90, fontface = "bold", size = (BASE + 3) / .pt, color = "#1F2937") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()
}
CT_SHORT <- local({ source(file.path(PROJ, "Scripts/Visualization/00_common.R"), local = TRUE); CELL_LAB })   # Bn / CD4 T cell / CD8 T cell
with_sep <- function(p) {
  g <- grid::grobTree(patchwork::patchworkGrob(p), grid::segmentsGrob(x0 = 0, x1 = 0, y0 = 0.02, y1 = 0.98, gp = grid::gpar(col = SEP_COL, lwd = 0.9)))
  wrap_elements(full = g)
}
build_row <- function(ct, show_x, tags) {   # tags = 这一行两个 pair 的字母 (左 = A 差别小, 右 = B 差别大; 只有第一行有, 其余为 NULL)
  pr <- function(grp, tg = NULL) { k <- SEL[cell_type == ct & group == grp, paste(cell_type, gene_id, peak_id, sep = "|")]; if (is.null(OBJ[[k]])) plot_spacer() else build_pair(OBJ[[k]], show_x, tg) }
  wrap_plots(row_label(ct), wrap_elements(full = patchwork::patchworkGrob(pr("small", tags[1]))), with_sep(pr("large", tags[2])), nrow = 1, widths = c(0.035, 1, 1))
}

legend_grob <- function() {
  p <- mr_core(OBJ[[1]]$fits[[2]], list(x = c(0, 1), y = c(0, 1))) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  wrap_elements(full = cowplot_free_legend(p))
}
cowplot_free_legend <- function(p) {   # ggplot2 >= 3.5: 图例 grob 名为 guide-box-<位置>, 未用的位置是 zeroGrob
  g <- ggplotGrob(p); idx <- grep("^guide-box", g$layout$name)
  g$grobs[[idx[!vapply(g$grobs[idx], inherits, logical(1), "zeroGrob")][1]]]
}

build_page <- function() {
  cts <- intersect(CELLTYPES, SEL$cell_type)
  rows <- lapply(seq_along(cts), function(i) wrap_elements(full = patchwork::patchworkGrob(build_row(cts[i], show_x = (i == length(cts)), tags = if (i == 1) c("A", "B") else NULL))))
  wrap_plots(c(rows, list(legend_grob())), ncol = 1, heights = c(rep(1, length(cts)), 0.07))
}

PAGE_W <- 13; PAGE_H <- 9.9
pdf_file <- file.path(FIG_DIR, "Figure_S2.pdf")
cairo_pdf(pdf_file, width = PAGE_W, height = PAGE_H)
print(build_page(), newpage = FALSE)   # newpage = FALSE: 构图时 ggplotGrob 会先在设备上开一页, 否则多出空白页
invisible(dev.off())
msg("Saved: ", pdf_file)

unlink(LD_DIR, recursive = TRUE)
msg("All done.")
