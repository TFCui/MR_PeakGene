#!/usr/bin/env Rscript
# ==============================================================================
# 2026.10.07  Fig2 v2: 示例 peak-gene pairs, 每个 pair 一页, 每个细胞类型一个 PDF (Fig2_examples_<CT>_v2.pdf)
# 由 v1 (Plot2_scripts.R) 升级: 左 2/3 = 2x2 证据面板, 右 1/3 = 位点 link 图 (来自 Plot8, 只突出当前 pair)
# 版面 (2 行 x 3 列):
#   a (左上) = 两条轨: caQTL | eQTL (共用 xlim = peak/TSS 中点 ±500 kb); 红色三角 = MR 用到的 SNP; 金色三角 = SMR lead SNP
#   b (中上) = MR (flip 后散点 + Causal SMR / IVW / Egger 线, 图例在图内左上角)
#   c (左下) = Marginal Multiomics: WNN metacell Spearman
#   d (中下) = Marginal Pseudobulk: 个体层 Spearman
#   e (右侧, 跨两行) = 位点图: Gene / Peaks (当前 pair 的 peak 红色, 其他显著 peak 黑色) / MR SNPs (当前 pair) /
#       6 个方法各一行的显著 link (FDR < 0.05; 颜色 = 方向, 线宽 = -log10 FDR; 无图例, 说明见页面 caption)
# 数值 (p / FDR / beta) 全部取自 MR_results_formal; SNP 层数据由 run_one_pair 重跑得到并与全量结果对账
# panel b-d 的 subtitle 只显示 p; FDR < 0.1 的结果仍用该方法的颜色加粗显示 (需 ggtext; 未安装时不高亮)
#
# Usage:
#   Rscript Plot2_scripts_v2.R                                 # 默认: Bn_TCL1A, 前 3 个 example (试跑)
#   Rscript Plot2_scripts_v2.R --celltype all --n 0            # 三个细胞类型, 全部 example
#   可选: --outdir DIR
# ==============================================================================
# ------------------------------------------------------------------------------
# Claude_project 移植版 (2026-10-09): 由 20261006_paper/Scripts/Source_visual/Plot2_scripts_v2.R 复制, 只改:
#   输入 = Results/MR/<CT>_AllPairs_results.tsv.gz (已过滤, 不再做 sigma2 过滤); LD 参考 = data/ref_chr/hg38
#   example pairs 改为 source 01_example_pairs.R; 输出 = Results/Figures; PLINK / 临时文件 = tmp/plink_fig2/
#   新 MR 表无 egger_int 列 -> 用重跑得到的 Egger 截距画 Egger 线, 并核对重跑 Egger 斜率与表一致
# 运行: sbatch Scripts/Visualization/Plot2_examples_slurm.sh (array 0-2, 每个 CT 依次跑 v1 和 v2)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(GenomicRanges); library(EnsDb.Hsapiens.v86); library(Signac)
})
options(width = 300)

# ---------------------------- 参数 (试跑 / 全量 只改这里或命令行) -------------
args    <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) { i <- match(paste0("--", name), args); if (is.na(i) || i == length(args)) default else args[i + 1] }
CT_ARG     <- get_arg("celltype", "Bn_TCL1A")      # "Bn_TCL1A" / "CD4_Tn_CCR7" / "CD8_Tn_CCR7" / "all"
N_EXAMPLES <- as.integer(get_arg("n", "3"))        # 每个细胞类型画前 N 个 example; 0 = 全部
FIG_DIR    <- get_arg("outdir", "/scratch/user/tfcui2025/CIMA/Claude_project/Results/Figures")

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ---------------------------- 路径 ---------------------------------------------
CIMA_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(CIMA_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
PLINK    <- file.path(CIMA_DIR, "R_package/my_plink/plink")
REF_DIR  <- file.path(PROJ, "data/ref_chr/hg38")
LD_DIR   <- file.path(PROJ, "tmp/plink_fig2", paste0("ld_fig2_", Sys.getpid()))   # PLINK / 临时文件放项目 tmp/, 结束自动删除
OMIX     <- file.path(CIMA_DIR, "Data/OMIX009786")
PB_DIR   <- file.path(OMIX, "cell_rna/QC_pseudobulk")
WNN_FILE    <- file.path(CIMA_DIR, "Data/Valid_Multiome/1_pbmc_multiomics.reprocessed.rda")
WNN_CT_FILE <- file.path(CIMA_DIR, "Data/Valid_Multiome/ATAC_RNA_WNN.rda")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LD_DIR,  recursive = TRUE, showWarnings = FALSE)

eqtl_file        <- function(ct) file.path(OMIX, "cis_eQTL_Summary_Statics/cis_eQTL", paste0(ct, ".txt.gz"))
caqtl_clump_file <- function(ct) file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL/LD_clump", paste0(ct, "_Clumped_r2_0.99_kb_20000.txt.gz"))
caqtl_full_file  <- function(ct) file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL", paste0(ct, ".txt.gz"))

CELLTYPES <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
SIGMA2_MAX <- 10; FDR_CUT <- 0.1

# MR 固定参数 (与 MR_results_formal 一致)
MIN_SNP <- 3; DEDUP_TOL <- 1e-3; KEEP_ALLELE <- TRUE; FLIP_POINTS <- TRUE
EIG_TRUNCATE <- TRUE; EIG_VAR_KEEP <- 0.999; EIG_MIN_KEEP <- 5; EIG_TOL <- 1e-10; P_FLOOR <- 1e-300; SIGMA_MODE <- "over"

# ---------------------------- Example pairs ------------------------------------
source(file.path(PROJ, "Scripts/Visualization/01_example_pairs.R"))   # CATS / EXAMPLE_PAIRS / prep_examples
EXAMPLES <- lapply(EXAMPLE_PAIRS, prep_examples)

# ---------------------------- 配色 / 主题 (方法名 / 配色全文统一, 见 CLAUDE.md) ---
pal_method <- c(`Causal IVW` = "#E64B35", `Causal Egger` = "#3C5488", `Causal SMR` = "#F0A14A",
                `Causal Joint` = "#8491B4", `Marginal Multiomics` = "#4DBBD5", `Marginal Pseudobulk` = "#00A087")
MR_LINES   <- c("Causal SMR", "Causal IVW", "Causal Egger")
COL_PEAK  <- "#6F99AD"     # caQTL 点 / peak 方块
COL_TSS   <- "#B09C85"     # eQTL 点 / TSS 三角
COL_MRSNP <- "#DC0000"     # MR 用到的 SNP (红色三角)
LEAD_FILL <- "#FFD34E"; LEAD_BORDER <- "#5C4A00"
COL_ERR   <- "#C8D1D9"; COL_REF <- "grey60"; COL_NS <- "grey45"

WIN_HALF     <- 0.5        # caQTL / eQTL 轨: peak-TSS 中点 ± 0.5 Mb

# panel e (位点 link 图, 来自 Plot8)
LINK_FDR_CUT <- 0.05                 # link 显著阈值 (与 Plot8 一致; 与 subtitle 用的 FDR_CUT 分开)
SCORE_CAP    <- 20                   # -log10(FDR) 上限, 线宽标尺三个 CT 共用
LK_PAD_FRAC  <- 0.08; LK_PAD_MIN <- 1e4; LK_EMPTY_HALF <- 5e4
LW_RANGE     <- c(0.25, 1.0)         # link 线宽范围 (比 Plot8 细)
DIR_COL      <- c("Positive" = "#D7191C", "Negative" = "#2C7BB6", "No direction" = "grey15")
GENE_COL     <- "#1B1B1B"; PEAK_CUR <- "#E64B35"; PEAK_OTH <- "grey25"; SNP_COL <- "#3C5488"
LINK_METHODS <- data.table(method  = c("Causal SMR", "Causal IVW", "Causal Egger", "Causal Joint", "Marginal Multiomics", "Marginal Pseudobulk"),
                           fdr_col = c("smr_fdr", "ivw_fdr", "egger_fdr", "cct_fdr", "wnn_fdr", "pseudo_fdr"),
                           eff_col = c("smr_b", "ivw_beta", "egger_beta", NA, "wnn_rho", "pseudo_rho"))

PAGE_W <- 21; PAGE_H <- 12
BASE   <- 11
ANNOT  <- 3.4                                      # annotate 字号 (mm), 约 9.5 pt
PT_MR  <- 3.6; PT_LEAD <- 5.2; PT_QTL <- 1.5; PT_MARG <- 2.7; LWD <- 1.1
SUB_SIZE <- BASE + 1                               # panel b-d subtitle 字号 (想再大就改成 BASE + 2)

USE_MD <- requireNamespace("ggtext", quietly = TRUE)
if (!USE_MD) msg("提示: 未安装 ggtext, subtitle 不做颜色高亮")
sub_el <- if (USE_MD) ggtext::element_markdown(size = SUB_SIZE, color = "grey20", lineheight = 1.25, margin = margin(t = 1, b = 5)) else
                      element_text(size = SUB_SIZE, color = "grey20", lineheight = 1.25, margin = margin(t = 1, b = 5))

legend_inside <- function(pos = c(0.02, 0.98)) {
  base <- theme(legend.justification = c(0, 1), legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
                legend.key = element_blank(), legend.title = element_blank(), legend.text = element_text(size = BASE - 1),
                legend.key.width = unit(9, "mm"), legend.key.height = unit(4, "mm"), legend.margin = margin(2, 4, 2, 2))
  if (utils::packageVersion("ggplot2") >= "3.5.0") base + theme(legend.position = "inside", legend.position.inside = pos) else base + theme(legend.position = pos)
}
mb_labels <- function(w) { acc <- min(0.1, 10^floor(log10(w / 4))); label_number(accuracy = acc, big.mark = "") }

theme_fig <- theme_classic(base_size = BASE) + theme(
  text = element_text(color = "black"),
  plot.title = element_text(face = "bold", size = BASE + 1, hjust = 0, margin = margin(b = 3)), plot.subtitle = sub_el,
  axis.title = element_text(size = BASE), axis.text = element_text(size = BASE - 1, color = "black"),
  axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4), axis.ticks.length = unit(1.2, "mm"),
  legend.position = "none", plot.margin = margin(4, 8, 4, 4))

# ---------------------------- MR 函数 (与 MR_results_formal 一致) ---------------
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

run_plink_ld <- function(ids, bfile, pre, keep_allele_order = KEEP_ALLELE) {
  fs <- paste0(pre, c(".snps.txt", ".ld", ".snplist", ".log", ".nosex")); on.exit(unlink(fs), add = TRUE)
  fwrite(data.table(unique(ids[!is.na(ids)])), fs[1], col.names = FALSE)
  kao <- if (keep_allele_order) "--keep-allele-order" else ""
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
             alpha = if (intercept) est[1] else NA_real_)
}

run_smr_pair <- function(d) {
  i <- order(d$p_x, -abs(d$beta_x / d$se_x))[1]
  z_x <- qnorm(pmax(d$p_x[i], P_FLOOR) / 2, lower.tail = FALSE); z_y <- qnorm(pmax(d$p_y[i], P_FLOOR) / 2, lower.tail = FALSE)
  chi <- (z_x^2 * z_y^2) / (z_x^2 + z_y^2)
  data.table(smr_variant_id = d$variant_id[i], smr_b = d$beta_y[i] / d$beta_x[i], smr_p = pchisq(chi, df = 1, lower.tail = FALSE))
}

run_one_pair <- function(gene, peak, ca_data, eq_data) {
  pid <- paste(peak, gene, sep = "__")
  st  <- function(s) list(res = NULL, smr = NULL, dat = NULL, status = s)
  ca <- ca_data[.(peak, gene), nomatch = 0L]; eq <- eq_data[.(gene), nomatch = 0L]
  if (!nrow(ca) || !nrow(eq)) return(st("no_qtl"))
  ca[, snp_id := make_plink_id(variant_id, `A2(REF)`, `A1(ALT/effect allele)`)]
  ca <- ca[, .(variant_id, snp_id, ea_x = `A1(ALT/effect allele)`, oa_x = `A2(REF)`, beta_x = slope, se_x = slope_se, p_x = pval_nominal)]
  eq <- eq[, .(variant_id,          ea_y = `A1(ALT/effect allele)`, oa_y = `A2(REF)`, beta_y = slope, se_y = slope_se, p_y = pval_nominal)]
  dat <- merge(ca, eq, by = "variant_id")[ea_x == ea_y & oa_x == oa_y]
  dat <- dat[dedup_keep(beta_x, beta_y)]
  if (nrow(dat) < MIN_SNP) return(st("few_snp_qtl"))

  ld <- run_plink_ld(dat$snp_id, ref_bfile(peak), file.path(LD_DIR, gsub("[^A-Za-z0-9]", "_", pid)))
  if (is.null(ld)) return(st("plink_fail"))
  d <- dat[match(rownames(ld), snp_id)][!is.na(variant_id)]
  R <- ld[d$snp_id, d$snp_id, drop = FALSE]
  if (anyNA(R)) { bad <- which(apply(R, 1, function(x) any(!is.finite(x)))); d <- d[-bad]; R <- R[-bad, -bad, drop = FALSE] }
  if (nrow(d) < MIN_SNP) return(list(res = NULL, smr = NULL, dat = d, status = "few_snp_ld"))

  s_x <- if (FLIP_POINTS) ifelse(d$beta_x < 0, -1, 1) else rep(1, nrow(d))
  bx <- d$beta_x * s_x; by <- d$beta_y * s_x; R_use <- R * outer(s_x, s_x)
  eg <- tryCatch(eigen(outer(d$se_y, d$se_y) * R_use, symmetric = TRUE), error = function(e) NULL)
  res <- if (is.null(eg)) NULL else rbindlist(lapply(c(FALSE, TRUE), function(itc)
    tryCatch(run_gls_mr(bx, by, d$se_y, R_use, intercept = itc, eg = eg), error = function(e) NULL)))
  list(res = res, smr = run_smr_pair(d), dat = d, status = if (!is.null(res) && nrow(res) == 2L) "ok" else "gls_partial")
}

# ---------------------------- 读数据小工具 --------------------------------------
read_pb <- function(path) { dt <- fread(path); m <- as.matrix(dt[, -1, with = FALSE]); rownames(m) <- dt[[1]]; m }

read_filtered <- function(file, keys) {           # 只读第一列在 keys 里的行
  kf <- tempfile(tmpdir = LD_DIR, fileext = ".txt"); writeLines(unique(keys), kf); on.exit(unlink(kf), add = TRUE)
  fread(cmd = sprintf("zcat %s | awk 'NR==FNR{a[$1]=1; next} FNR==1 || ($1 in a)' %s -", shQuote(file), shQuote(kf)))
}
parse_pos <- function(v) as.numeric(sub("^[^_]*_([0-9]+).*$", "\\1", v))

# ---------------------------- 绘图小工具 ----------------------------------------
fmt_p <- function(p) if (length(p) == 0 || is.na(p)) "NA" else if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p)
getv  <- function(rf, col) if (is.null(rf) || !nrow(rf) || !col %in% names(rf)) NA else rf[[col]][1]
mspan <- function(txt, sig, col) if (!USE_MD) txt else if (isTRUE(sig)) sprintf("<span style='color:%s'>**%s**</span>", col, txt) else sprintf("<span style='color:%s'>%s</span>", COL_NS, txt)
join_lines <- function(x) paste(x, collapse = if (USE_MD) "<br>" else "\n")
# 只显示 p; 加粗 / 上色仍按 FDR < FDR_CUT 判断
pq_txt <- function(rf, method, key) {
  q <- getv(rf, paste0(key, "_fdr"))
  mspan(sprintf("%s: p=%s", method, fmt_p(getv(rf, paste0(key, "_p")))), isTRUE(q < FDR_CUT), pal_method[[method]])
}

add_corner <- function(p, labels, colors = rep("black", length(labels))) {
  for (k in seq_along(labels)) p <- p + annotate("text", x = -Inf, y = Inf, label = labels[k], hjust = -0.05, vjust = 1.5 + (k - 1) * 1.4, color = colors[k], size = ANNOT)
  p
}
empty_plot <- function(title_txt, txt = "Data not available") {
  ggplot() + theme_void() + labs(title = title_txt) + annotate("text", x = 0.5, y = 0.5, label = txt, color = "grey60", fontface = "italic", size = ANNOT + 0.4) +
    theme(plot.title = element_text(face = "bold", size = BASE + 1, hjust = 0))
}

# ---------------------------- Panel a: caQTL / eQTL 两轨 ------------------------
qtl_track <- function(d, title_txt, pt_col, xlim, pk_x, tss_x, corner_txt, chr, show_x = FALSE) {
  d_oth <- d[!d$is_mr & !d$is_lead, , drop = FALSE]; d_mr <- d[d$is_mr & !d$is_lead, , drop = FALSE]; d_lead <- d[d$is_lead, , drop = FALSE]
  p <- ggplot(d, aes(x, y)) +
    geom_vline(xintercept = pk_x,  color = COL_PEAK, linewidth = 0.6) +
    geom_vline(xintercept = tss_x, color = COL_TSS,  linewidth = 0.6, linetype = "dashed") +
    geom_point(data = d_oth, shape = 21, fill = pt_col, color = "white", size = PT_QTL, stroke = 0.15, alpha = 0.8)
  if (nrow(d_mr))   p <- p + geom_point(data = d_mr, shape = 24, fill = COL_MRSNP, color = "white", size = PT_QTL + 1.2, stroke = 0.3)
  if (nrow(d_lead)) p <- p + geom_point(data = d_lead, shape = 24, fill = LEAD_FILL, color = LEAD_BORDER, size = PT_LEAD - 1, stroke = 0.5)
  p <- p + scale_x_continuous(limits = xlim, expand = c(0, 0), oob = oob_keep, labels = mb_labels(diff(xlim))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.25))) +
    labs(title = title_txt, x = if (show_x) sprintf("%s (Mb)", chr) else NULL, y = expression(-log[10](italic(P)))) + theme_fig
  if (!show_x) p <- p + theme(axis.text.x = element_blank())
  add_corner(p, corner_txt$txt, corner_txt$col)
}

# peak 中点到 TSS 的距离标签 (写在页面 subtitle)
peak_tss_lab <- function(pk, tss) {
  dist <- abs(mean(pk) - tss) * 1e6
  if (dist < 1e3) sprintf("%.0f bp", dist) else sprintf("%.1f kb", dist / 1e3)
}

build_track_panel <- function(ca_d, eq_d, peak, gene, mr_ids, lead_id) {
  if (!nrow(ca_d) || !nrow(eq_d)) return(list(plot = empty_plot("caQTL / eQTL", "QTL data not available"), lab = "NA"))
  chr  <- sub(":.*$", "", peak)
  pk   <- as.numeric(strsplit(sub("^[^:]+:", "", peak), "-")[[1]]) / 1e6
  tss  <- median(eq_d$pos - eq_d$start_distance) / 1e6             # tensorQTL: start_distance = variant - phenotype start
  xlim <- (mean(pk) + tss) / 2 + c(-WIN_HALF, WIN_HALF)             # peak-TSS 中点 ± 500 kb

  mk <- function(dd) {
    o <- data.frame(x = dd$pos / 1e6, y = -log10(pmax(dd$pval_nominal, P_FLOOR)), variant_id = dd$variant_id, stringsAsFactors = FALSE)
    o$is_mr <- o$variant_id %in% mr_ids; o$is_lead <- !is.na(lead_id) & o$variant_id == lead_id
    o[o$x >= xlim[1] & o$x <= xlim[2], , drop = FALSE]
  }
  d_ca <- mk(ca_d); d_eq <- mk(eq_d)
  p_ca <- qtl_track(d_ca, "caQTL", COL_PEAK, xlim, mean(pk), tss, list(txt = sprintf("%d variants", nrow(d_ca)), col = "black"), chr)
  p_eq <- qtl_track(d_eq, "eQTL",  COL_TSS,  xlim, mean(pk), tss,
                    list(txt = sprintf("%d variants", nrow(d_eq)), col = "black"), chr, show_x = TRUE)
  list(plot = p_ca / p_eq + plot_layout(heights = c(1, 1)), lab = peak_tss_lab(pk, tss))
}

# ---------------------------- Panel b: MR ---------------------------------------
build_mr_panel <- function(obj, rf) {
  ttl <- "Mendelian randomization"
  d <- obj$dat
  if (is.null(d) || !nrow(d)) return(empty_plot(ttl, sprintf("No SNP data (%s)", obj$status)))
  sgn <- ifelse(d$beta_x < 0, -1, 1); xf <- d$beta_x * sgn; yf <- d$beta_y * sgn
  lead_id <- getv(rf, "smr_variant_id"); is_lead <- !is.na(lead_id) & d$variant_id == lead_id
  pts <- data.frame(x = xf, y = yf, xmin = xf - 1.96 * d$se_x, xmax = xf + 1.96 * d$se_x, ymin = yf - 1.96 * d$se_y, ymax = yf + 1.96 * d$se_y,
                    is_lead = is_lead, label = ifelse(is_lead, d$variant_id, NA_character_))
  lines <- data.table(series = MR_LINES,
                      slope = as.numeric(c(getv(rf, "smr_b"), getv(rf, "ivw_beta"), getv(rf, "egger_beta"))),
                      intercept = as.numeric(c(0, 0, getv(rf, "egger_int"))))[is.finite(slope) & is.finite(intercept)]
  lines[, series := factor(series, levels = MR_LINES)]
  sub_txt <- join_lines(c(paste0(pq_txt(rf, "Causal SMR", "smr"),   "  |  ", pq_txt(rf, "Causal IVW", "ivw")),
                          paste0(pq_txt(rf, "Causal Egger", "egger"), "  |  ", pq_txt(rf, "Causal Joint", "cct"))))

  p <- ggplot(pts, aes(x, y)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = COL_REF, linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = COL_REF, linewidth = 0.4) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), color = COL_ERR, width = 0, linewidth = 0.45) +
    geom_errorbar(aes(xmin = xmin, xmax = xmax), orientation = "y", color = COL_ERR, width = 0, linewidth = 0.45) +
    geom_abline(data = lines, aes(slope = slope, intercept = intercept, color = series, linetype = series), linewidth = LWD) +
    geom_point(data = pts[!pts$is_lead, ], shape = 24, fill = COL_MRSNP, color = "white", size = PT_MR, stroke = 0.4) +
    geom_point(data = pts[pts$is_lead, ], shape = 24, fill = LEAD_FILL, color = LEAD_BORDER, size = PT_LEAD, stroke = 0.5) +
    geom_text_repel(data = pts[pts$is_lead, ], aes(label = label), size = ANNOT, fontface = "bold", color = LEAD_BORDER,
                    box.padding = 0.6, min.segment.length = 0, segment.size = 0.3, seed = 1) +
    scale_color_manual(values = pal_method, breaks = levels(lines$series)) +
    scale_linetype_manual(values = c(`Causal SMR` = "dashed", `Causal IVW` = "solid", `Causal Egger` = "solid"), breaks = levels(lines$series)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) + expand_limits(x = 0) +
    labs(title = ttl, subtitle = sub_txt, x = "caQTL effect size (aligned)", y = "eQTL effect size (aligned)") +
    theme_fig + legend_inside(c(0.02, 0.98))
  p
}

# ---------------------------- Panel c / d: 边际相关 -----------------------------
marg_subtitle <- function(rf, method, key) {
  q <- getv(rf, paste0(key, "_fdr")); rho <- getv(rf, paste0(key, "_rho"))
  txt <- join_lines(c(sprintf("Spearman rho = %s", if (is.na(rho)) "NA" else sprintf("%.3f", rho)),
                      sprintf("p=%s", fmt_p(getv(rf, paste0(key, "_p"))))))
  mspan(txt, isTRUE(q < FDR_CUT), pal_method[[method]])
}

plot_wnn <- function(gene, peak_seurat, rf) {
  ttl <- "Marginal Multiomics (WNN metacells)"
  if (!(gene %chin% colnames(wnn_rna_log)) || !(peak_seurat %chin% colnames(wnn_peak_log))) return(empty_plot(ttl))
  df <- data.frame(X = as.numeric(wnn_peak_log[, peak_seurat]), Y = as.numeric(wnn_rna_log[, gene]))
  p <- ggplot(df, aes(X, Y)) + geom_smooth(method = "lm", formula = y ~ x, color = "grey30", fill = "grey85", linetype = "longdash", linewidth = LWD, alpha = 0.5)
  if (HAS_CT) {
    df$CT <- as.character(wnn.celltype); df$Col <- as.character(wnn.celltype.col)
    lab <- as.data.table(df)[, .(X = median(X), Y = median(Y)), by = .(CT, Col)]
    p <- p + geom_point(data = df, aes(fill = Col), shape = 21, color = "white", size = PT_MARG + 0.4, stroke = 0.3, alpha = 0.95) +
      geom_text_repel(data = lab, aes(X, Y, label = CT, color = Col), inherit.aes = FALSE, size = ANNOT - 0.3, fontface = "bold",
                      bg.color = "white", bg.r = 0.12, box.padding = 0.4, min.segment.length = 0.2, segment.size = 0.25, max.overlaps = 25, seed = 1) +
      scale_fill_identity() + scale_color_identity()
  } else p <- p + geom_point(shape = 21, fill = pal_method[["Marginal Multiomics"]], color = "white", size = PT_MARG + 0.4, stroke = 0.3)
  p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = ttl, subtitle = marg_subtitle(rf, "Marginal Multiomics", "wnn"), x = "ATAC log1p(CPM)", y = "RNA log1p(CPM)") + theme_fig
  add_corner(p, sprintf("%d metacells", nrow(df)))
}

plot_pseudo <- function(gene, peak_seurat, rf) {
  ttl <- "Marginal Pseudobulk (individuals)"
  if (!(gene %chin% rownames(subject_gene_logcpm)) || !(peak_seurat %chin% rownames(subject_peak_logcpm))) return(empty_plot(ttl))
  df <- data.frame(X = as.numeric(subject_peak_logcpm[peak_seurat, ]), Y = as.numeric(subject_gene_logcpm[gene, ]))
  p <- ggplot(df, aes(X, Y)) +
    geom_point(shape = 21, fill = pal_method[["Marginal Pseudobulk"]], color = "white", size = PT_MARG, stroke = 0.25, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#00624F", fill = pal_method[["Marginal Pseudobulk"]], linewidth = LWD, alpha = 0.18) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title = ttl, subtitle = marg_subtitle(rf, "Marginal Pseudobulk", "pseudo"), x = "ATAC log1p(CPM)", y = "RNA log1p(CPM)") + theme_fig
  add_corner(p, sprintf("%d individuals", nrow(df)))
}

# ---------------------------- Panel e: 位点 link 图 (来自 Plot8, 只突出当前 pair) ----
# 基因注释 (EnsDb.Hsapiens.v86, hg38). 不用 seqlevelsStyle<- "UCSC" (会联网, 计算节点报错), 直接加 "chr" 前缀
to_ucsc <- function(x) sub("^chrMT$", "chrM", paste0("chr", as.character(x)))
gene_gr <- genes(EnsDb.Hsapiens.v86, filter = SeqNameFilter(c(1:22, "X", "Y")))
gene_dt <- data.table(chr = to_ucsc(seqnames(gene_gr)), start = start(gene_gr), end = end(gene_gr),
                      strand = as.character(strand(gene_gr)), gene_name = gene_gr$gene_name)
gene_dt[, tss := fifelse(strand == "+", start, end)]
exon_gr <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
exon_dt <- data.table(chr = to_ucsc(seqnames(exon_gr)), start = start(exon_gr), end = end(exon_gr),
                      gene_name = exon_gr$gene_name, type = as.character(exon_gr$type))[type %chin% c("exon", "cds", "utr")]
rm(gene_gr, exon_gr); invisible(gc())
msg(sprintf("EnsDb v86: %d genes | %d exon/cds/utr 片段", nrow(gene_dt), nrow(exon_dt)))

parse_peak <- function(dt) {
  dt[, c("pk_chr", "pk_start", "pk_end") := tstrsplit(peak_id, "[:-]")]
  dt[, `:=`(pk_start = as.integer(pk_start), pk_end = as.integer(pk_end))]
  dt[, pk_mid := (pk_start + pk_end) / 2]
  dt[]
}

# 每个方法 FDR < cut 的 links (长表: 一行 = 一个 gene-peak-method)
get_sig_links <- function(res, genes, fdr_cut = LINK_FDR_CUT) {
  sub <- res[gene_id %chin% genes]
  lk <- rbindlist(lapply(seq_len(nrow(LINK_METHODS)), function(k) {
    m <- LINK_METHODS[k]; f <- sub[[m$fdr_col]]
    eff <- if (is.na(m$eff_col)) rep(NA_real_, nrow(sub)) else sub[[m$eff_col]]
    idx <- which(!is.na(f) & f < fdr_cut)
    data.table(gene_id = sub$gene_id[idx], peak_id = sub$peak_id[idx], method = m$method, fdr = f[idx], effect = eff[idx])
  }))
  if (!nrow(lk)) lk <- data.table(gene_id = character(), peak_id = character(), method = character(), fdr = numeric(), effect = numeric())
  parse_peak(lk)
  lk[, direction := fcase(is.na(effect), "No direction", effect > 0, "Positive", default = "Negative")]
  lk[, score := -log10(pmax(fdr, 1e-300))]
  lk[]
}

plot_locus_pair <- function(g, peak, lk_all, mr_dat, lead_id, score_max) {
  ttl <- "Peak-gene links"
  lk <- copy(lk_all[gene_id == g])
  ex <- parse_peak(data.table(gene_id = g, peak_id = peak))
  cand <- gene_dt[gene_name == g & chr == ex$pk_chr][order(abs(tss - ex$pk_mid))]
  if (!nrow(cand)) return(empty_plot(ttl, sprintf("%s not found in EnsDb v86", g)))
  tg <- cand[1]
  lk <- lk[pk_chr == tg$chr]
  sd <- if (is.null(mr_dat) || !nrow(mr_dat)) data.table(variant_id = character(), pos = numeric()) else
          data.table(variant_id = mr_dat$variant_id, pos = parse_pos(mr_dat$variant_id))
  sd[, is_lead := !is.na(lead_id) & variant_id == lead_id]

  # 窗口: TSS + 显著 peak + 当前 peak (不含 MR SNP, 否则 SNP 分布太宽会把 link 压扁; 窗口外的 SNP 不画, 只计数)
  xs <- c(tg$tss, lk$pk_start, lk$pk_end, ex$pk_start, ex$pk_end)
  if (length(xs) > 1) { pad <- max(LK_PAD_MIN, diff(range(xs)) * LK_PAD_FRAC); xlim <- c(min(xs) - pad, max(xs) + pad) } else xlim <- tg$tss + c(-LK_EMPTY_HALF, LK_EMPTY_HALF)
  wid <- diff(xlim)
  n_snp <- nrow(sd); sd <- sd[pos >= xlim[1] & pos <= xlim[2]]

  ytitle   <- element_text(angle = 90, vjust = 0.5, hjust = 0.5, face = "bold", size = BASE - 2)
  th_track <- theme_classic(base_size = BASE) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank(),
          axis.title.y = ytitle, legend.position = "none", plot.margin = margin(1, 4, 1, 4))
  tss_line <- geom_vline(xintercept = tg$tss, linetype = "dashed", color = "grey60", linewidth = 0.4)

  ## Gene 轨 (只画目标基因)
  gx  <- exon_dt[chr == tg$chr & gene_name == g]
  arr <- tg$tss + ifelse(tg$strand == "+", 1, -1) * wid * 0.04
  pG <- ggplot() + tss_line +
    annotate("segment", x = tg$start, xend = tg$end, y = 0, yend = 0, color = GENE_COL, linewidth = 0.7) +
    geom_rect(data = gx[type != "cds"], aes(xmin = start, xmax = end, ymin = -0.15, ymax = 0.15), fill = GENE_COL, color = NA) +
    geom_rect(data = gx[type == "cds"], aes(xmin = start, xmax = end, ymin = -0.28, ymax = 0.28), fill = GENE_COL, color = NA) +
    annotate("segment", x = tg$tss, xend = tg$tss, y = 0, yend = 0.45, color = GENE_COL, linewidth = 0.7) +
    annotate("segment", x = tg$tss, xend = arr, y = 0.45, yend = 0.45, color = GENE_COL, linewidth = 0.7, arrow = arrow(length = unit(0.18, "cm"))) +
    annotate("text", x = (max(tg$start, xlim[1]) + min(tg$end, xlim[2])) / 2, y = 0.8, label = g, fontface = "bold.italic", size = ANNOT + 0.4) +
    coord_cartesian(xlim = xlim, ylim = c(-0.4, 1.05), expand = FALSE) +
    labs(title = ttl, y = "Gene") + th_track +
    theme(plot.title = element_text(face = "bold", size = BASE + 1, hjust = 0, margin = margin(b = 3)))

  ## Peak 轨: 该基因所有显著 peak (黑色) + 当前 pair 的 peak (红色, 画在最上层)
  pk <- unique(lk[, .(peak_id, pk_start, pk_end, pk_mid)], by = "peak_id")[peak_id != peak]
  hw_min <- wid * 0.003
  if (nrow(pk)) pk[, hw := pmax((pk_end - pk_start) / 2, hw_min)]
  ex[, hw := pmax((pk_end - pk_start) / 2, hw_min)]
  pP <- ggplot() + tss_line
  if (nrow(pk)) pP <- pP + geom_rect(data = pk, aes(xmin = pk_mid - hw, xmax = pk_mid + hw, ymin = 0.15, ymax = 0.85), fill = PEAK_OTH, color = NA)
  pP <- pP + geom_rect(data = ex, aes(xmin = pk_mid - hw, xmax = pk_mid + hw, ymin = 0.05, ymax = 0.95), fill = PEAK_CUR, color = NA) +
    coord_cartesian(xlim = xlim, ylim = c(0, 1), expand = FALSE) + labs(y = "Peaks") + th_track

  ## MR SNP 轨 (当前 pair): 竖线 = MR instruments, 金色三角 = SMR lead SNP
  pS <- ggplot() + tss_line +
    annotate("rect", xmin = ex$pk_mid - ex$hw, xmax = ex$pk_mid + ex$hw, ymin = -0.45, ymax = 0.45, fill = PEAK_CUR, alpha = 0.25)
  if (nrow(sd)) {
    pS <- pS + geom_segment(data = sd, aes(x = pos, xend = pos, y = -0.32, yend = 0.32), color = SNP_COL, linewidth = 0.45)
    if (sd[is_lead == TRUE, .N]) pS <- pS +
      geom_point(data = sd[is_lead == TRUE], aes(pos, 0.42), shape = 25, fill = LEAD_FILL, color = LEAD_BORDER, size = 3, stroke = 0.5) +
      geom_text_repel(data = sd[is_lead == TRUE], aes(pos, 0.42, label = variant_id), size = ANNOT - 0.6, fontface = "bold", color = LEAD_BORDER,
                      nudge_y = 0.35, direction = "x", min.segment.length = 0, segment.size = 0.3, seed = 1)
  } else pS <- pS + annotate("text", x = mean(xlim), y = 0, label = "no MR SNPs", color = "grey40", fontface = "italic", size = ANNOT - 0.4)
  pS <- pS + annotate("text", x = xlim[1] + wid * 0.01, y = -0.62, label = if (nrow(sd) < n_snp) sprintf("n = %d (%d in view)", n_snp, nrow(sd)) else sprintf("n = %d", n_snp), hjust = 0, size = ANNOT - 0.6, color = "grey20") +
    coord_cartesian(xlim = xlim, ylim = c(-0.8, 1.05), expand = FALSE, clip = "off") + labs(y = "MR SNPs") + th_track

  ## Link 轨: 6 个方法各一行; 当前 pair 的 link 画在最上层
  lev <- LINK_METHODS$method
  cnt <- lk[, .N, by = method][data.table(method = lev), on = "method"][is.na(N), N := 0L]
  cnt[, method := factor(method, levels = lev)]
  s_lo <- -log10(LINK_FDR_CUT)
  if (nrow(lk)) {
    lk[, `:=`(x1 = pk_mid, x2 = tg$tss, link_id = .I, is_cur = peak_id == peak, score_plot = pmin(pmax(score, s_lo), score_max))]
    lk[, h := pmax(0.2, abs(x2 - x1) / max(1, max(abs(x2 - x1))))]
    tt <- seq(0, 1, length.out = 100)
    arcs <- lk[, .(x = x1 + (x2 - x1) * tt, y = -0.92 * h * sin(pi * tt)), by = .(link_id, method, direction, score_plot, is_cur)]
    setorder(arcs, is_cur, score_plot)
  } else arcs <- data.table(link_id = integer(), method = character(), direction = character(), score_plot = numeric(), is_cur = logical(), x = numeric(), y = numeric())
  arcs[, `:=`(method = factor(method, levels = lev), direction = factor(direction, levels = names(DIR_COL)))]

  pL <- ggplot() + tss_line +
    geom_blank(data = data.table(method = factor(lev, levels = lev), x = tg$tss, y = 0), aes(x, y)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_path(data = arcs, aes(x, y, group = link_id, color = direction, linewidth = score_plot), alpha = 0.95, lineend = "round") +
    geom_text(data = cnt, aes(x = xlim[1] + wid * 0.01, y = -0.9, label = sprintf("%s  (n = %d)", method, N), color = I(pal_method[as.character(method)])),
              hjust = 0, size = ANNOT - 0.4, fontface = "bold") +
    facet_grid(method ~ ., drop = FALSE) +
    scale_color_manual(values = DIR_COL, drop = FALSE, guide = "none") +
    scale_linewidth_continuous(range = LW_RANGE, limits = c(s_lo, score_max), guide = "none") +
    scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6), breaks = scales::breaks_pretty(n = 4)) +
    coord_cartesian(xlim = xlim, ylim = c(-1, 0.04), expand = FALSE) +
    labs(x = sprintf("%s position (Mb)", tg$chr), y = "Links") + theme_classic(base_size = BASE) +
    theme(strip.background = element_blank(), strip.text = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.title.y = element_text(angle = 90, face = "bold", size = BASE - 2),
          axis.text.x = element_text(size = BASE - 1, color = "black"), axis.title.x = element_text(size = BASE),
          panel.spacing.y = unit(0.6, "mm"), legend.position = "none", plot.margin = margin(1, 4, 2, 4))

  (pG / pP / pS / pL) + plot_layout(heights = c(0.8, 0.35, 0.55, 6.6))
}

# ---------------------------- WNN (全局一份) ------------------------------------
e1 <- new.env(); load(WNN_FILE, envir = e1)
wnn_rna_log <- log1p(e1$wnn.rna); wnn_peak_log <- log1p(e1$wnn.peak); rm(e1)
e2 <- new.env(); load(WNN_CT_FILE, envir = e2)
if (exists("wnn.celltype", envir = e2, inherits = FALSE))     wnn.celltype     <- get("wnn.celltype", envir = e2)
if (exists("wnn.celltype.col", envir = e2, inherits = FALSE)) wnn.celltype.col <- get("wnn.celltype.col", envir = e2)
rm(e2); invisible(gc())
HAS_CT <- exists("wnn.celltype") && exists("wnn.celltype.col") && length(wnn.celltype) == nrow(wnn_rna_log)
msg("WNN: ", nrow(wnn_rna_log), " metacells | cell-type colours: ", HAS_CT)

# ---------------------------- 主循环 ---------------------------------------------
cts <- if (CT_ARG == "all") CELLTYPES else CT_ARG
stopifnot(all(cts %in% CELLTYPES))

# MR 结果 (三个 CT; 统一 sigma2 规则). panel e 的线宽标尺在三个 CT 的 example 基因间共用
RES <- setNames(lapply(CELLTYPES, function(ct) fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")))), CELLTYPES)
SCORE_MAX <- min(SCORE_CAP, max(c(2, unlist(lapply(CELLTYPES, function(ct) get_sig_links(RES[[ct]], unique(EXAMPLES[[ct]]$gene_id))$score)))))
msg(sprintf("panel e 线宽标尺: -log10(FDR) in [%.2f, %.2f]", -log10(LINK_FDR_CUT), SCORE_MAX))
LAYOUT <- c(patchwork::area(1, 1), patchwork::area(1, 2), patchwork::area(2, 1), patchwork::area(2, 2), patchwork::area(1, 3, 2, 3))
CAPTION <- sprintf(paste0("e: significant peak-gene links (FDR < %.2f) per method; arc colour = effect direction (red positive, blue negative, black = no direction for Causal Joint); ",
                          "arc width ∝ -log10 FDR (shared scale, %.1f-%.0f); red peak = this pair's peak, black = other significant peaks; ticks = MR instruments within the window, gold triangle = SMR lead SNP."),
                   LINK_FDR_CUT, -log10(LINK_FDR_CUT), SCORE_MAX)
open_pdf <- function(f) if (capabilities("cairo")) cairo_pdf(f, width = PAGE_W, height = PAGE_H, onefile = TRUE) else pdf(f, width = PAGE_W, height = PAGE_H, useDingbats = FALSE)

for (CT in cts) {
  msg("==================== ", CT, " ====================")
  mr_all <- RES[[CT]]
  ex <- EXAMPLES[[CT]]
  pc <- mr_all[ex, on = .(gene_id, peak_id), nomatch = NULL]
  miss <- ex[!pc, on = .(gene_id, peak_id)]
  if (nrow(miss)) { msg("警告: 以下 example 不在 MR 结果中 (或被 sigma2 过滤), 跳过:"); print(miss) }
  setorder(pc, ord)
  if (N_EXAMPLES > 0) pc <- head(pc, N_EXAMPLES)
  pc[, peak_seurat := gsub("[:-]", "-", peak_id)]
  msg("作图 pairs: ", nrow(pc))

  # QTL 数据 (只读需要的 gene / peak)
  eq_full <- read_filtered(eqtl_file(CT), unique(pc$gene_id));        eq_full[, pos := parse_pos(variant_id)]
  ca_full <- read_filtered(caqtl_full_file(CT), unique(pc$peak_id));  ca_full[, pos := parse_pos(variant_id)]
  eq_sub  <- eq_full[, .(phenotype_id, variant_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se, pval_nominal)]
  setkey(eq_sub, phenotype_id, variant_id)
  ca_clump <- fread(caqtl_clump_file(CT))
  ca_sub   <- ca_clump[unique(pc[, .(phenotype_id = peak_id, gene_id)]), on = .(phenotype_id, gene_id), nomatch = 0L,
                       .(phenotype_id, variant_id, gene_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se, pval_nominal)]
  rm(ca_clump); invisible(gc()); setkey(ca_sub, phenotype_id, gene_id)

  # 重跑 MR 取 SNP 层数据, 并与全量结果对账
  objs <- setNames(lapply(seq_len(nrow(pc)), function(i) run_one_pair(pc$gene_id[i], pc$peak_id[i], ca_sub, eq_sub)), paste(pc$peak_id, pc$gene_id, sep = "__"))
  chk <- rbindlist(lapply(seq_len(nrow(pc)), function(i) {
    o <- objs[[i]]; b <- if (is.null(o$res) || !nrow(o$res)) NA_real_ else o$res[model == "IVW", beta][1]
    data.table(pair = names(objs)[i], status = o$status, ivw_beta_rerun = b, ivw_beta_formal = pc$ivw_beta[i],
               smr_snp_match = !is.null(o$smr) && identical(o$smr$smr_variant_id, pc$smr_variant_id[i]))
  }))
  print(chk)
  # 新 MR 表没有 Egger 截距: 用重跑的 Egger 截距 (alpha) 画 Egger 线; 重跑 Egger 斜率应与表中 egger_beta 一致
  egg_get <- function(o, col) if (is.null(o$res) || !nrow(o$res[model == "Egger"])) NA_real_ else o$res[model == "Egger"][[col]][1]
  pc[, `:=`(egger_int = sapply(objs, egg_get, col = "alpha"), egger_beta_rerun = sapply(objs, egg_get, col = "beta"))]
  print(pc[, .(gene_id, peak_id, egger_beta, egger_beta_rerun, egger_int, abs_diff = abs(egger_beta - egger_beta_rerun))])

  # pseudobulk
  pk_cnt <- read_pb(file.path(PB_DIR, paste0("pseudobulk_atac_", CT, ".tsv.gz")))
  gn_cnt <- read_pb(file.path(PB_DIR, paste0("pseudobulk_rna_",  CT, ".tsv.gz")))
  smp <- intersect(colnames(pk_cnt), colnames(gn_cnt))
  lib <- fread(file.path(PB_DIR, paste0("libsize_", CT, ".csv")))
  subject_gene_logcpm <- log1p(sweep(gn_cnt[, smp, drop = FALSE], 2, setNames(lib$rna_lib_size,  lib$sample)[smp], "/") * 1e6)
  subject_peak_logcpm <- log1p(sweep(pk_cnt[, smp, drop = FALSE], 2, setNames(lib$atac_lib_size, lib$sample)[smp], "/") * 1e6)
  rm(pk_cnt, gn_cnt); invisible(gc())

  # panel e: 该 CT 中 example 基因的全部显著 links
  lk_all <- get_sig_links(mr_all, unique(pc$gene_id))

  # 作图
  suffix   <- if (N_EXAMPLES > 0) sprintf("_top%d", N_EXAMPLES) else ""
  pdf_file <- file.path(FIG_DIR, sprintf("Fig2_examples_%s%s_v2.pdf", CT, suffix))
  open_pdf(pdf_file)
  for (i in seq_len(nrow(pc))) {
    rf <- pc[i]; gene <- rf$gene_id; peak <- rf$peak_id; pid <- paste(peak, gene, sep = "__")
    msg(sprintf("  [%d/%d] %s", i, nrow(pc), pid))
    tryCatch({
      obj <- objs[[pid]]; mr_ids <- if (is.null(obj$dat)) character(0) else obj$dat$variant_id
      tr  <- build_track_panel(ca_full[phenotype_id == peak], eq_full[phenotype_id == gene], peak, gene, mr_ids, getv(rf, "smr_variant_id"))
      pA  <- wrap_elements(full = tr$plot)
      pB  <- build_mr_panel(obj, rf)
      pC  <- plot_wnn(gene, rf$peak_seurat, rf)
      pD  <- plot_pseudo(gene, rf$peak_seurat, rf)
      pE  <- wrap_elements(full = plot_locus_pair(gene, peak, lk_all, obj$dat, getv(rf, "smr_variant_id"), SCORE_MAX))
      fig <- wrap_plots(pA, pB, pC, pD, pE, design = LAYOUT) +
        plot_annotation(tag_levels = "a", caption = CAPTION,
          title    = sprintf("%s  |  %s  (%s)", gene, peak, CT),
          subtitle = sprintf("Example category: %s   |   peak-TSS distance: %s   |   %d SNPs in MR", rf$category, tr$lab, if (is.null(obj$dat)) 0L else nrow(obj$dat)),
          theme = theme(plot.title = element_text(size = BASE + 4, face = "bold", hjust = 0),
                        plot.subtitle = element_text(size = BASE, color = "grey30", hjust = 0, margin = margin(b = 6)),
                        plot.caption  = element_text(size = BASE - 2, color = "grey35", hjust = 0))) &
        theme(plot.tag = element_text(face = "bold", size = BASE + 3))
      print(fig)
    }, error = function(e) msg("    作图失败: ", conditionMessage(e)))
  }
  dev.off()
  msg("Saved: ", pdf_file)
  rm(eq_full, ca_full, eq_sub, ca_sub, objs, subject_gene_logcpm, subject_peak_logcpm, lk_all); invisible(gc())
}

unlink(LD_DIR, recursive = TRUE)
msg("All done.")



