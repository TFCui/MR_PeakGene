# QC2_1_FlipScan: 全基因组 MR-Egger flip 前后对比 (流程见 Scripts/Pipeline.md 第 4.2 节)
# 用法: Rscript QC2_1_FlipScan.R <CT>
# pair: Results/MR/<CT>_AllPairs_results.tsv.gz 中的全部 pair (已过 NA / ivw_sigma2 <= 10)
# 每个 pair: 与 Scripts/Model/MR_AllPairs.R 相同的对齐 -> 去重 -> PLINK LD (只做一次), 然后 GLS Egger 拟合两次:
#   flip = TRUE (与正式流程相同: bx, by, R 同乘 sign(beta_x)) / flip = FALSE (原始等位编码)
#   函数拷自 MR_AllPairs.R; 分 144 个 chunk, mclapply 并行, chunk 已存在则跳过 (断点续跑)
# 输出: Results/QC/2_FlipPoints/<CT>_Egger_FlipScan.tsv.gz
#   两版 Egger 的 beta / se / p / 截距 alpha / p_alpha; 不 flip 版按 gene 内 BH 得 egger_fdr_noflip; flip 版 FDR 用正式表 egger_fdr
# 环境变量 TEST_N > 0: 只跑前 N 个 pair (冒烟测试), 文件名带 _TEST<N>

library(data.table)
library(parallel)
options(width = 250)

args      <- commandArgs(trailingOnly = TRUE)
CELL_TYPE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "Bn_TCL1A"

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
MR_FILE  <- file.path(PROJ, "Results/MR", paste0(CELL_TYPE, "_AllPairs_results.tsv.gz"))
OUT_DIR  <- file.path(PROJ, "Results/QC/2_FlipPoints")
PLINK    <- file.path(BASE_DIR, "R_package/my_plink/plink")
REF_DIR  <- file.path(PROJ, "data/ref_chr/hg38")                   # 软链接, 只读
EQTL_FILE  <- file.path(BASE_DIR, "Data/OMIX009786/cis_eQTL_Summary_Statics/cis_eQTL", paste0(CELL_TYPE, ".txt.gz"))
CAQTL_FILE <- file.path(BASE_DIR, "Data/OMIX009786/cis_caQTL_Summary_Statics/cis_caQTL/LD_clump", paste0(CELL_TYPE, "_Clumped_r2_0.99_kb_20000.txt.gz"))

# 固定参数 (与 MR_AllPairs.R 一致)
MIN_SNP <- 3; DEDUP_TOL <- 1e-3; KEEP_ALLELE <- TRUE
EIG_TRUNCATE <- TRUE; EIG_VAR_KEEP <- 0.999; EIG_MIN_KEEP <- 5; EIG_TOL <- 1e-10; SIGMA_MODE <- "over"

NCORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "36"))
N_CHUNK <- 144
TEST_N  <- as.integer(Sys.getenv("TEST_N", "0"))
TAG     <- if (TEST_N > 0) paste0(CELL_TYPE, "_TEST", TEST_N) else CELL_TYPE

CHUNK_DIR <- file.path(PROJ, "tmp/QC2_FlipScan", TAG)
PLINK_TMP <- file.path(PROJ, "tmp/QC2_FlipScan/plink")
for (d in c(CHUNK_DIR, PLINK_TMP, OUT_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 函数 (拷自 Scripts/Model/MR_AllPairs.R)

dedup_keep <- function(x, y, tol = DEDUP_TOL) {
  keep <- logical(length(x)); ref_x <- ref_y <- numeric(0)
  for (j in seq_along(x)) {
    if (length(ref_x) && any(abs(ref_x - x[j]) <= tol & abs(ref_y - y[j]) <= tol)) next
    keep[j] <- TRUE; ref_x <- c(ref_x, x[j]); ref_y <- c(ref_y, y[j])
  }
  keep
}

make_plink_id <- function(variant_id, ref, alt) {
  chr <- sub("^chr", "", sub("_.*$", "", variant_id))
  bp  <- sub("^[^_]*_", "", variant_id)
  chr <- fcase(chr == "X", "23", chr == "Y", "24", chr %chin% c("M", "MT"), "26", default = chr)
  paste(chr, bp, ref, alt, sep = ":")
}

ref_bfile <- function(peak) file.path(REF_DIR, sub(":.*$", "", peak))

run_plink_ld <- function(ids, bfile, pre, keep_allele_order = KEEP_ALLELE) {
  fs <- paste0(pre, c(".snps.txt", ".ld", ".snplist", ".log", ".nosex"))
  on.exit(unlink(fs), add = TRUE)
  fwrite(data.table(unique(ids[!is.na(ids)])), fs[1], col.names = FALSE)
  kao <- if (keep_allele_order) "--keep-allele-order" else ""
  cmd <- paste(PLINK, "--bfile", bfile, "--extract", fs[1], kao, "--r square --write-snplist --threads 1 --memory 2000 --out", pre, "--silent")
  st  <- suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))
  if (st != 0 || !file.exists(fs[2]) || !file.exists(fs[3])) return(NULL)
  snps <- fread(fs[3], header = FALSE)[[1]]
  R    <- as.matrix(fread(fs[2], header = FALSE, na.strings = c("nan", "-nan", "NaN", "NA")))
  storage.mode(R) <- "numeric"
  dimnames(R) <- list(snps, snps)
  R
}

# 对称白化 Sigma^(-1/2) = V diag(lam^-1/2) V', 截断方向置 0; df = n - p
run_gls_mr <- function(bx, by, se_y, R, intercept = FALSE, eg = NULL) {
  n <- length(bx)
  if (is.null(eg)) eg <- eigen(outer(se_y, se_y) * R, symmetric = TRUE)
  lam <- eg$values; V <- eg$vectors
  pos  <- lam > EIG_TOL * lam[1]
  cum  <- cumsum(pmax(lam, 0)) / sum(pmax(lam, 0))
  keep <- if (EIG_TRUNCATE) (c(TRUE, cum[-n] < EIG_VAR_KEEP) | seq_len(n) <= EIG_MIN_KEEP) & pos else pos
  k <- sum(keep); p <- if (intercept) 2L else 1L; df <- n - p
  if (k < p)  stop(sprintf("k = %d, p = %d", k, p))
  if (df < 1) stop(sprintf("df = %d", df))
  d_half <- numeric(n); d_half[keep] <- lam[keep]^(-0.5)
  W  <- V %*% diag(d_half, n) %*% t(V)
  X  <- if (intercept) cbind(alpha = 1, beta = bx) else cbind(beta = bx)
  ys <- as.numeric(W %*% by); Xs <- W %*% X
  fit <- lm(ys ~ Xs - 1); est <- unname(coef(fit))
  se_fixed <- sqrt(diag(solve(crossprod(Xs))))
  sigma <- sqrt(sum(residuals(fit)^2) / df)
  infl  <- switch(SIGMA_MODE, over = max(1, sigma), under = 1 / min(1, sigma), none = 1)
  se <- se_fixed * infl; pv <- 2 * pt(abs(est / se), df = df, lower.tail = FALSE)
  c(beta = est[p], se = se[p], p = pv[p], alpha = if (intercept) est[1] else NA_real_, p_alpha = if (intercept) pv[1] else NA_real_)
}

# 一个 pair: 对齐 -> 去重 -> LD 一次; Egger flip / 不 flip 各拟合一次
run_pair_flip <- function(gene, peak, ca_data, eq_data, tmp_dir) {
  na_row <- function(status, n = NA_integer_) data.table(gene_id = gene, peak_id = peak, n_snp = n, status = status)
  ca <- ca_data[.(peak, gene), nomatch = 0L]; eq <- eq_data[.(gene), nomatch = 0L]
  if (!nrow(ca) || !nrow(eq)) return(na_row("no_qtl"))
  ca[, snp_id := make_plink_id(variant_id, `A2(REF)`, `A1(ALT/effect allele)`)]
  ca <- ca[, .(variant_id, snp_id, ea_x = `A1(ALT/effect allele)`, oa_x = `A2(REF)`, beta_x = slope)]
  eq <- eq[, .(variant_id,          ea_y = `A1(ALT/effect allele)`, oa_y = `A2(REF)`, beta_y = slope, se_y = slope_se)]
  dat <- merge(ca, eq, by = "variant_id")[ea_x == ea_y & oa_x == oa_y]
  dat <- dat[dedup_keep(beta_x, beta_y, tol = DEDUP_TOL)]
  if (nrow(dat) < MIN_SNP) return(na_row("few_snp_qtl", nrow(dat)))

  ld <- run_plink_ld(dat$snp_id, ref_bfile(peak), file.path(tmp_dir, gsub("[^A-Za-z0-9]", "_", paste(peak, gene, sep = "__"))))
  if (is.null(ld)) return(na_row("plink_fail"))
  d <- dat[match(rownames(ld), snp_id)][!is.na(variant_id)]
  R <- ld[d$snp_id, d$snp_id, drop = FALSE]
  if (nrow(d) < MIN_SNP) return(na_row("few_snp_ld", nrow(d)))
  if (anyNA(R)) {
    bad <- which(apply(R, 1, function(x) any(!is.finite(x))))
    d <- d[-bad]; R <- R[-bad, -bad, drop = FALSE]
    if (nrow(d) < MIN_SNP) return(na_row("few_snp_nan", nrow(d)))
  }

  fit <- function(flip) {
    s_x <- if (flip) ifelse(d$beta_x < 0, -1, 1) else rep(1, nrow(d))
    bx <- d$beta_x * s_x; by <- d$beta_y * s_x; R_use <- R * outer(s_x, s_x)
    tryCatch(run_gls_mr(bx, by, d$se_y, R_use, intercept = TRUE), error = function(e) rep(NA_real_, 5))
  }
  f1 <- fit(TRUE); f0 <- fit(FALSE)
  data.table(gene_id = gene, peak_id = peak, n_snp = nrow(d), n_neg_x = sum(d$beta_x < 0), status = "ok",
             egger_beta_flip = f1[1], egger_se_flip = f1[2], egger_p_flip = f1[3], egger_alpha_flip = f1[4], egger_p_alpha_flip = f1[5],
             egger_beta_noflip = f0[1], egger_se_noflip = f0[2], egger_p_noflip = f0[3], egger_alpha_noflip = f0[4], egger_p_alpha_noflip = f0[5])
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 读数据: pair = 正式 MR 结果; QTL 表 fork 前缩到这些 pair

t_start <- Sys.time()
cat(sprintf("===== %s | cores = %d | chunks = %d | TEST_N = %d =====\n", CELL_TYPE, NCORES, N_CHUNK, TEST_N))
mr <- fread(MR_FILE, select = c("gene_id", "peak_id", "ivw_n_snp", "ivw_beta", "ivw_fdr", "egger_beta", "egger_p", "egger_fdr"))
pairs <- mr[, .(gene_id, peak_id)]
setorder(pairs, peak_id, gene_id)
if (TEST_N > 0) { pairs <- head(pairs, TEST_N); cat(sprintf("*** TEST_N = %d ***\n", TEST_N)) }
cat(sprintf("pairs: %d\n", nrow(pairs)))

ca_all  <- fread(CAQTL_FILE)[, .(phenotype_id, variant_id, gene_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se)]
ca_data <- ca_all[unique(pairs[, .(phenotype_id = peak_id, gene_id)]), on = .(phenotype_id, gene_id), nomatch = 0L]
rm(ca_all); gc()
eq_all  <- fread(EQTL_FILE, select = c("phenotype_id", "variant_id", "A1(ALT/effect allele)", "A2(REF)", "slope", "slope_se"))
setkey(eq_all, phenotype_id, variant_id)
u       <- unique(merge(unique(ca_data[, .(phenotype_id, gene_id, variant_id)]), pairs, by.x = c("phenotype_id", "gene_id"), by.y = c("peak_id", "gene_id"))[, .(gene_id, variant_id)])
eq_data <- eq_all[.(u$gene_id, u$variant_id), nomatch = 0L]        # 多等位位点的多条记录全部保留
rm(eq_all, u); gc()
setkey(ca_data, phenotype_id, gene_id); setkey(eq_data, phenotype_id, variant_id)
cat(sprintf("caQTL 子集: %d 行 | eQTL 子集: %d 行\n", nrow(ca_data), nrow(eq_data)))

setDTthreads(1)
pairs[, chunk := (seq_len(.N) - 1L) %% N_CHUNK + 1L]

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 分块运行 (chunk 已存在则跳过)

chunk_file <- function(ch) file.path(CHUNK_DIR, sprintf("chunk_%03d.rds", ch))
run_chunk <- function(ch) {
  t0  <- Sys.time()
  tmp <- file.path(PLINK_TMP, sprintf("%s_chunk%03d", TAG, ch)); dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  pp  <- pairs[chunk == ch]
  out <- rbindlist(lapply(seq_len(nrow(pp)), function(i)
    tryCatch(run_pair_flip(pp$gene_id[i], pp$peak_id[i], ca_data, eq_data, tmp),
             error = function(e) data.table(gene_id = pp$gene_id[i], peak_id = pp$peak_id[i], status = paste("error:", conditionMessage(e))))), fill = TRUE)
  unlink(tmp, recursive = TRUE)
  f <- chunk_file(ch); saveRDS(out, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  cat(sprintf("  chunk %03d: %d pairs, %.1f min\n", ch, nrow(pp), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(ch)
}
todo <- which(!file.exists(chunk_file(seq_len(N_CHUNK))))
cat(sprintf("待跑 chunks: %d / %d\n", length(todo), N_CHUNK))
if (length(todo)) {
  ret <- mclapply(todo, run_chunk, mc.cores = NCORES, mc.preschedule = FALSE)
  bad <- vapply(ret, inherits, logical(1), "try-error")
  if (any(bad)) { cat("失败的 chunk:\n"); print(todo[bad]); print(unique(unlist(ret[bad]))) }
}
miss <- which(!file.exists(chunk_file(seq_len(N_CHUNK))))
if (length(miss)) { cat(sprintf("还有 %d 个 chunk 未完成, 重投即可续跑\n", length(miss))); quit(status = 1) }

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 合并 -> 对账 -> 不 flip 版 gene 内 BH -> 落盘

setDTthreads(NCORES)
res <- rbindlist(lapply(chunk_file(seq_len(N_CHUNK)), readRDS), fill = TRUE)
cat("pair status 分布:\n"); print(res[, .N, by = status][order(-N)])
res <- merge(res, mr, by = c("gene_id", "peak_id"), all.y = (TEST_N == 0))
stopifnot(nrow(res) == nrow(pairs))

cat("\n===== 对账: flip 版 Egger 与正式表\n")
print(res[status == "ok", .(n = .N, n_snp_match = sum(n_snp == ivw_n_snp), max_abs_d_beta = max(abs(egger_beta_flip - egger_beta), na.rm = TRUE),
                            max_abs_d_log10p = max(abs(log10(egger_p_flip) - log10(egger_p)), na.rm = TRUE))])

res[, egger_fdr_noflip := p.adjust(egger_p_noflip, method = "BH"), by = gene_id]
res[, egger_fdr_flip := egger_fdr]
out_cols <- c("gene_id", "peak_id", "status", "n_snp", "n_neg_x", "ivw_beta", "ivw_fdr",
              "egger_beta_flip", "egger_se_flip", "egger_p_flip", "egger_fdr_flip", "egger_alpha_flip", "egger_p_alpha_flip",
              "egger_beta_noflip", "egger_se_noflip", "egger_p_noflip", "egger_fdr_noflip", "egger_alpha_noflip", "egger_p_alpha_noflip")
res <- res[, ..out_cols]

cat("\n===== Egger 显著 (p.adj < 0.05): 不 flip (行) x flip (列)\n")
print(res[status == "ok", table(noflip_sig = egger_fdr_noflip < 0.05, flip_sig = egger_fdr_flip < 0.05, useNA = "ifany")])

out_file <- file.path(OUT_DIR, paste0(TAG, "_Egger_FlipScan.tsv.gz"))
fwrite(res, out_file, sep = "\t")
if (TEST_N == 0) unlink(CHUNK_DIR, recursive = TRUE)
cat(sprintf("\nSaved: %s (%d pairs) | %.1f min\n", out_file, nrow(res), as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
