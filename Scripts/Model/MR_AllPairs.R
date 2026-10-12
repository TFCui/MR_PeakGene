# 2026-10-09
# 全量 peak-gene pairs 的 MR (GLS-IVW / GLS-Egger / SMR / CCT) + 边际相关 (WNN / pseudobulk Spearman)
# 流程说明见 Scripts/Pipeline.md 第 2 节; 移植自 Scripts_reference/Server5_wgen6_MR9.sh, 核心统计函数不变
# 用法: Rscript MR_AllPairs.R <CELL_TYPE>
#   pairs 分 N_CHUNK 块, mclapply 并行; chunk 文件已存在则跳过 -> 重投同一作业即可断点续跑
#   PLINK 临时文件放 Claude_project/tmp/plink/, 每个 chunk 独立子目录, 用完即删
#   全部 chunk 完成后: 合并 -> 删 NA -> ivw_sigma2 <= 10 -> CCT -> 按 gene 做 BH -> 落盘
# 环境变量 TEST_N > 0 时只跑前 N 个 pair (冒烟测试), 结果文件名带 _TEST<N>, 不覆盖正式结果

library(data.table)
library(parallel)
options(width = 300)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 配置

args      <- commandArgs(trailingOnly = TRUE)
CELL_TYPE <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "Bn_TCL1A"

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
OUT_DIR  <- file.path(PROJ, "Results/MR")
TMP_DIR  <- file.path(PROJ, "tmp")
PLINK    <- file.path(BASE_DIR, "R_package/my_plink/plink")
REF_DIR  <- file.path(PROJ, "data/ref_chr/hg38")                   # 软链接, 只读

EQTL_FILE      <- file.path(BASE_DIR, "Data/OMIX009786/cis_eQTL_Summary_Statics/cis_eQTL", paste0(CELL_TYPE, ".txt.gz"))
CAQTL_FILE     <- file.path(BASE_DIR, "Data/OMIX009786/cis_caQTL_Summary_Statics/cis_caQTL/LD_clump", paste0(CELL_TYPE, "_Clumped_r2_0.99_kb_20000.txt.gz"))
GENE_PEAK_FILE <- file.path(BASE_DIR, "Data/OMIX009786/gene_peak", paste0("All_GenePeak_", CELL_TYPE, ".rds"))
PB_DIR         <- file.path(BASE_DIR, "Data/OMIX009786/cell_rna/QC_pseudobulk")
WNN_FILE       <- file.path(BASE_DIR, "Data/Valid_Multiome/1_pbmc_multiomics.reprocessed.rda")

# 固定参数
MIN_SNP      <- 3
DEDUP_TOL    <- 1e-3
KEEP_ALLELE  <- TRUE
FLIP_POINTS  <- TRUE
EIG_TRUNCATE <- TRUE
EIG_VAR_KEEP <- 0.999
EIG_MIN_KEEP <- 5
EIG_TOL      <- 1e-10
P_FLOOR      <- 1e-300
SIGMA_MODE   <- "over"                                            # se * max(1, sigma)
SIGMA2_MAX   <- 10                                                # 只保留 ivw_sigma2 <= 10
BASE_P_COLS  <- c("ivw_p", "egger_p", "smr_p", "wnn_p", "pseudo_p")
CCT_COLS     <- c("smr_p", "ivw_p", "egger_p")

NCORES  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "36"))
N_CHUNK <- 144
TEST_N  <- as.integer(Sys.getenv("TEST_N", "0"))
TAG     <- if (TEST_N > 0) paste0(CELL_TYPE, "_TEST", TEST_N) else CELL_TYPE

CHUNK_DIR <- file.path(TMP_DIR, "MR_chunks", TAG)
PLINK_TMP <- file.path(TMP_DIR, "plink")
dir.create(CHUNK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLINK_TMP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 函数 (与 Scripts_reference/Server5_wgen6_MR9.sh 相同)

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
run_gls_mr <- function(bx, by, se_y, R, intercept = FALSE, truncate = EIG_TRUNCATE, var_keep = EIG_VAR_KEEP,
                       min_keep = EIG_MIN_KEEP, tol = EIG_TOL, sigma_mode = SIGMA_MODE, eg = NULL, verbose = FALSE) {
  n <- length(bx)
  stopifnot(length(by) == n, length(se_y) == n, nrow(R) == n, ncol(R) == n, !anyNA(R), all(se_y > 0))

  if (is.null(eg)) eg <- eigen(outer(se_y, se_y) * R, symmetric = TRUE)
  lam <- eg$values; V <- eg$vectors
  pos  <- lam > tol * lam[1]
  cum  <- cumsum(pmax(lam, 0)) / sum(pmax(lam, 0))
  keep <- if (truncate) (c(TRUE, cum[-n] < var_keep) | seq_len(n) <= min_keep) & pos else pos
  k <- sum(keep); p <- if (intercept) 2L else 1L
  df <- n - p
  if (verbose) cat(sprintf("    eigen(Sigma): n = %d, 保留 k = %d, df = %d\n", n, k, df))
  if (k < p)  stop(sprintf("保留方向数不足: k = %d, p = %d", k, p))
  if (df < 1) stop(sprintf("自由度不足: df = %d, p = %d", df, p))

  d_half <- numeric(n); d_half[keep] <- lam[keep]^(-0.5)
  W  <- V %*% diag(d_half, n) %*% t(V)
  X  <- if (intercept) cbind(alpha = 1, beta = bx) else cbind(beta = bx)
  ys <- as.numeric(W %*% by); Xs <- W %*% X

  fit      <- lm(ys ~ Xs - 1)
  est      <- unname(coef(fit))
  se_fixed <- sqrt(diag(solve(crossprod(Xs))))
  Q        <- sum(residuals(fit)^2)
  sigma    <- sqrt(Q / df)
  infl     <- switch(sigma_mode, over = max(1, sigma), under = 1 / min(1, sigma), none = 1)
  se       <- se_fixed * infl
  pv       <- 2 * pt(abs(est / se), df = df, lower.tail = FALSE)

  data.table(model = if (intercept) "Egger" else "IVW", n_snp = n, n_eff = k, df = df, sigma2 = sigma^2, se_infl = infl,
             alpha = if (intercept) est[1] else NA_real_, se_alpha = if (intercept) se[1] else NA_real_, p_alpha = if (intercept) pv[1] else NA_real_,
             beta = est[p], se_beta = se[p], p_beta = pv[p])
}

# SMR: top SNP = caQTL nominal p 最小 (|Z| tie-break); z 由 p 反推
run_smr_pair <- function(d) {
  if (nrow(d) < 1) return(NULL)
  i    <- order(d$p_x, -abs(d$beta_x / d$se_x))[1]
  p_x  <- pmax(d$p_x[i], P_FLOOR); p_y <- pmax(d$p_y[i], P_FLOOR)
  z_x  <- qnorm(p_x / 2, lower.tail = FALSE); z_y <- qnorm(p_y / 2, lower.tail = FALSE)
  b    <- d$beta_y[i] / d$beta_x[i]
  seb  <- abs(b) * sqrt(1 / z_y^2 + 1 / z_x^2)
  chi  <- (z_x^2 * z_y^2) / (z_x^2 + z_y^2)
  data.table(smr_variant_id = d$variant_id[i], smr_snp_id = d$snp_id[i], smr_b = b, smr_se = seb,
             smr_chisq = chi, smr_p = pchisq(chi, df = 1, lower.tail = FALSE),
             smr_z_caqtl = z_x, smr_z_eqtl = z_y, smr_se_caqtl = abs(d$beta_x[i]) / z_x, smr_se_eqtl = abs(d$beta_y[i]) / z_y)
}

cauchy_combine_matrix <- function(P, weights = NULL, tail_cutoff = 1e-15) {
  P <- as.matrix(P); storage.mode(P) <- "double"
  n <- nrow(P); k <- ncol(P)
  if (n == 0L) return(numeric(0))
  if (is.null(weights)) weights <- rep(1, k)
  W <- matrix(rep(weights, each = n), n, k)
  valid <- !is.na(P) & P >= 0 & P <= 1
  W[!valid] <- 0
  wsum <- rowSums(W); W <- W / ifelse(wsum > 0, wsum, 1)

  small  <- valid & P < tail_cutoff
  large  <- valid & P > 1 - tail_cutoff
  middle <- valid & !small & !large
  Z <- matrix(0, n, k)
  Z[small]  <-  1 / (pi * pmax(P[small], .Machine$double.xmin))
  Z[large]  <- -1 / (pi * pmax(1 - P[large], .Machine$double.xmin))
  Z[middle] <- tan((0.5 - P[middle]) * pi)

  Tstat <- rowSums(W * Z); out <- rep(NA_real_, n)
  ok <- wsum > 0; pos <- ok & Tstat > 0; neg <- ok & !(Tstat > 0)
  out[pos] <- atan(1 / Tstat[pos]) / pi
  out[neg] <- 0.5 - atan(Tstat[neg]) / pi
  pmin(1, pmax(0, out))
}

run_one_pair <- function(gene, peak, ca_data, eq_data, tmp_dir, verbose = FALSE) {
  pid <- paste(peak, gene, sep = "__")
  mk_st <- function(n_qtl, n_ld, status) data.table(gene_id = gene, peak_id = peak, n_qtl = n_qtl, n_ld = n_ld, status = status)

  ### Step 1: QTL 对齐
  ca <- ca_data[.(peak, gene), nomatch = 0L]
  eq <- eq_data[.(gene), nomatch = 0L]
  if (!nrow(ca) || !nrow(eq)) return(list(res = NULL, smr = NULL, status = mk_st(0L, NA_integer_, "no_qtl")))

  ca[, snp_id := make_plink_id(variant_id, `A2(REF)`, `A1(ALT/effect allele)`)]
  ca <- ca[, .(variant_id, snp_id, ea_x = `A1(ALT/effect allele)`, oa_x = `A2(REF)`, beta_x = slope, se_x = slope_se, p_x = pval_nominal)]
  eq <- eq[, .(variant_id,          ea_y = `A1(ALT/effect allele)`, oa_y = `A2(REF)`, beta_y = slope, se_y = slope_se, p_y = pval_nominal)]

  dat <- merge(ca, eq, by = "variant_id")
  dat <- dat[ea_x == ea_y & oa_x == oa_y]
  dat <- dat[dedup_keep(beta_x, beta_y, tol = DEDUP_TOL)]
  n_qtl <- nrow(dat)
  if (n_qtl < MIN_SNP) return(list(res = NULL, smr = NULL, status = mk_st(n_qtl, NA_integer_, "few_snp_qtl")))

  ### Step 2: LD
  pre <- file.path(tmp_dir, gsub("[^A-Za-z0-9]", "_", pid))
  ld  <- run_plink_ld(dat$snp_id, ref_bfile(peak), pre)
  if (is.null(ld)) return(list(res = NULL, smr = NULL, status = mk_st(n_qtl, NA_integer_, "plink_fail")))

  d <- dat[match(rownames(ld), snp_id)]
  d <- d[!is.na(variant_id)]
  R <- ld[d$snp_id, d$snp_id, drop = FALSE]
  if (nrow(d) < MIN_SNP) return(list(res = NULL, smr = NULL, status = mk_st(n_qtl, nrow(d), "few_snp_ld")))
  if (anyNA(R)) {
    bad <- which(apply(R, 1, function(x) any(!is.finite(x))))
    d <- d[-bad]; R <- R[-bad, -bad, drop = FALSE]
    if (nrow(d) < MIN_SNP) return(list(res = NULL, smr = NULL, status = mk_st(n_qtl, nrow(d), "few_snp_nan")))
  }

  ### Step 3: flip
  s_x   <- if (FLIP_POINTS) ifelse(d$beta_x < 0, -1, 1) else rep(1, nrow(d))
  bx    <- d$beta_x * s_x; by <- d$beta_y * s_x; se_y <- d$se_y
  R_use <- R * outer(s_x, s_x)

  ### Step 4: GLS IVW / Egger + SMR
  eg <- tryCatch(eigen(outer(se_y, se_y) * R_use, symmetric = TRUE), error = function(e) NULL)
  smr <- run_smr_pair(d)
  smr <- if (is.null(smr)) NULL else cbind(data.table(gene_id = gene, peak_id = peak), smr)
  if (is.null(eg)) return(list(res = NULL, smr = smr, status = mk_st(n_qtl, nrow(d), "eigen_fail")))

  res <- list()
  for (mod in c("IVW", "Egger")) {
    r <- tryCatch(run_gls_mr(bx, by, se_y, R_use, intercept = (mod == "Egger"), eg = eg, verbose = verbose), error = function(e) NULL)
    if (!is.null(r)) res[[mod]] <- cbind(data.table(gene_id = gene, peak_id = peak), r)
  }

  list(res = rbindlist(res), smr = smr,
       status = mk_st(n_qtl, nrow(d), if (length(res) == 2L) "ok" else if (length(res)) "partial" else "gls_fail"))
}

read_pb <- function(path) {
  dt <- fread(path); m <- as.matrix(dt[, -1, with = FALSE]); rownames(m) <- dt[[1]]; m
}

# 边际相关: Spearman, 不 residualize, 无 covariate; 返回 n / rho / p
spearman_pairs <- function(pairs, mat_gene, mat_peak, gene_in_row = TRUE, prefix, ncores = 1L) {
  gs <- if (gene_in_row) rownames(mat_gene) else colnames(mat_gene)
  ps <- if (gene_in_row) rownames(mat_peak) else colnames(mat_peak)
  ok <- pairs$gene_id %chin% gs & pairs$peak_seurat %chin% ps

  f <- function(i) {
    if (!ok[i]) return(c(NA_real_, NA_real_, NA_real_))
    vg <- if (gene_in_row) mat_gene[pairs$gene_id[i], ]    else mat_gene[, pairs$gene_id[i]]
    vp <- if (gene_in_row) mat_peak[pairs$peak_seurat[i], ] else mat_peak[, pairs$peak_seurat[i]]
    sp <- cor.test(vg, vp, method = "spearman", exact = FALSE)
    c(length(vg), unname(sp$estimate), sp$p.value)
  }
  m <- if (ncores > 1L) simplify2array(mclapply(seq_len(nrow(pairs)), f, mc.cores = ncores))
       else vapply(seq_len(nrow(pairs)), f, numeric(3))

  out <- data.table(gene_id = pairs$gene_id, peak_id = pairs$peak_id, n = as.integer(m[1, ]), rho = m[2, ], p = m[3, ])
  setnames(out, c("n", "rho", "p"), paste0(prefix, c("_n", "_rho", "_p")))
  out
}

# 日志: 当前表的 pair / gene / peak 数
report_n <- function(dt, label) cat(sprintf("  %-28s %8d pairs | gene %6d | peak %6d\n", label, nrow(dt), uniqueN(dt$gene_id), uniqueN(dt$peak_id)))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 读数据

t_start <- Sys.time()
cat(sprintf("===== %s | cores = %d | chunks = %d | TEST_N = %d =====\n", CELL_TYPE, NCORES, N_CHUNK, TEST_N))

eq_all <- fread(EQTL_FILE, select = c("phenotype_id", "variant_id", "A1(ALT/effect allele)", "A2(REF)", "slope", "slope_se", "pval_nominal"))
setkey(eq_all, phenotype_id, variant_id)

ca_all <- fread(CAQTL_FILE)
ca_all <- ca_all[, intersect(c("phenotype_id", "variant_id", "gene_id", "A1(ALT/effect allele)", "A2(REF)", "slope", "slope_se", "pval_nominal"), names(ca_all)), with = FALSE]
stopifnot("gene_id" %chin% names(ca_all))
setkey(ca_all, phenotype_id, gene_id)

gp <- unique(readRDS(GENE_PEAK_FILE)[, .(gene_id, phenotype_id = sub("-", ":", peak_id))])
cat(sprintf("eQTL: %d 行 | caQTL: %d 行 | all_gene_peak: %d pairs\n", nrow(eq_all), nrow(ca_all), nrow(gp)))

subject_peak_counts <- read_pb(file.path(PB_DIR, paste0("pseudobulk_atac_", CELL_TYPE, ".tsv.gz")))
subject_gene_counts <- read_pb(file.path(PB_DIR, paste0("pseudobulk_rna_",  CELL_TYPE, ".tsv.gz")))
common_samples      <- intersect(colnames(subject_peak_counts), colnames(subject_gene_counts))
subject_peak_counts <- subject_peak_counts[, common_samples, drop = FALSE]
subject_gene_counts <- subject_gene_counts[, common_samples, drop = FALSE]

lib <- fread(file.path(PB_DIR, paste0("libsize_", CELL_TYPE, ".csv")))
rna_lib  <- setNames(lib$rna_lib_size,  lib$sample)[common_samples]
atac_lib <- setNames(lib$atac_lib_size, lib$sample)[common_samples]
stopifnot(!anyNA(rna_lib), !anyNA(atac_lib))
subject_gene_logcpm <- log1p(sweep(subject_gene_counts, 2, rna_lib,  "/") * 1e6)
subject_peak_logcpm <- log1p(sweep(subject_peak_counts, 2, atac_lib, "/") * 1e6)
rm(subject_peak_counts, subject_gene_counts); gc()

load(WNN_FILE)                                                    # wnn.rna, wnn.peak (metacell x feature)
wnn_rna_log  <- log1p(wnn.rna)
wnn_peak_log <- log1p(wnn.peak)
cat(sprintf("pseudobulk: %d 个体 | WNN: %d metacell\n", length(common_samples), nrow(wnn_rna_log)))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 候选 pairs: all_gene_peak x caQTL x eQTL 共享 variants >= MIN_SNP, 再过滤到 WNN / pseudobulk 都有的

cand <- unique(ca_all[gp, on = .(phenotype_id, gene_id), nomatch = 0L, .(phenotype_id, gene_id, variant_id)])
hit  <- eq_all[.(cand$gene_id, cand$variant_id), which = TRUE, mult = "first"]
ov   <- cand[!is.na(hit)]

pairs_all <- ov[, .(n_variant = uniqueN(variant_id)), by = .(gene_id, peak_id = phenotype_id)][n_variant >= MIN_SNP]
pairs_all[, peak_seurat := gsub("[:-]", "-", peak_id)]

pairs <- pairs_all[gene_id %chin% rownames(subject_gene_logcpm) & peak_seurat %chin% rownames(subject_peak_logcpm) &
                   gene_id %chin% colnames(wnn_rna_log)        & peak_seurat %chin% colnames(wnn_peak_log)]
setorder(pairs, peak_id, gene_id)
cat("候选 pairs:\n")
report_n(pairs_all, "QTL 共享 variant >= 3")
report_n(pairs,     "WNN + pseudobulk 均可用")

if (TEST_N > 0) { pairs <- head(pairs, TEST_N); cat(sprintf("*** TEST_N = %d, 只跑前 %d 个 pair ***\n", TEST_N, nrow(pairs))) }
stopifnot(nrow(pairs) > 0)

# fork 前把 QTL 表缩到只含这些 pairs, 并建好 key (子进程共享, 不复制)
ca_data <- ca_all[unique(pairs[, .(phenotype_id = peak_id, gene_id)]), on = .(phenotype_id, gene_id), nomatch = 0L]
ov_sub  <- ov[pairs, on = .(gene_id, phenotype_id = peak_id), nomatch = 0L]
u       <- unique(ov_sub[, .(gene_id, variant_id)])
eq_data <- eq_all[.(u$gene_id, u$variant_id), nomatch = 0L]        # 多等位位点的多条记录全部保留
setkey(ca_data, phenotype_id, gene_id); setkey(eq_data, phenotype_id, variant_id)
rm(ca_all, eq_all, cand, hit, ov, ov_sub, u); gc()
cat(sprintf("caQTL 子集: %d 行 | eQTL 子集: %d 行\n", nrow(ca_data), nrow(eq_data)))

setDTthreads(1)
pairs[, chunk := (seq_len(.N) - 1L) %% N_CHUNK + 1L]               # 轮流分配, 平衡各块负载

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 分块跑 MR + SMR (chunk 已存在则跳过, 可断点续跑)

chunk_file <- function(ch) file.path(CHUNK_DIR, sprintf("chunk_%03d.rds", ch))

run_chunk <- function(ch) {
  t0  <- Sys.time()
  tmp <- file.path(PLINK_TMP, sprintf("%s_chunk%03d", TAG, ch))
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  pp  <- pairs[chunk == ch]
  out <- lapply(seq_len(nrow(pp)), function(i)
    tryCatch(run_one_pair(pp$gene_id[i], pp$peak_id[i], ca_data, eq_data, tmp),
             error = function(e) list(res = NULL, smr = NULL,
               status = data.table(gene_id = pp$gene_id[i], peak_id = pp$peak_id[i], n_qtl = NA_integer_, n_ld = NA_integer_,
                                   status = paste("error:", conditionMessage(e))))))
  unlink(tmp, recursive = TRUE)
  f <- chunk_file(ch)
  saveRDS(list(res    = rbindlist(lapply(out, `[[`, "res"),    fill = TRUE),
               smr    = rbindlist(lapply(out, `[[`, "smr"),    fill = TRUE),
               status = rbindlist(lapply(out, `[[`, "status"), fill = TRUE)), paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)                               # 写完再改名, 避免半成品被当成已完成
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
if (length(miss)) {
  cat(sprintf("还有 %d 个 chunk 未完成, 重投同一个作业即可续跑: %s\n", length(miss), paste(miss, collapse = ",")))
  quit(status = 1)
}

all_chunks <- lapply(chunk_file(seq_len(N_CHUNK)), readRDS)
res_mr  <- rbindlist(lapply(all_chunks, `[[`, "res"),    fill = TRUE)
res_smr <- rbindlist(lapply(all_chunks, `[[`, "smr"),    fill = TRUE)
status  <- rbindlist(lapply(all_chunks, `[[`, "status"), fill = TRUE)
rm(all_chunks); gc()
cat(sprintf("\nMR 结果: %d 行 | SMR: %d 行\n", nrow(res_mr), nrow(res_smr)))
cat("pair status 分布:\n"); print(status[, .N, by = status][order(-N)])

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### WNN / pseudobulk 边际相关

setDTthreads(NCORES)
cat("计算 pseudobulk / WNN Spearman ...\n")
res_pseudo <- spearman_pairs(pairs, subject_gene_logcpm, subject_peak_logcpm, gene_in_row = TRUE,  prefix = "pseudo", ncores = NCORES)
res_wnn    <- spearman_pairs(pairs, wnn_rna_log,         wnn_peak_log,        gene_in_row = FALSE, prefix = "wnn",    ncores = NCORES)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 合并 -> 删 NA -> ivw_sigma2 <= 10 -> CCT -> 按 gene 做 BH

ivw   <- res_mr[model == "IVW",   .(gene_id, peak_id, ivw_n_snp = n_snp, ivw_n_eff = n_eff, ivw_df = df, ivw_sigma2 = sigma2, ivw_beta = beta, ivw_se = se_beta, ivw_p = p_beta)]
egger <- res_mr[model == "Egger", .(gene_id, peak_id, egger_beta = beta, egger_se = se_beta, egger_p = p_beta)]
smr   <- res_smr[, .(gene_id, peak_id, smr_variant_id, smr_b, smr_se, smr_chisq, smr_p)]

result_all <- Reduce(function(x, y) merge(x, y, by = c("gene_id", "peak_id"), all.x = TRUE),
                     list(pairs[, .(gene_id, peak_id, n_variant)], ivw, egger, smr, res_wnn, res_pseudo))

cat("\n合并与过滤:\n")
report_n(result_all, "合并后")
cat("  各基础方法 p 为 NA 的 pair 数:\n"); print(result_all[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = BASE_P_COLS])
result_all <- result_all[complete.cases(result_all[, ..BASE_P_COLS])]
report_n(result_all, "删除 NA 后")
result_all <- result_all[ivw_sigma2 <= SIGMA2_MAX]
report_n(result_all, sprintf("ivw_sigma2 <= %g 后", SIGMA2_MAX))

result_all[, cct_n := rowSums(!is.na(as.matrix(.SD))), .SDcols = CCT_COLS]
result_all[, cct_p := cauchy_combine_matrix(as.matrix(.SD)), .SDcols = CCT_COLS]

p_cols   <- c("ivw_p", "egger_p", "smr_p", "wnn_p", "pseudo_p", "cct_p")
fdr_cols <- sub("_p$", "_fdr", p_cols)
result_all[, (fdr_cols) := lapply(.SD, p.adjust, method = "BH"), by = gene_id, .SDcols = p_cols]

out_cols <- c("gene_id", "peak_id", "n_variant",
              "ivw_n_snp", "ivw_n_eff", "ivw_df", "ivw_sigma2", "ivw_beta", "ivw_se", "ivw_p",
              "egger_beta", "egger_se", "egger_p",
              "smr_variant_id", "smr_b", "smr_se", "smr_chisq", "smr_p",
              "wnn_n", "wnn_rho", "wnn_p", "pseudo_n", "pseudo_rho", "pseudo_p",
              "cct_n", "cct_p", fdr_cols)
result_all <- result_all[, ..out_cols]
setkey(result_all, gene_id, peak_id)
stopifnot(ncol(result_all) == 32L, !anyNA(result_all[, c(p_cols, fdr_cols), with = FALSE]), all(result_all$cct_n == 3L))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 落盘 + 汇总

out_file <- file.path(OUT_DIR, paste0(TAG, "_AllPairs_results.tsv.gz"))
fwrite(result_all, out_file, sep = "\t")
unlink(CHUNK_DIR, recursive = TRUE)                               # chunk 合并后删除, 控制文件数

cat(sprintf("\n最终结果: %d pairs (gene %d, peak %d)\n", nrow(result_all), uniqueN(result_all$gene_id), uniqueN(result_all$peak_id)))
print(data.table(
  method    = sub("_p$", "", p_cols),
  sig_p05   = sapply(result_all[, ..p_cols],   function(x) sum(x < 0.05)),
  sig_fdr05 = sapply(result_all[, ..fdr_cols], function(x) sum(x < 0.05)),
  sig_fdr10 = sapply(result_all[, ..fdr_cols], function(x) sum(x < 0.10))))
cat(sprintf("\n全部完成, 用时 %.1f min\n输出: %s\n", as.numeric(difftime(Sys.time(), t_start, units = "mins")), out_file))
