# QC_OneSampleMR 第 2 步: null SNP 上 caQTL z 与 eQTL z 的相关性 (流程见 Scripts/Pipeline.md 第 4.1 节)
# 用法: Rscript QC_OneSampleMR.R <CT>
# 输入: caQTL / eQTL 全量; All_GenePeak_<CT>.rds; tmp/QC_OneSampleMR/<CT>_independent_variants.txt
# 输出: Results/QC/<CT>_pair_corr.tsv.gz, <CT>_global_rho.tsv

suppressPackageStartupMessages(library(data.table))
options(width = 200)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

CT   <- commandArgs(TRUE)[1]
BASE <- "/scratch/user/tfcui2025/CIMA"
PROJ <- file.path(BASE, "Claude_project")
TMP  <- file.path(PROJ, "tmp/QC_OneSampleMR")
OUT  <- file.path(PROJ, "Results/QC"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
D    <- file.path(BASE, "Data/OMIX009786")
NTHR <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")); setDTthreads(NTHR)
P_NULL  <- 0.05     # null SNP: caQTL 和 eQTL 的 p 都大于它
DUP_THR <- 1e-3     # 去重阈值 (正式结果)
OLD_THR <- 0.01     # 旧脚本阈值, 只用于对照

# ---- 1. 读 pruned 且 p > P_NULL 的 SNP (awk 流式过滤, 不把全量读入内存) ----
keep_file <- file.path(TMP, paste0(CT, "_independent_variants.txt"))
stopifnot(file.exists(keep_file))
read_null <- function(f, idname) {
  cmd <- sprintf("zcat %s | awk -F'\\t' -v P=%s 'NR==FNR{k[$1]=1; next} FNR>1 && ($7+0)>P && ($2 in k)' %s - | cut -f1,2,8,9,11,12",
                 f, P_NULL, keep_file)
  x <- fread(cmd = cmd, header = FALSE, sep = "\t", col.names = c(idname, "variant_id", "slope", "se", "A1", "A2"))
  x
}
msg(CT, ": 读 caQTL / eQTL (pruned, p > ", P_NULL, ")")
ca <- read_null(file.path(D, "cis_caQTL_Summary_Statics/cis_caQTL", paste0(CT, ".txt.gz")), "peak")
eq <- read_null(file.path(D, "cis_eQTL_Summary_Statics/cis_eQTL", paste0(CT, ".txt.gz")), "gene")
msg("caQTL 行 = ", nrow(ca), " | eQTL 行 = ", nrow(eq))

# ---- 2. 候选 pair 配对, 按 variant_id 求共享 SNP ----
gp <- as.data.table(readRDS(file.path(D, "gene_peak", paste0("All_GenePeak_", CT, ".rds"))))
pairs <- unique(gp[, .(gene = gene_id, peak = sub("-", ":", peak_id))])   # chr1-9928-10429 -> chr1:9928-10429
setnames(ca, c("slope", "se", "A1", "A2"), c("ca_slope", "ca_se", "ca_A1", "ca_A2"))
setnames(eq, c("slope", "se", "A1", "A2"), c("eq_slope", "eq_se", "eq_A1", "eq_A2"))
x <- merge(pairs, ca, by = "peak", allow.cartesian = TRUE)
x <- merge(x, eq, by = c("gene", "variant_id"))
msg("pair x SNP 行 = ", nrow(x), " | pair = ", uniqueN(x[, .(gene, peak)]))

# 等位一致 (同向直接用; 反向则翻转 eQTL slope; 其他丢掉)
same <- x$ca_A1 == x$eq_A1 & x$ca_A2 == x$eq_A2
flip <- x$ca_A1 == x$eq_A2 & x$ca_A2 == x$eq_A1
msg(sprintf("等位: 一致 %d | 反向 %d | 不一致丢弃 %d", sum(same), sum(flip), sum(!same & !flip)))
x[flip, eq_slope := -eq_slope]
x <- x[same | flip]
x[, `:=`(ca_z = ca_slope / ca_se, eq_z = eq_slope / eq_se)]
x <- x[is.finite(ca_z) & is.finite(eq_z)]

# ---- 3. 去重 + 相关性 ----
dedup <- function(x, thr) {
  y <- unique(x, by = c("peak", "gene", "ca_slope", "eq_slope"))
  setorder(y, peak, gene, ca_slope, eq_slope)
  y[, `:=`(d_ca = abs(ca_slope - shift(ca_slope)), d_eq = abs(eq_slope - shift(eq_slope))), by = .(peak, gene)]
  y <- y[is.na(d_ca) | d_ca > thr | d_eq > thr]
  y[, c("d_ca", "d_eq") := NULL]
  y
}
tp <- function(r, n) {   # 相关系数 -> t 检验的 p 与 -log10 p
  t <- r * sqrt((n - 2) / (1 - r^2)); lp <- pt(-abs(t), df = n - 2, log.p = TRUE) + log(2)
  list(t_stat = t, p_value = exp(lp), neglog10p = -lp / log(10))
}
pair_cor <- function(y) {
  r <- y[, .(snp_count = .N, rho_hat = if (.N > 2) cor(ca_z, eq_z) else NA_real_), by = .(peak, gene)]
  r[, c("t_stat", "p_value", "neglog10p") := tp(rho_hat, snp_count)]
  r
}
global_cor <- function(y) {
  r <- y[, cor(ca_z, eq_z)]; n <- nrow(y); c(rho = r, n_snp = n, tp(r, n)[c("p_value", "neglog10p")] |> unlist())
}

y  <- dedup(x, DUP_THR); y_old <- dedup(x, OLD_THR)
msg(sprintf("去重后 SNP 行: %d (阈值 %g) | %d (阈值 %g)", nrow(y), DUP_THR, nrow(y_old), OLD_THR))
pc <- pair_cor(y)
g  <- global_cor(y); g_old <- global_cor(y_old)
gl <- data.table(cell_type = CT,
                 rho = g["rho"], p_value = g["p_value"], neglog10p = g["neglog10p"], n_snp = g["n_snp"],
                 n_pair = nrow(y[, .N, by = .(peak, gene)]), n_pair_cor = sum(!is.na(pc$rho_hat)),
                 median_pair_rho = median(pc$rho_hat, na.rm = TRUE), dup_thr = DUP_THR,
                 rho_dup0.01 = g_old["rho"], p_dup0.01 = g_old["p_value"], n_snp_dup0.01 = g_old["n_snp"])
print(gl)
fwrite(pc[order(peak, gene)], file.path(OUT, paste0(CT, "_pair_corr.tsv.gz")), sep = "\t")
fwrite(gl, file.path(OUT, paste0(CT, "_global_rho.tsv")), sep = "\t")
msg("完成")
