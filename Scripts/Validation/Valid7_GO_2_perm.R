# 2026-10-15
# Valid7_GO 第 2 步: 随机抽样 GSEA (经验 P 的零分布)
# 对一组 gene set (Bn_TCL1A / CD4_Tn_CCR7 / CD8_Tn_CCR7 / Combined):
#   第 b 次 (b = 1..N_PERM): 用种子 SEED + b 从随机池 (tested gene) 中抽与 MR-Only gene 相同个数的 gene,
#   score = 第 1 步存的随机池 score (-log10(全部 link 中最小的 MR p)), 做与第 1 步相同的 GSEA, 记录每个 term 的 NES
# 提速: 不调用 gseGO, GO BP 注释 (与 gseGO 内部同一份, clusterProfiler:::get_GO_data) 只取一次, 每次直接调 fgsea
#   fgsea 参数与 gseGO 一致: minSize 10, maxSize 500, gseaParam = 1 (exponent), scoreType = "pos", eps = 0
# 开始前核对: 观测 geneList 直接调 fgsea 的 NES 与第 1 步 gseGO 的 NES 一致
# 输入: tmp/Valid7_GO/genelists.rds, Results/Validation/Valid7_GO/GO_GSEA_MROnly.tsv.gz (第 1 步)
# 输出: tmp/Valid7_GO/perm/<set>_perm_<from>-<to>.rds (perm, ID, NES); 每 BATCH 次一个文件, 已存在跳过 (重投即续跑)
# 用法: Rscript Valid7_GO_2_perm.R <set>; 环境变量 N_PERM (默认 5000)

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(fgsea)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(BiocParallel)
})
options(width = 300)
register(SerialParam())

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
TMP_DIR  <- file.path(PROJ, "tmp/Valid7_GO")
PERM_DIR <- file.path(TMP_DIR, "perm")
OBS_FILE <- file.path(PROJ, "Results/Validation/Valid7_GO/GO_GSEA_MROnly.tsv.gz")
dir.create(PERM_DIR, recursive = TRUE, showWarnings = FALSE)

EXPONENT <- 1
SEED     <- 2026
N_PERM   <- as.integer(Sys.getenv("N_PERM", "5000"))
BATCH    <- 250
n_cores  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))

s <- commandArgs(trailingOnly = TRUE)[1]
genelists <- readRDS(file.path(TMP_DIR, "genelists.rds"))
stopifnot(s %in% names(genelists))
gl   <- genelists[[s]]$obs
pool <- genelists[[s]]$pool
n    <- length(gl)
cat(sprintf("%s: MR-Only gene %d, 随机池 %d, N_PERM = %d, BATCH = %d, cores = %d\n", s, n, length(pool), N_PERM, BATCH, n_cores))

# GO BP 注释 (term -> SYMBOL), 与 gseGO 内部相同
pathways <- DOSE:::getGeneSet(clusterProfiler:::get_GO_data(org.Hs.eg.db, "BP", "SYMBOL"))
cat(sprintf("GO BP term: %d\n", length(pathways)))

run_fgsea <- function(stats, seed) {
  set.seed(seed)
  r <- suppressWarnings(fgsea(pathways = pathways, stats = stats, minSize = 10, maxSize = 500, eps = 0,
                              gseaParam = EXPONENT, scoreType = "pos", BPPARAM = SerialParam(progressbar = FALSE)))
  r[!is.na(NES), .(ID = pathway, NES)]
}

## 核对: 观测 geneList 直接调 fgsea vs 第 1 步 gseGO
obs <- fread(OBS_FILE)[gene_set == s, .(ID, NES_gseGO = NES)]
chk <- merge(obs, run_fgsea(gl, SEED), by = "ID", all = TRUE)
cat(sprintf("核对 fgsea vs gseGO: gseGO term %d, fgsea term %d, 共同 %d, NES 最大差 %.2e\n",
            sum(!is.na(chk$NES_gseGO)), sum(!is.na(chk$NES)), sum(!is.na(chk$NES_gseGO) & !is.na(chk$NES)),
            max(abs(chk$NES_gseGO - chk$NES), na.rm = TRUE)))
stopifnot(!anyNA(chk$NES_gseGO), !anyNA(chk$NES), max(abs(chk$NES_gseGO - chk$NES)) < 1e-6)

## 随机
t0 <- Sys.time()
for (from in seq(1, N_PERM, by = BATCH)) {
  to <- min(from + BATCH - 1, N_PERM)
  f  <- file.path(PERM_DIR, sprintf("%s_perm_%05d-%05d.rds", s, from, to))
  if (file.exists(f)) { cat("跳过已存在:", basename(f), "\n"); next }
  t1  <- Sys.time()
  res <- mclapply(from:to, function(b) {
    set.seed(SEED + b)
    g <- sample(names(pool), n)
    r <- run_fgsea(sort(pool[g], decreasing = TRUE), SEED + b)
    r[, perm := b]
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  err <- vapply(res, inherits, logical(1), "try-error")
  if (any(err)) { print(res[err][1]); stop("随机 GSEA 出错") }
  res <- rbindlist(res)
  setcolorder(res, "perm")
  saveRDS(res, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  cat(sprintf("perm %d-%d 完成: %.1f min (每次 %.1f 核·秒)\n", from, to, as.numeric(difftime(Sys.time(), t1, units = "mins")),
              as.numeric(difftime(Sys.time(), t1, units = "secs")) * n_cores / (to - from + 1)))
}
cat(sprintf("%s 完成, 共 %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
