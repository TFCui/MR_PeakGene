# 2026-10-15
# Valid7_GO 第 1 步: MR-Only gene 的 GO BP GSEA (观测)
# MR-Only link: IVW / Egger / SMR 任一 FDR < 0.1, 且 WNN、Pseudo 的 FDR 都 >= 0.1
# 4 组 gene set: Bn_TCL1A / CD4_Tn_CCR7 / CD8_Tn_CCR7 的 MR-Only gene, 以及 3 个 CT 的并集 (Combined)
# geneList = MR-Only gene, score = -log10(该 gene 的 MR-Only link 中最小的 MR p), MR p = min(ivw_p, egger_p, smr_p); Combined 取 3 个 CT 中的最小 p
# gseGO: ont = BP, exponent = 1 (score 不截断), scoreType = "pos", minGSSize 10, maxGSSize 500, eps = 0, 保存全部 term
# 另存第 2 步 (随机) 要用的 geneList 和随机池 (tested gene, score = -log10(全部 link 中最小的 MR p); Combined 取 3 个 CT 中的最小)
# 输出: Results/Validation/Valid7_GO/GO_GSEA_MROnly.tsv.gz (第 3 步加上经验 P 后覆盖); tmp/Valid7_GO/genelists.rds

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(BiocParallel)
})
options(width = 300)
register(SerialParam())

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
TMP_DIR  <- file.path(PROJ, "tmp/Valid7_GO")
OUT_DIR  <- file.path(PROJ, "Results/Validation/Valid7_GO")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

celltypes <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
gsets     <- c(celltypes, "Combined")
FDR_CUT   <- 0.1
EXPONENT  <- 1
SEED      <- 2026
n_cores   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 1. geneList (MR-Only gene) 和随机池 (tested gene)

gene_scores <- function(ct) {
  mr <- fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")))
  mr[, mr_p := pmin(ivw_p, egger_p, smr_p)]
  mr[, mr_only := (ivw_fdr < FDR_CUT | egger_fdr < FDR_CUT | smr_fdr < FDR_CUT) & wnn_fdr >= FDR_CUT & pseudo_fdr >= FDR_CUT]
  list(obs    = mr[mr_only == TRUE, .(p = min(mr_p)), by = gene_id],
       pool   = mr[, .(p = min(mr_p)), by = gene_id],
       n_link = sum(mr$mr_only))
}
to_list <- function(dt) sort(setNames(-log10(pmax(dt$p, 1e-300)), dt$gene_id), decreasing = TRUE)

l <- setNames(lapply(celltypes, gene_scores), celltypes)
l$Combined <- list(obs    = rbindlist(lapply(l[celltypes], `[[`, "obs"))[,  .(p = min(p)), by = gene_id],
                   pool   = rbindlist(lapply(l[celltypes], `[[`, "pool"))[, .(p = min(p)), by = gene_id],
                   n_link = sum(vapply(l[celltypes], `[[`, numeric(1), "n_link")))
genelists <- lapply(l, function(x) list(obs = to_list(x$obs), pool = to_list(x$pool), n_link = x$n_link))
for (s in gsets) {
  g <- genelists[[s]]
  stopifnot(all(names(g$obs) %chin% names(g$pool)))
  cat(sprintf("%-12s MR-Only link %6d | MR-Only gene %5d | tested gene (随机池) %6d | 观测 score 中位 %.2f [%.2f, %.2f] | 随机池 score 中位 %.2f [%.2f, %.2f]\n",
              s, g$n_link, length(g$obs), length(g$pool), median(g$obs), min(g$obs), max(g$obs), median(g$pool), min(g$pool), max(g$pool)))
}
saveRDS(genelists, file.path(TMP_DIR, "genelists.rds"))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 2. gseGO: 4 组

run_one <- function(s) {
  gl <- genelists[[s]]$obs
  set.seed(SEED)
  gse <- suppressWarnings(gseGO(geneList = gl, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", exponent = EXPONENT,
                                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0,
                                seed = FALSE, verbose = FALSE, scoreType = "pos"))
  res <- as.data.table(as.data.frame(gse))
  res[, `:=`(gene_set = s, n_genes = length(gl))]
  res
}

t0 <- Sys.time()
out <- vector("list", length(gsets))
out[1]  <- list(run_one(gsets[1]))   # 主进程先跑一个, 缓存 GO 注释
out[-1] <- mclapply(gsets[-1], run_one, mc.cores = n_cores, mc.preschedule = FALSE)
err <- vapply(out, inherits, logical(1), "try-error")
if (any(err)) { print(out[err]); stop("gseGO 出错") }
cat(sprintf("\ngseGO 完成: %d 组, %.1f min\n", length(out), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

res <- rbindlist(out, use.names = TRUE)
setcolorder(res, c("gene_set", "n_genes"))
res <- res[order(match(gene_set, gsets), pvalue)]

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 3. 汇总 + 输出

cat("\n每组: term 数 / 最小 p / p < 0.05 的 term 数 (期望 5%) / FDR < 0.1 / FDR < 0.05\n")
print(res[, .(n_genes = n_genes[1], n_term = .N, min_p = signif(min(pvalue), 2), p_lt_05 = sum(pvalue < 0.05), exp_05 = round(0.05 * .N),
              fdr_lt_01 = sum(p.adjust < 0.1), fdr_lt_005 = sum(p.adjust < 0.05), min_fdr = signif(min(p.adjust), 2)), by = gene_set])

cat("\n每组 p 最小的 10 个 term\n")
print(res[, head(.SD, 10), by = gene_set, .SDcols = c("Description", "setSize", "NES", "pvalue", "p.adjust")][
  , `:=`(Description = substr(Description, 1, 60), NES = round(NES, 2), pvalue = signif(pvalue, 2), p.adjust = signif(p.adjust, 2))][], nrows = 200)

fwrite(res, file.path(OUT_DIR, "GO_GSEA_MROnly.tsv.gz"), sep = "\t")
cat("\n写出:", file.path(OUT_DIR, "GO_GSEA_MROnly.tsv.gz"), "|", nrow(res), "行\n")
