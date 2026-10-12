# 2026-10-09
# Valid4_AUC: 每个 gene 内的 ROC AUC (流程见 Scripts/Pipeline.md 3.4)
# 移植自 Scripts_reference/AUC2_all.sh: pROC::roc 默认 direction, score = -log10(p), 正负例都有的 gene 才算
# 金标准读 Valid1_Hyper 保存的 <CT>_gold_standard.rds (必须先跑 Valid1); GTEx 额外要求与 validation multiome 的 peak 重叠
# 用法: Rscript Valid4_AUC.R <CT>   (array 0-2, 每个 CT 一个 task)
# 输出: Results/Validation/Valid4_AUC/<CT>_AUC_long.tsv.gz

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(IRanges)
  library(pROC)
  library(parallel)
})
options(width = 300)

args <- commandArgs(trailingOnly = TRUE)
CT   <- if (length(args) >= 1) args[1] else stop("需要提供细胞类型, 例如: Rscript Valid4_AUC.R Bn_TCL1A")
NC   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat(sprintf("[start] %s | cores = %d | %s\n", CT, NC, format(Sys.time())))

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
MR_FILE  <- file.path(PROJ, "Results/MR", paste0(CT, "_AllPairs_results.tsv.gz"))
GOLD     <- file.path(PROJ, "Results/Validation/Valid1_Hyper", paste0(CT, "_gold_standard.rds"))
OUT_DIR  <- file.path(PROJ, "Results/Validation/Valid4_AUC")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(GOLD)) stop("找不到 Valid1_Hyper 的金标准, 先跑 Valid1: ", GOLD)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 金标准: 读 Valid1 的结果, 展开成带 gene 列的 GRanges

flatten_db <- function(x) {
  gl <- if (is(x, "GRangesList")) x else GRangesList(lapply(x, granges))
  u  <- unlist(gl, use.names = FALSE)
  mcols(u) <- NULL
  u$gene <- rep(names(gl), lengths(gl))
  u
}

DB <- readRDS(GOLD)
gold.list <- lapply(DB, flatten_db)

# GTEx 额外过滤: 与 validation multiome 的 peak 重叠 (沿用旧 AUC 脚本)
load(file.path(BASE_DIR, "Data/Valid_Multiome/1_pbmc_multiomics.reprocessed.rda"))
n0 <- length(gold.list[["GTEx"]])
gold.list[["GTEx"]] <- gold.list[["GTEx"]][countOverlaps(gold.list[["GTEx"]], pbmc@assays$ATAC@ranges) > 0]
cat(sprintf("[GTEx] 与 validation multiome peak 重叠: %s / %s\n", format(length(gold.list[["GTEx"]]), big.mark = ","), format(n0, big.mark = ",")))
rm(list = setdiff(ls(), c("CT", "NC", "MR_FILE", "OUT_DIR", "gold.list"))); invisible(gc())
print(sapply(gold.list, length))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### MR 结果

res <- fread(MR_FILE)
cat(sprintf("[MR] pairs = %s, genes = %s\n", format(nrow(res), big.mark = ","), format(uniqueN(res$gene_id), big.mark = ",")))

res.gr <- StringToGRanges(res$peak_id, sep = c(":", "-"))
res.gr$gene            <- res$gene_id
res.gr$ivw.pval        <- res$ivw_p
res.gr$egger.pval      <- res$egger_p
res.gr$smr.pval        <- res$smr_p
res.gr$cct.pval        <- res$cct_p
res.gr$multi.pval      <- res$wnn_p
res.gr$pseudobulk.pval <- res$pseudo_p

methods <- c(ivw.pval = "ivw", egger.pval = "egger", smr.pval = "smr", cct.pval = "cct", multi.pval = "wnn", pseudobulk.pval = "pseudo")

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 打标签 + 逐 gene AUC

# 同一个 gene 内与该 gene 的金标准区域重叠 -> TRUE
# (旧代码逐 gene 循环 findOverlaps; 这里一次 findOverlaps 再要求 gene 相同, 结果相同、更快)
label_by_gene_overlap <- function(res.gr, truth.gr) {
  hits <- findOverlaps(res.gr, truth.gr, ignore.strand = TRUE)
  same <- res.gr$gene[queryHits(hits)] == truth.gr$gene[subjectHits(hits)]
  res.gr$truth <- seq_along(res.gr) %in% queryHits(hits)[same]
  res.gr[res.gr$gene %in% unique(truth.gr$gene)]                   # 只保留金标准中有记录的 gene (同旧代码)
}

gene_wise_auc <- function(res.gr, nc = 1L) {
  idx_by_gene <- split(seq_along(res.gr), res.gr$gene)
  out <- mclapply(names(idx_by_gene), function(g) {
    idx <- idx_by_gene[[g]]
    y   <- res.gr$truth[idx]
    if (sum(y) == 0 || sum(!y) == 0) return(NULL)
    rbindlist(lapply(names(methods), function(m) {
      pval <- mcols(res.gr)[[m]][idx]
      keep <- !is.na(pval) & is.finite(pval) & pval > 0
      y2 <- y[keep]
      if (length(y2) < 2 || length(unique(y2)) < 2) return(NULL)
      roc.obj <- roc(response = y2, predictor = -log10(pval[keep]), quiet = TRUE)
      data.table(gene = g, method = methods[[m]], auc = as.numeric(auc(roc.obj)),
                 n_peaks = length(y2), n_pos = sum(y2), n_neg = sum(!y2))
    }))
  }, mc.cores = nc)
  bad <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("mclapply 出错: ", paste(unique(unlist(out[bad])), collapse = " | "))
  rbindlist(out)
}

auc_all <- list()
for (gold.name in names(gold.list)) {
  t0 <- Sys.time()
  bench <- label_by_gene_overlap(res.gr, gold.list[[gold.name]])
  ga <- gene_wise_auc(bench, nc = NC)
  if (nrow(ga)) ga <- cbind(data.table(cell_type = CT, gold_standard = gold.name), ga)
  auc_all[[gold.name]] <- ga
  cat(sprintf("[%s] genes in gold = %s | genes with AUC = %s | %.1f min\n", gold.name,
              format(uniqueN(bench$gene), big.mark = ","), format(uniqueN(ga$gene), big.mark = ","),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

auc_long <- rbindlist(auc_all)
f_out <- file.path(OUT_DIR, paste0(CT, "_AUC_long.tsv.gz"))
fwrite(auc_long, f_out, sep = "\t")
cat(sprintf("\n[done] %s 行 -> %s\n", format(nrow(auc_long), big.mark = ","), f_out))
print(dcast(auc_long[, .(med = round(median(auc), 3)), by = .(gold_standard, method)], gold_standard ~ method, value.var = "med"))
