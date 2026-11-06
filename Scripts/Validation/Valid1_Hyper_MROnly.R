# 2026-11-03
# Valid1_Hyper_MROnly: MR-Only link 的 gene 内金标准命中 (3 个细胞类型合并; Plot13 / 论文当前 S5 Fig 用)
#   MR-Only link (每个 CT 内): IVW / Egger / SMR 任一 *_fdr < 0.1, 且 wnn_fdr、pseudo_fdr 都 >= 0.1 (与 Valid2_SLDSC 相同; CCT 不参与)
#   3 个 CT 合并, 每个 gene 一行:
#     背景 N   = 3 个 CT All_GenePeak 中该 gene 候选 peak 的并集
#     n_sig    = 3 个 CT 中该 gene MR-Only peak 的并集
#     金标准   = Valid1 保存的 <CT>_gold_standard.rds 按 gene 取 3 个 CT 的并集 (HiChIP / PCHiC = 3 个 CT loop 的并集; 其余 4 个库本来不分 CT)
#     K / k    = 背景 / MR-Only peak 中与该 gene 金标准区域重叠的数; pval = 超几何 (与 Valid1 同一公式)
#   与 Valid1 的 <CT>_hyper_long.tsv 不同: 那里按 CT 分开, HiChIP / PCHiC 用各 CT 自己的 loop
# 输入: Results/MR/<CT>_AllPairs_results.tsv.gz; All_GenePeak_<CT>.rds; Results/Validation/Valid1_Hyper/<CT>_gold_standard.rds
# 输出: Results/Validation/Valid1_Hyper/MROnly_gene_overlap_combined.tsv (gene_id, database, N, K, n_sig, k, pval, n_ct)
# 运行: sbatch Scripts/Validation/Valid1_Hyper_MROnly_slurm.sh
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})
options(width = 300)

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
OUT_DIR  <- file.path(PROJ, "Results/Validation/Valid1_Hyper")
GP_DIR   <- file.path(BASE_DIR, "Data/OMIX009786/gene_peak")
CELLS    <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
FDR_CUT  <- 0.1
DB_NAMES <- c("enhancerAtlas", "FANTOM5", "4DGenome", "GTEx", "HiChIP", "PCHiC")
fmt <- function(x) format(x, big.mark = ",")

phyper_p <- function(N, K, n_sig, k) if (K == 0L || n_sig == 0L || N == K) 1 else phyper(k - 1L, m = K, n = N - K, k = n_sig, lower.tail = FALSE)

# 金标准对象 (GRangesList 或 list of GRanges) -> data.table(gene, chr, start, end)
db_to_dt <- function(obj) {
  gr <- if (is(obj, "GRangesList")) unlist(obj) else unlist(GRangesList(lapply(obj, granges)))
  data.table(gene = names(gr), chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
}

###------------------------------------------------------------------------------------------------
### 1. 3 个 CT: MR-Only link / 背景 / 金标准
MRO <- list(); BG <- list(); DBT <- setNames(vector("list", length(DB_NAMES)), DB_NAMES)
for (ct in CELLS) {
  res <- fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")))
  mro <- res[(ivw_fdr < FDR_CUT | egger_fdr < FDR_CUT | smr_fdr < FDR_CUT) & wnn_fdr >= FDR_CUT & pseudo_fdr >= FDR_CUT,
             .(gene_id, peak_id, cell_type = ct)]
  cat(sprintf("[%s] MR pairs %s; MR-Only links %s (%s genes)\n", ct, fmt(nrow(res)), fmt(nrow(mro)), fmt(uniqueN(mro$gene_id))))
  MRO[[ct]] <- mro
  gp <- as.data.table(readRDS(file.path(GP_DIR, paste0("All_GenePeak_", ct, ".rds"))))[, .(gene_id, peak_id = sub("-", ":", peak_id))]
  BG[[ct]] <- unique(gp)
  DB <- readRDS(file.path(OUT_DIR, paste0(ct, "_gold_standard.rds")))
  stopifnot(identical(names(DB), DB_NAMES))
  for (d in DB_NAMES) DBT[[d]] <- rbind(DBT[[d]], db_to_dt(DB[[d]]))
  rm(res, gp, DB); invisible(gc())
}
MRO <- rbindlist(MRO)
n_ct <- MRO[, .(n_ct = uniqueN(cell_type)), by = gene_id]
mro_pairs <- unique(MRO[, .(gene_id, peak_id)])
cat(sprintf("\n[combined] MR-Only links (按 gene + peak 去重) %s, genes %s; 来自 1 / 2 / 3 个 CT 的 gene: %s\n",
            fmt(nrow(mro_pairs)), fmt(uniqueN(mro_pairs$gene_id)), paste(fmt(tabulate(n_ct$n_ct, 3)), collapse = " / ")))

genes <- sort(unique(mro_pairs$gene_id))
BG <- unique(rbindlist(BG)[gene_id %chin% genes])
bg_by_gene  <- split(BG$peak_id, BG$gene_id)
mro_by_gene <- split(mro_pairs$peak_id, mro_pairs$gene_id)
in_bg <- sum(vapply(genes, function(g) sum(mro_by_gene[[g]] %chin% bg_by_gene[[g]]), integer(1)))
cat(sprintf("[check] MR-Only peak 落在自身 gene 的背景 (3 CT 并集) 内: %s / %s (%.2f%%)\n",
            fmt(in_bg), fmt(nrow(mro_pairs)), 100 * in_bg / nrow(mro_pairs)))

# 金标准: 3 个 CT 取并集, 只保留 MR-Only gene
DBG <- lapply(DB_NAMES, function(d) {
  x <- unique(DBT[[d]][gene %chin% genes])
  split(GRanges(x$chr, IRanges(x$start, x$end)), x$gene)
})
names(DBG) <- DB_NAMES
rm(DBT); invisible(gc())
cat("\n[db] 合并后 (只算 MR-Only gene):\n")
for (d in DB_NAMES) cat(sprintf("  %-14s %5s genes, 中位区间数 %4.0f\n", d, fmt(length(DBG[[d]])), median(lengths(DBG[[d]]))))

###------------------------------------------------------------------------------------------------
### 2. 每个 gene × 库: N / K / n_sig / k
ALL_PEAK_IDS <- unique(BG$peak_id)
.mat <- do.call(rbind, strsplit(ALL_PEAK_IDS, "[:-]"))
ALL_PEAK_GR <- GRanges(.mat[, 1], IRanges(start = as.integer(.mat[, 2]), end = as.integer(.mat[, 3])))
rm(.mat); invisible(gc())
gr_for <- function(peak_ids) ALL_PEAK_GR[chmatch(peak_ids, ALL_PEAK_IDS)]

t0 <- Sys.time()
out <- vector("list", length(genes))
for (i in seq_along(genes)) {
  g      <- genes[i]
  bg_pk  <- unique(bg_by_gene[[g]])
  bg_gr  <- gr_for(bg_pk)
  is_sig <- bg_pk %chin% mro_by_gene[[g]]
  rows <- list()
  for (d in DB_NAMES) {
    if (!g %in% names(DBG[[d]])) next
    h <- suppressWarnings(countOverlaps(bg_gr, DBG[[d]][[g]])) > 0L
    N <- length(bg_pk); K <- sum(h); n_sig <- sum(is_sig); k <- sum(h & is_sig)
    rows[[length(rows) + 1L]] <- list(gene_id = g, database = d, N = N, K = K, n_sig = n_sig, k = k, pval = phyper_p(N, K, n_sig, k))
  }
  if (length(rows)) out[[i]] <- rbindlist(rows)
}
tab <- rbindlist(out)
tab <- merge(tab, n_ct, by = "gene_id", sort = FALSE)
tab[, database := factor(database, levels = DB_NAMES)]
setorder(tab, database, gene_id)
f_out <- file.path(OUT_DIR, "MROnly_gene_overlap_combined.tsv")
fwrite(tab, f_out, sep = "\t")
cat(sprintf("\n[done] %s 行 -> %s (%.1f min)\n", fmt(nrow(tab)), f_out, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(tab[, .(n_gene = .N, n_K0 = sum(K == 0), n_used = sum(K > 0 & N > K), sum_N = sum(N), sum_K = sum(K), sum_nsig = sum(n_sig), sum_k = sum(k)),
          keyby = database])
cat("\n全部完成\n")
