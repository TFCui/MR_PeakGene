# 2026-10-09
# Valid3_Motif 步骤 1: 各方法显著 peak 的 motif 富集 (流程见 Scripts/Pipeline.md 3.3)
# 移植自 Scripts_reference/run_motif_enrichment.R; 改动: 读新表 (不再做 sigma2 过滤), 不画 MotifPlot
#   前景: 每个方法 *_fdr < 0.1 的 peak; universe: All_GenePeak peak ∩ ATAC 对象
#   背景: MatchRegionStats 按 GC.percent + log1p(全部细胞 peak count) 从 universe 匹配 5 万个, 种子 2026 + 方法序号
#   Motif: JASPAR2020 CORE vertebrates
# 用法: Rscript Valid3_Motif_1_enrich.R <CT>   (array 0-2)
# 输出: Results/Validation/Valid3_Motif/<CT>_motif_enrichment.tsv.gz

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
  library(BiocParallel)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript Valid3_Motif_1_enrich.R <celltype>")
celltype <- args[1]

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
FDR_CUT <- 0.1
N_BG    <- 50000
SEED    <- 2026
MATCH_FEATURES <- c("GC.percent", "log_count")

methods <- c(IVW = "ivw_fdr", MR_Egger = "egger_fdr", SMR = "smr_fdr", CCT = "cct_fdr", WNN = "wnn_fdr", Pseudo = "pseudo_fdr")

BASE_DIR  <- "/scratch/user/tfcui2025/CIMA"
PROJ      <- file.path(BASE_DIR, "Claude_project")
mr_file   <- file.path(PROJ, "Results/MR", paste0(celltype, "_AllPairs_results.tsv.gz"))
link_file <- file.path(BASE_DIR, "Data/OMIX009786/gene_peak", paste0("All_GenePeak_", celltype, ".rds"))
atac_file <- file.path(BASE_DIR, "Data/OMIX009786/cell_rna/QC_filtered", paste0("atac_", celltype, "_qc.rds"))
out_dir   <- file.path(PROJ, "Results/Validation/Valid3_Motif")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_tsv   <- file.path(out_dir, paste0(celltype, "_motif_enrichment.tsv.gz"))

msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
msg("Cell type: ", celltype, " | cores: ", n_cores)
for (f in c(mr_file, link_file, atac_file)) if (!file.exists(f)) stop("File not found: ", f)

# ---------------------------- 1. MR results ----------------------------------
mr <- fread(mr_file, select = c("gene_id", "peak_id", unname(methods)))
msg("MR pairs: ", nrow(mr), " | unique peaks: ", uniqueN(mr$peak_id))
mr[, peak_id := sub(":", "-", peak_id, fixed = TRUE)]             # chr:start-end -> chr-start-end

fg_list <- lapply(methods, function(col) unique(mr[get(col) < FDR_CUT, peak_id]))
msg("Significant peaks (FDR < ", FDR_CUT, ") per method:")
print(sapply(fg_list, length))

# ---------------------------- 2. Background universe --------------------------
universe <- unique(readRDS(link_file)$peak_id)
msg("AllLinks unique peaks: ", length(universe))

# ---------------------------- 3. snATAC object + motifs -----------------------
msg("Loading snATAC object ...")
atac <- readRDS(atac_file)
assay_name <- DefaultAssay(atac)
peaks_obj  <- rownames(atac)
msg("Assay: ", assay_name, " | peaks in object: ", length(peaks_obj), " | cells: ", ncol(atac))

# accessibility per peak, computed across ALL cells (taken before AddMotifs)
mf0 <- GetAssayData(atac, assay = assay_name, layer = "meta.features")
if ("count" %in% colnames(mf0) && !anyNA(mf0[peaks_obj, "count"])) {
  peak_count <- setNames(mf0[peaks_obj, "count"], peaks_obj)
  msg("Accessibility: using 'count' column from meta.features")
} else {
  msg("Accessibility: 'count' not found in meta.features; computing rowSums of counts")
  peak_count <- Matrix::rowSums(GetAssayData(atac, assay = assay_name, layer = "counts"))[peaks_obj]
}
rm(mf0); invisible(gc())

n_univ_out <- sum(!universe %in% peaks_obj)
if (n_univ_out > 0) msg("WARNING: ", n_univ_out, " AllLinks peaks not in snATAC object; removed from universe")
universe <- intersect(universe, peaks_obj)
for (m in names(fg_list)) {
  n_out <- sum(!fg_list[[m]] %in% universe)
  if (n_out > 0) msg("WARNING: ", m, ": ", n_out, " significant peaks not in universe; removed")
  fg_list[[m]] <- intersect(fg_list[[m]], universe)
}

pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
msg("JASPAR2020 motifs: ", length(pfm))

register(MulticoreParam(workers = n_cores))
msg("Running AddMotifs ...")
atac <- AddMotifs(object = atac, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm, assay = assay_name)
register(SerialParam())                                           # avoid nested forking below
msg("AddMotifs done")

meta.feature <- GetAssayData(atac, assay = assay_name, layer = "meta.features")
if (!"GC.percent" %in% colnames(meta.feature)) stop("GC.percent not found in meta.features")
meta.feature$count_all <- as.numeric(peak_count[rownames(meta.feature)])
meta.feature$log_count <- log1p(meta.feature$count_all)
if (anyNA(meta.feature[universe, MATCH_FEATURES])) stop("NA in matching features for universe peaks")
mf_universe <- meta.feature[universe, , drop = FALSE]

# ---------------------------- 4. Enrichment per method (parallel) ------------
run_one <- function(m) {
  fg <- fg_list[[m]]
  set.seed(SEED + match(m, names(methods)))
  bg <- MatchRegionStats(meta.feature = mf_universe, query.feature = meta.feature[fg, , drop = FALSE],
                         features.match = MATCH_FEATURES, n = N_BG, verbose = FALSE)
  res <- as.data.table(FindMotifs(object = atac, features = fg, background = bg, assay = assay_name, verbose = FALSE))
  res[, p.adjust_BH := p.adjust(pvalue, method = "BH")]
  setorder(res, pvalue, -fold.enrichment)
  res[, rank := .I]
  res[, `:=`(celltype = celltype, method = m, n_fg = length(fg), n_bg = length(bg),
             GC_median_fg = median(meta.feature[fg, "GC.percent"]), GC_median_bg = median(meta.feature[bg, "GC.percent"]),
             logcount_median_fg = median(meta.feature[fg, "log_count"]), logcount_median_bg = median(meta.feature[bg, "log_count"]))]
  res
}

msg("Running FindMotifs for ", length(methods), " methods on ", min(n_cores, length(methods)), " cores ...")
res_list <- mclapply(names(methods), function(m) tryCatch(run_one(m), error = function(e) {
  message("ERROR in method ", m, ": ", conditionMessage(e)); NULL
}), mc.cores = min(n_cores, length(methods)))
names(res_list) <- names(methods)
failed <- names(res_list)[vapply(res_list, function(x) is.null(x) || inherits(x, "try-error"), logical(1))]
if (length(failed)) msg("WARNING: failed methods: ", paste(failed, collapse = ", "))

res_all <- rbindlist(res_list[setdiff(names(res_list), failed)], fill = TRUE)
setcolorder(res_all, c("celltype", "method", "rank", "motif", "motif.name", "observed", "background",
                       "percent.observed", "percent.background", "fold.enrichment", "pvalue", "p.adjust_BH"))
res_all[, method := factor(method, levels = names(methods))]
setorder(res_all, method, rank)
fwrite(res_all, out_tsv, sep = "\t")
msg("Saved: ", out_tsv, " (", nrow(res_all), " rows)")

print(res_all[, .(n_fg = n_fg[1], n_bg = n_bg[1], GC_fg = round(GC_median_fg[1], 2), GC_bg = round(GC_median_bg[1], 2),
                  logcnt_fg = round(logcount_median_fg[1], 2), logcnt_bg = round(logcount_median_bg[1], 2),
                  n_sig = sum(p.adjust < 0.05 & fold.enrichment > 1), top5 = paste(head(motif.name, 5), collapse = ",")), by = method])
msg("Done.")
if (length(failed)) quit(status = 1)
