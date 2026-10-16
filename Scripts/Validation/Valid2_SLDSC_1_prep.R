# 2026-10-09
# 2026-10-14: 集合改为 MROnly + pool (6 个方法各自的 FDR 0.1 结果暂存为 SLDSC_6methods_FDR0.1.tsv / set_stats_6methods_FDR0.1.tsv)
# Valid2_SLDSC 步骤 1: 准备 annotation (流程见 Scripts/Pipeline.md 3.2)
# 移植自 Scripts_reference/D1_prep.R; 改动: 读新表的 *_fdr (不再过滤 / 重新 BH), 只做两侧各扩 500 bp, 文件名去掉 0911 / ext 后缀
# 每个 CT 2 个集合 -> hg19 -> 两侧各扩 500 bp -> reduce -> thin-annot
#   MROnly: IVW / Egger / SMR 任一 FDR < 0.1, 且 WNN、Pseudo 的 FDR 都 >= 0.1 的 pair 的 peak (link 层)
#   pool  : All_GenePeak 的全部 peak (与 FDR 无关)
# 2026-10-14: 另加 Combined (第 4 个 "CT") = 3 个 CT 的 MROnly 并集 / pool 并集 (hg19 id 取并集, 扩展后 reduce 合并跨 CT 的重叠)
# 一个 job 处理 3 个 CT + Combined; 22 条染色体 mclapply 并行
# 输出: tmp/Valid2_SLDSC/annot_ld/<CT>/<CT>_<set>.<chr>.annot.gz
#       Results/Validation/Valid2_SLDSC/set_stats.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(parallel)
})
options(width = 300)

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
OUT_DIR  <- file.path(PROJ, "Results/Validation/Valid2_SLDSC")
ANNOT    <- file.path(PROJ, "tmp/Valid2_SLDSC/annot_ld")
POOL_DIR <- file.path(BASE_DIR, "Data/OMIX009786/gene_peak")
CHAIN    <- file.path(BASE_DIR, "Data/PBMC_validation/hg38ToHg19.over.chain")
PLINK    <- file.path(BASE_DIR, "R_package/S_LDSC/reference_genome/1000G_EUR_Phase3_plink/1000G.EUR.QC.")

CELLS   <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
CELLS_ALL <- c(CELLS, "Combined")
MR_M    <- c("ivw", "egger", "smr")
MG_M    <- c("wnn", "pseudo")
SETS    <- c("MROnly", "pool")
FDR_CUT <- 0.1
EXT     <- 500L
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "22"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (ct in CELLS_ALL) dir.create(file.path(ANNOT, ct), recursive = TRUE, showWarnings = FALSE)
chain <- import.chain(CHAIN)

# pool 格式 "chr1-9928-10429" -> "chr1:9928-10429"
norm_id <- function(x) { sp <- tstrsplit(x, "[:-]"); paste0(sp[[1]], ":", sp[[2]], "-", sp[[3]]) }

# peak id (0-based) -> 两侧各扩 ext -> 合并重叠 (GRanges, 1-based)
id2gr <- function(ids, ext = 0L) {
  sp <- tstrsplit(ids, "[:-]", type.convert = TRUE)
  reduce(GRanges(sp[[1]], IRanges(sp[[2]] + 1L - ext, sp[[3]] + ext)))
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Step 1: 每个 CT: MROnly + pool -> liftOver hg19

peaksets  <- list()
set_stats <- list()
for (CT in CELLS) {
  pairs <- fread(file.path(MR_DIR, paste0(CT, "_AllPairs_results.tsv.gz")), select = c("gene_id", "peak_id", paste0(c(MR_M, MG_M), "_fdr")))
  stopifnot(!anyNA(pairs))
  pairs[, mr := Reduce(`|`, lapply(paste0(MR_M, "_fdr"), function(k) get(k) < FDR_CUT))]
  pairs[, mg := Reduce(`|`, lapply(paste0(MG_M, "_fdr"), function(k) get(k) < FDR_CUT))]
  pool_peaks38 <- unique(norm_id(as.data.table(readRDS(file.path(POOL_DIR, paste0("All_GenePeak_", CT, ".rds"))))$peak_id))

  # liftOver hg38 -> hg19 (只保留 1-to-1 + chr1-22)
  all_peaks <- unique(c(pool_peaks38, pairs$peak_id))
  sp    <- tstrsplit(all_peaks, "[:-]", type.convert = TRUE)
  gl    <- liftOver(GRanges(sp[[1]], IRanges(sp[[2]] + 1L, sp[[3]]), peak_id = all_peaks), chain)
  n_hit <- lengths(gl)
  gr19  <- unlist(gl[n_hit == 1])
  peak_map <- data.table(peak_id = gr19$peak_id, chr19 = as.character(seqnames(gr19)), start19 = start(gr19) - 1L, end19 = end(gr19))
  peak_map <- peak_map[chr19 %in% paste0("chr", 1:22)]
  peak_map[, peak_id_hg19 := paste0(chr19, ":", start19, "-", end19)]

  pairs19 <- merge(pairs, peak_map[, .(peak_id, peak_id_hg19)], by = "peak_id")
  pool19  <- unique(peak_map[peak_id %in% pool_peaks38, peak_id_hg19])
  cat(sprintf("[%s] liftOver 1-to-1: %d / %d | MR pairs: %d -> %d | pool peaks (hg19): %d\n",
              CT, sum(n_hit == 1), length(n_hit), nrow(pairs), nrow(pairs19), length(pool19)))

  peaksets[[CT]][["MROnly"]] <- unique(pairs19[mr & !mg, peak_id_hg19])
  peaksets[[CT]][["pool"]] <- pool19

  set_stats[[CT]] <- data.table(cell_type = CT, set = SETS,
                                n_peak  = sapply(SETS, function(a) length(peaksets[[CT]][[a]])),
                                in_pool = sapply(SETS, function(a) mean(peaksets[[CT]][[a]] %in% pool19)))
}

# Combined: 3 个 CT 的并集
for (s in SETS) peaksets[["Combined"]][[s]] <- unique(unlist(lapply(CELLS, function(ct) peaksets[[ct]][[s]])))
set_stats[["Combined"]] <- data.table(cell_type = "Combined", set = SETS,
                                      n_peak  = sapply(SETS, function(a) length(peaksets[["Combined"]][[a]])),
                                      in_pool = sapply(SETS, function(a) mean(peaksets[["Combined"]][[a]] %in% peaksets[["Combined"]][["pool"]])))
cat(sprintf("[Combined] MROnly peaks (hg19): %d | pool peaks (hg19): %d\n", length(peaksets$Combined$MROnly), length(peaksets$Combined$pool)))
set_stats <- rbindlist(set_stats)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Step 2: 写 thin-annot (22 条染色体并行)

jobs <- CJ(cell_type = CELLS_ALL, set = SETS, sorted = FALSE)
jobs[, name := paste0(cell_type, "_", set)]
cat("annotation sets:", nrow(jobs), "(应为 4 x 2 = 8)\n")
gr_list <- lapply(seq_len(nrow(jobs)), function(j) id2gr(peaksets[[jobs$cell_type[j]]][[jobs$set[j]]], EXT))

write_chr <- function(chr) {
  bim <- fread(paste0(PLINK, chr, ".bim"), header = FALSE, select = c(1, 4), col.names = c("chr", "bp"))
  snp <- GRanges(paste0("chr", bim$chr), IRanges(bim$bp, width = 1))
  rbindlist(lapply(seq_len(nrow(jobs)), function(j) {
    hit <- overlapsAny(snp, gr_list[[j]])
    fwrite(data.table(ANNOT = as.integer(hit)), file.path(ANNOT, jobs$cell_type[j], sprintf("%s.%d.annot.gz", jobs$name[j], chr)))
    data.table(name = jobs$name[j], chr = chr, n_snp = length(snp), n_in = sum(hit))
  }))
}
snp_stats <- mclapply(1:22, write_chr, mc.cores = min(N_CORES, 22))
bad <- vapply(snp_stats, inherits, logical(1), "try-error")
if (any(bad)) stop("写 annot 出错: ", paste(unique(unlist(snp_stats[bad])), collapse = " | "))
snp_stats <- rbindlist(snp_stats)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Step 3: 诊断 —— peak 数 / bp / SNP 覆盖

snp_stats <- snp_stats[, .(n_snp = sum(n_snp), n_in = sum(n_in)), by = name]
jobs[, bp := sapply(gr_list, function(g) sum(as.numeric(width(g))))]
stats_all <- merge(merge(jobs, snp_stats, by = "name"), set_stats, by = c("cell_type", "set"))
stats_all[, prop_snp := n_in / n_snp]
stats_all[, `:=`(cell_type = factor(cell_type, levels = CELLS_ALL), set = factor(set, levels = SETS))]
setorder(stats_all, cell_type, set)
stats_all[, name := NULL]
fwrite(stats_all, file.path(OUT_DIR, "set_stats.tsv"), sep = "\t")
print(stats_all, digits = 3)

n_files <- length(list.files(ANNOT, pattern = "annot.gz$", recursive = TRUE))
cat("annot files:", n_files, "(应为 8 x 22 = 176)\n")
stopifnot(n_files == 176L)
