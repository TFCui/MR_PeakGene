# QC_OneSampleMR 第 1 步: 对某个 CT 的全部 caQTL variant 做 LD pruning (流程见 Scripts/Pipeline.md 第 4.1 节)
# 用法: Rscript QC_OneSampleMR_prune.R <CT>
# 输入: caQTL 全量 summary (只取 unique variant_id + A1/A2); hg38 LD 面板 (ID = CHR:BP:REF:ALT)
# 输出: tmp/QC_OneSampleMR/<CT>_independent_variants.txt (pruned 后的 variant_id, 一行一个)
#       tmp/QC_OneSampleMR/<CT>_variant_map.tsv.gz (variant_id <-> 面板 ID)

suppressPackageStartupMessages(library(data.table))
options(width = 200)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

CT   <- commandArgs(TRUE)[1]
BASE <- "/scratch/user/tfcui2025/CIMA"
PROJ <- file.path(BASE, "Claude_project")
OUT  <- file.path(PROJ, "tmp/QC_OneSampleMR"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PLINK <- file.path(BASE, "R_package/my_plink/plink")
REF   <- file.path(PROJ, "data/ref_chr/hg38")
CAQTL <- file.path(BASE, "Data/OMIX009786/cis_caQTL_Summary_Statics/cis_caQTL", paste0(CT, ".txt.gz"))
NTHR  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
stopifnot(file.exists(CAQTL))
setDTthreads(NTHR)

# 1. unique variant (流式去重, 避免把全部行读进内存)
msg(CT, ": 读取 unique variant")
v <- fread(cmd = sprintf("zcat %s | cut -f2,11,12 | awk -F'\\t' 'NR>1 && !s[$1]++'", CAQTL), header = FALSE,
           col.names = c("variant_id", "A1", "A2"), sep = "\t")
msg("unique variant: ", nrow(v))
v[, c("chr", "pos") := tstrsplit(variant_id, "_", fixed = TRUE)[1:2]]
v[, chr := sub("^chr", "", chr)]
v <- v[chr %in% as.character(1:22)]
v[, id_a := paste0(chr, ":", pos, ":", A2, ":", A1)]   # REF=A2, ALT=A1
v[, id_b := paste0(chr, ":", pos, ":", A1, ":", A2)]   # 等位顺序相反

# 2. 匹配面板 ID
bim <- rbindlist(lapply(1:22, function(i) fread(file.path(REF, paste0("chr", i, ".bim")), select = 2, col.names = "id")))
panel <- unique(bim$id)
v[, panel_id := fifelse(id_a %chin% panel, id_a, fifelse(id_b %chin% panel, id_b, NA_character_))]
msg(sprintf("匹配到面板: %d / %d (%.1f%%)", sum(!is.na(v$panel_id)), nrow(v), 100 * mean(!is.na(v$panel_id))))
v <- v[!is.na(panel_id)]
v <- v[!duplicated(panel_id)]
fwrite(v[, .(variant_id, panel_id, chr)], file.path(OUT, paste0(CT, "_variant_map.tsv.gz")), sep = "\t")

# 3. 逐染色体 PLINK pruning
keep <- character(0)
for (i in 1:22) {
  ext <- file.path(OUT, sprintf("%s_chr%d_extract.txt", CT, i))
  writeLines(v[chr == as.character(i), panel_id], ext)
  pre <- file.path(OUT, sprintf("%s_prune_chr%d", CT, i))
  msg("chr", i, " | variants = ", sum(v$chr == as.character(i)))
  rc <- system2(PLINK, c("--bfile", file.path(REF, paste0("chr", i)), "--extract", ext, "--maf", "0.01",
                         "--indep-pairwise", "500kb", "50", "0.1", "--threads", NTHR, "--out", pre),
                stdout = FALSE, stderr = FALSE)
  stopifnot(rc == 0, file.exists(paste0(pre, ".prune.in")))
  keep <- c(keep, readLines(paste0(pre, ".prune.in")))
  file.remove(ext)
}
res <- v[panel_id %chin% keep, variant_id]
writeLines(res, file.path(OUT, paste0(CT, "_independent_variants.txt")))
msg(sprintf("%s: pruned 后保留 %d 个 variant (%d 个进入 pruning)", CT, length(res), nrow(v)))
