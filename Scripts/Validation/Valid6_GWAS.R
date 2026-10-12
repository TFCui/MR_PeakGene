# 2026-10-09
# Valid6_GWAS: 只被 MR 支持的强 peak-gene link, 其 peak 内是否有 GWAS 全基因组显著 SNP (流程见 Scripts/Pipeline.md 3.6)
# 移植自 20260525_QC1/Scripts/source/G1_find_gwas.R (旧输出 Valid8_matchGWAS/G1_Focus_MRonly_GWAS.tsv); 逻辑和阈值不变, 只改:
#   输入 = Results/MR/<CT>_AllPairs_results.tsv.gz (已过滤, 不再做 sigma2 过滤); GWAS 只用项目的 9 个 trait (旧版 13 个)
# 筛选 (三个都要满足):
#   1. MR 强: IVW / Egger / SMR / CCT 任一 FDR < 0.001
#   2. WNN 与 Pseudo FDR 都 >= 0.05 (MR_only; 否则记为 Both)
#   3. peak 内有 GWAS SNP P < 5e-8
# GWAS: selected_sumstats (munged, SNP A1 A2 Z N), P = 2*pnorm(-|Z|); 缓存按 P < 1e-5 生成
#       rsID -> hg19 (1000G EUR bim) -> liftOver hg38 (1-to-1, chr1-22) -> EAS_hg38.bim 取 CHR:BP:REF:ALT
# 重合: SNP hg38 位置落在 peak [start+1, end] 内 (EXT = 0)
# 输出 Results/Validation/Valid6_GWAS/:
#   GWAS_MRonly_overlap.tsv  每行 = MR_only pair x 重合 SNP x trait (列与旧 G1_Focus_MRonly_GWAS.tsv 相同); Fig8 用
#   GWAS_pair_counts.tsv     每个 CT x category: MR 强的 pair 数 / 其中 peak 含 GWAS SNP 的 pair 数
# 缓存: tmp/Valid6_GWAS/GWAS_SNPs_P1e-5_hg38.tsv.gz (存在即跳过; 改了 trait 或 P_CACHE 要删)

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})
options(width = 300)

## ---- 0. 参数和路径 ----
BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
OUT_DIR  <- file.path(PROJ, "Results/Validation/Valid6_GWAS")
TMP_DIR  <- file.path(PROJ, "tmp/Valid6_GWAS")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_DIR   <- file.path(BASE_DIR, "R_package/S_LDSC/GWAS_data/selected_sumstats")
BIM19      <- file.path(BASE_DIR, "R_package/S_LDSC/reference_genome/1000G_EUR_Phase3_plink/1000G.EUR.QC.")
BIM38      <- file.path(BASE_DIR, "R_package/my_plink/plink_hg38/EAS_hg38.bim")
CHAIN      <- file.path(BASE_DIR, "Data/PBMC_validation/hg19ToHg38.over.chain")
GWAS_CACHE <- file.path(TMP_DIR, "GWAS_SNPs_P1e-5_hg38.tsv.gz")
OUT_FILE   <- file.path(OUT_DIR, "GWAS_MRonly_overlap.tsv")
CNT_FILE   <- file.path(OUT_DIR, "GWAS_pair_counts.tsv")
OLD_FILE   <- file.path(BASE_DIR, "Results/RunMR/Materials_runMR/20260430_LinkPeakGene_more/MRresult_20260911/Valid8_matchGWAS/G1_Focus_MRonly_GWAS.tsv")

CTS        <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
MR_M       <- c("ivw", "egger", "smr", "cct")
ASSOC_M    <- c("wnn", "pseudo")
METHODS    <- c(MR_M, ASSOC_M)
MR_LAB     <- c(ivw = "IVW", egger = "Egger", smr = "SMR", cct = "CCT")
MR_FDR     <- 0.001                 # MR 强: MR_M 中任一方法 FDR < MR_FDR
ASSOC_FDR  <- 0.05                  # association 不显著: ASSOC_M 两个方法 FDR 都 >= ASSOC_FDR
P_GWAS     <- 5e-8                  # peak 内 GWAS SNP 阈值
P_CACHE    <- 1e-5                  # 缓存按 P < P_CACHE 生成
EXT        <- 0L

# 项目的 9 个 GWAS (标签与 Valid2_SLDSC 一致)
TRAITS <- c(AS     = "CUSTOM.AS",
            CD     = "CUSTOM.CD",
            IBD    = "PASS.Inflammatory_Bowel_Disease.deLange2017",
            RA     = "PASS.Rheumatoid_Arthritis.Saevarsdottir2022",
            SLE    = "PASS.Lupus.Bentham2015",
            MS     = "PASS_Multiple_sclerosis",
            T1D    = "PASS.Type_1_Diabetes.Chiou2021",
            Asthma = "GBMI.Asthma",
            PBC    = "PASS.Primary_Biliary_Cirrhosis.Cordell2015")

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")


####################################################################################################
#### Step 1: MR 强的 pair + 分类
####################################################################################################

keep_cols <- c("gene_id", "peak_id", "smr_variant_id",
               paste0(METHODS, "_p"), paste0(METHODS, "_fdr"),
               "ivw_beta", "egger_beta", "smr_b", "wnn_rho", "pseudo_rho")

pairs <- rbindlist(lapply(CTS, function(ct) {
  dt <- fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")), select = keep_cols)
  dt[, cell_type := ct]
}))
setcolorder(pairs, "cell_type")
stopifnot(!anyDuplicated(pairs, by = c("cell_type", "gene_id", "peak_id")))

for (m in METHODS) {
  f   <- pairs[[paste0(m, "_fdr")]]
  cut <- if (m %in% MR_M) MR_FDR else ASSOC_FDR
  set(pairs, j = paste0("sig_", m), value = !is.na(f) & f < cut)
}
pairs[, mr_any    := Reduce(`|`, .SD), .SDcols = paste0("sig_", MR_M)]
pairs[, assoc_any := Reduce(`|`, .SD), .SDcols = paste0("sig_", ASSOC_M)]
pairs <- pairs[mr_any == TRUE]
pairs[, category := fifelse(assoc_any, "Both", "MR_only")]
pairs[, mr_methods_sig := apply(.SD, 1, function(x) paste(MR_LAB[MR_M[x]], collapse = ",")), .SDcols = paste0("sig_", MR_M)]
pairs[, min_mr_fdr := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = paste0(MR_M, "_fdr")]

msg("MR 强的 pairs (", paste(MR_LAB, collapse = "/"), " 任一 FDR < ", MR_FDR, "); Both = WNN 或 Pseudo FDR < ", ASSOC_FDR, ":")
print(dcast(pairs[, .N, by = .(cell_type, category)], cell_type ~ category, value.var = "N", fill = 0))


####################################################################################################
#### Step 2: GWAS (P < P_CACHE) -> hg19 -> hg38 (缓存)
####################################################################################################

if (file.exists(GWAS_CACHE)) {
  msg("读取缓存: ", GWAS_CACHE)
  gwas <- fread(GWAS_CACHE)
  stopifnot(setequal(unique(gwas$trait), names(TRAITS)))
} else {
  ## ---- 2.1 9 个 trait, 只留 P < P_CACHE ----
  gwas <- rbindlist(lapply(names(TRAITS), function(tr) {
    g <- fread(file.path(GWAS_DIR, paste0(TRAITS[[tr]], ".sumstats.gz")), select = c("SNP", "A1", "A2", "Z"))
    g <- g[!is.na(Z)]
    g[, P := 2 * pnorm(-abs(Z))]
    msg(sprintf("  %-7s SNPs: %8d | P<1e-5: %6d | P<5e-8: %5d", tr, nrow(g), sum(g$P < P_CACHE), sum(g$P < P_GWAS)))
    g[P < P_CACHE][, trait := tr]
  }))
  setnames(gwas, "SNP", "rsid")

  ## ---- 2.2 rsID -> hg19 (1000G EUR bim) ----
  bim19 <- rbindlist(lapply(1:22, function(chr)
    fread(paste0(BIM19, chr, ".bim"), header = FALSE, select = c(1, 2, 4), col.names = c("chr", "rsid", "pos_hg19"))))
  bim19 <- bim19[rsid %chin% gwas$rsid]
  bim19 <- bim19[!duplicated(rsid)]
  snp <- data.table(rsid = unique(gwas$rsid))
  snp <- bim19[snp, on = "rsid"]
  msg(sprintf("rsID -> hg19: %d / %d (%.2f%%)", sum(!is.na(snp$pos_hg19)), nrow(snp), 100 * mean(!is.na(snp$pos_hg19))))
  snp <- snp[!is.na(pos_hg19)]

  ## ---- 2.3 liftOver hg19 -> hg38 (只保留 1-to-1 + chr1-22) ----
  chain <- import.chain(CHAIN)
  gr19  <- GRanges(paste0("chr", snp$chr), IRanges(snp$pos_hg19, width = 1), rsid = snp$rsid)
  gl    <- liftOver(gr19, chain)
  n_hit <- lengths(gl)
  gr38  <- unlist(gl[n_hit == 1])
  lift  <- data.table(rsid = gr38$rsid, chr38 = as.character(seqnames(gr38)), pos_hg38 = start(gr38))
  lift  <- lift[chr38 %in% paste0("chr", 1:22)]
  msg(sprintf("liftOver 1-to-1 (chr1-22): %d / %d (%.2f%%)", nrow(lift), nrow(snp), 100 * nrow(lift) / nrow(snp)))
  snp <- snp[lift, on = "rsid"]
  snp[, chr := chr38][, chr38 := NULL]

  ## ---- 2.4 hg38 variant ID (EAS_hg38.bim: CHR:BP:REF:ALT, 等位基因两种方向都试) ----
  gwas <- snp[gwas, on = "rsid", nomatch = 0L]
  bim38 <- fread(BIM38, header = FALSE, select = c(1, 2, 4, 5, 6), col.names = c("chr", "varid_hg38", "pos_hg38", "b1", "b2"))
  bim38[, chr := paste0("chr", chr)]
  bim38 <- bim38[unique(gwas[, .(chr, pos_hg38)]), on = .(chr, pos_hg38), nomatch = 0L]

  allele_key <- function(a, b) fifelse(a < b, paste(a, b, sep = "/"), paste(b, a, sep = "/"))
  bim38[, akey := allele_key(b1, b2)]
  gwas[,  akey := allele_key(A1, A2)]
  bim38 <- unique(bim38[, .(chr, pos_hg38, akey, varid_hg38)], by = c("chr", "pos_hg38", "akey"))
  gwas  <- bim38[gwas, on = .(chr, pos_hg38, akey)]
  gwas[, in_EAS_ref := !is.na(varid_hg38)]
  gwas[is.na(varid_hg38), varid_hg38 := paste0(sub("^chr", "", chr), ":", pos_hg38)]   # 与 EAS_hg38.bim 一致, 不带 chr
  gwas[, akey := NULL]
  msg(sprintf("在 EAS_hg38.bim 中匹配到 (位置 + 等位基因): %.2f%%", 100 * mean(gwas$in_EAS_ref)))

  setcolorder(gwas, c("trait", "rsid", "varid_hg38", "chr", "pos_hg38", "pos_hg19", "A1", "A2", "Z", "P", "in_EAS_ref"))
  setorder(gwas, trait, chr, pos_hg38)
  fwrite(gwas, GWAS_CACHE, sep = "\t")
  msg("保存缓存: ", GWAS_CACHE)
}

## ---- 抽查: rs2476601 (PTPN22) hg38 应为 chr1:113834946 ----
chk <- unique(gwas[rsid == "rs2476601", .(rsid, chr, pos_hg19, pos_hg38, varid_hg38)])
if (nrow(chk)) { msg("抽查 rs2476601 (期望 chr1:113834946):"); print(chk) } else msg("抽查: rs2476601 不在 P<1e-5 集合中")

## ---- 只留 P < P_GWAS; 列名改为 gwas_p / gwas_z, 免得和 MR 的 p 混淆 ----
gwas <- gwas[P < P_GWAS]
setnames(gwas, c("P", "Z"), c("gwas_p", "gwas_z"))
msg("GWAS SNP (P < ", P_GWAS, ") per trait:")
print(gwas[, .N, by = trait])


####################################################################################################
#### Step 3: peak x SNP 重合
####################################################################################################

peaks <- unique(pairs$peak_id)
sp    <- tstrsplit(peaks, "[:-]", type.convert = TRUE)
pk_gr <- GRanges(sp[[1]], IRanges(sp[[2]] + 1L - EXT, sp[[3]] + EXT))

snp_u  <- unique(gwas[, .(varid_hg38, chr, pos_hg38)])
snp_gr <- GRanges(snp_u$chr, IRanges(snp_u$pos_hg38, width = 1))

hits <- findOverlaps(pk_gr, snp_gr)
ov   <- data.table(peak_id = peaks[queryHits(hits)], varid_hg38 = snp_u$varid_hg38[subjectHits(hits)])
ov   <- gwas[ov, on = "varid_hg38", allow.cartesian = TRUE]               # 加上 trait 层信息
msg(sprintf("重合: peaks %d / %d 含 P<5e-8 SNP; peak x SNP x trait 行数 %d", uniqueN(ov$peak_id), length(peaks), nrow(ov)))


####################################################################################################
#### Step 4: MR_only pairs x 重合 SNP
####################################################################################################

link_cols <- c("cell_type", "gene_id", "peak_id", "mr_methods_sig")
p_cols    <- paste0(METHODS, "_p")
gwas_cols <- c("trait", "rsid", "varid_hg38", "gwas_p")
ref_cols  <- c(paste0(METHODS, "_fdr"), "ivw_beta", "egger_beta", "smr_b", "wnn_rho", "pseudo_rho", "smr_variant_id",
               "chr", "pos_hg38", "pos_hg19", "A1", "A2", "gwas_z", "in_EAS_ref")

pair_cols <- c("cell_type", "gene_id", "peak_id", "category", "mr_methods_sig", "min_mr_fdr", "smr_variant_id",
               p_cols, paste0(METHODS, "_fdr"), "ivw_beta", "egger_beta", "smr_b", "wnn_rho", "pseudo_rho")
snp_cols  <- c("trait", "rsid", "varid_hg38", "chr", "pos_hg38", "pos_hg19", "A1", "A2", "gwas_z", "gwas_p", "in_EAS_ref")

all_ov <- pairs[, ..pair_cols][ov[, c("peak_id", snp_cols), with = FALSE], on = "peak_id", nomatch = 0L, allow.cartesian = TRUE]

# 计数: MR 强的 pair 数, 其中 peak 含 P<5e-8 SNP 的 pair 数
cnt <- merge(pairs[, .(n_pairs_MR = .N), by = .(cell_type, category)],
             all_ov[, .(n_pairs_GWAS = uniqueN(paste(gene_id, peak_id))), by = .(cell_type, category)],
             by = c("cell_type", "category"), all.x = TRUE)
cnt[is.na(n_pairs_GWAS), n_pairs_GWAS := 0L]
msg("每个 CT x category: MR 强的 pairs / 其中 peak 含 GWAS SNP (P < ", P_GWAS, ") 的 pairs:")
print(cnt)
fwrite(cnt, CNT_FILE, sep = "\t")

focus <- all_ov[category == "MR_only"]
setorder(focus, min_mr_fdr, cell_type, gene_id, peak_id, gwas_p)
focus <- focus[, c(link_cols, p_cols, gwas_cols, ref_cols), with = FALSE]
fwrite(focus, OUT_FILE, sep = "\t")

msg(sprintf("Focus (MR_only + GWAS P<5e-8): %d 行 | pairs %d | genes %d | peaks %d",
            nrow(focus), uniqueN(focus[, .(cell_type, gene_id, peak_id)]), uniqueN(focus$gene_id), uniqueN(focus$peak_id)))
print(focus[, .N, by = .(cell_type, gene_id, peak_id)][, .(n_pairs = .N), by = cell_type])
print(focus[, .(n_rows = .N, n_pairs = uniqueN(paste(cell_type, gene_id, peak_id))), by = trait])

# 每个 pair 一行: 6 个方法的 P 值 + 最强 GWAS SNP
pair_sum <- focus[, c(.SD[1, c(p_cols), with = FALSE],
                      .(lead_rsid = rsid[which.min(gwas_p)], lead_gwas_p = min(gwas_p),
                        traits = paste(unique(trait), collapse = ","))),
                  by = c(link_cols)]
num <- c(p_cols, "lead_gwas_p")
pair_sum[, (num) := lapply(.SD, signif, 2), .SDcols = num]
msg("每个 pair 的摘要 (6 个方法 P 值 + 最强 GWAS SNP):")
print(pair_sum)


####################################################################################################
#### Step 5: 与旧结果 (Valid8_matchGWAS, 13 traits, 旧 FDR) 对照, pair 层
####################################################################################################

old <- unique(fread(OLD_FILE)[, .(cell_type, gene_id, peak_id, old_traits = trait)])[
  , .(old_traits = paste(sort(unique(old_traits)), collapse = ",")), by = .(cell_type, gene_id, peak_id)]
new <- pair_sum[, .(cell_type, gene_id, peak_id, new_traits = traits)]
cmp <- merge(old, new, by = c("cell_type", "gene_id", "peak_id"), all = TRUE)
cmp[, status := fcase(!is.na(old_traits) & !is.na(new_traits), "both",
                      is.na(new_traits), "old_only", default = "new_only")]
msg(sprintf("与旧 Valid8 对照 (pair 层): 旧 %d | 新 %d | 共有 %d | 只在旧 %d | 只在新 %d",
            nrow(old), nrow(new), cmp[status == "both", .N], cmp[status == "old_only", .N], cmp[status == "new_only", .N]))
print(cmp[status != "both"][order(status, cell_type, gene_id)])

msg("Saved: ", OUT_FILE)
msg("Saved: ", CNT_FILE)
