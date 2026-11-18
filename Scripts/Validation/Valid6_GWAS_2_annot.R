# 2026-11-09
# Valid6_GWAS_2_annot: Fig 4B 选例子用的候选表 (MR-only link 的 peak 内有 GWAS SNP) + 注释
#   候选池 (放宽版, 只有一个): IVW / Egger / SMR / CCT 任一 FDR < 0.05, 且 WNN、Pseudo FDR 都 >= 0.05, peak 内有 GWAS P < 5e-8 SNP
#     GWAS 只用 Fig 4A 的 7 个性状 (去掉 Asthma、PBC); SNP 取 Valid6_GWAS.R 的缓存, peak x SNP 重叠逻辑同 Valid6_GWAS.R (EXT = 0)
#     原 Valid6_GWAS (MR FDR < 0.001) 的 pair 是子集, 列 in_old_strict 标出
#   注释 (每个 pair 一行):
#     1 link: 显著 MR 方法, min MR FDR, IVW / Egger / SMR 效应同号, WNN / Pseudo p, n_ct_mronly (同一 pair 在几个 CT 中是 MR-only)
#     2 位置: peak 中心到 TSS 距离 (All_GenePeak 的 distance), promoter (|distance| <= 2 kb) / genic / intergenic, GWAS lead SNP 最近的 protein-coding 基因
#     3 金标准: Valid1_Hyper/<CT>_gold_standard.rds, countOverlaps(peak, DB[[d]][[gene]]) > 0 (同 Valid1_Hyper.R);
#              HiChIP / PCHiC 另给 3 个 CT loop 并集的结果 (*_anyCT); 命中区间坐标写入 gold_regions
#     4 GWAS: 性状, SNP, P; lead = 7 个性状中 P 最小的 SNP
#     5 GWAS Catalog (tmp/Valid6_GWAS/gwas_catalog/, P <= 5e-8): gene 是否为同一性状 / 任一 7 个性状的 reported / mapped gene; peak ±250 kb 内是否有同一性状的关联
#     6 QTL: peak 内 GWAS SNP 是否为 MR 用到的 SNP (gwas_is_iv) 及与 MR SNP 的最大 r^2 (max_r2_gwas_iv);
#            lead GWAS SNP 与 SMR lead caQTL SNP 的 r^2 (PLINK, EAS hg38); lead GWAS SNP 对该 peak 的 caQTL p / 对该 gene 的 eQTL p (CIMA 完整 summary statistics)
#     7 标记: MHC (hg38 chr6:28,510,120-33,480,577); multi_gene_peak (同一 CT 中该 peak 连到 >= 2 个候选 gene)
#   排序: 非 MHC -> n_gold -> CT 匹配的 HiChIP / PCHiC -> gwascat_trait_gene -> r2 >= 0.8 -> n_ct_mronly -> MR 显著方法数 -> min MR FDR
# 输出 tmp/Valid6_GWAS/annot/ (中间表; 整理后的例子清单由 Valid6_GWAS_4_select.py 写到 Results):
#   GWAS_example_candidates.tsv       每个 pair 一行 (含排序 rank)
#   GWAS_example_candidates_snps.tsv  pair x SNP x trait (列同 GWAS_MRonly_overlap.tsv + QTL 列); Plot10 画 GWAS SNP 轨道用
# 运行: sbatch Scripts/Validation/Valid6_GWAS_2_annot_slurm.sh

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)
  library(parallel)
})
options(width = 300)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ---- 0. 参数和路径 ----
BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
OUT_DIR  <- file.path(PROJ, "Results/Validation/Valid6_GWAS")
GS_DIR   <- file.path(PROJ, "Results/Validation/Valid1_Hyper")
TMP_DIR  <- file.path(PROJ, "tmp/Valid6_GWAS")
LD_DIR   <- file.path(TMP_DIR, paste0("ld_annot_", Sys.getpid()))
GP_DIR   <- file.path(BASE_DIR, "Data/OMIX009786/gene_peak")
EQTL_DIR <- file.path(BASE_DIR, "Data/OMIX009786/cis_eQTL_Summary_Statics/cis_eQTL")
CAQTL_DIR <- file.path(BASE_DIR, "Data/OMIX009786/cis_caQTL_Summary_Statics/cis_caQTL")
REF_DIR  <- file.path(PROJ, "data/ref_chr/hg38")
PLINK    <- file.path(BASE_DIR, "R_package/my_plink/plink")
GWAS_CACHE <- file.path(TMP_DIR, "GWAS_SNPs_P1e-5_hg38.tsv.gz")
GWASCAT  <- file.path(TMP_DIR, "gwas_catalog/gwas-catalog-download-associations-alt-full.tsv")
OLD_FILE <- file.path(OUT_DIR, "GWAS_MRonly_overlap.tsv")
ANNOT_DIR <- file.path(TMP_DIR, "annot")                             # 中间表 (Fig 4B 选例子用), 不放 Results
OUT_PAIR <- file.path(ANNOT_DIR, "GWAS_example_candidates.tsv")
OUT_SNP  <- file.path(ANNOT_DIR, "GWAS_example_candidates_snps.tsv")
dir.create(ANNOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LD_DIR, recursive = TRUE, showWarnings = FALSE)

CTS       <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
MR_M      <- c("ivw", "egger", "smr", "cct")
ASSOC_M   <- c("wnn", "pseudo")
METHODS   <- c(MR_M, ASSOC_M)
MR_LAB    <- c(ivw = "IVW", egger = "Egger", smr = "SMR", cct = "CCT")
MR_FDR    <- 0.05
OLD_MR_FDR <- 0.001
ASSOC_FDR <- 0.05
P_GWAS    <- 5e-8
TRAITS7   <- c("AS", "CD", "IBD", "RA", "SLE", "MS", "T1D")
PROM_BP   <- 2000L
LOCUS_BP  <- 250000L
MHC_GR    <- GRanges("chr6", IRanges(28510120, 33480577))
DB_NAMES  <- c("enhancerAtlas", "FANTOM5", "4DGenome", "GTEx", "HiChIP", "PCHiC")
NCORE     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# GWAS Catalog MAPPED_TRAIT (按 ", " 拆开后整词匹配) -> 7 个性状
CAT_TRAIT <- list(AS  = "ankylosing spondylitis",
                  CD  = "Crohn disease",
                  IBD = c("inflammatory bowel disease", "Crohn disease", "ulcerative colitis"),
                  RA  = c("rheumatoid arthritis", "ACPA-positive rheumatoid arthritis", "ACPA-negative rheumatoid arthritis"),
                  SLE = c("systemic lupus erythematosus", "lupus nephritis"),
                  MS  = "multiple sclerosis",
                  T1D = "type 1 diabetes mellitus")

pk_gr <- function(ids) { sp <- tstrsplit(ids, "[:-]", type.convert = TRUE); GRanges(sp[[1]], IRanges(sp[[2]] + 1L, sp[[3]])) }


####################################################################################################
#### Step 1: MR-only pair (放宽版) + GWAS SNP 重叠
####################################################################################################

keep_cols <- c("gene_id", "peak_id", "smr_variant_id", paste0(METHODS, "_p"), paste0(METHODS, "_fdr"),
               "ivw_beta", "egger_beta", "smr_b", "wnn_rho", "pseudo_rho")
pairs <- rbindlist(lapply(CTS, function(ct)
  fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")), select = keep_cols)[, cell_type := ct]))
setcolorder(pairs, "cell_type")
for (m in METHODS) {
  f <- pairs[[paste0(m, "_fdr")]]
  set(pairs, j = paste0("sig_", m), value = !is.na(f) & f < (if (m %in% MR_M) MR_FDR else ASSOC_FDR))
}
pairs[, mr_any    := Reduce(`|`, .SD), .SDcols = paste0("sig_", MR_M)]
pairs[, assoc_any := Reduce(`|`, .SD), .SDcols = paste0("sig_", ASSOC_M)]
pairs <- pairs[mr_any & !assoc_any]
pairs[, mr_methods_sig := apply(.SD, 1, function(x) paste(MR_LAB[MR_M[x]], collapse = ",")), .SDcols = paste0("sig_", MR_M)]
pairs[, n_mr_sig := rowSums(.SD), .SDcols = paste0("sig_", MR_M)]
pairs[, min_mr_fdr := do.call(pmin, c(.SD, na.rm = TRUE)), .SDcols = paste0(MR_M, "_fdr")]
pairs[, in_old_strict := min_mr_fdr < OLD_MR_FDR]
pairs[, n_ct_mronly := .N, by = .(gene_id, peak_id)]
msg("MR-only pairs (MR 任一 FDR < ", MR_FDR, ", WNN / Pseudo FDR >= ", ASSOC_FDR, "):")
print(pairs[, .N, by = cell_type])

gwas <- fread(GWAS_CACHE)[P < P_GWAS & trait %in% TRAITS7]
setnames(gwas, c("P", "Z"), c("gwas_p", "gwas_z"))
snp_u  <- unique(gwas[, .(varid_hg38, chr, pos_hg38)])
peaks  <- unique(pairs$peak_id)
hits   <- findOverlaps(pk_gr(peaks), GRanges(snp_u$chr, IRanges(snp_u$pos_hg38, width = 1)))
ov     <- data.table(peak_id = peaks[queryHits(hits)], varid_hg38 = snp_u$varid_hg38[subjectHits(hits)])
ov     <- gwas[ov, on = "varid_hg38", allow.cartesian = TRUE]

snp_cols <- c("trait", "rsid", "varid_hg38", "chr", "pos_hg38", "pos_hg19", "A1", "A2", "gwas_z", "gwas_p", "in_EAS_ref")
long <- pairs[ov[, c("peak_id", snp_cols), with = FALSE], on = "peak_id", nomatch = 0L, allow.cartesian = TRUE]
cand <- unique(long[, .(cell_type, gene_id, peak_id)])
pairs <- pairs[cand, on = .(cell_type, gene_id, peak_id)]
msg(sprintf("候选: %d pairs | %d genes | %d peaks", nrow(pairs), uniqueN(pairs$gene_id), uniqueN(pairs$peak_id)))
print(pairs[, .(n_pairs = .N, n_old_strict = sum(in_old_strict)), by = cell_type])

# 核对: 旧 Valid6 (FDR < 0.001) 去掉只有 Asthma / PBC 的 pair 后, 全部在新池且被标为 in_old_strict
old <- unique(fread(OLD_FILE)[trait %in% TRAITS7, .(cell_type, gene_id, peak_id)])
stopifnot(nrow(old[!pairs, on = .(cell_type, gene_id, peak_id)]) == 0L,
          nrow(old) == sum(pairs$in_old_strict))
msg("核对通过: 旧 Valid6 的 ", nrow(old), " 个 pair (7 个性状) 全部在新池中, 且 in_old_strict 数一致")


####################################################################################################
#### Step 2: 位置 (TSS 距离 / promoter-genic-intergenic / SNP 最近基因)
####################################################################################################

gp <- rbindlist(lapply(CTS, function(ct) {
  d <- as.data.table(readRDS(file.path(GP_DIR, paste0("All_GenePeak_", ct, ".rds"))))
  d[, .(cell_type = ct, gene_id, peak_id = sub("-", ":", peak_id), TSS, peak_center, distance)]
}))
pairs <- gp[pairs, on = .(cell_type, gene_id, peak_id)]
stopifnot(!anyNA(pairs$distance))

ucsc <- function(gr) { sl <- seqlevels(gr); genome(gr) <- NA
  renameSeqlevels(gr, setNames(ifelse(sl == "MT", "chrM", paste0("chr", sl)), sl)) }
genes_all <- ucsc(genes(EnsDb.Hsapiens.v86, columns = c("gene_name", "gene_biotype")))
genes_all <- genes_all[seqnames(genes_all) %in% paste0("chr", 1:22)]
pc_tss    <- resize(genes_all[genes_all$gene_biotype == "protein_coding"], 1, fix = "start")

pg <- pk_gr(pairs$peak_id)
pairs[, genic := vapply(seq_len(.N), function(i) {
  g <- genes_all[genes_all$gene_name == gene_id[i] & as.character(seqnames(genes_all)) == as.character(seqnames(pg[i]))]
  length(g) > 0 && overlapsAny(pg[i], g, ignore.strand = TRUE)
}, logical(1))]
pairs[, location := fifelse(abs(distance) <= PROM_BP, "promoter", fifelse(genic, "genic", "intergenic"))]

# lead GWAS SNP (每个 pair 7 个性状中 P 最小)
lead <- long[order(gwas_p), .SD[1], by = .(cell_type, gene_id, peak_id)][
  , .(cell_type, gene_id, peak_id, lead_rsid = rsid, lead_varid = varid_hg38, lead_chr = chr, lead_pos = pos_hg38,
      lead_trait = trait, lead_gwas_p = gwas_p, lead_in_EAS_ref = in_EAS_ref)]
pairs <- lead[pairs, on = .(cell_type, gene_id, peak_id)]
lg <- GRanges(pairs$lead_chr, IRanges(pairs$lead_pos, width = 1))
nn <- nearest(lg, pc_tss, ignore.strand = TRUE)
pairs[, snp_nearest_gene := pc_tss$gene_name[nn]]
pairs[, snp_nearest_dist := distance(lg, pc_tss[nn], ignore.strand = TRUE)]
pairs[, nearest_is_target := snp_nearest_gene == gene_id]

trait_sum <- long[, .(traits = paste(sort(unique(trait)), collapse = ","),
                      gwas_snps = paste(unique(sprintf("%s(%s:%.1e)", rsid, trait, gwas_p)), collapse = ";"),
                      n_gwas_snp = uniqueN(rsid)), by = .(cell_type, gene_id, peak_id)]
pairs <- trait_sum[pairs, on = .(cell_type, gene_id, peak_id)]


####################################################################################################
#### Step 3: 金标准 (CT 匹配; HiChIP / PCHiC 另给 3 CT 并集)
####################################################################################################

GS <- setNames(lapply(CTS, function(ct) readRDS(file.path(GS_DIR, paste0(ct, "_gold_standard.rds")))), CTS)
stopifnot(all(vapply(GS, function(x) identical(names(x), DB_NAMES), logical(1))))
get_gr <- function(db_obj, g) if (g %in% names(db_obj)) granges(db_obj[[g]]) else GRanges()
reg_str <- function(gr) paste(unique(sprintf("%s:%d-%d", seqnames(gr), start(gr), end(gr))), collapse = "|")

gs_rows <- lapply(seq_len(nrow(pairs)), function(i) {
  ct <- pairs$cell_type[i]; g <- pairs$gene_id[i]; p <- pk_gr(pairs$peak_id[i])
  out <- list(); regs <- character(0)
  for (d in DB_NAMES) {
    r  <- get_gr(GS[[ct]][[d]], g)
    ov <- suppressWarnings(subsetByOverlaps(r, p))
    out[[paste0("gs_", d)]] <- length(ov) > 0L
    if (length(ov)) regs <- c(regs, paste0(d, "=", reg_str(ov)))
  }
  for (d in c("HiChIP", "PCHiC")) {
    r <- do.call(c, lapply(CTS, function(c2) get_gr(GS[[c2]][[d]], g)))
    out[[paste0("gs_", d, "_anyCT")]] <- length(r) > 0L && suppressWarnings(overlapsAny(p, r))
  }
  out$gold_regions <- paste(regs, collapse = "; ")
  out
})
pairs <- cbind(pairs, rbindlist(gs_rows))
pairs[, n_gold := rowSums(.SD), .SDcols = paste0("gs_", DB_NAMES)]
pairs[, ct_loop := gs_HiChIP | gs_PCHiC]
msg("金标准命中 (pair 数):")
print(pairs[, lapply(.SD, sum), .SDcols = c(paste0("gs_", DB_NAMES), "gs_HiChIP_anyCT", "gs_PCHiC_anyCT")])

# 抽查: 第一个 n_gold > 0 的 pair, 单独手算
chk_i <- which(pairs$n_gold > 0)[1]
if (!is.na(chk_i)) {
  ct <- pairs$cell_type[chk_i]; g <- pairs$gene_id[chk_i]; p <- pk_gr(pairs$peak_id[chk_i])
  manual <- vapply(DB_NAMES, function(d) g %in% names(GS[[ct]][[d]]) && suppressWarnings(countOverlaps(p, GS[[ct]][[d]][[g]])) > 0L, logical(1))
  stopifnot(identical(unname(manual), unname(unlist(pairs[chk_i, paste0("gs_", DB_NAMES), with = FALSE]))))
  msg("抽查金标准通过: ", ct, " ", g, " ", pairs$peak_id[chk_i])
}


####################################################################################################
#### Step 4: GWAS Catalog (P <= 5e-8)
####################################################################################################

cat_dt <- fread(GWASCAT, sep = "\t", quote = "", select = c("PUBMEDID", "CHR_ID", "CHR_POS", "REPORTED GENE(S)", "MAPPED_GENE",
                                                             "SNPS", "PVALUE_MLOG", "MAPPED_TRAIT"), colClasses = "character")
setnames(cat_dt, c("pmid", "chr", "pos", "reported", "mapped", "snps", "mlog", "mtrait"))
cat_dt <- cat_dt[suppressWarnings(as.numeric(mlog)) >= -log10(P_GWAS)]
cat_dt[, tr_list := strsplit(mtrait, ", ", fixed = TRUE)]
for (tr in TRAITS7) set(cat_dt, j = paste0("is_", tr), value = vapply(cat_dt$tr_list, function(x) any(x %in% CAT_TRAIT[[tr]]), logical(1)))
cat_dt[, any7 := Reduce(`|`, .SD), .SDcols = paste0("is_", TRAITS7)]
cat_dt <- cat_dt[any7 == TRUE]
cat_dt[, genes := lapply(paste(reported, mapped, sep = ","), function(s) unique(trimws(unlist(strsplit(s, "[,;]| - | x ")))))]
msg(sprintf("GWAS Catalog: %d 条 P <= 5e-8 关联属于 7 个性状", nrow(cat_dt)))

cat_pos <- cat_dt[grepl("^[0-9]+$", chr) & grepl("^[0-9]+$", pos)]
cat_pos[, `:=`(chr = paste0("chr", chr), pos = as.integer(pos))]

cat_rows <- lapply(seq_len(nrow(pairs)), function(i) {
  g <- pairs$gene_id[i]; trs <- strsplit(pairs$traits[i], ",")[[1]]
  has_g <- vapply(cat_dt$genes, function(x) g %in% x, logical(1))
  tr_hit <- Reduce(`|`, lapply(trs, function(tr) cat_dt[[paste0("is_", tr)]]))
  same <- has_g & tr_hit
  p <- pk_gr(pairs$peak_id[i])
  near <- cat_pos[chr == as.character(seqnames(p)) & abs(pos - (start(p) + end(p)) / 2) <= LOCUS_BP]
  near_tr <- near[Reduce(`|`, lapply(trs, function(tr) near[[paste0("is_", tr)]]))]
  list(gwascat_trait_gene = any(same),
       gwascat_any_immune = any(has_g),
       gwascat_gene_pmid  = paste(head(unique(cat_dt$pmid[same]), 5), collapse = ","),
       gwascat_gene_traits = paste(unique(unlist(cat_dt$tr_list[has_g])[unlist(cat_dt$tr_list[has_g]) %in% unlist(CAT_TRAIT)]), collapse = ","),
       gwascat_trait_locus = nrow(near_tr) > 0L,
       gwascat_locus_snps  = paste(head(unique(near_tr$snps), 5), collapse = ","))
})
pairs <- cbind(pairs, rbindlist(cat_rows))
msg(sprintf("GWAS Catalog: gene 为同一性状的 reported/mapped gene: %d / %d; peak ±250 kb 有同一性状关联: %d / %d",
            sum(pairs$gwascat_trait_gene), nrow(pairs), sum(pairs$gwascat_trait_locus), nrow(pairs)))


####################################################################################################
#### Step 5: QTL (lead GWAS SNP 对 peak 的 caQTL / 对 gene 的 eQTL; r^2 与 SMR lead caQTL SNP)
####################################################################################################

long[, qtl_var := paste0(chr, "_", pos_hg38)]
pairs[, lead_qtl_var := paste0(lead_chr, "_", lead_pos)]

lookup <- function(file, keys, tag) {
  kf <- file.path(LD_DIR, paste0("keys_", tag, ".tsv"))
  fwrite(unique(keys), kf, sep = "\t", col.names = FALSE)
  cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==FNR{k[$1\"\\t\"$2]=1;next} FNR>1 && (($1\"\\t\"$2) in k){print $1\"\\t\"$2\"\\t\"$7\"\\t\"$8\"\\t\"$11\"\\t\"$12}' %s -",
                 shQuote(file), shQuote(kf))
  res <- fread(cmd = cmd, header = FALSE, sep = "\t", colClasses = c("character", "character", "numeric", "numeric", "character", "character"))
  if (nrow(res) == 0L) res <- data.table(V1 = character(), V2 = character(), V3 = numeric(), V4 = numeric(), V5 = character(), V6 = character())
  setnames(res, c("phen", "qtl_var", "p", "slope", "alt", "ref"))
  res[, tag := tag]
}
jobs <- c(lapply(CTS, function(ct) list(file = file.path(CAQTL_DIR, paste0(ct, ".txt.gz")),
                                         keys = unique(long[cell_type == ct, .(peak_id, qtl_var)]), tag = paste0("ca_", ct))),
          lapply(CTS, function(ct) list(file = file.path(EQTL_DIR, paste0(ct, ".txt.gz")),
                                         keys = unique(long[cell_type == ct, .(gene_id, qtl_var)]), tag = paste0("e_", ct))))
msg("查 QTL (6 个文件, ", NCORE, " 核) ...")
qtl <- rbindlist(mclapply(jobs, function(j) lookup(j$file, j$keys, j$tag), mc.cores = min(NCORE, 6L)))
qtl[, `:=`(type = sub("_.*", "", tag), cell_type = sub("^(ca|e)_", "", tag))]
msg("QTL 命中行数:"); print(qtl[, .N, by = tag])
qtl <- unique(qtl, by = c("type", "cell_type", "phen", "qtl_var"))

ca <- qtl[type == "ca", .(cell_type, peak_id = phen, qtl_var, caqtl_p = p, caqtl_slope = slope)]
eq <- qtl[type == "e",  .(cell_type, gene_id = phen, qtl_var, eqtl_p = p, eqtl_slope = slope)]
long <- ca[long, on = .(cell_type, peak_id, qtl_var)]
long <- eq[long, on = .(cell_type, gene_id, qtl_var)]
pairs <- ca[, .(cell_type, peak_id, lead_qtl_var = qtl_var, lead_caqtl_p = caqtl_p)][pairs, on = .(cell_type, peak_id, lead_qtl_var)]
pairs <- eq[, .(cell_type, gene_id, lead_qtl_var = qtl_var, lead_eqtl_p = eqtl_p)][pairs, on = .(cell_type, gene_id, lead_qtl_var)]

## r^2: lead GWAS SNP vs SMR lead caQTL SNP (按位置在 EAS hg38 bim 中找 ID; 一个位置多个 ID 时取第一个)
pairs[, `:=`(smr_chr = sub("_.*", "", smr_variant_id), smr_pos = as.integer(sub(".*_", "", smr_variant_id)))]
r2_res <- rbindlist(lapply(split(seq_len(nrow(pairs)), pairs$lead_chr), function(idx) {
  ch  <- sub("chr", "", pairs$lead_chr[idx[1]])
  bim <- fread(file.path(REF_DIR, paste0("chr", ch, ".bim")), header = FALSE, select = c(2, 4), col.names = c("id", "pos"))
  bim <- bim[!duplicated(pos)]
  a <- bim$id[match(pairs$lead_pos[idx], bim$pos)]
  b <- bim$id[match(pairs$smr_pos[idx], bim$pos)]
  ids <- unique(na.omit(c(a, b)))
  r2 <- rep(NA_real_, length(idx))
  if (length(ids) >= 1L) {
    pre <- file.path(LD_DIR, paste0("chr", ch)); writeLines(ids, paste0(pre, ".ids"))
    system(paste(PLINK, "--bfile", file.path(REF_DIR, paste0("chr", ch)), "--extract", paste0(pre, ".ids"),
                 "--r2 square --write-snplist --threads 1 --memory 4000 --out", pre, "--silent"))
    sl <- readLines(paste0(pre, ".snplist")); m <- as.matrix(fread(paste0(pre, ".ld"), header = FALSE))
    ia <- match(a, sl); ib <- match(b, sl); ok <- !is.na(ia) & !is.na(ib)
    r2[ok] <- m[cbind(ia[ok], ib[ok])]
  }
  data.table(row = idx, r2_gwas_caqtl = r2)
}))
pairs[r2_res$row, r2_gwas_caqtl := r2_res$r2_gwas_caqtl]
pairs[lead_pos == smr_pos & is.na(r2_gwas_caqtl), r2_gwas_caqtl := 1]
msg(sprintf("r2(lead GWAS SNP, SMR lead caQTL SNP): NA %d | >= 0.8 %d / %d", sum(is.na(pairs$r2_gwas_caqtl)),
            sum(pairs$r2_gwas_caqtl >= 0.8, na.rm = TRUE), nrow(pairs)))

## ---- 5b. peak 内的 GWAS SNP 是否是 MR 用到的 SNP (工具变量) ----
# MR 的 SNP = 该 pair 在 clumped caQTL (r2 0.99) 中的变异 ∩ 该 gene 的 eQTL 变异 (同 MR_AllPairs.R; 等位一致 / 去重步骤这里不重做)
#   gwas_is_iv     : peak 内任一 GWAS SNP 本身就是 MR 的 SNP
#   max_r2_gwas_iv : peak 内 GWAS SNP 与 MR SNP 的最大 r2 (EAS hg38; 被 clumping 合并掉的 SNP 会以高 r2 的代理出现)
CLUMP_DIR <- file.path(CAQTL_DIR, "LD_clump")
iv <- rbindlist(lapply(CTS, function(ct) {
  cl <- fread(file.path(CLUMP_DIR, paste0(ct, "_Clumped_r2_0.99_kb_20000.txt.gz")), select = c("phenotype_id", "variant_id", "gene_id"))
  setnames(cl, c("peak_id", "qtl_var", "gene_id"))
  cl[pairs[cell_type == ct, .(gene_id, peak_id)], on = .(peak_id, gene_id), nomatch = 0L][, cell_type := ct]
}))
eq_iv <- rbindlist(mclapply(CTS, function(ct)
  lookup(file.path(EQTL_DIR, paste0(ct, ".txt.gz")), unique(iv[cell_type == ct, .(gene_id, qtl_var)]), paste0("eiv_", ct))[, cell_type := ct],
  mc.cores = min(NCORE, 3L)))
iv <- unique(iv[unique(eq_iv[, .(cell_type, gene_id = phen, qtl_var)]), on = .(cell_type, gene_id, qtl_var), nomatch = 0L])
iv_n <- iv[, .(n_iv = .N), by = .(cell_type, gene_id, peak_id)]
pairs <- iv_n[pairs, on = .(cell_type, gene_id, peak_id)]
stopifnot(all(pairs$n_iv >= 3L))                                   # MR 要求共享变异 >= 3
long[, is_iv := FALSE][iv, is_iv := TRUE, on = .(cell_type, gene_id, peak_id, qtl_var)]
pairs <- long[, .(gwas_is_iv = any(is_iv)), by = .(cell_type, gene_id, peak_id)][pairs, on = .(cell_type, gene_id, peak_id)]

pos_of <- function(v) as.integer(sub(".*_", "", v))
key_cols <- c("cell_type", "gene_id", "peak_id")
ivr <- rbindlist(lapply(split(seq_len(nrow(pairs)), pairs$lead_chr), function(idx) {
  ch  <- sub("chr", "", pairs$lead_chr[idx[1]])
  bim <- fread(file.path(REF_DIR, paste0("chr", ch, ".bim")), header = FALSE, select = c(2, 4), col.names = c("id", "pos"))
  bim <- bim[!duplicated(pos)]
  sp  <- pairs[idx, ..key_cols]
  g_pos <- long[sp, on = key_cols][, c(key_cols, "pos_hg38"), with = FALSE]
  i_pos <- iv[sp, on = key_cols, nomatch = 0L][, pos := pos_of(qtl_var)]
  ids <- unique(na.omit(bim$id[match(unique(c(g_pos$pos_hg38, i_pos$pos)), bim$pos)]))
  pre <- file.path(LD_DIR, paste0("iv_chr", ch)); writeLines(ids, paste0(pre, ".ids"))
  system(paste(PLINK, "--bfile", file.path(REF_DIR, paste0("chr", ch)), "--extract", paste0(pre, ".ids"),
               "--r2 square --write-snplist --threads 1 --memory 4000 --out", pre, "--silent"))
  sl <- readLines(paste0(pre, ".snplist")); m <- as.matrix(fread(paste0(pre, ".ld"), header = FALSE))
  rbindlist(lapply(idx, function(i) {
    k  <- pairs[i, ..key_cols]
    ga <- match(bim$id[match(unique(g_pos[k, on = key_cols]$pos_hg38), bim$pos)], sl)
    ia <- match(bim$id[match(unique(i_pos[k, on = key_cols]$pos), bim$pos)], sl)
    ga <- ga[!is.na(ga)]; ia <- ia[!is.na(ia)]
    data.table(row = i, max_r2_gwas_iv = if (length(ga) && length(ia)) max(m[ga, ia, drop = FALSE], na.rm = TRUE) else NA_real_)
  }))
}))
pairs[ivr$row, max_r2_gwas_iv := ivr$max_r2_gwas_iv]
pairs[gwas_is_iv == TRUE, max_r2_gwas_iv := 1]
msg(sprintf("GWAS SNP 是 MR 的 SNP: %d / %d pairs; 与 MR SNP r2 >= 0.8 (含本身): %d / %d", sum(pairs$gwas_is_iv), nrow(pairs),
            sum(pairs$max_r2_gwas_iv >= 0.8, na.rm = TRUE), nrow(pairs)))


####################################################################################################
#### Step 6: 标记、方向一致、排序、输出
####################################################################################################

pairs[, MHC := overlapsAny(pk_gr(peak_id), MHC_GR)]
pairs[, multi_gene_peak := .N >= 2L, by = .(cell_type, peak_id)]
pairs[, multi_gene_peak_genes := paste(sort(unique(gene_id)), collapse = ","), by = .(cell_type, peak_id)]
pairs[, mr_dir_consistent := sign(ivw_beta) == sign(egger_beta) & sign(ivw_beta) == sign(smr_b)]

pairs[, r2_hi := !is.na(r2_gwas_caqtl) & r2_gwas_caqtl >= 0.8]
setorderv(pairs, c("MHC", "n_gold", "ct_loop", "gwascat_trait_gene", "r2_hi", "n_ct_mronly", "n_mr_sig", "min_mr_fdr"),
          order = c(1, -1, -1, -1, -1, -1, -1, 1))
pairs[, rank := .I]

out_cols <- c("rank", "cell_type", "gene_id", "peak_id", "in_old_strict", "mr_methods_sig", "n_mr_sig", "min_mr_fdr",
              paste0(MR_M, "_fdr"), "wnn_p", "pseudo_p", "ivw_beta", "egger_beta", "smr_b", "mr_dir_consistent", "n_ct_mronly",
              "distance", "location", "traits", "n_gwas_snp", "gwas_snps", "lead_rsid", "lead_trait", "lead_gwas_p",
              "snp_nearest_gene", "snp_nearest_dist", "nearest_is_target",
              "n_gold", paste0("gs_", DB_NAMES), "gs_HiChIP_anyCT", "gs_PCHiC_anyCT", "gold_regions",
              "gwascat_trait_gene", "gwascat_any_immune", "gwascat_gene_traits", "gwascat_gene_pmid", "gwascat_trait_locus", "gwascat_locus_snps",
              "smr_variant_id", "r2_gwas_caqtl", "lead_caqtl_p", "lead_eqtl_p", "n_iv", "gwas_is_iv", "max_r2_gwas_iv",
              "MHC", "multi_gene_peak", "multi_gene_peak_genes")
fwrite(pairs[, ..out_cols], OUT_PAIR, sep = "\t")

long_cols <- c("cell_type", "gene_id", "peak_id", "mr_methods_sig", paste0(METHODS, "_p"), "trait", "rsid", "varid_hg38", "gwas_p",
               paste0(METHODS, "_fdr"), "ivw_beta", "egger_beta", "smr_b", "wnn_rho", "pseudo_rho", "smr_variant_id",
               "chr", "pos_hg38", "pos_hg19", "A1", "A2", "gwas_z", "in_EAS_ref", "caqtl_p", "caqtl_slope", "eqtl_p", "eqtl_slope", "is_iv")
setorder(long, cell_type, gene_id, peak_id, gwas_p)
fwrite(long[, ..long_cols], OUT_SNP, sep = "\t")

## ---- 汇总 ----
msg("候选按 CT / 性状 (pair 数):")
print(pairs[, .N, by = cell_type])
print(long[, .(n_pairs = uniqueN(paste(cell_type, gene_id, peak_id))), by = trait][order(-n_pairs)])
msg("n_gold 分布 (非 MHC):"); print(pairs[MHC == FALSE, .N, by = n_gold][order(-n_gold)])
show <- c("rank", "cell_type", "gene_id", "peak_id", "mr_methods_sig", "min_mr_fdr", "n_ct_mronly", "distance", "location", "traits",
          "lead_rsid", "lead_gwas_p", "n_gold", paste0("gs_", DB_NAMES), "gwascat_trait_gene", "r2_gwas_caqtl", "lead_caqtl_p", "lead_eqtl_p", "MHC")
msg("排序前 40:")
print(copy(pairs[1:min(40, .N), ..show])[, `:=`(min_mr_fdr = signif(min_mr_fdr, 2), lead_gwas_p = signif(lead_gwas_p, 2),
                                          r2_gwas_caqtl = round(r2_gwas_caqtl, 2), lead_caqtl_p = signif(lead_caqtl_p, 2),
                                          lead_eqtl_p = signif(lead_eqtl_p, 2))][])

unlink(LD_DIR, recursive = TRUE)
msg("Saved: ", OUT_PAIR)
msg("Saved: ", OUT_SNP)
