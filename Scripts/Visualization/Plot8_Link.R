# ==============================================================================
# Claude_project 移植版 (2026-10-09): 由 20261006_paper/Scripts/Plot8_val_Link.sh 复制, 作图代码不动, 只改:
#   GWAS pairs = Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv (替代旧 Valid8_matchGWAS/G1_Focus_MRonly_GWAS.tsv)
#   TRAIT_FILES = 项目的 9 个 GWAS (与 Valid6 / Valid2 一致)
#   MR = Results/MR/<CT>_AllPairs_results.tsv.gz (已过滤, 不再做 sigma2 过滤); LD 参考 = data/ref_chr/hg38
#   example pairs 改为 source 01_example_pairs.R (只在 PAIR_SET = "example" 时用; 本项目只做 gwas 模式)
#   输出 = Results/Figures/Fig8_Link_GWAS.pdf; PLINK / 临时文件 = tmp/plink_fig8/
# 运行: sbatch Scripts/Visualization/Plot8_Link_slurm.sh
# ==============================================================================

# 2026.10.06
# plot Peak-Gene link
# Figures address: /scratch/user/tfcui2025/CIMA/Results/RunMR/Materials_runMR/20261006_paper/Materials
# Script: /scratch/user/tfcui2025/CIMA/Results/RunMR/Materials_runMR/20261006_paper/Scripts

# 2026.10.07
# Peak-gene link locus plot + MR SNP track
# 轨道 (上 -> 下): Gene / Peaks / MR SNPs / 6 种方法各自的 links
#   links: 直接来自 MR_results_formal (不调用 LinkPeaks); 只画 FDR < 0.05; 颜色 = 方向, 线宽 = -log10(FDR)
#   MR SNPs: 对每个 example pair 按 MR_results_formal 的流程重跑 QTL 对齐 -> 去重 -> PLINK LD 过滤, 得到 MR 实际使用的 SNP;
#            一行 = 该基因的一个 example peak; 竖线 = MR instrument, 金色三角 = SMR lead SNP, 浅红底 = example peak 位置
# 三个细胞类型画在同一个 PDF 的不同页上, 每个基因一页
# 输出: 只输出一个 PDF (Fig8_Link.pdf), 不再输出 CSV
# 运行: 交互式 R 中粘贴本文件代码 (Scripts/Plot8_val_Link.sh)

# 2026.10.08
# GWAS 模式: 画 Valid8 文件 (G1_Focus_MRonly_GWAS.tsv, 已严格筛选) 中的全部 pairs, 输出 Fig8_Link_GWAS.pdf
#   PAIR_SET = "gwas" | "example" (example = 原 24x3 个 example pairs, 输出 Fig8_Link.pdf)
#   GWAS SNP (紫色) 标在 Peaks 轨 / MR SNP 轨 / 新增的 caQTL-eQTL 轨 (-log10 P, caQTL 用全量 summary stats)
#   example 版式 (上 -> 下): Gene / Peaks / MR SNPs / caQTL-eQTL / 6 种方法各自的 links
#   gwas 版式: 一页 = 一个 link x 一个 trait; 左 = GWAS (该 trait, selected_sumstats) / caQTL / eQTL 的 -log10 P,
#              右 = Gene / Peaks / links; 两边 x 窗口相同; 不画图例; peak 位置加紫色竖带 (alpha 0.5)
# matched GWAS files: /scratch/user/tfcui2025/CIMA/Results/RunMR/Materials_runMR/20260430_LinkPeakGene_more/MRresult_20260911/Valid8_matchGWAS/G1_Focus_MRonly_GWAS.tsv

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
  library(GenomicRanges); library(EnsDb.Hsapiens.v86); library(Signac)
})
options(width = 300)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
 
###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 路径与参数

CIMA_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(CIMA_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
PLINK    <- file.path(CIMA_DIR, "R_package/my_plink/plink")
REF_DIR  <- file.path(PROJ, "data/ref_chr/hg38")
OMIX     <- file.path(CIMA_DIR, "Data/OMIX009786")
LD_DIR   <- file.path(PROJ, "tmp/plink_fig8", paste0("ld_fig8_", Sys.getpid()))   # PLINK / 临时文件放项目 tmp/, 结束自动删除
FIG_DIR  <- file.path(PROJ, "Results/Figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LD_DIR,  recursive = TRUE, showWarnings = FALSE)

PAIR_SET   <- "gwas"                 # "gwas": Valid6 GWAS 重合 pairs (本项目只做这个); "example": 原 example pairs
GWAS_FILE  <- file.path(PROJ, "Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv")
PDF_FILE   <- file.path(FIG_DIR, if (PAIR_SET == "gwas") "Fig8_Link_GWAS.pdf" else "Fig8_Link.pdf")

# GWAS 轨 (gwas 版式): munged sumstats 只有 rsID + Z, 位置用 1000G EUR hg19 bim + liftOver 到 hg38 (同 Valid8 G1_find_gwas.R)
GWAS_SS_DIR <- file.path(CIMA_DIR, "R_package/S_LDSC/GWAS_data/selected_sumstats")
BIM19       <- file.path(CIMA_DIR, "R_package/S_LDSC/reference_genome/1000G_EUR_Phase3_plink/1000G.EUR.QC.")
CHAIN19to38 <- file.path(CIMA_DIR, "Data/PBMC_validation/hg19ToHg38.over.chain")
GWS_P       <- 5e-8
TRAIT_FILES <- c(AS     = "CUSTOM.AS",                                       # 项目的 9 个 GWAS, 与 Valid6_GWAS.R 一致
                 CD     = "CUSTOM.CD",
                 IBD    = "PASS.Inflammatory_Bowel_Disease.deLange2017",
                 RA     = "PASS.Rheumatoid_Arthritis.Saevarsdottir2022",
                 SLE    = "PASS.Lupus.Bentham2015",
                 MS     = "PASS_Multiple_sclerosis",
                 T1D    = "PASS.Type_1_Diabetes.Chiou2021",
                 Asthma = "GBMI.Asthma",
                 PBC    = "PASS.Primary_Biliary_Cirrhosis.Cordell2015")

eqtl_file        <- function(ct) file.path(OMIX, "cis_eQTL_Summary_Statics/cis_eQTL", paste0(ct, ".txt.gz"))
caqtl_clump_file <- function(ct) file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL/LD_clump", paste0(ct, "_Clumped_r2_0.99_kb_20000.txt.gz"))
caqtl_full_file  <- function(ct) file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL", paste0(ct, ".txt.gz"))   # GWAS SNP 常不在 clump 文件里

# 作图参数
FDR_CUT    <- 0.05
SCORE_CAP  <- 20                     # -log10(FDR) 上限, 所有页面共用同一线宽标尺
PAD_FRAC   <- 0.08; PAD_MIN <- 1e4   # 窗口两侧留白
EMPTY_HALF <- 5e4                    # 无显著 link 且无 SNP 时窗口 = TSS +/- 50kb
LW_RANGE   <- c(1.2, 3.2)            # link 线宽范围
BASE_SIZE  <- 17                     # 全局字号
PAGE_W     <- 14; PAGE_H <- 18       # example 版式 (竖排)
PAGE_W_G   <- 22; PAGE_H_G <- 11     # gwas 版式 (左右两栏)
BAND_FRAC  <- 0.0025; BAND_ALPHA <- 0.5   # 紫色竖带: 最小半宽 = 窗口宽 x BAND_FRAC

# MR 固定参数 (与 MR_results_formal 一致)
MIN_SNP <- 3; DEDUP_TOL <- 1e-3; KEEP_ALLELE <- TRUE; P_FLOOR <- 1e-300

# 方法名 / 顺序: 全文统一 (见 CLAUDE.md); Causal Joint (CCT) 无方向
METHODS <- data.table(method  = c("Causal SMR", "Causal IVW", "Causal Egger", "Causal Joint", "Marginal Multiomics", "Marginal Pseudobulk"),
                      fdr_col = c("smr_fdr", "ivw_fdr", "egger_fdr", "cct_fdr", "wnn_fdr", "pseudo_fdr"),
                      eff_col = c("smr_b", "ivw_beta", "egger_beta", NA, "wnn_rho", "pseudo_rho"),
                      p_col   = c("smr_p", "ivw_p", "egger_p", "cct_p", "wnn_p", "pseudo_p"))
DIR_COL  <- c("Positive" = "#D7191C", "Negative" = "#2C7BB6", "No direction" = "grey15")
GENE_COL <- "#1B1B1B"
PEAK_COL <- c("TRUE" = "#E64B35", "FALSE" = "grey25")
SNP_COL  <- "#3C5488"
LEAD_FILL <- "#FFD34E"; LEAD_BORDER <- "#5C4A00"
GWAS_COL  <- "#7E1F86"; QTL_BG <- "grey70"

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Example pairs

source(file.path(PROJ, "Scripts/Visualization/01_example_pairs.R"))   # EXAMPLE_PAIRS (gene_id, peak_id)
CELLTYPES <- EXAMPLE_PAIRS

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### GWAS 重合 SNP (Valid8): 一行 = 一个 (cell_type, gene, peak, rsid), 多个 trait 合并

GWAS_RAW <- fread(GWAS_FILE)
GWAS <- GWAS_RAW[, .(variant_id = paste0("chr", sub("^chr", "", chr[1]), "_", pos_hg38[1]), pos = as.numeric(pos_hg38[1]), P = min(gwas_p),
                     traits = paste(unique(trait[order(gwas_p)]), collapse = ",")),
                 by = .(cell_type, gene_id, peak_id, rsid)]
GWAS[, label := sprintf("%s (%s; P=%.1e)", rsid, traits, P)]
msg(sprintf("GWAS: %d SNP-pair 行 | %d pairs", nrow(GWAS), uniqueN(GWAS[, .(cell_type, gene_id, peak_id)])))

# gwas 版式按 trait 分页: 一行 = 一个 (cell_type, gene, peak, trait, rsid)
GWAS_T <- unique(GWAS_RAW[, .(cell_type, gene_id, peak_id, trait, rsid, pos = as.numeric(pos_hg38), P = gwas_p,
                              variant_id = paste0("chr", sub("^chr", "", chr), "_", pos_hg38))])
msg(sprintf("GWAS 分页: %d 页 (link x trait)", uniqueN(GWAS_T[, .(cell_type, gene_id, peak_id, trait)])))

if (PAIR_SET == "gwas") {   # 文件中的 pairs 全部画 (已在上游严格筛选)
  gp <- unique(GWAS_RAW[, .(cell_type, gene_id, peak_id)])
  CELLTYPES <- lapply(setNames(names(CELLTYPES), names(CELLTYPES)), function(ct) gp[cell_type == ct, .(gene_id, peak_id)])
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 基因注释 (EnsDb.Hsapiens.v86, hg38, UCSC 风格染色体名)

# 不用 seqlevelsStyle<- "UCSC": 它会联网下载 UCSC chromInfo, 计算节点无网络会报错; 直接加 "chr" 前缀
to_ucsc <- function(x) sub("^chrMT$", "chrM", paste0("chr", as.character(x)))
gene_gr <- genes(EnsDb.Hsapiens.v86, filter = SeqNameFilter(c(1:22, "X", "Y")))
gene_dt <- data.table(chr = to_ucsc(seqnames(gene_gr)), start = start(gene_gr), end = end(gene_gr),
                      strand = as.character(strand(gene_gr)), gene_name = gene_gr$gene_name)
gene_dt[, tss := fifelse(strand == "+", start, end)]

exon_gr <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
exon_dt <- data.table(chr = to_ucsc(seqnames(exon_gr)), start = start(exon_gr), end = end(exon_gr),
                      gene_name = exon_gr$gene_name, type = as.character(exon_gr$type))[type %chin% c("exon", "cds", "utr")]
msg(sprintf("EnsDb v86: %d genes | %d exon/cds/utr 片段", nrow(gene_dt), nrow(exon_dt)))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### MR SNP 函数 (QTL 对齐 -> 去重 -> LD 过滤, 与 MR_results_formal 一致; 不重跑 GLS)

dedup_keep <- function(x, y, tol = DEDUP_TOL) {
  keep <- logical(length(x)); ref_x <- ref_y <- numeric(0)
  for (j in seq_along(x)) {
    if (length(ref_x) && any(abs(ref_x - x[j]) <= tol & abs(ref_y - y[j]) <= tol)) next
    keep[j] <- TRUE; ref_x <- c(ref_x, x[j]); ref_y <- c(ref_y, y[j])
  }
  keep
}

make_plink_id <- function(variant_id, ref, alt) {
  chr <- sub("^chr", "", sub("_.*$", "", variant_id)); bp <- sub("^[^_]*_", "", variant_id)
  chr <- fcase(chr == "X", "23", chr == "Y", "24", chr %chin% c("M", "MT"), "26", default = chr)
  paste(chr, bp, ref, alt, sep = ":")
}
ref_bfile <- function(peak) file.path(REF_DIR, sub(":.*$", "", peak))

run_plink_ld <- function(ids, bfile, pre) {
  fs <- paste0(pre, c(".snps.txt", ".ld", ".snplist", ".log", ".nosex")); on.exit(unlink(fs), add = TRUE)
  fwrite(data.table(unique(ids[!is.na(ids)])), fs[1], col.names = FALSE)
  kao <- if (KEEP_ALLELE) "--keep-allele-order" else ""
  cmd <- paste(PLINK, "--bfile", bfile, "--extract", fs[1], kao, "--r square --write-snplist --threads 1 --memory 4000 --out", pre, "--silent")
  st  <- suppressWarnings(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE))
  if (st != 0 || !file.exists(fs[2]) || !file.exists(fs[3])) return(NULL)
  snps <- fread(fs[3], header = FALSE)[[1]]
  R <- as.matrix(fread(fs[2], header = FALSE, na.strings = c("nan", "-nan", "NaN", "NA"))); storage.mode(R) <- "numeric"
  dimnames(R) <- list(snps, snps); R
}

get_mr_snps <- function(gene, peak, ca_data, eq_data) {
  out <- function(d, s) list(dat = d, status = s)
  ca <- ca_data[.(peak, gene), nomatch = 0L]; eq <- eq_data[.(gene), nomatch = 0L]
  if (!nrow(ca) || !nrow(eq)) return(out(NULL, "no_qtl"))
  ca[, snp_id := make_plink_id(variant_id, `A2(REF)`, `A1(ALT/effect allele)`)]
  ca <- ca[, .(variant_id, snp_id, ea_x = `A1(ALT/effect allele)`, oa_x = `A2(REF)`, beta_x = slope, se_x = slope_se, p_x = pval_nominal)]
  eq <- eq[, .(variant_id,          ea_y = `A1(ALT/effect allele)`, oa_y = `A2(REF)`, beta_y = slope, se_y = slope_se, p_y = pval_nominal)]
  dat <- merge(ca, eq, by = "variant_id")[ea_x == ea_y & oa_x == oa_y]
  dat <- dat[dedup_keep(beta_x, beta_y)]
  if (nrow(dat) < MIN_SNP) return(out(dat, "few_snp_qtl"))

  ld <- run_plink_ld(dat$snp_id, ref_bfile(peak), file.path(LD_DIR, gsub("[^A-Za-z0-9]", "_", paste(peak, gene, sep = "__"))))
  if (is.null(ld)) return(out(NULL, "plink_fail"))
  d <- dat[match(rownames(ld), snp_id)][!is.na(variant_id)]
  R <- ld[d$snp_id, d$snp_id, drop = FALSE]
  if (anyNA(R)) { bad <- which(apply(R, 1, function(x) any(!is.finite(x)))); if (length(bad)) d <- d[-bad] }
  out(d, if (nrow(d) >= MIN_SNP) "ok" else "few_snp_ld")
}

# 只读第一列在 keys 里的行
read_filtered <- function(file, keys) {
  kf <- tempfile(tmpdir = LD_DIR, fileext = ".txt"); writeLines(unique(keys), kf); on.exit(unlink(kf), add = TRUE)
  fread(cmd = sprintf("zcat %s | awk 'NR==FNR{a[$1]=1; next} FNR==1 || ($1 in a)' %s -", shQuote(file), shQuote(kf)))
}

# 一个细胞类型: 读 QTL -> 每个 example pair 取 MR SNP -> 与 formal 结果对账
collect_mr_snps <- function(ct, pairs, res) {
  pairs <- unique(pairs[, .(gene_id, peak_id)])
  eq <- read_filtered(eqtl_file(ct), unique(pairs$gene_id))[, .(phenotype_id, variant_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se, pval_nominal)]
  setkey(eq, phenotype_id, variant_id)
  ca <- fread(caqtl_clump_file(ct))[unique(pairs[, .(phenotype_id = peak_id, gene_id)]), on = .(phenotype_id, gene_id), nomatch = 0L,
                                    .(phenotype_id, variant_id, gene_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se, pval_nominal)]
  setkey(ca, phenotype_id, gene_id)

  # caQTL / eQTL 轨用的全量 summary stats (caQTL 全量文件第一列 = peak)
  ca_full <- read_filtered(caqtl_full_file(ct), unique(pairs$peak_id))[, .(phenotype_id, variant_id, pval_nominal)]
  qtl <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    g <- pairs$gene_id[i]; pk <- pairs$peak_id[i]
    rbind(ca_full[phenotype_id == pk, .(layer = "caQTL", variant_id, p = pval_nominal)],
          eq[.(g), nomatch = 0L, .(layer = "eQTL", variant_id, p = pval_nominal)])[, `:=`(gene_id = g, peak_id = pk)]
  }))
  qtl[, `:=`(pos = as.numeric(sub("^[^_]*_([0-9]+).*$", "\\1", variant_id)), mlp = -log10(pmax(p, P_FLOOR)))]
  qtl <- unique(qtl, by = c("gene_id", "peak_id", "layer", "variant_id"))
  rm(ca_full); invisible(gc())

  snps <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    g <- pairs$gene_id[i]; pk <- pairs$peak_id[i]
    o <- tryCatch(get_mr_snps(g, pk, ca, eq), error = function(e) list(dat = NULL, status = paste("error:", conditionMessage(e))))
    msg(sprintf("  [%s] %s__%s: %s, %d SNPs", ct, pk, g, o$status, if (is.null(o$dat)) 0L else nrow(o$dat)))
    if (is.null(o$dat) || !nrow(o$dat)) return(NULL)
    o$dat[, .(celltype = ct, gene_id = g, peak_id = pk, variant_id, pos = as.numeric(sub("^[^_]*_([0-9]+).*$", "\\1", variant_id)),
              beta_x, p_x, beta_y, p_y, status = o$status)]
  }), fill = TRUE)
  if (!nrow(snps)) return(list(snps = snps, qtl = qtl))

  # 对账: SNP 数 vs ivw_n_snp, lead SNP vs smr_variant_id
  fm <- res[pairs, on = .(gene_id, peak_id), nomatch = NULL, .(gene_id, peak_id, ivw_n_snp, smr_variant_id)]
  snps[fm, smr_lead := i.smr_variant_id, on = .(gene_id, peak_id)]
  snps[, is_lead := !is.na(smr_lead) & variant_id == smr_lead][, smr_lead := NULL]
  chk <- snps[, .(n_snp_rerun = .N, lead_found = any(is_lead)), by = .(gene_id, peak_id)][fm, on = .(gene_id, peak_id)]
  msg(sprintf("[%s] 对账: SNP 数一致 %d / %d, lead SNP 找到 %d / %d", ct,
              chk[, sum(n_snp_rerun == ivw_n_snp, na.rm = TRUE)], nrow(chk), chk[, sum(lead_found, na.rm = TRUE)], nrow(chk)))
  if (chk[!(n_snp_rerun %in% ivw_n_snp) | !lead_found, .N]) print(chk[!(n_snp_rerun %in% ivw_n_snp) | !lead_found])
  list(snps = snps[], qtl = qtl)
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Link 函数

parse_peak <- function(dt) {
  dt[, c("pk_chr", "pk_start", "pk_end") := tstrsplit(peak_id, "[:-]")]
  dt[, `:=`(pk_start = as.integer(pk_start), pk_end = as.integer(pk_end))]
  dt[, pk_mid := (pk_start + pk_end) / 2]
  dt[]
}

# 每个方法 FDR < cut 的 links (长表: 一行 = 一个 gene-peak-method)
get_sig_links <- function(res, genes, fdr_cut = FDR_CUT) {
  sub <- res[gene_id %chin% genes]
  lk <- rbindlist(lapply(seq_len(nrow(METHODS)), function(k) {
    m <- METHODS[k]; f <- sub[[m$fdr_col]]
    eff <- if (is.na(m$eff_col)) rep(NA_real_, nrow(sub)) else sub[[m$eff_col]]
    idx <- which(!is.na(f) & f < fdr_cut)
    data.table(gene_id = sub$gene_id[idx], peak_id = sub$peak_id[idx], method = m$method, fdr = f[idx], effect = eff[idx])
  }))
  if (!nrow(lk)) lk <- data.table(gene_id = character(), peak_id = character(), method = character(), fdr = numeric(), effect = numeric())
  parse_peak(lk)
  lk[, direction := fcase(is.na(effect), "No direction", effect > 0, "Positive", default = "Negative")]
  lk[, score := -log10(pmax(fdr, 1e-300))]
  lk[]
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### GWAS 轨数据: munged sumstats (SNP, Z) -> P = 2*pnorm(-|Z|); rsID -> hg19 (1000G EUR bim) -> liftOver hg38 (1-to-1); 同 Valid8 G1_find_gwas.R

.gw_cache <- new.env()
get_gwas_track <- function(trait, chr) {
  key <- paste(trait, chr, sep = "__")
  if (!is.null(.gw_cache[[key]])) return(.gw_cache[[key]])
  ss_key <- paste0("ss__", trait)
  if (is.null(.gw_cache[[ss_key]])) {
    ss <- fread(file.path(GWAS_SS_DIR, paste0(TRAIT_FILES[[trait]], ".sumstats.gz")), select = c("SNP", "Z"))[!is.na(Z)]
    .gw_cache[[ss_key]] <- ss[, .(rsid = SNP, P = 2 * pnorm(-abs(Z)))]
  }
  if (is.null(.gw_cache$chain)) .gw_cache$chain <- rtracklayer::import.chain(CHAIN19to38)
  bim <- fread(paste0(BIM19, sub("^chr", "", chr), ".bim"), header = FALSE, select = c(2, 4), col.names = c("rsid", "pos19"))
  d   <- .gw_cache[[ss_key]][bim, on = "rsid", nomatch = 0L][!duplicated(rsid)]
  gl  <- rtracklayer::liftOver(GRanges(chr, IRanges(d$pos19, width = 1)), .gw_cache$chain)
  ok  <- lengths(gl) == 1L
  g38 <- unlist(gl[ok])
  out <- data.table(rsid = d$rsid[ok], chr38 = as.character(seqnames(g38)), pos = start(g38), P = d$P[ok])[chr38 == chr]
  out[, `:=`(mlp = -log10(pmax(P, P_FLOOR)), chr38 = NULL)]
  msg(sprintf("    GWAS 轨 %s %s: %d SNPs (sumstats %d 个在 bim 中, liftOver 1-to-1 %d)", trait, chr, nrow(out), nrow(d), sum(ok)))
  .gw_cache[[key]] <- out
  out
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 作图: 共用的 locus 计算 + 各条轨

x_axis_off  <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank())
x_axis_on   <- theme(axis.text.x = element_text(size = BASE_SIZE - 2, color = "black"), axis.title.x = element_text(face = "bold"))
strip_theme <- theme(strip.placement = "outside", strip.background = element_blank(),
                     strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = BASE_SIZE), panel.spacing.y = unit(0.25, "cm"))
track_theme <- function() theme_classic(base_size = BASE_SIZE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, face = "bold", size = BASE_SIZE))
mb_scale <- scale_x_continuous(labels = function(x) sprintf("%.2f", x / 1e6))

# Gene / Peaks 轨用与 links 相同的 facet strip 标签 (右栏标签列对齐); 返回 geom_blank + facet + strip 主题
strip_facet <- function(lab, S) list(geom_blank(data = data.table(lab = lab, x = S$tg$tss, y = 0), aes(x, y)),
                                     facet_grid(lab ~ ., switch = "y"), strip_theme)

# 目标基因 TSS / 显著 links / example peaks / MR SNP / GWAS SNP / 窗口
locus_setup <- function(g, lk_all, ex_pairs, snp_dt, qtl_dt, gw_dt) {
  lk <- copy(lk_all[gene_id == g])
  ex <- parse_peak(unique(ex_pairs[gene_id == g, .(gene_id, peak_id)]))
  sd <- snp_dt[gene_id == g]
  qd <- qtl_dt[gene_id == g]
  gw <- gw_dt[gene_id == g & peak_id %chin% ex$peak_id]

  # 目标基因 TSS: 同名多条时取与 example peak 同染色体、离 peak 最近的
  hint_chr <- names(which.max(table(c(ex$pk_chr, lk$pk_chr))))
  hint_pos <- median(c(ex$pk_mid, lk$pk_mid))
  cand <- gene_dt[gene_name == g & chr == hint_chr][order(abs(tss - hint_pos))]
  if (!nrow(cand)) { msg(sprintf("    %s: EnsDb v86 中找不到 (或染色体不符), 跳过", g)); return(NULL) }
  tg <- cand[1]
  lk <- lk[pk_chr == tg$chr]; ex <- ex[pk_chr == tg$chr]
  ex[, row := .I]
  if (nrow(sd)) sd <- sd[ex[, .(peak_id, row)], on = "peak_id", nomatch = 0L]
  if (nrow(gw)) gw <- gw[ex[, .(peak_id, row)], on = "peak_id", nomatch = 0L]
  ex[, has_gwas := peak_id %chin% gw$peak_id]

  # 窗口: TSS + 显著 peak + example peak + MR SNP + GWAS SNP
  xs <- c(tg$tss, lk$pk_start, lk$pk_end, ex$pk_start, ex$pk_end, sd$pos, gw$pos)
  if (length(xs) > 1) { pad <- max(PAD_MIN, diff(range(xs)) * PAD_FRAC); xlim <- c(min(xs) - pad, max(xs) + pad) } else xlim <- tg$tss + c(-EMPTY_HALF, EMPTY_HALF)
  wid <- diff(xlim)
  ex[, hw := pmax((pk_end - pk_start) / 2, wid * 0.003)]
  list(g = g, tg = tg, lk = lk, ex = ex, sd = sd, qd = qd, gw = gw, xlim = xlim, wid = wid,
       tss_line = geom_vline(xintercept = tg$tss, linetype = "dashed", color = "grey60", linewidth = 0.5))
}

# 紫色细竖带: example peak (内含 GWAS SNP) 的位置, 贯穿整条轨 / 所有 facet
peak_band <- function(S) {
  b <- S$ex[has_gwas == TRUE]; if (!nrow(b)) b <- S$ex
  hw <- pmax((b$pk_end - b$pk_start) / 2, S$wid * BAND_FRAC)
  annotate("rect", xmin = b$pk_mid - hw, xmax = b$pk_mid + hw, ymin = -Inf, ymax = Inf, fill = GWAS_COL, alpha = BAND_ALPHA)
}

## A. Gene track (只画目标基因)
track_gene <- function(S, band = FALSE, strip = FALSE) {
  tg <- S$tg; xlim <- S$xlim
  gx  <- exon_dt[chr == tg$chr & gene_name == S$g]
  arr <- tg$tss + ifelse(tg$strand == "+", 1, -1) * S$wid * 0.03
  ggplot() + { if (band) peak_band(S) } + S$tss_line +
    annotate("segment", x = tg$start, xend = tg$end, y = 0, yend = 0, color = GENE_COL, linewidth = 0.9) +
    geom_rect(data = gx[type != "cds"], aes(xmin = start, xmax = end, ymin = -0.15, ymax = 0.15), fill = GENE_COL, color = NA) +
    geom_rect(data = gx[type == "cds"], aes(xmin = start, xmax = end, ymin = -0.28, ymax = 0.28), fill = GENE_COL, color = NA) +
    annotate("segment", x = tg$tss, xend = tg$tss, y = 0, yend = 0.45, color = GENE_COL, linewidth = 1) +
    annotate("segment", x = tg$tss, xend = arr, y = 0.45, yend = 0.45, color = GENE_COL, linewidth = 1, arrow = arrow(length = unit(0.25, "cm"))) +
    annotate("text", x = (max(tg$start, xlim[1]) + min(tg$end, xlim[2])) / 2, y = 0.75, label = S$g, fontface = "bold.italic", size = 6.5) +
    coord_cartesian(xlim = xlim, ylim = c(-0.45, 1), expand = FALSE) +
    labs(y = if (strip) NULL else "Gene") + track_theme() + x_axis_off + (if (strip) strip_facet("Gene", S) else NULL)
}

## B. Peak track (只显示至少一个方法显著的 peak; example peak 标红; GWAS SNP = 紫色菱形)
track_peaks <- function(S, band = FALSE, strip = FALSE) {
  pk <- unique(S$lk, by = "peak_id"); gw <- S$gw
  p <- ggplot() + { if (band) peak_band(S) } + S$tss_line +
    coord_cartesian(xlim = S$xlim, ylim = c(0, 1), expand = FALSE) + labs(y = if (strip) NULL else "Peaks") + track_theme() + x_axis_off +
    (if (strip) strip_facet("Peaks", S) else NULL)
  if (nrow(pk)) {
    pk[, hw := pmax((pk_end - pk_start) / 2, S$wid * 0.003)]
    p <- p + geom_rect(data = pk, aes(xmin = pk_mid - hw, xmax = pk_mid + hw, ymin = 0.1, ymax = 0.9, fill = is_example), color = NA) +
      scale_fill_manual(values = PEAK_COL, breaks = c("TRUE", "FALSE"), labels = c("TRUE" = "Example peak", "FALSE" = "Other sig. peak"), name = "Peak")
  } else p <- p + annotate("text", x = mean(S$xlim), y = 0.5, label = "no significant peaks", color = "grey40", fontface = "italic", size = 5.5)
  if (nrow(gw)) p <- p + geom_point(data = unique(gw, by = "pos"), aes(pos, 0.5), shape = 23, fill = GWAS_COL, color = "white", size = 4.2, stroke = 0.6)
  p
}

## C. MR SNP track (一行 = 一个 example peak; 只用于 example 版式)
track_mrsnp <- function(S) {
  ex <- S$ex; sd <- S$sd; gw <- S$gw; xlim <- S$xlim
  n_ex <- max(1L, nrow(ex))
  n_by <- if (nrow(sd)) sd[, .N, by = row] else data.table(row = integer(), N = integer())
  ex[n_by, n_snp := i.N, on = "row"][is.na(n_snp), n_snp := 0L]
  p <- ggplot() + S$tss_line +
    geom_rect(data = ex, aes(xmin = pk_mid - hw, xmax = pk_mid + hw, ymin = -row - 0.42, ymax = -row + 0.42), fill = PEAK_COL[["TRUE"]], alpha = 0.3, color = NA)
  if (nrow(sd)) {
    p <- p + geom_segment(data = sd, aes(x = pos, xend = pos, y = -row - 0.3, yend = -row + 0.3), color = SNP_COL, linewidth = 0.9)
    if (sd[is_lead == TRUE, .N]) p <- p +
      geom_point(data = sd[is_lead == TRUE], aes(pos, -row + 0.3), shape = 25, fill = LEAD_FILL, color = LEAD_BORDER, size = 5, stroke = 0.7) +
      geom_text_repel(data = sd[is_lead == TRUE], aes(pos, -row + 0.3, label = variant_id), size = 4.6, fontface = "bold", color = LEAD_BORDER,
                      nudge_y = 0.3, direction = "x", min.segment.length = 0, segment.size = 0.4, seed = 1)
  } else p <- p + annotate("text", x = mean(xlim), y = -1, label = "no MR SNPs", color = "grey40", fontface = "italic", size = 5.5)
  # GWAS SNP: 紫色长竖线 + 菱形 (若同时是 MR instrument / SMR lead, 原标记仍保留); rsID 标签放在行下方 (trait / P 见 subtitle)
  if (nrow(gw)) p <- p +
    geom_segment(data = gw, aes(x = pos, xend = pos, y = -row - 0.42, yend = -row + 0.42), color = GWAS_COL, linewidth = 1.3) +
    geom_point(data = gw, aes(pos, -row), shape = 23, fill = GWAS_COL, color = "white", size = 4.2, stroke = 0.6) +
    geom_text_repel(data = gw, aes(pos, -row - 0.42, label = rsid), size = 4.2, fontface = "bold", color = GWAS_COL,
                    nudge_y = -0.25, min.segment.length = 0, segment.size = 0.4, box.padding = 0.4, seed = 1)
  p + scale_y_continuous(breaks = -ex$row, labels = sprintf("P%d (n=%d)%s", ex$row, ex$n_snp, ifelse(ex$has_gwas, " *GWAS", ""))) +
    coord_cartesian(xlim = xlim, ylim = c(-n_ex - if (nrow(gw)) 0.9 else 0.5, 0.1), expand = FALSE) +
    labs(y = "MR\nSNPs") + theme_classic(base_size = BASE_SIZE) + x_axis_off +
    theme(axis.text.y = element_text(size = BASE_SIZE - 4, color = "black"), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, face = "bold", size = BASE_SIZE))
}

## D. -log10 P 轨: [GWAS] / caQTL / eQTL (第一个 example peak; gwas 版式每个基因只有一个)
#  灰点 = 窗口内全部 SNP, 蓝点 = MR instrument, 金色三角 = SMR lead, 紫色菱形 = GWAS SNP
#  gt / gl: GWAS 轨数据 (variant_id, pos, mlp) 与 facet 名; NULL = 只画 caQTL / eQTL
track_assoc <- function(S, gt = NULL, gl = NULL, band = FALSE, x_axis = FALSE) {
  xlim <- S$xlim; pk1 <- S$ex$peak_id[1]
  lev <- c(gl, "caQTL", "eQTL")
  q <- S$qd[peak_id == pk1 & pos >= xlim[1] & pos <= xlim[2], .(layer, variant_id, pos, mlp)]
  if (!is.null(gt)) q <- rbind(data.table(layer = gl, gt[pos >= xlim[1] & pos <= xlim[2], .(variant_id, pos, mlp)]), q)
  s1 <- S$sd[peak_id == pk1]
  g1 <- S$gw[peak_id == pk1]
  qm <- q[layer %chin% c("caQTL", "eQTL") & variant_id %chin% s1$variant_id]
  ql <- q[layer %chin% c("caQTL", "eQTL") & variant_id %chin% s1[is_lead == TRUE, variant_id]]
  qg <- q[g1[, .(variant_id, rsid)], on = "variant_id", nomatch = 0L]
  if (!is.null(gt) && nrow(g1)) {   # GWAS SNP 不在 GWAS 轨 (liftOver / bim 缺失) 时, 用 Valid8 文件中的 P 补点
    miss <- g1[!variant_id %chin% qg[layer == gl, variant_id], .(layer = gl, variant_id, pos, mlp = -log10(pmax(P, P_FLOOR)), rsid)]
    if (nrow(miss)) { q <- rbind(q, miss[, !"rsid"]); qg <- rbind(qg, miss) }
  }
  q[, layer := factor(layer, levels = lev)]
  if (nrow(qg)) qg[, layer := factor(layer, levels = lev)]
  if (nrow(qm)) qm[, layer := factor(layer, levels = lev)]
  if (nrow(ql)) ql[, layer := factor(layer, levels = lev)]

  # 每个 layer 各自 y 轴 (free_y), 用 geom_blank 固定 0 与顶部留白; GWAS 轨至少显示到 5e-8 线
  y_rng <- data.table(layer = factor(lev, levels = lev))
  y_rng[q[, .(m = max(mlp)), by = layer], m := i.m, on = "layer"][is.na(m), m := 1]
  if (!is.null(gl)) y_rng[layer == gl, m := pmax(m, -log10(GWS_P))]
  y_blank <- rbind(y_rng[, .(layer, y = 0)], y_rng[, .(layer, y = pmax(1, m) * 1.22)])[, x := S$tg$tss]

  p <- ggplot() + { if (band) peak_band(S) } + S$tss_line +
    geom_blank(data = y_blank, aes(x, y)) +
    { if (!is.null(gl)) geom_hline(data = data.table(layer = factor(gl, levels = lev), y = -log10(GWS_P)), aes(yintercept = y),
                                   linetype = "dashed", color = "grey40", linewidth = 0.4) } +
    { if (!band && nrow(g1)) geom_vline(xintercept = unique(g1$pos), color = GWAS_COL, linetype = "dashed", linewidth = 0.5) } +
    geom_point(data = q, aes(pos, mlp), color = "grey62", size = 1.5, alpha = 0.45, shape = 16) +
    geom_point(data = qm, aes(pos, mlp), fill = SNP_COL, color = "white", size = 3, stroke = 0.4, shape = 21) +
    geom_point(data = ql, aes(pos, mlp), shape = 25, fill = LEAD_FILL, color = LEAD_BORDER, size = 5, stroke = 0.8) +
    geom_point(data = qg, aes(pos, mlp), shape = 23, fill = GWAS_COL, color = "white", size = 5.5, stroke = 0.7) +
    geom_text_repel(data = qg, aes(pos, mlp, label = rsid), size = 6.5, fontface = "bold", color = GWAS_COL,
                    min.segment.length = 0, segment.size = 0.5, box.padding = 0.6, point.padding = 0.4, seed = 1) +
    facet_grid(layer ~ ., switch = "y", scales = "free_y") +
    scale_y_continuous(breaks = function(l) pretty(c(0, l[2]), n = 3)) + mb_scale +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    labs(y = NULL, x = sprintf("%s position (Mb)", S$tg$chr)) + theme_classic(base_size = BASE_SIZE) + strip_theme +
    theme(axis.text.y = element_text(size = BASE_SIZE - 6, color = "black")) + (if (x_axis) x_axis_on else x_axis_off)
  if (!nrow(q)) p <- p + annotate("text", x = mean(xlim), y = 0.5, label = "no QTL data", color = "grey40", fontface = "italic", size = 5)
  p
}

## E. Link tracks (一个方法一行, facet 保证 6 行都在)
link_counts <- function(lk) {
  lev <- METHODS$method
  cnt <- lk[, .N, by = method][data.table(method = lev), on = "method"][is.na(N), N := 0L]
  cnt[, method := factor(method, levels = lev)]
  cnt[]
}

track_links <- function(S, score_max, band = FALSE, show_n = TRUE) {
  lk <- copy(S$lk); tg <- S$tg; xlim <- S$xlim
  lev <- METHODS$method; cnt <- link_counts(lk)
  s_lo <- -log10(FDR_CUT)
  if (nrow(lk)) {
    lk[, `:=`(x1 = pk_mid, x2 = tg$tss, link_id = .I, score_plot = pmin(pmax(score, s_lo), score_max))]   # 手动 squish (scale_linewidth 不支持 oob)
    lk[, h := pmax(0.2, abs(x2 - x1) / max(1, max(abs(x2 - x1))))]
    tt <- seq(0, 1, length.out = 100)
    arcs <- lk[, .(x = x1 + (x2 - x1) * tt, y = -h * sin(pi * tt)), by = .(link_id, method, direction, score_plot)]
    setorder(arcs, score_plot)                                    # 强的 link 画在上层
  } else arcs <- data.table(link_id = integer(), method = character(), direction = character(), score_plot = numeric(), x = numeric(), y = numeric())
  arcs[, `:=`(method = factor(method, levels = lev), direction = factor(direction, levels = names(DIR_COL)))]

  ggplot() + { if (band) peak_band(S) } + S$tss_line +
    geom_blank(data = data.table(method = factor(lev, levels = lev), x = tg$tss, y = 0), aes(x, y)) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.5) +
    geom_path(data = arcs, aes(x, y, group = link_id, color = direction, linewidth = score_plot), alpha = 0.95, lineend = "round") +
    { if (show_n) geom_text(data = cnt, aes(x = xlim[1] + S$wid * 0.006, y = -0.88, label = paste0("n = ", N)), hjust = 0, size = 5, color = "grey20") } +
    facet_grid(method ~ ., switch = "y", drop = FALSE) +
    scale_color_manual(values = DIR_COL, drop = FALSE, name = "Direction") +
    scale_linewidth_continuous(range = LW_RANGE, limits = c(s_lo, score_max), name = expression(-log[10]~FDR)) + mb_scale +
    coord_cartesian(xlim = xlim, ylim = c(-1.08, 0.08), expand = FALSE) +
    guides(color = guide_legend(order = 1, override.aes = list(linewidth = 2.5)), linewidth = guide_legend(order = 2)) +
    labs(x = sprintf("%s position (Mb)", tg$chr), y = NULL) + theme_classic(base_size = BASE_SIZE) + strip_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) + x_axis_on
}

title_theme <- theme(plot.title = element_text(face = "bold", size = BASE_SIZE + 6),
                     plot.subtitle = element_text(size = BASE_SIZE - 3, color = "grey20", lineheight = 1.2),
                     plot.caption = element_text(size = BASE_SIZE - 5, color = "grey35", hjust = 0))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 版式 1 (example): 竖排 Gene / Peaks / MR SNPs / caQTL-eQTL / links

plot_locus <- function(g, lk_all, ex_pairs, snp_dt, qtl_dt, gw_dt, ct, score_max) {
  S <- locus_setup(g, lk_all, ex_pairs, snp_dt, qtl_dt, gw_dt); if (is.null(S)) return(NULL)
  pG <- track_gene(S); pP <- track_peaks(S); pS <- track_mrsnp(S); pQ <- track_assoc(S); pL <- track_links(S, score_max)
  cnt <- link_counts(S$lk); ex <- S$ex; gw <- S$gw; tg <- S$tg

  sub1 <- sprintf("Significant links (FDR < %.2f):  %s", FDR_CUT, paste(sprintf("%s %d", cnt$method, cnt$N), collapse = "  |  "))
  sub2 <- sprintf("Example peak(s): %s", paste(sprintf("P%d = %s", ex$row, ex$peak_id), collapse = ";  "))
  sub3 <- if (nrow(gw)) sprintf("GWAS SNP(s) in peak: %s", paste(unique(gw$label), collapse = ";  ")) else "GWAS SNP(s) in peak: none"
  cap  <- sprintf(paste0("TSS %s:%s (%s strand)   |   MR SNPs: ticks = instruments used in MR (after LD filtering), gold triangle = SMR lead SNP, shaded = example peak\n",
                         "Purple diamond = GWAS SNP (Valid8)   |   caQTL (P1 peak) / eQTL (%s): -log10 P, grey = all SNPs, blue = MR instruments"),
                  tg$chr, format(tg$tss, big.mark = ","), tg$strand, g)
  (pG / pP / pS / pQ / pL) + plot_layout(heights = c(0.9, 0.45, 0.45 + 0.45 * max(1L, nrow(ex)) + if (nrow(gw)) 0.4 else 0, 3, 7), guides = "collect") +
    plot_annotation(title = sprintf("%s  -  %s", g, ct), subtitle = paste(sub1, sub2, sub3, sep = "\n"), caption = cap, theme = title_theme) &
    theme(legend.title = element_text(face = "bold", size = BASE_SIZE - 1), legend.text = element_text(size = BASE_SIZE - 3))
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 版式 2 (gwas): 一页 = 一个 link x 一个 trait; 左 = GWAS / caQTL / eQTL (-log10 P), 右 = Gene / Peaks / links; 两边 x 窗口相同, 无图例

plot_locus_gwas <- function(g, tr, lk_all, ex_pairs, snp_dt, qtl_dt, gw_tr, ct, score_max, res) {
  S <- locus_setup(g, lk_all, ex_pairs, snp_dt, qtl_dt, gw_tr); if (is.null(S)) return(NULL)
  gl <- "GWAS"
  gt <- get_gwas_track(tr, S$tg$chr)[, .(variant_id = paste0(S$tg$chr, "_", pos), pos, mlp)]

  pQ <- track_assoc(S, gt, gl, band = TRUE, x_axis = TRUE)
  pG <- track_gene(S, band = TRUE, strip = TRUE); pP <- track_peaks(S, band = TRUE, strip = TRUE); pL <- track_links(S, score_max, band = TRUE, show_n = FALSE)
  right <- (pG / pP / pL) + plot_layout(heights = c(0.9, 0.45, 7))

  pk <- S$ex$peak_id[1]
  pv <- get_pair_p(res, g, pk)
  sub <- paste("P-value:  ", paste(sprintf("%s %s", METHODS$method, ifelse(pv$found, sprintf("%.1e", pv$p), "1 (not found)")), collapse = "   |   "))
  (pQ | right) + plot_layout(widths = c(1, 1)) +
    plot_annotation(title = sprintf("%s - %s - %s - %s", ct, g, pk, tr), subtitle = sub, theme = title_theme) &
    theme(legend.position = "none")
}

# 该 (gene, peak) 6 个方法的 P (formal 结果, 已过 sigma2); NA / 找不到 -> 1
get_pair_p <- function(res, g, pk) {
  r <- res[gene_id == g & peak_id == pk]
  p <- vapply(METHODS$p_col, function(cn) if (nrow(r) && cn %in% names(r) && !is.na(r[[cn]][1])) as.numeric(r[[cn]][1]) else NA_real_, numeric(1))
  data.table(method = METHODS$method, p = fifelse(is.na(p), 1, p), found = !is.na(p))
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 运行: 三个细胞类型 -> 同一个 PDF

SIGMA2_MAX <- 10   # 与其他图统一: 只保留 ivw_sigma2 <= 10 (NA 丢掉)
RES <- lapply(setNames(names(CELLTYPES), names(CELLTYPES)), function(ct) fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz"))))

# 线宽标尺在三个细胞类型间共用
SCORE_MAX <- min(SCORE_CAP, max(c(2, unlist(lapply(names(CELLTYPES), function(ct) get_sig_links(RES[[ct]], unique(CELLTYPES[[ct]]$gene_id))$score)))))
msg(sprintf("线宽标尺: -log10(FDR) in [%.2f, %.2f]", -log10(FDR_CUT), SCORE_MAX))

if (PAIR_SET == "gwas") pdf(PDF_FILE, width = PAGE_W_G, height = PAGE_H_G, useDingbats = FALSE) else pdf(PDF_FILE, width = PAGE_W, height = PAGE_H, useDingbats = FALSE)
for (ct in names(CELLTYPES)) {
  msg("==================== ", ct, " ====================")
  pairs <- CELLTYPES[[ct]]; res <- RES[[ct]]; genes <- unique(pairs$gene_id)

  if (!nrow(pairs)) { msg("  无 pairs, 跳过"); next }
  cs <- collect_mr_snps(ct, pairs, res); snp_dt <- cs$snps; qtl_dt <- cs$qtl
  if (!nrow(snp_dt)) snp_dt <- data.table(gene_id = character(), peak_id = character(), variant_id = character(), pos = numeric(), is_lead = logical())
  gw_dt <- GWAS[cell_type == ct]
  gw_chk <- gw_dt[pairs, on = .(gene_id, peak_id), nomatch = 0L]
  msg(sprintf("[%s] GWAS SNP 在 QTL 中找到: caQTL %d / %d, eQTL %d / %d", ct,
              gw_chk[qtl_dt[layer == "caQTL"], on = .(gene_id, peak_id, variant_id), nomatch = 0L, .N], nrow(gw_chk),
              gw_chk[qtl_dt[layer == "eQTL"],  on = .(gene_id, peak_id, variant_id), nomatch = 0L, .N], nrow(gw_chk)))

  pm <- pairs[, get_pair_p(res, gene_id, peak_id)[found == FALSE, .N], by = .(gene_id, peak_id)][, sum(V1)]
  msg(sprintf("[%s] 方法 P 值缺失 (按 1 处理) 的 pair x method: %d", ct, pm))
  lk_all <- get_sig_links(res, genes)
  lk_all[, is_example := FALSE][unique(pairs[, .(gene_id, peak_id)]), is_example := TRUE, on = .(gene_id, peak_id)]
  msg(sprintf("[%s] %d genes, %d 显著 links", ct, length(genes), nrow(lk_all)))
  if (nrow(lk_all)) print(dcast(lk_all[, .(gene_id, method = factor(method, levels = METHODS$method))], gene_id ~ method, fun.aggregate = length, drop = FALSE))

  for (g in genes) {
    if (PAIR_SET == "gwas") {   # 一页 = 一个 link x 一个 trait (按该 trait 最小 GWAS P 排序)
      gtr <- GWAS_T[cell_type == ct & gene_id == g]
      for (tr in gtr[, .(P = min(P)), by = trait][order(P), trait]) {
        msg("  plotting ", g, " | ", tr)
        p <- tryCatch(plot_locus_gwas(g, tr, lk_all, pairs, snp_dt, qtl_dt, gtr[trait == tr], ct, SCORE_MAX, res),
                      error = function(e) { msg("    失败: ", conditionMessage(e)); NULL })
        if (!is.null(p)) print(p)
      }
    } else {
      msg("  plotting ", g)
      p <- tryCatch(plot_locus(g, lk_all, pairs, snp_dt, qtl_dt, gw_dt, ct, SCORE_MAX), error = function(e) { msg("    失败: ", conditionMessage(e)); NULL })
      if (!is.null(p)) print(p)
    }
  }
  invisible(gc())
}
dev.off()

unlink(LD_DIR, recursive = TRUE)
msg("Saved:\n  ", PDF_FILE)