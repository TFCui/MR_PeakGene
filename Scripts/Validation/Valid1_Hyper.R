# 2026-10-09
# Valid1_Hyper: 金标准准备 + 每个 gene 内的超几何富集 (流程见 Scripts/Pipeline.md 3.1)
# 移植自 Scripts_reference/Hyper1and2_v1.sh, 逻辑不变; 改读新的 MR 结果表, 不画图
# 一个 job 依次处理 3 个 CT; 输出:
#   Results/Validation/Valid1_Hyper/<CT>_hyper_long.tsv
#   Results/Validation/Valid1_Hyper/<CT>_gold_standard.rds   (给 Valid4_AUC 用)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(readxl)
  library(EnsDb.Hsapiens.v86)
})
options(width = 300)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 配置

BASE_DIR  <- "/scratch/user/tfcui2025/CIMA"
PROJ      <- file.path(BASE_DIR, "Claude_project")
MR_DIR    <- file.path(PROJ, "Results/MR")
OUT_DIR   <- file.path(PROJ, "Results/Validation/Valid1_Hyper")
GP_DIR    <- file.path(BASE_DIR, "Data/OMIX009786/gene_peak")
DIR_VALID <- file.path(BASE_DIR, "Data/PBMC_validation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELLS       <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
FDR_CUT     <- 0.1     # 6 个方法共用
HICHIP_QCUT <- 0.10    # FitHiChIP FDR
PCHIC_SCORE <- 5

cell_cfg <- list(Bn_TCL1A    = list(hichip = "Bcell", pchic = "nB"),
                 CD4_Tn_CCR7 = list(hichip = "CD4",   pchic = "nCD4"),
                 CD8_Tn_CCR7 = list(hichip = "CD8",   pchic = "nCD8"))

method_fdr <- c(IVW = "ivw_fdr", Egger = "egger_fdr", SMR = "smr_fdr", CCT = "cct_fdr", WNN = "wnn_fdr", Pseudobulk = "pseudo_fdr")

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 金标准原始数据 (只读一次)

load(file.path(DIR_VALID, "gtex.rda"))                            # 对象 gtex
enhanceratlas_grl <- readRDS(file.path(DIR_VALID, "enhanceratlas.combined.grl.rds"))
fantom_grl        <- readRDS(file.path(DIR_VALID, "fantom.grl.rds"))
four_d_genome_grl <- readRDS(file.path(DIR_VALID, "four.d.genome.grl.rds"))
transcripts.gr    <- transcripts(EnsDb.Hsapiens.v86, columns = c("gene_name", "gene_id"))
# 等价于 seqlevelsStyle(transcripts.gr) <- "UCSC" (1 -> chr1, MT -> chrM); 那个写法要联网取 UCSC chromInfo, 计算节点没有网络
sl <- seqlevels(transcripts.gr)
genome(transcripts.gr) <- NA
transcripts.gr <- renameSeqlevels(transcripts.gr, setNames(ifelse(sl == "MT", "chrM", paste0("chr", sl)), sl))

hichip_cols <- c("Gene_Chr", "TSSs", "GeneID", "GeneName", "GeneType",
                 "Promoter_Chr", "Promoter_Start", "Promoter_End",
                 "Interact_Chr", "Interact_Start", "Interact_End", "Annotation",
                 "Mono_Count", "Mono_Q", "NK_Count", "NK_Q", "Bcell_Count", "Bcell_Q",
                 "CD4_Count", "CD4_Q", "CD8_Count", "CD8_Q",
                 "Mono_TPM", "NK_TPM", "Bcell_TPM", "CD4_TPM", "CD8_TPM")
loops_3 <- read_excel(file.path(DIR_VALID, "Validation_Data_scMutiMap/41588_2020_745_MOESM3_ESM.xlsx"),
                      sheet = "Supplementary Table 3", skip = 3, col_names = hichip_cols)
setDT(loops_3)
loops_3 <- loops_3[-1]                                            # 去掉残留表头行

pchic_promoter <- fread(file.path(DIR_VALID, "Validation_Data_scMutiMap/PCHiC/PCHiC_peak_matrix_cutoff5.txt"))

chain <- import.chain(file.path(DIR_VALID, "hg19ToHg38.over.chain"))
lift_to_hg38 <- function(gr19, label) {
  n_all  <- length(gr19)
  lifted <- liftOver(gr19, chain)
  idx  <- which(lengths(lifted) > 0L)
  ok   <- lengths(unique(seqnames(lifted[idx]))) == 1L
  keep <- idx[ok]
  out  <- unlist(range(lifted[keep]), use.names = FALSE)         # 多段取首尾 range
  stopifnot(length(out) == length(keep))
  mcols(out) <- mcols(gr19)[keep, , drop = FALSE]
  cat(sprintf("    [%s] liftover: %s / %s 保留 (%.1f%%); 未映射 %s, 跨染色体 %s; 中位宽度 %s -> %s bp\n",
              label, format(length(out), big.mark = ","), format(n_all, big.mark = ","), 100 * length(out) / n_all,
              format(n_all - length(idx), big.mark = ","), format(sum(!ok), big.mark = ","),
              median(width(gr19)), median(width(out))))
  out
}

phyper_test <- function(N, K, n_sig, k, N_qtl, K_qtl) {
  pval     <- if (K == 0L || n_sig == 0L || N == K) 1 else phyper(k - 1L, m = K, n = N - K, k = n_sig, lower.tail = FALSE)
  pval_qtl <- if (K_qtl == 0L || n_sig == 0L || N_qtl == K_qtl) 1 else phyper(k - 1L, m = K_qtl, n = N_qtl - K_qtl, k = n_sig, lower.tail = FALSE)
  list(pval = pval, pval_qtl = pval_qtl, N = N, K = K, n_sig = n_sig, k = k, N_qtl = N_qtl, K_qtl = K_qtl)
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 逐个 CT

for (CT in CELLS) {
  t0  <- Sys.time()
  cfg <- cell_cfg[[CT]]
  cat(sprintf("\n==================== %s (HiChIP %s_Q, PCHiC %s) ====================\n", CT, cfg$hichip, cfg$pchic))

  res <- fread(file.path(MR_DIR, paste0(CT, "_AllPairs_results.tsv.gz")))
  gp  <- as.data.table(readRDS(file.path(GP_DIR, paste0("All_GenePeak_", CT, ".rds"))))
  gp[, peak_id_std := sub("-", ":", peak_id)]
  bg_by_gene <- lapply(split(gp$peak_id_std, gp$gene_id), unique)
  cat(sprintf("[bg] %s pairs, %s genes, %s unique peaks; 每基因 peak 数中位 %.0f\n",
              format(nrow(gp), big.mark = ","), format(length(bg_by_gene), big.mark = ","),
              format(uniqueN(gp$peak_id_std), big.mark = ","), median(lengths(bg_by_gene))))
  rm(gp); invisible(gc())

  # 六个方法各自的显著 peak
  SIG <- lapply(method_fdr, function(fc) {
    sr <- res[get(fc) < FDR_CUT]
    if (nrow(sr) == 0) return(list())
    lapply(split(sr$peak_id, sr$gene_id), unique)
  })
  cat(sprintf("[sig] FDR < %.2f, MR pairs = %s\n", FDR_CUT, format(nrow(res), big.mark = ",")))
  for (m in names(SIG)) cat(sprintf("  %-10s sig = %5s pairs, %4s genes\n", m,
                                    format(sum(lengths(SIG[[m]])), big.mark = ","), format(length(SIG[[m]]), big.mark = ",")))
  gene_universe <- intersect(unique(res$gene_id), names(bg_by_gene))

  ## ---- 步骤 1: 金标准 ----
  DB <- list()
  DB[["enhancerAtlas"]] <- enhanceratlas_grl[names(enhanceratlas_grl) %in% gene_universe]
  DB[["FANTOM5"]]       <- fantom_grl[names(fantom_grl) %in% gene_universe]
  DB[["4DGenome"]]      <- four_d_genome_grl[names(four_d_genome_grl) %in% gene_universe]

  gtex_sub <- gtex[names(gtex) %in% gene_universe]
  DB[["GTEx"]] <- lapply(gtex_sub, function(gt) {
    gt <- gt[!is.na(gt$meta.pval) & gt$meta.pval < 0.05]
    if (length(gt) == 0) return(gt)
    gt[countOverlaps(gt, transcripts.gr) == 0]
  })
  DB[["GTEx"]] <- DB[["GTEx"]][lengths(DB[["GTEx"]]) > 0]

  q_col <- paste0(cfg$hichip, "_Q")
  stopifnot(q_col %in% names(loops_3), cfg$pchic %in% names(pchic_promoter))
  hi <- loops_3[, .(Gene  = GeneName,
                    chr   = Interact_Chr,
                    start = suppressWarnings(as.integer(Interact_Start)),
                    end   = suppressWarnings(as.integer(Interact_End)),
                    q     = suppressWarnings(as.numeric(get(q_col))))]
  hi <- unique(hi[!is.na(q) & q < HICHIP_QCUT & !is.na(start) & !is.na(end) & Gene %in% gene_universe])
  cat(sprintf("    [HiChIP] 过滤后 %s 条 loop, %s 基因\n", format(nrow(hi), big.mark = ","), format(uniqueN(hi$Gene), big.mark = ",")))
  hi_gr38 <- lift_to_hg38(GRanges(hi$chr, IRanges(hi$start, hi$end), Gene = hi$Gene), "HiChIP")
  DB[["HiChIP"]] <- split(granges(hi_gr38), hi_gr38$Gene)

  hc <- pchic_promoter[, c("baitName", "oeChr", "oeStart", "oeEnd", cfg$pchic), with = FALSE]
  setnames(hc, cfg$pchic, "score")
  hc <- hc[!is.na(score) & score >= PCHIC_SCORE, .(Gene = baitName, chr = paste0("chr", oeChr), start = oeStart, end = oeEnd)]
  hc <- hc[, .(Gene = unlist(strsplit(Gene, ";", fixed = TRUE))), by = .(chr, start, end)]
  hc <- unique(hc[Gene != "." & Gene %in% gene_universe])
  cat(sprintf("    [PCHiC] 过滤后 %s 条 link, %s 基因\n", format(nrow(hc), big.mark = ","), format(uniqueN(hc$Gene), big.mark = ",")))
  hc_gr38 <- lift_to_hg38(GRanges(hc$chr, IRanges(hc$start, hc$end), Gene = hc$Gene), "PCHiC")
  DB[["PCHiC"]] <- split(granges(hc_gr38), hc_gr38$Gene)

  cat(sprintf("[db] gene_universe = %s\n", format(length(gene_universe), big.mark = ",")))
  for (d in names(DB)) {
    obj <- DB[[d]]
    all_ranges <- if (is(obj, "GRangesList")) unlist(obj, use.names = FALSE) else unlist(GRangesList(lapply(obj, granges)), use.names = FALSE)
    cat(sprintf("  %-14s %5s genes (%.1f%% of universe), 中位区间数 %3.0f, 中位宽度 %6.0f bp\n",
                d, format(length(obj), big.mark = ","), 100 * length(obj) / length(gene_universe),
                median(lengths(obj)), median(width(all_ranges))))
  }
  saveRDS(DB, file.path(OUT_DIR, paste0(CT, "_gold_standard.rds")))

  ## ---- 步骤 2: 超几何检验 ----
  ALL_PEAK_IDS <- unique(unlist(bg_by_gene, use.names = FALSE))
  .mat <- do.call(rbind, strsplit(ALL_PEAK_IDS, "[:-]"))
  ALL_PEAK_GR <- GRanges(.mat[, 1], IRanges(start = as.integer(.mat[, 2]), end = as.integer(.mat[, 3])))
  rm(.mat); invisible(gc())
  gr_for <- function(peak_ids) ALL_PEAK_GR[chmatch(peak_ids, ALL_PEAK_IDS)]

  mr_by_gene <- lapply(split(res$peak_id, res$gene_id), unique)
  chk <- vapply(names(mr_by_gene), function(g) {
    bg <- bg_by_gene[[g]]
    if (is.null(bg)) return(c(0L, length(mr_by_gene[[g]])))
    c(length(intersect(mr_by_gene[[g]], bg)), length(mr_by_gene[[g]]))
  }, integer(2))
  cat(sprintf("[qtl] MR peaks 落在自身背景内: %s / %s (%.1f%%)\n", format(sum(chk[1, ]), big.mark = ","),
              format(sum(chk[2, ]), big.mark = ","), 100 * sum(chk[1, ]) / sum(chk[2, ])))

  db_names   <- names(DB)
  genes_test <- intersect(names(bg_by_gene), unique(unlist(lapply(SIG, names), use.names = FALSE)))
  cat(sprintf("[run] 待检验基因: %s\n", format(length(genes_test), big.mark = ",")))

  out <- vector("list", length(genes_test))
  for (i in seq_along(genes_test)) {
    g      <- genes_test[i]
    bg_pk  <- bg_by_gene[[g]]
    bg_gr  <- gr_for(bg_pk)
    N      <- length(bg_pk)
    is_qtl <- bg_pk %chin% mr_by_gene[[g]]
    N_qtl  <- sum(is_qtl)

    hit <- list()
    for (d in db_names) {
      if (!g %in% names(DB[[d]])) next
      hit[[d]] <- suppressWarnings(countOverlaps(bg_gr, DB[[d]][[g]])) > 0L
    }
    if (length(hit) == 0L) next

    rows <- list()
    for (m in names(SIG)) {
      sig_g <- SIG[[m]][[g]]
      if (is.null(sig_g)) next
      sig_pk <- intersect(sig_g, bg_pk)
      if (length(sig_pk) == 0L) next
      pos <- chmatch(sig_pk, bg_pk)
      for (d in names(hit)) {
        h <- hit[[d]]
        r <- phyper_test(N = N, K = sum(h), n_sig = length(sig_pk), k = sum(h[pos]), N_qtl = N_qtl, K_qtl = sum(h & is_qtl))
        rows[[length(rows) + 1L]] <- c(list(gene_id = g, method = m, database = d), r)
      }
    }
    if (length(rows)) out[[i]] <- rbindlist(rows)
  }

  hyper_long <- rbindlist(out)
  hyper_long[, database := factor(database, levels = db_names)]
  hyper_long[, method   := factor(method,   levels = names(SIG))]
  setorder(hyper_long, method, database, gene_id)
  f_out <- file.path(OUT_DIR, sprintf("%s_hyper_long.tsv", CT))
  fwrite(hyper_long, f_out, sep = "\t")

  cat(sprintf("[done] %s 行 -> %s (%.1f min)\n", format(nrow(hyper_long), big.mark = ","), f_out,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(hyper_long[, .(n_gene = .N, n_p05 = sum(pval < 0.05), n_p05_qtl = sum(pval_qtl < 0.05)), keyby = .(method, database)])
}

cat("\n全部完成\n")
