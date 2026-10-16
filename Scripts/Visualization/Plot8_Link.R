# ==============================================================================
# Fig8: MR-only 强 link 与 GWAS SNP 的位点图 (Figures_20261019/Fig8_Link_GWAS.pdf)
# 2026-10-13 重画: 一页 = 一个 (cell type, gene, peak); 只画研究的这一条 peak-gene link, 不画其他 peak
#   pairs / GWAS SNP = Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv (一行 = pair x trait x GWAS SNP; 按文件顺序给 pair 编号)
#   轨道 (上 -> 下): Gene / SNPs / Peak / Links (Causal SMR, IVW, Egger, Marginal Multiomics, Pseudobulk 各一行; 不画 Causal Joint)
#   X 轴分段缩放 (牺牲比例, 无灰底, 靠刻度标位置): Peak 区 (peak 中点 ± ZOOM_HALF) 和 Gene 区 (整个基因体, 固定宽度) 各自放大,
#     中间 gap 和两侧 flank 大幅压缩; 刻度只在 Peak / Gene 区, 标签 = 真实位置 (Mb); gap 处画 // 断裂标记; 区域重叠时合并
#   Links: FDR < 0.05 = 实线 (颜色 = 方向, 线宽 ∝ -log10 FDR, 所有页面共用标尺); 不显著 = 灰色细虚线; 行标签写 FDR; 无图例
#   SNPs: 蓝色竖线 = MR instrument (按 MR 流程重跑 QTL 对齐 -> 去重 -> PLINK LD, 与 ivw_n_snp 对账; 只挑 ~N_SNP_TARGET 个, 离 peak 近优先, 纯为美观);
#     紫色菱形 = GWAS SNP (P < 5e-8, 全部画)
# 旧版 (左栏 GWAS / caQTL / eQTL + 全部显著 link, 一页 = link x trait) 见 Results/Figures/Figures_20261012/Fig8_Link_GWAS.pdf
# 运行: sbatch Scripts/Visualization/Plot8_Link_slurm.sh   (PLINK / 临时文件 = tmp/plink_fig8/)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
  library(GenomicRanges); library(EnsDb.Hsapiens.v86); library(Signac); library(grid)
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
FIG_DIR  <- file.path(PROJ, "Results/Figures/Figures_20261019")
PAIR_FILE <- file.path(PROJ, "Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv")
PDF_FILE  <- file.path(FIG_DIR, "Fig8_Link_GWAS.pdf")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LD_DIR,  recursive = TRUE, showWarnings = FALSE)

eqtl_file        <- function(ct) file.path(OMIX, "cis_eQTL_Summary_Statics/cis_eQTL", paste0(ct, ".txt.gz"))
caqtl_clump_file <- function(ct) file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL/LD_clump", paste0(ct, "_Clumped_r2_0.99_kb_20000.txt.gz"))

# 作图参数
FDR_CUT    <- 0.05
SCORE_CAP  <- 20                      # -log10(FDR) 上限, 所有页面共用同一线宽标尺
ZOOM_HALF  <- 1000                    # Peak 区 = peak 中点 ± 1 kb
W_PEAK <- 0.30; W_GENE <- 0.10        # x 轴宽度权重 (归一化): Peak 区 / Gene 区 (整个基因体)
W_GAP  <- 0.20; GAP_FULL_BP <- 20000  # Peak 区与 Gene 区之间的 gap (短于 GAP_FULL_BP 时按比例缩小, 避免被拉伸)
W_FLANK <- 0.12; FLANK_BP <- 20000    # 两侧 flank: 宽度权重 W_FLANK, 至少覆盖 FLANK_BP; MR SNP 不够 N_SNP_TARGET 时向外延伸 (最多 FLANK_MAX)
FLANK_MAX  <- 300000
TICK_SEP   <- 0.14                    # 刻度最小间距 (x 轴比例)
N_SNP_TARGET <- 8; SNP_MIN_SEP <- 0.025 # SNP 轨: 每页目标 SNP 数 (含 GWAS SNP) / 最小间距 (x 轴比例)
LW_RANGE   <- c(1, 3.2); LW_NS <- 0.6 # 显著 link 线宽范围; 不显著 link 线宽
BASE_SIZE  <- 13
PAGE_W     <- 11; PAGE_H <- 6.5
FIG_FRAC   <- 0.70                    # 图 (含黑色边框) 占页宽比例, 靠左放置, 右侧留白

# MR 固定参数 (与 MR 结果一致)
MIN_SNP <- 3; DEDUP_TOL <- 1e-3; KEEP_ALLELE <- TRUE; P_FLOOR <- 1e-300

# 方法名 / 顺序 / 配色: 全文统一 (见 00_common.R); 不画 Causal Joint
METHODS <- data.table(method  = c("Causal SMR", "Causal IVW", "Causal Egger", "Marginal Multiomics", "Marginal Pseudobulk"),
                      fdr_col = c("smr_fdr", "ivw_fdr", "egger_fdr", "wnn_fdr", "pseudo_fdr"),
                      eff_col = c("smr_b", "ivw_beta", "egger_beta", "wnn_rho", "pseudo_rho"),
                      col     = c("#F0A14A", "#E64B35", "#3C5488", "#4DBBD5", "#00A087"))
DIR_COL   <- c("Positive" = "#D7191C", "Negative" = "#2C7BB6")
NS_COL    <- "grey65"
GENE_COL  <- "#4E79A7"; PEAK_FILL <- "#E1B13F"
SNP_COL   <- "#3C5488"; GWAS_COL <- "#7E1F86"

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Pairs 与 GWAS SNP

EX <- fread(PAIR_FILE)[, ord := rleid(cell_type, gene_id, peak_id)]
PAIRS <- unique(EX[, .(ord, cell_type, gene_id, peak_id, mr_methods_sig)], by = "ord")[order(ord)]
GW <- EX[, .(pos = as.numeric(pos_hg38[1]), P = min(gwas_p), traits = paste(sprintf("%s P=%.1e", trait, gwas_p)[order(gwas_p)], collapse = "; ")),
         by = .(cell_type, gene_id, peak_id, rsid)]
msg(sprintf("Fig8 pairs: %d | GWAS SNP-pair 行: %d", nrow(PAIRS), nrow(GW)))

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
rm(gene_gr, exon_gr); invisible(gc())
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

# 一个细胞类型: 读 QTL -> 每个 pair 取 MR SNP -> 与 MR 结果对账 (SNP 数 vs ivw_n_snp, lead SNP vs smr_variant_id)
collect_mr_snps <- function(ct, pairs, res) {
  pairs <- unique(pairs[, .(gene_id, peak_id)])
  eq <- read_filtered(eqtl_file(ct), unique(pairs$gene_id))[, .(phenotype_id, variant_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se, pval_nominal)]
  setkey(eq, phenotype_id, variant_id)
  ca <- fread(caqtl_clump_file(ct))[unique(pairs[, .(phenotype_id = peak_id, gene_id)]), on = .(phenotype_id, gene_id), nomatch = 0L,
                                    .(phenotype_id, variant_id, gene_id, `A1(ALT/effect allele)`, `A2(REF)`, slope, slope_se, pval_nominal)]
  setkey(ca, phenotype_id, gene_id)
  snps <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    g <- pairs$gene_id[i]; pk <- pairs$peak_id[i]
    o <- tryCatch(get_mr_snps(g, pk, ca, eq), error = function(e) list(dat = NULL, status = paste("error:", conditionMessage(e))))
    msg(sprintf("  [%s] %s__%s: %s, %d SNPs", ct, pk, g, o$status, if (is.null(o$dat)) 0L else nrow(o$dat)))
    if (is.null(o$dat) || !nrow(o$dat)) return(NULL)
    o$dat[, .(gene_id = g, peak_id = pk, variant_id, pos = as.numeric(sub("^[^_]*_([0-9]+).*$", "\\1", variant_id)))]
  }))
  if (!nrow(snps)) return(data.table(gene_id = character(), peak_id = character(), variant_id = character(), pos = numeric(), is_lead = logical()))
  fm <- res[pairs, on = .(gene_id, peak_id), nomatch = NULL, .(gene_id, peak_id, ivw_n_snp, smr_variant_id)]
  snps[fm, smr_lead := i.smr_variant_id, on = .(gene_id, peak_id)]
  snps[, is_lead := !is.na(smr_lead) & variant_id == smr_lead][, smr_lead := NULL]
  chk <- snps[, .(n_snp_rerun = .N, lead_found = any(is_lead)), by = .(gene_id, peak_id)][fm, on = .(gene_id, peak_id)]
  msg(sprintf("[%s] 对账: SNP 数一致 %d / %d, lead SNP 找到 %d / %d", ct,
              chk[, sum(n_snp_rerun == ivw_n_snp, na.rm = TRUE)], nrow(chk), chk[, sum(lead_found, na.rm = TRUE)], nrow(chk)))
  if (chk[!(n_snp_rerun %in% ivw_n_snp) | !lead_found, .N]) print(chk[!(n_snp_rerun %in% ivw_n_snp) | !lead_found])
  snps[]
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 分段缩放 X 轴: 真实坐标 (bp) -> [0, 1]
# 切成 flank | Peak 区 | gap | Gene 区 | flank (顺序按真实位置; 区域重叠时自动合并), 每段线性; 宽度权重见 W_*
#   Peak 区 = peak 中点 ± ZOOM_HALF; Gene 区 = 整个基因体 (大基因压缩、小基因拉伸, 固定宽度)

make_axis <- function(pk_start, pk_end, tss, g_start, g_end, snp_pos = numeric(), n_gwas = 0) {
  pk_mid <- (pk_start + pk_end) / 2
  P <- c(min(pk_start, pk_mid - ZOOM_HALF), max(pk_end, pk_mid + ZOOM_HALF)); G <- c(g_start, g_end)
  zlo <- min(P[1], G[1]); zhi <- max(P[2], G[2])
  # flank 长度: 默认 FLANK_BP; 窗口内 MR SNP 不够 (N_SNP_TARGET - n_gwas) 个时, 向外延伸到够为止 (<= FLANK_MAX)
  need <- max(N_SNP_TARGET - n_gwas, 0) - sum(snp_pos >= zlo & snp_pos <= zhi)
  dout <- sort(pmax(zlo - snp_pos, snp_pos - zhi)[snp_pos < zlo | snp_pos > zhi])
  flank <- if (need > 0 && length(dout)) min(FLANK_MAX, max(FLANK_BP, dout[min(need, length(dout))] + 500)) else FLANK_BP
  win <- c(max(0, zlo - flank), zhi + flank)
  bp  <- sort(unique(c(win, P, G)))
  sg  <- data.table(lo = head(bp, -1), hi = tail(bp, -1))[, `:=`(mid = (lo + hi) / 2, len = hi - lo)]
  sg[, kind := fcase(mid >= P[1] & mid <= P[2], "peak", mid >= G[1] & mid <= G[2], "gene", mid < zlo | mid > zhi, "flank", default = "gap")]
  gap_len <- max(sg[kind == "gap", sum(len)], 1)
  sg[, w := fcase(kind == "peak", W_PEAK * len / diff(P),
                  kind == "gene", W_GENE * len / diff(G),
                  kind == "gap",  W_GAP * min(1, gap_len / GAP_FULL_BP) * len / gap_len,
                  default = W_FLANK * pmin(1, len / flank))]
  sg[, x1 := cumsum(w) / sum(w)][, x0 := data.table::shift(x1, fill = 0)]
  list(fx = approxfun(c(sg$lo, sg$hi[nrow(sg)]), c(0, sg$x1), rule = 2), win = win, P = P, G = G, tss = tss, sg = sg)
}

# 刻度: 只在 Peak 区 / Gene 区 (含 TSS), 标签 = 真实位置 (Mb, 小数位按刻度间距); 优先级 Peak > TSS > Gene, 挤在一起的去掉
axis_breaks <- function(A) {
  tk <- function(s, n, prio) {
    b <- pretty(s, n = n); b <- b[b >= s[1] & b <= s[2]]
    if (!length(b)) return(NULL)
    step <- if (length(b) > 1) min(diff(b)) else diff(s)
    data.table(pos = b, digits = max(2, ceiling(-log10(step / 1e6))), prio = prio)
  }
  g  <- tk(A$G, 3, 3)
  br <- rbind(tk(A$P, 3, 1), data.table(pos = A$tss, digits = max(3, if (is.null(g)) 3 else g$digits[1]), prio = 2), g)
  br <- unique(br[order(prio, pos)], by = "pos")
  br[, x := A$fx(pos)]
  keep <- rep(FALSE, nrow(br)); xs <- numeric(0)
  for (i in seq_len(nrow(br))) if (all(abs(br$x[i] - xs) >= TICK_SEP)) { keep[i] <- TRUE; xs <- c(xs, br$x[i]) }
  br <- br[keep][order(pos)]
  br[, label := mapply(function(p, d) formatC(p / 1e6, format = "f", digits = d), pos, digits)]
  br[]
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 作图

# 每条轨: TSS 虚线 (贯穿各轨, 对齐 link 终点); y 轴文字作为轨道标签 (保证 patchwork 左侧对齐)
track_base <- function(S, ylim) {
  ggplot() +
    geom_vline(xintercept = S$A$fx(S$tg$tss), linetype = "dashed", color = "grey60", linewidth = 0.4) +
    coord_cartesian(xlim = c(0, 1), ylim = ylim, expand = FALSE, clip = "off") +
    theme_classic(base_size = BASE_SIZE) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(),
          axis.text.y = element_text(face = "bold", size = BASE_SIZE - 1, color = "black"), plot.margin = margin(1, 6, 1, 4))
}
ylab_axis <- function(at, lab) scale_y_continuous(breaks = at, labels = lab)

## Gene 轨 (沿 Gene 区的缩放坐标画)
track_gene <- function(S) {
  tg <- S$tg; fx <- S$A$fx; win <- S$A$win
  gx <- exon_dt[chr == tg$chr & gene_name == S$g & end >= win[1] & start <= win[2]]
  gx[, `:=`(x1 = fx(pmax(start, win[1])), x2 = fx(pmin(end, win[2])))]
  g1 <- fx(max(tg$start, win[1])); g2 <- fx(min(tg$end, win[2])); xt <- fx(tg$tss)
  arr <- xt + ifelse(tg$strand == "+", 1, -1) * 0.035
  track_base(S, c(-0.5, 1.1)) +
    annotate("segment", x = g1, xend = g2, y = 0, yend = 0, color = GENE_COL, linewidth = 0.8) +
    geom_rect(data = gx[type != "cds"], aes(xmin = x1, xmax = x2, ymin = -0.16, ymax = 0.16), fill = GENE_COL, color = NA) +
    geom_rect(data = gx[type == "cds"], aes(xmin = x1, xmax = x2, ymin = -0.3, ymax = 0.3), fill = GENE_COL, color = NA) +
    annotate("segment", x = xt, xend = xt, y = 0, yend = 0.5, color = GENE_COL, linewidth = 0.9) +
    annotate("segment", x = xt, xend = arr, y = 0.5, yend = 0.5, color = GENE_COL, linewidth = 0.9, arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("text", x = (g1 + g2) / 2, y = 0.85, label = S$g, fontface = "bold.italic", size = 5.5, color = GENE_COL) +
    ylab_axis(0, "Gene")
}

## 挑 SNP (纯为美观): GWAS SNP 全留; MR instrument 在整个窗口内按离 peak 近优先挑, 凑到 N_SNP_TARGET 个, x 轴上间隔 >= SNP_MIN_SEP
pick_mr_snps <- function(sd, gw, A, pk_mid) {
  cand <- sd[pos >= A$win[1] & pos <= A$win[2] & !(pos %in% gw$pos)][, `:=`(x = A$fx(pos), d = abs(pos - pk_mid))][order(d)]
  n_want <- max(N_SNP_TARGET, nrow(gw)) - nrow(gw); xs <- A$fx(gw$pos); sel <- integer()
  for (i in seq_len(nrow(cand))) {
    if (length(sel) >= n_want) break
    if (all(abs(cand$x[i] - xs) >= SNP_MIN_SEP)) { sel <- c(sel, i); xs <- c(xs, cand$x[i]) }
  }
  cand[sel]
}

## SNP 轨: MR instrument (蓝色竖线) / GWAS SNP (紫色菱形 + rsID)
track_snps <- function(S) {
  fx <- S$A$fx; sd <- S$mr; gw <- S$gw
  p <- track_base(S, c(-0.15, 1.5))
  if (nrow(sd)) p <- p + geom_segment(data = sd, aes(x = fx(pos), xend = fx(pos), y = 0, yend = 0.55), color = SNP_COL, linewidth = 0.7)
  if (nrow(gw)) p <- p +
    geom_point(data = gw, aes(fx(pos), 0.28), shape = 23, fill = GWAS_COL, color = "white", size = 4.6, stroke = 0.6) +
    geom_text_repel(data = gw, aes(fx(pos), 0.28, label = rsid), size = 4, fontface = "bold", color = GWAS_COL,
                    nudge_y = 0.8, direction = "both", min.segment.length = 0, segment.size = 0.3, box.padding = 0.4, seed = 2)
  p + ylab_axis(0.3, "SNPs")
}

## Peak 轨
track_peak <- function(S) {
  fx <- S$A$fx; pk <- S$pk
  track_base(S, c(-0.1, 1)) +
    annotate("rect", xmin = fx(pk$pk_start), xmax = fx(pk$pk_end), ymin = 0.1, ymax = 0.8, fill = PEAK_FILL, color = NA) +
    ylab_axis(0.45, "Peak")
}

## Link 轨: 5 个方法各一行 (行内 y: 0 = 基线, 弧线向下); 无图例
ROW_H <- 1.0
track_links <- function(S, score_max) {
  fx <- S$A$fx; x1 <- fx(S$pk$pk_mid); x2 <- fx(S$tg$tss)
  lk <- copy(S$lk)[, `:=`(k = .I, base = -(.I - 1) * ROW_H)]
  s_lo <- -log10(FDR_CUT)
  lk[, `:=`(sig = !is.na(fdr) & fdr < FDR_CUT, score_plot = pmin(pmax(-log10(pmax(fdr, 1e-300)), s_lo), score_max))]
  lk[, direction := factor(fifelse(effect > 0, "Positive", "Negative"), levels = names(DIR_COL))]
  tt <- seq(0, 1, length.out = 120)
  arcs <- lk[, .(x = x1 + (x2 - x1) * tt, y = base - 0.8 * sin(pi * tt)), by = .(k, method, sig, direction, score_plot)]
  # y 轴标签 = 方法名 + FDR; ggplot 按 y 从小到大取标签样式, 所以全部反转
  lab <- rev(lk[, sprintf("%s\nFDR = %s", method, ifelse(is.na(fdr), "NA", formatC(fdr, format = "e", digits = 1)))])
  ybr <- rev(lk$base - 0.5)
  suppressWarnings(
  track_base(S, c(min(lk$base) - 0.95, 0.08)) +
    geom_hline(yintercept = lk$base, color = "grey70", linewidth = 0.4) +
    geom_path(data = arcs[sig == FALSE], aes(x, y, group = k), color = NS_COL, linetype = "22", linewidth = LW_NS) +
    geom_path(data = arcs[sig == TRUE], aes(x, y, group = k, color = direction, linewidth = score_plot), lineend = "round") +
    scale_color_manual(values = DIR_COL, drop = FALSE, guide = "none") +
    scale_linewidth_continuous(range = LW_RANGE, limits = c(s_lo, score_max), guide = "none") +
    ylab_axis(ybr, lab) +
    theme(axis.text.y = element_text(color = rev(lk$col), face = rev(ifelse(lk$sig, "bold", "plain")), size = BASE_SIZE - 2, lineheight = 0.95, hjust = 1),
          legend.position = "none"))
}

## x 轴 (单独一条): 刻度只在 Peak / Gene 区, 标签 = 真实位置; gap 处画 // 断裂 (轴线在 gap 处断开)
track_axis <- function(S) {
  br <- axis_breaks(S$A); sg <- S$A$sg
  gm <- sg[kind == "gap", (x0 + x1) / 2]; d <- 0.014
  cuts <- sort(c(0, 1, gm - d, gm + d))
  ln <- data.table(x = cuts[seq(1, length(cuts), 2)], xend = cuts[seq(2, length(cuts), 2)])
  track_base(S, c(0, 1)) +
    geom_segment(data = ln, aes(x = x, xend = xend, y = 0, yend = 0), linewidth = 0.4, color = "black") +
    { if (length(gm)) annotate("text", x = gm, y = 0, label = "//", size = 4.5, fontface = "bold") } +
    scale_x_continuous(breaks = br$x, labels = br$label) +
    labs(x = sprintf("%s position (Mb)  -  x-axis rescaled around peak and gene", S$tg$chr)) +
    theme(axis.ticks.x = element_line(linewidth = 0.4), axis.ticks.length.x = unit(1.5, "mm"),
          axis.text.x = element_text(size = BASE_SIZE - 3, color = "black"), axis.title.x = element_text(size = BASE_SIZE - 1, face = "bold"),
          axis.text.y = element_blank(), panel.background = element_blank())
}

dist_lab <- function(d) if (d < 1e3) sprintf("%.0f bp", d) else sprintf("%.1f kb", d / 1e3)

plot_pair <- function(pr, res, snp_dt, score_max) {
  g <- pr$gene_id; peak <- pr$peak_id
  pk <- data.table(peak_id = peak)[, c("pk_chr", "pk_start", "pk_end") := tstrsplit(peak_id, "[:-]")]
  pk[, `:=`(pk_start = as.numeric(pk_start), pk_end = as.numeric(pk_end))][, pk_mid := (pk_start + pk_end) / 2]
  cand <- gene_dt[gene_name == g & chr == pk$pk_chr][order(abs(tss - pk$pk_mid))]
  if (!nrow(cand)) { msg(sprintf("    %s: EnsDb v86 中找不到, 跳过", g)); return(NULL) }
  tg <- cand[1]
  r <- res[gene_id == g & peak_id == peak]
  if (!nrow(r)) { msg(sprintf("    %s__%s 不在 MR 结果中, 跳过", peak, g)); return(NULL) }
  lk <- METHODS[, .(method, col, fdr = vapply(fdr_col, function(cn) as.numeric(r[[cn]][1]), 0), effect = vapply(eff_col, function(cn) as.numeric(r[[cn]][1]), 0))]
  gw <- GW[cell_type == pr$cell_type & gene_id == g & peak_id == peak]
  sd <- snp_dt[gene_id == g & peak_id == peak]
  A  <- make_axis(pk$pk_start, pk$pk_end, tg$tss, tg$start, tg$end, snp_pos = setdiff(sd$pos, gw$pos), n_gwas = nrow(gw))
  S  <- list(g = g, tg = tg, pk = pk, lk = lk, A = A, gw = gw, mr = pick_mr_snps(sd, gw, A, pk$pk_mid))
  gw_out <- gw[pos < A$P[1] | pos > A$P[2]]
  if (nrow(gw_out)) msg("    警告: GWAS SNP 不在 Peak 区内: ", paste(gw_out$rsid, collapse = ", "))
  msg(sprintf("    SNP 轨: %d MR + %d GWAS (MR instrument 共 %d)", nrow(S$mr), nrow(gw), nrow(sd)))

  # 副标题: Gene TSS / Peak-Gene distance / peak 内 GWAS signal (trait, rsID, P)
  sub1 <- sprintf("Gene TSS: %s:%s (%s strand)   |   Peak-Gene distance: %s", tg$chr, format(tg$tss, big.mark = ","), tg$strand, dist_lab(abs(pk$pk_mid - tg$tss)))
  gs   <- unique(EX[cell_type == pr$cell_type & gene_id == g & peak_id == peak, .(trait, rsid, gwas_p)])[order(trait, gwas_p)]
  ent  <- sprintf("%s, %s, P = %s", gs$trait, gs$rsid, formatC(gs$gwas_p, format = "e", digits = 1))
  sub2 <- paste0("GWAS signal in peak:  ", paste(vapply(split(ent, ceiling(seq_along(ent) / 2)), paste, "", collapse = "   |   "), collapse = "\n                                       "))
  (track_gene(S) / track_snps(S) / track_peak(S) / track_links(S, score_max) / track_axis(S)) +
    plot_layout(heights = c(1.0, 1.2, 0.4, 4.4, 0.01)) +
    plot_annotation(title = sprintf("%s  |  %s  (%s)", g, peak, pr$cell_type), subtitle = paste(sub1, sub2, sep = "\n"),
                    theme = theme(plot.title = element_text(face = "bold", size = BASE_SIZE + 4),
                                  plot.subtitle = element_text(size = BASE_SIZE - 2, color = "grey20", lineheight = 1.15)))
}

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 运行

cts <- intersect(c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7"), PAIRS$cell_type)
RES <- setNames(lapply(cts, function(ct) fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")))), cts)
SNPS <- setNames(lapply(cts, function(ct) collect_mr_snps(ct, PAIRS[cell_type == ct], RES[[ct]])), cts)

# 线宽标尺在所有页面共用
all_fdr <- unlist(lapply(cts, function(ct) as.matrix(RES[[ct]][PAIRS[cell_type == ct], on = .(gene_id, peak_id), nomatch = NULL, METHODS$fdr_col, with = FALSE])))
SCORE_MAX <- min(SCORE_CAP, max(2, -log10(pmax(all_fdr[!is.na(all_fdr) & all_fdr < FDR_CUT], 1e-300))))
msg(sprintf("线宽标尺: -log10(FDR) in [%.2f, %.2f]", -log10(FDR_CUT), SCORE_MAX))

cairo_pdf(PDF_FILE, width = PAGE_W, height = PAGE_H, onefile = TRUE)
n_ok <- 0L
for (i in seq_len(nrow(PAIRS))) {
  pr <- PAIRS[i]
  msg(sprintf("  [%d/%d] %s | %s | %s", i, nrow(PAIRS), pr$cell_type, pr$gene_id, pr$peak_id))
  p <- tryCatch(plot_pair(pr, RES[[pr$cell_type]], SNPS[[pr$cell_type]], SCORE_MAX), error = function(e) { msg("    失败: ", conditionMessage(e)); NULL })
  if (!is.null(p)) {
    grid.newpage(); pushViewport(viewport(x = 0.01, y = 0.5, width = FIG_FRAC, height = 0.98, just = "left"))
    print(p, newpage = FALSE)
    grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1.2))
    popViewport(); n_ok <- n_ok + 1L
  }
}
invisible(dev.off())
unlink(LD_DIR, recursive = TRUE)
msg(sprintf("Saved (%d pages): %s", n_ok, PDF_FILE))
