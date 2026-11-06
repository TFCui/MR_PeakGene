# Plot15 (论文暂定 S7 Fig + 正文 Figure 3A): 全基因组 Manhattan / Miami 图 (全部候选 peak–gene pair)
# Figure S7: 3 页 (8.5 x 11 in), 每个 CT 一页, 6 条轨道 (Causal SMR / IVW / Egger / Joint, Marginal Multiomics / Pseudobulk)
# Figure 3A: CD4 T cell (EX_CT) 的 Miami 图: 上半 Causal Joint (向上), 下半 Marginal Multiomics (镜像向下, 高度为上半的 0.55);
#   2026-11-04 起从 S7 第 1 页移入正文 Figure 3; 这里只构建并 saveRDS 到 tmp/Plot15_Manhattan/miami_CD4.rds, 由 Plot9_GO_GRN.R 读入拼图
#   (改 Miami 图后要先跑本脚本, 再跑 Plot9)
#   Figure 2 的例子 pair (Example_pairs.tsv 中 pair_set = Figure_2; 改例子只需改这张表) 在其所在 CT 的页面上每条轨道都标,
#   不显著也标在实际位置 (体现 MR 强、marginal 弱的对比): 洋红菱形点 + 浅洋红底色的粗体标签 (与 top 基因标签同一图层, 互不遮挡)
# 口径:
#   x = peak 中心 ((start + end) / 2), 用 hg38 chr1–22 长度累加成全基因组坐标 (结果中只有常染色体)
#   y = −log10 P (各方法的 <m>_p), 线性等距 (2026-11-04 起不再分段压缩); 每条轨道用自己的范围 (各方法量级差很多)
#   阈值线 = Bonferroni (0.05 / 该 CT 的 pair 数); FDR 是 gene 内 BH, 没有统一阈值线, 只写进 tsv
#   图中不写 pair 数 / 过线数等数字, 这些只存在 Figure_S7_data.tsv 中 (写 Results 时查)
#   配色 (2026-11-04): 不按方法分色; 过线的点 causal 用亮橙 COL_CAUSAL, marginal 用亮绿 COL_MARG;
#   其余按奇偶染色体交替浅灰 / 中灰; 点用 ggrastr 栅格化
#   top 基因: 过线 pair 中每个 gene 取最强 pair, 按 P 贪心挑选, 同一染色体 5 Mb 内只标一个 (Miami 每半边 10 个, S7 每轨道 8 个)
# 输入: Results/MR/<CT>_AllPairs_results.tsv.gz; Scripts/Visualization/Example_pairs.tsv
# 输出: Results/Figures/Figures_20261102/Figure_S7.pdf (只出 PDF) + Figure_S7_data.tsv (type = count: 每个 CT × 方法的计数; type = label: 标注的基因, panel = Figure_3A / S7_p1–p3)
#       tmp/Plot15_Manhattan/miami_CD4.rds (Figure 3A 的 ggplot 对象)
# 运行: Rscript Scripts/Visualization/Plot15_Manhattan.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressPackageStartupMessages({ library(patchwork); library(ggrepel); library(ggrastr) })

M6     <- METHODS                         # Causal SMR / IVW / Egger / Joint, Marginal Multiomics / Pseudobulk
IS_CAUSAL <- setNames(grepl("^Causal", M6), M6)
WIN    <- 5e6                             # 标签去重窗口 (同一染色体)
COL_CAUSAL <- "#FF9F1C"; COL_MARG <- "#3BB273"   # 亮橙 / 亮绿
COL_HL  <- "#C2185B"; COL_HL_TXT <- "#8A0F3F"; COL_HL_BG <- "#FBE1EC"   # 例子 pair 的高亮色 (洋红)、标签文字色、标签底色
EX_CT   <- "CD4_Tn_CCR7"                         # 第 1 页 Miami 图的 CT
COL_ODD <- "#D4D4D4"; COL_EVEN <- "#ADADAD"
BASEF <- 10
FONT  <- "sans"
col_of <- function(mm) if (IS_CAUSAL[[mm]]) COL_CAUSAL else COL_MARG

# hg38 chr1–22 长度
CHR_LEN <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622,
             133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468)
GAP <- 1.2e7
CHR <- data.table(chr = paste0("chr", 1:22), len = CHR_LEN)
CHR[, offset := c(0, cumsum(len + GAP)[-22])]
CHR[, mid := offset + len / 2]
X_MAX <- CHR[22, offset + len]

###------------------------------------------------------------------------------------------------
### 1. 读数据: 每个 CT 一张长表 (pair × 6 方法)
load_ct <- function(ct) {
  d <- read_mr(ct)
  n0 <- nrow(d)
  d[, c("chr", "s", "e") := tstrsplit(peak_id, "[:-]")]
  d[, pos := (as.numeric(s) + as.numeric(e)) / 2]
  d <- merge(d, CHR[, .(chr, offset)], by = "chr")
  d[, `:=`(x = offset + pos, chr_i = as.integer(sub("chr", "", chr)))]
  stopifnot(nrow(d) == n0)
  rbindlist(lapply(M6, function(mm) {
    k <- METHOD_KEY[[mm]]
    d[, .(cell_type = ct, method = mm, gene_id, peak_id, chr, chr_i, pos, x, p = get(paste0(k, "_p")), fdr = get(paste0(k, "_fdr")))]
  }))[, mlp := -log10(p)]
}
D <- setNames(lapply(CELLS, load_ct), CELLS)
N_PAIR <- sapply(D, function(d) nrow(d) / length(M6))
BONF   <- -log10(0.05 / N_PAIR)

EX <- fread(file.path(PROJ, "Scripts/Visualization/Example_pairs.tsv"))[pair_set == "Figure_2", .(cell_type, gene_id, peak_id)]

### counts
CNT <- rbindlist(lapply(CELLS, function(ct) D[[ct]][, .(n_pair = .N, bonf_mlp = round(BONF[[ct]], 3),
  n_above = sum(mlp > BONF[[ct]]), n_gene_above = uniqueN(gene_id[mlp > BONF[[ct]]]),
  n_fdr0.1 = sum(fdr < 0.1), max_mlp = round(max(mlp), 1)), by = .(cell_type, method)]))
cat("\n===== 每个 CT × 方法: pair 数 / Bonferroni / 过线 pair 与 gene 数 / FDR < 0.1 / 最大 −log10 P\n"); print(CNT)
ex_rank <- rbindlist(lapply(seq_len(nrow(EX)), function(i) {
  D[[EX$cell_type[i]]][, rk := frank(p, ties.method = "min"), by = method][peak_id == EX$peak_id[i] & gene_id == EX$gene_id[i],
    .(cell_type, gene_id, method, mlp = round(mlp, 1), rank = rk)]
}))
cat("\n===== Figure 2 例子 pair 在各方法中的 −log10 P 和名次\n"); print(dcast(ex_rank, cell_type + gene_id ~ method, value.var = "mlp"))

###------------------------------------------------------------------------------------------------
### 2. 标签挑选
pick_labels <- function(dd, bonf, n_lab, pre = NULL) {
  cand <- dd[mlp > bonf][order(p)][, .SD[1], by = gene_id][order(p)]
  kept <- pre[, .(chr, pos)]
  out <- list()
  for (i in seq_len(nrow(cand))) {
    if (length(out) >= n_lab) break
    r <- cand[i]
    if (!is.null(kept) && nrow(kept) && any(kept$chr == r$chr & abs(kept$pos - r$pos) < WIN)) next
    out[[length(out) + 1]] <- r
    kept <- rbind(kept, r[, .(chr, pos)])
  }
  if (length(out)) rbindlist(out) else dd[0]
}

###------------------------------------------------------------------------------------------------
### 3. 一条轨道 (flip = TRUE: Miami 图下半, 画 y = −(−log10 P), 刻度写绝对值)
#   top 基因和例子的标签放在同一个 geom_label_repel 图层里, 互相避让; 例子用浅洋红底色 + 粗体洋红字, top 基因无底色
# cols: 非显著点 / 例子高亮的配色 (默认 = S7; Figure 3A 传薄荷色板); strip: 非 NULL 时在轨道上方加标签条 (Figure 3A 用)
COLS_S7 <- list(odd = COL_ODD, even = COL_EVEN, hl = COL_HL, hl_txt = COL_HL_TXT, hl_bg = COL_HL_BG)
manhattan_track <- function(dd, name, col, bonf, n_lab, highlight = NULL, show_x = FALSE, flip = FALSE, size_pt = 0.55,
                            cols = COLS_S7, strip = NULL) {
  sgn <- if (flip) -1 else 1
  dd <- copy(dd)[, `:=`(y = sgn * mlp, strip_lab = "s")]
  dd[, grp := fifelse(mlp > bonf, "sig", fifelse(chr_i %% 2 == 1, "odd", "even"))]
  ymax <- max(dd$mlp); lim <- ymax * 1.16
  hl  <- if (!is.null(highlight) && nrow(highlight)) merge(dd, highlight, by = c("gene_id", "peak_id")) else dd[0]
  lab <- pick_labels(dd[!gene_id %in% hl$gene_id], bonf, n_lab, pre = hl)
  L <- rbind(hl[, .(gene_id, x, y, key = "ex")], lab[, .(gene_id, x, y, key = "top")])
  L[, `:=`(ny = sgn * fifelse(key == "ex", 0.18, 0.08) * ymax, nx = fifelse(key == "ex", 4e7, 0),
           ff = fifelse(key == "ex", "bold.italic", "italic"), sz = fifelse(key == "ex", 2.9, 2.5),
           segc = fifelse(key == "ex", cols$hl, "grey45"))]
  bk <- pretty(c(0, ymax), 4); bk <- bk[bk <= ymax * 1.1]
  g <- ggplot(dd, aes(x, y)) +
    rasterise(geom_point(data = dd[grp != "sig"], aes(color = grp), size = size_pt, stroke = 0, shape = 16), dpi = 300) +
    rasterise(geom_point(data = dd[grp == "sig"], color = col, size = size_pt * 1.8, stroke = 0, shape = 16, alpha = 0.9), dpi = 300) +
    geom_hline(yintercept = sgn * bonf, color = "grey35", linewidth = 0.35, linetype = "22") +
    geom_point(data = lab, color = "grey20", fill = col, shape = 21, size = size_pt * 3.2, stroke = 0.35) +
    geom_point(data = hl, shape = 23, fill = cols$hl, color = "grey10", size = 3.2, stroke = 0.6) +
    geom_label_repel(data = L, aes(label = gene_id, color = key, fill = key, fontface = ff, size = sz, segment.colour = segc),
                     family = FONT, label.size = NA, label.padding = 0.12, label.r = 0.12, segment.size = 0.35,
                     min.segment.length = 0, box.padding = 0.3, point.padding = 0.15, nudge_x = L$nx, nudge_y = L$ny,
                     max.overlaps = Inf, seed = 1, show.legend = FALSE) +
    scale_color_manual(values = c(odd = cols$odd, even = cols$even, ex = cols$hl_txt, top = "grey10"), guide = "none") +
    scale_fill_manual(values = c(ex = cols$hl_bg, top = NA), na.value = NA, guide = "none") +
    scale_size_identity() +
    scale_x_continuous(breaks = CHR$mid, labels = if (show_x) 1:22 else NULL, expand = expansion(add = GAP / 2), limits = c(0, X_MAX)) +
    scale_y_continuous(breaks = sgn * bk, labels = bk, expand = expansion(mult = 0.01)) +
    coord_cartesian(ylim = sort(c(0, sgn * lim)), clip = "off") +
    labs(x = NULL, y = bquote(atop(bold(.(name)), -log[10] ~ italic(P)))) +
    theme_classic(base_size = BASEF, base_family = FONT) +
    theme(axis.text = element_text(color = "black", size = BASEF - 2), axis.text.x = element_text(size = BASEF - 2.5),
          axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.3),
          axis.ticks.x = if (show_x) element_line(linewidth = 0.3) else element_blank(),
          axis.line.x = if (flip) element_blank() else element_line(linewidth = 0.35),
          axis.title.y = element_text(size = BASEF - 1.5, lineheight = 0.9, color = "grey15"), plot.margin = margin(3, 6, 3, 4))
  # strip = "" : 也做 facet 但不显示标签条. Miami 两条轨道必须同为 facet 图, 否则 patchwork 加宽共用的 y 轴列时,
  #   非 facet 图的 y 轴线会画在该列最左侧, 与上半的轴线错开 (2026-11-04 第 2 轮修正)
  if (!is.null(strip)) g <- g + facet_wrap(~ strip_lab, labeller = as_labeller(function(x) strip)) +
    (if (nzchar(strip)) theme_strip_mint(10.5) else theme(strip.text = element_blank(), strip.background = element_blank()))
  list(plot = g, labels = rbind(hl[, .(gene_id, peak_id, mlp, fdr, is_example = TRUE)],
                                lab[, .(gene_id, peak_id, mlp, fdr, is_example = FALSE)]))
}

###------------------------------------------------------------------------------------------------
### 4. 拼页
LABS <- list()
add_lab <- function(r, panel, ct, mm) LABS[[length(LABS) + 1]] <<- r$labels[, `:=`(panel = panel, cell_type = ct, method = mm)]
page_title <- function(txt) plot_annotation(title = txt, theme = theme(plot.title = element_text(face = "bold", size = BASEF + 3, family = FONT, hjust = 0)))
ex_of <- function(ct) EX[cell_type == ct, .(gene_id, peak_id)]

# Figure 3A: Miami 图 (EX_CT; 上 Causal Joint, 下 Marginal Multiomics); 上半加标签条 "CD4 T cell" (2026-11-04 起 Plot9 不再加标题)
MI_TOP <- "Causal Joint"; MI_BOT <- "Marginal Multiomics"
dct <- D[[EX_CT]]
# Figure 3A 配色 (2026-11-04, 薄荷色板 PAL_MINT): causal 钢蓝 / marginal 薄荷绿, 非显著点浅冷灰; 例子 pair 用柔和珊瑚色 (全图唯一暖色点缀)
COLS_MI <- list(odd = "#DADFE3", even = "#BCC4CA", hl = "#E07A5F", hl_txt = "#9C3D25", hl_bg = "#FBE5DE")
COL_MI  <- c(top = PAL_MINT[["m6"]], bot = PAL_MINT[["m4"]])
r_top <- manhattan_track(dct[method == MI_TOP], MI_TOP, COL_MI[["top"]], BONF[[EX_CT]], 10, ex_of(EX_CT), show_x = TRUE, size_pt = 0.75,
                         cols = COLS_MI, strip = sprintf("Genome-wide peak\u2013gene associations (%s)", CELL_LAB[[EX_CT]]))   # 标签条: 简洁标题 + 括号内细胞类型
r_bot <- manhattan_track(dct[method == MI_BOT], MI_BOT, COL_MI[["bot"]], BONF[[EX_CT]], 10, ex_of(EX_CT), flip = TRUE, size_pt = 0.75, cols = COLS_MI, strip = "")
add_lab(r_top, "Figure_3A", EX_CT, MI_TOP); add_lab(r_bot, "Figure_3A", EX_CT, MI_BOT)
# 压缩 y: 下半 (Multiomics 点多但无基因名) 只占上半高度的 MI_BOT_H; 两半之间去掉空白, 只留中间一条染色体编号
MI_BOT_H <- 0.75   # 2026-11-04 第 2 轮: 0.55 → 0.75, 下半要放得下竖排的 "Marginal multiomics" 标题
# y 轴标题 (2026-11-04): 不用 ggplot 的 y 标题 (它紧贴各自的刻度数字, 上半 "200" 与下半 "30" 宽度不同, 两行标题会错位),
#   改为左边单独一列文字格, 两个标题在同一竖线上; 方法名用论文写法 (lab_paper: Causal joint / Marginal multiomics)
mi_title <- function(mm, y = 0.5) wrap_elements(full = grid::textGrob(bquote(atop(bold(.(lab_paper(mm))), -log[10] ~ italic(P))),
  rot = 90, y = y, gp = grid::gpar(fontsize = 11, fontfamily = FONT, col = "grey10"))) + theme(plot.margin = margin(0, 0, 0, 0))
no_ytitle <- theme(axis.title.y = element_blank(), axis.text.y = element_text(hjust = 1))
miami <- wrap_plots(mi_title(MI_TOP, y = 0.46), r_top$plot + no_ytitle + theme(plot.margin = margin(3, 6, 0, 2)),   # y = 0.46: 避开上方标签条
                    mi_title(MI_BOT), r_bot$plot + no_ytitle + theme(plot.margin = margin(0, 6, 3, 2)), ncol = 2) +
  plot_layout(widths = unit(c(0.5, 1), c("in", "null")), heights = c(1, MI_BOT_H))
TMP15 <- file.path(TMP_DIR, "Plot15_Manhattan"); dir.create(TMP15, recursive = TRUE, showWarnings = FALSE)
saveRDS(miami, file.path(TMP15, "miami_CD4.rds"))
msg("Saved: ", file.path(TMP15, "miami_CD4.rds"))

# S7 第 1–3 页: 每个 CT × 6 个方法
pages_ct <- lapply(seq_along(CELLS), function(j) {
  ct <- CELLS[j]; d <- D[[ct]]
  pl <- lapply(seq_along(M6), function(i) {
    mm <- M6[i]
    r <- manhattan_track(d[method == mm], mm, col_of(mm), BONF[[ct]], 8, ex_of(ct), show_x = i == length(M6))
    add_lab(r, paste0("S7_p", j), ct, mm)
    r$plot
  })
  wrap_plots(pl, ncol = 1) + page_title(CELL_LAB[[ct]])
})

FIG_W <- 8.5; FIG_H <- 11
OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261102")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- file.path(OUT_DIR, "Figure_S7.pdf")
cairo_pdf(out, width = FIG_W, height = FIG_H, family = FONT, onefile = TRUE); for (pg in pages_ct) print(pg); invisible(dev.off())
msg("Saved: ", out)

###------------------------------------------------------------------------------------------------
### 5. tsv (一个文件)
LAB <- rbindlist(LABS)
LAB[, rank := seq_len(.N), by = .(panel, cell_type, method)]
LAB <- LAB[, .(panel, cell_type, method, rank, gene_id, peak_id, mlp = round(mlp, 2), fdr = signif(fdr, 3), is_example)]
OUT_TSV <- file.path(OUT_DIR, "Figure_S7_data.tsv")
fwrite(rbind(CNT[, type := "count"], LAB[, type := "label"], fill = TRUE)[, c("type", setdiff(names(CNT), "type"), setdiff(names(LAB), c("type", names(CNT)))), with = FALSE],
       OUT_TSV, sep = "\t", na = "")
cat("\n===== 标注的基因\n"); print(LAB[, .(genes = paste(gene_id, collapse = ", ")), by = .(panel, cell_type, method)])
msg("Saved: ", OUT_TSV)
