# Fig3: gene 内超几何富集的 global test 表 (6 方法 x 6 金标准), 3 个细胞类型合并在同一个 sheet; 不再画 p 值直方图
# 移植自 20261006_paper/Scripts/Plot3_val_Hyper.sh, global P 的算法不变 (各 gene 的 N/K/n_sig/k 求和后做一次超几何检验)
# 新增: 合并 fold enrichment; global P 用 log.p 计算, 极小的 P 也不下溢
# 输入: Results/Validation/Valid1_Hyper/<CT>_hyper_long.tsv
# 输出: Results/Figures/Figures_20261027/Table_1_2_global.tsv (长表), Table_1_2.xlsx (单个 sheet: 上块 Global P = 论文 S1 Table, 下块 pooled fold = S2 Table; 每块 3 个 CT x 6 方法)
# 运行: Rscript Scripts/Visualization/Plot3_Hyper.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressMessages(library(openxlsx))

in_dir <- file.path(VAL_DIR, "Valid1_Hyper")
hyper_list <- setNames(lapply(CELLS, function(ct) fread(file.path(in_dir, paste0(ct, "_hyper_long.tsv")))), CELLS)

m_lbls   <- c(SMR = "Causal SMR", IVW = "Causal IVW", Egger = "Causal Egger",
              CCT = "Causal Joint", WNN = "Marginal Multiomics", Pseudobulk = "Marginal Pseudobulk")   # 数据中的 method -> 显示名
db_names <- c("enhancerAtlas", "FANTOM5", "4DGenome", "GTEx", "HiChIP", "PCHiC")
db_lbls  <- c(enhancerAtlas = "EnhancerAtlas", FANTOM5 = "FANTOM5", `4DGenome` = "4DGenome", GTEx = "GTEx", HiChIP = "HiChIP", PCHiC = "PCHiC")

# ---------- 单个细胞类型: 各格子 (方法 x 金标准) 的 N/K/n/k 求和后算一个 global P 和 pooled fold ----------
global_p <- function(hyper_long, ct) {
  p_dt <- hyper_long[!is.na(pval) & N > 0 & K > 0 & n_sig > 0 & N > K][
    , `:=`(method_lbl = factor(m_lbls[as.character(method)], levels = METHODS),
           database   = factor(database, levels = db_names))]
  stopifnot(!anyNA(p_dt$method_lbl), !anyNA(p_dt$database))
  p_dt[, .(n_gene = .N, sum_N = sum(N), sum_K = sum(K), sum_nsig = sum(n_sig), sum_k = sum(k)),
       by = .(method_lbl, database)][
    , `:=`(celltype  = ct,
           gp        = phyper(sum_k - 1L, sum_K, sum_N - sum_K, sum_nsig, lower.tail = FALSE),
           log10_gp  = phyper(sum_k - 1L, sum_K, sum_N - sum_K, sum_nsig, lower.tail = FALSE, log.p = TRUE) / log(10),
           fold      = (sum_k / sum_nsig) / (sum_K / sum_N))]
}

res <- rbindlist(lapply(CELLS, function(ct) global_p(hyper_list[[ct]], ct)))
res[, celltype := factor(celltype, levels = CELLS)]
setorder(res, celltype, method_lbl, database)
stopifnot(nrow(res) == length(CELLS) * length(METHODS) * length(db_names), all(is.finite(res$log10_gp)))

# log10 P -> "x.xe-NN" 文本 (不下溢)
fmt_p <- function(l) {
  e <- floor(l); m <- 10^(l - e)
  up <- round(m, 1) >= 10; m[up] <- 1; e[up] <- e[up] + 1
  sprintf("%.1fe%d", m, as.integer(e))
}
res[, p_txt := fmt_p(log10_gp)]

# 日志: 每个格子的 gene 数、global P、fold
for (ct in CELLS) {
  cat("\n=====", ct, "\n")
  x <- res[celltype == ct]
  print(dcast(x, method_lbl ~ database, value.var = "n_gene"))
  print(dcast(x, method_lbl ~ database, value.var = "p_txt"))
  print(dcast(x, method_lbl ~ database, value.var = "fold", fun.aggregate = function(v) round(v, 2)))
}

# ---------- tsv ----------
OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261027")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(OUT_DIR, "Table_1_2_global.tsv")
fwrite(res[, .(celltype, method = method_lbl, database, n_gene, sum_N, sum_K, sum_nsig, sum_k, fold, log10_gp, gp, p_txt)],
       out_tsv, sep = "\t", quote = FALSE)
msg("Saved: ", out_tsv)

# ---------- xlsx: 每个 CT 一个 sheet, 上 = Global P, 下 = pooled fold enrichment (每列最高 fold 高亮) ----------
st_title <- createStyle(fontSize = 16, textDecoration = "bold")
st_sect  <- createStyle(fontSize = 13, textDecoration = "bold")
st_head  <- createStyle(fontSize = 12, textDecoration = "bold", fontColour = "#FFFFFF", fgFill = "#1F3A5F", halign = "center", valign = "center",
                        border = "TopBottomLeftRight", borderColour = "#BBBBBB", wrapText = TRUE)
st_row   <- createStyle(fontSize = 12, textDecoration = "bold", fgFill = "#DCE6F1", halign = "left", valign = "center",
                        border = "TopBottomLeftRight", borderColour = "#BBBBBB")
st_cell  <- createStyle(fontSize = 12, halign = "center", valign = "center", border = "TopBottomLeftRight", borderColour = "#BBBBBB")
st_fold  <- createStyle(fontSize = 12, halign = "center", valign = "center", border = "TopBottomLeftRight", borderColour = "#BBBBBB", numFmt = "0.00")
st_best  <- createStyle(fontSize = 12, textDecoration = "bold", fgFill = "#FFF2CC", halign = "center", valign = "center",
                        border = "TopBottomLeftRight", borderColour = "#BBBBBB", numFmt = "0.00")
foot <- c("Global P: for each method x gold standard, N, K, n_sig and k are summed over genes (only genes with N > 0, K > 0, n_sig > 0 and N > K) and one hypergeometric test (upper tail) is run on the pooled counts.",
          "Fold enrichment (pooled) = (sum k / sum n_sig) / (sum K / sum N), i.e. the fraction of significant peaks overlapping the gold standard divided by the fraction of background peaks overlapping it. Shaded bold = highest fold in that column.",
          "N = background peaks of the gene; K = background peaks in the gold standard; n_sig = significant peaks; k = significant peaks in the gold standard.")

# 一个 sheet: 上块 = Global P (论文 S1 Table), 下块 = pooled fold (论文 S2 Table); 每块 18 行 = 3 个 CT x 6 方法, CT 单元格纵向合并; 每个 CT 内每列最高 fold 高亮
wb <- createWorkbook(); sh <- "Table_1_2"; addWorksheet(wb, sh)
CT_LBL  <- c(Bn_TCL1A = "Bn", CD4_Tn_CCR7 = "CD4", CD8_Tn_CCR7 = "CD8")
M_SHOW  <- setNames(METHODS, METHODS); M_SHOW[c("Causal Joint", "Marginal Multiomics")] <- c("Causal Joint (CCT)", "Marginal Multiomics (WNN)")
NM <- length(METHODS); ND <- length(db_names); NC <- 2 + ND
mk_long <- function(val) {   # 18 行 x (Cell type, Method, 6 个数据库)
  rbindlist(lapply(CELLS, function(ct) {
    m <- dcast(res[celltype == ct], method_lbl ~ database, value.var = val)[match(METHODS, method_lbl)]
    data.table(`Cell type` = CT_LBL[[ct]], Method = unname(M_SHOW[as.character(m$method_lbl)]), m[, db_names, with = FALSE])
  }))
}
tp <- mk_long("p_txt"); tf <- mk_long("fold"); setnames(tp, db_names, db_lbls[db_names]); setnames(tf, db_names, db_lbls[db_names])
stopifnot(nrow(tp) == length(CELLS) * NM, nrow(tf) == nrow(tp))

writeData(wb, sh, "Hypergeometric enrichment of method-specific significant peaks in gold-standard links (Bn, CD4, CD8)", startRow = 1, startCol = 1)
addStyle(wb, sh, st_title, rows = 1, cols = 1)
blk <- function(top, ttl, tb, st_val, best = FALSE) {   # top = 标题行; 表头 top + 1; 数据 top + 2 起
  writeData(wb, sh, ttl, startRow = top, startCol = 1); addStyle(wb, sh, st_sect, rows = top, cols = 1)
  writeData(wb, sh, as.data.frame(tb), startRow = top + 1, startCol = 1, headerStyle = st_head)
  r <- top + 1 + seq_len(nrow(tb))
  addStyle(wb, sh, st_row, rows = r, cols = 1:2, gridExpand = TRUE)
  addStyle(wb, sh, st_val, rows = r, cols = 3:NC, gridExpand = TRUE)
  for (k in seq_along(CELLS)) {   # CT 单元格纵向合并
    rr <- top + 1 + (k - 1) * NM + seq_len(NM)
    mergeCells(wb, sh, cols = 1, rows = rr)
    addStyle(wb, sh, createStyle(fontSize = 12, textDecoration = "bold", fgFill = "#DCE6F1", halign = "center", valign = "center",
                                 border = "TopBottomLeftRight", borderColour = "#BBBBBB"), rows = rr, cols = 1, gridExpand = TRUE)
    if (best) for (j in seq_len(ND)) addStyle(wb, sh, st_best, rows = rr[which.max(tb[[2 + j]][(k - 1) * NM + seq_len(NM)])], cols = 2 + j)
  }
  setRowHeights(wb, sh, rows = top + 1, heights = 24)
}
blk(3, "Global hypergeometric P (S1 Table)", tp, st_cell)
blk(3 + NM * length(CELLS) + 3, "Fold enrichment, pooled (S2 Table)", tf, st_fold, best = TRUE)
writeData(wb, sh, foot, startRow = 3 + 2 * (NM * length(CELLS) + 3) + 1, startCol = 1)
setColWidths(wb, sh, cols = 1:NC, widths = c(11, 28, rep(16, ND)))
freezePane(wb, sh, firstCol = TRUE)
out_xlsx <- file.path(OUT_DIR, "Table_1_2.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
msg("Saved: ", out_xlsx)
