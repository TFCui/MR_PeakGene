# Table 4 (论文 S4 Table): Bn / CD4 / CD8 三个细胞类型的数据基本信息 (只读汇总, 登录节点, 几分钟)
# 块: B 单细胞 / pseudobulk (QC 后), C cis-QTL 汇总统计, D 候选 pair 漏斗与 MR 分析集, E 边际相关方法的样本
# 口径 (见 Scripts/Pipeline.md 第 2 节、CLAUDE.md 第 6 节):
#   - 个体数 / 细胞数: Data/OMIX009786/cell_rna/QC_pseudobulk/libsize_<CT>.csv (同时通过 RNA 和 ATAC 过滤的个体)
#   - QTL 样本量: 汇总统计 sample_sizes 列 (众数; 同时给范围)
#   - "peaks with >= 1 caQTL at P < 1e-8": caQTL_sig/<CT>_Significant.txt.gz 里的 peak 数 (名义 P < 1e-8, 即选工具变量的阈值; 不是 CIMA 论文 FDR 定义的 caPeak)
#   - 漏斗: Scripts/Model/logs/MR_AllPairs_20074030_{0,1,2}.out (3 个 CT 的正式运行), 最后一行与 Results/MR 表核对
# 输出: Results/Figures/Figures_20261027/Table_4.xlsx (单个 sheet) + Table_4.tsv (长表, 用于和论文 tex 表核对)
# 运行: Rscript Scripts/Visualization/Table4_DataOverview.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressMessages({ library(openxlsx); library(scales) })

OMIX <- file.path(BASE, "Data/OMIX009786")
CT_LBL <- c(Bn_TCL1A = "Bn", CD4_Tn_CCR7 = "CD4", CD8_Tn_CCR7 = "CD8")
LOG_ID <- c(Bn_TCL1A = 0, CD4_Tn_CCR7 = 1, CD8_Tn_CCR7 = 2)
fmt <- function(x) comma(x)
medq <- function(x) sprintf("%s (%s\u2013%s)", fmt(round(median(x))), fmt(round(quantile(x, 0.25))), fmt(round(quantile(x, 0.75))))

shell_num <- function(cmd) as.numeric(system(cmd, intern = TRUE))
mode_range <- function(file, col, header = TRUE) {   # 全文件扫描某一列: 众数 + 范围
  tb <- fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NR>1{c[$%d]++} END{for(k in c) print k\"\\t\"c[k]}'", shQuote(file), col), header = FALSE, col.names = c("v", "n"))
  list(mode = tb[which.max(n), v], min = min(tb$v), max = max(tb$v), n_distinct = nrow(tb))
}

rows <- list(); add <- function(block, item, vals, note = "") rows[[length(rows) + 1]] <<- data.table(block = block, item = item, Bn = vals[1], CD4 = vals[2], CD8 = vals[3], note = note)

# ---- B. 单细胞 / pseudobulk ----
ls_dt <- lapply(CELLS, function(ct) fread(file.path(OMIX, "cell_rna/QC_pseudobulk", paste0("libsize_", ct, ".csv"))))
add("B", "Individuals (passing RNA and ATAC cell-count filters)", sapply(ls_dt, function(d) fmt(nrow(d))))
add("B", "scRNA-seq cells", sapply(ls_dt, function(d) fmt(sum(d$rna_cells))))
add("B", "scATAC-seq cells", sapply(ls_dt, function(d) fmt(sum(d$atac_cells))))
add("B", "scRNA-seq cells per individual, median (IQR)", sapply(ls_dt, function(d) medq(d$rna_cells)))
add("B", "scATAC-seq cells per individual, median (IQR)", sapply(ls_dt, function(d) medq(d$atac_cells)))
n_gene <- sapply(CELLS, function(ct) shell_num(sprintf("zcat %s | wc -l", shQuote(file.path(OMIX, "cell_rna/QC_pseudobulk", paste0("pseudobulk_rna_", ct, ".tsv.gz"))))) - 1)
n_peak <- sapply(CELLS, function(ct) shell_num(sprintf("zcat %s | wc -l", shQuote(file.path(OMIX, "cell_rna/QC_pseudobulk", paste0("pseudobulk_atac_", ct, ".tsv.gz"))))) - 1)
n_col  <- sapply(CELLS, function(ct) shell_num(sprintf("zcat %s | head -1 | tr '\\t' '\\n' | wc -l", shQuote(file.path(OMIX, "cell_rna/QC_pseudobulk", paste0("pseudobulk_rna_", ct, ".tsv.gz")))))) # 含行名列
stopifnot(all(n_col - 1 == sapply(ls_dt, nrow)))   # pseudobulk 矩阵列数 = 个体数
add("B", "Genes in the pseudobulk RNA matrix", fmt(n_gene))
add("B", "Peaks in the pseudobulk ATAC matrix", fmt(n_peak))

# ---- C. cis-QTL 汇总统计 ----
eq <- lapply(CELLS, function(ct) mode_range(file.path(OMIX, "cis_eQTL_Summary_Statics/cis_eQTL", paste0(ct, ".txt.gz")), 10))
ca <- lapply(CELLS, function(ct) mode_range(file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL", paste0(ct, ".txt.gz")), 10))
rng <- function(l) if (l$min == l$max) fmt(l$mode) else sprintf("%s (%s\u2013%s)", fmt(l$mode), fmt(l$min), fmt(l$max))
add("C", "eQTL sample size", sapply(eq, rng), "sample_sizes column of the cis-eQTL summary statistics")
add("C", "caQTL sample size", sapply(ca, rng), "sample_sizes column of the cis-caQTL summary statistics")
sig_peaks <- sapply(CELLS, function(ct) shell_num(sprintf("zcat %s | awk -F'\\t' 'NR>1{p[$1]=1} END{k=0; for(x in p) k++; print k}'",
                    shQuote(file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL/caQTL_sig", paste0(ct, "_Significant.txt.gz"))))))
sig_max <- sapply(CELLS, function(ct) shell_num(sprintf("zcat %s | awk -F'\\t' 'NR>1{ if($7>m) m=$7 } END{print m}'",
                    shQuote(file.path(OMIX, "cis_caQTL_Summary_Statics/cis_caQTL/caQTL_sig", paste0(ct, "_Significant.txt.gz"))))))
stopifnot(all(sig_max < 1e-8))
add("C", "Peaks with at least one caQTL at P < 1e-8", fmt(sig_peaks), "nominal P < 1e-8, the threshold used to select instruments (not an FDR-defined caPeak)")

# ---- D. 候选 pair 漏斗 ----
parse_log <- function(id) {
  l <- readLines(file.path(PROJ, "Scripts/Model/logs", sprintf("MR_AllPairs_20074030_%d.out", id)))
  num <- function(pat, k = 1) as.numeric(gsub(",", "", regmatches(l[grep(pat, l)[1]], gregexpr("[0-9]+", l[grep(pat, l)[1]]))[[1]][k]))
  list(all_gp = num("all_gene_peak", 3), shared = num("QTL 共享 variant", 2), wnn = num("WNN \\+ pseudobulk 均可用"),   # shared: 第 1 个数字是 ">= 3" 里的 3
       nona = num("删除 NA 后"), final = as.numeric(sub(".*后 +([0-9]+) pairs.*", "\\1", l[grep("ivw_sigma2 <= 10 后", l)[1]])),   # 行内含 ivw_sigma2 / 10, 不能取第一个数字
       final_gene = as.numeric(sub(".*gene +([0-9]+).*", "\\1", l[grep("ivw_sigma2 <= 10 后", l)[1]])),
       final_peak = as.numeric(sub(".*peak +([0-9]+).*", "\\1", l[grep("ivw_sigma2 <= 10 后", l)[1]])))
}
lg <- lapply(LOG_ID, parse_log)
stopifnot(!anyNA(unlist(lg)), lg[[1]]$shared > 3)   # 防止日志解析取错 / 取空
add("D", "Candidate gene\u2013peak pairs (All_GenePeak, TSS\u2013peak distance < 1 Mb)", sapply(lg, function(x) fmt(x$all_gp)))
add("D", "With caQTL and eQTL variants in common (at least 3 shared variants)", sapply(lg, function(x) fmt(x$shared)))
add("D", "Gene and peak present in both pseudobulk and WNN data", sapply(lg, function(x) fmt(x$wnn)))
add("D", "After removing pairs with a missing p-value in any base method", sapply(lg, function(x) fmt(x$nona)))
add("D", "Final analysis set (IVW residual variance \u2264 10)", sapply(lg, function(x) fmt(x$final)))
mr <- lapply(CELLS, read_mr)
n_snp <- lapply(mr, function(d) d$ivw_n_snp)
add("D", "Genes / peaks in the final analysis set", sapply(mr, function(d) sprintf("%s / %s", fmt(uniqueN(d$gene_id)), fmt(uniqueN(d$peak_id)))))
add("D", "Variants (instruments) per pair in the IVW fit, median (IQR)", sapply(n_snp, medq), "number of variants used by the IVW fit")
# 核对: 日志的最后一行 = MR 结果表
for (i in seq_along(CELLS)) stopifnot(nrow(mr[[i]]) == lg[[i]]$final, uniqueN(mr[[i]]$gene_id) == lg[[i]]$final_gene, uniqueN(mr[[i]]$peak_id) == lg[[i]]$final_peak)

# ---- E. 边际相关方法的样本 ----
add("E", "Pseudobulk correlation: individuals", sapply(ls_dt, function(d) fmt(nrow(d))), "same individuals as in block B")
add("E", "WNN correlation: metacells (external PBMC multiome)", rep("80", 3), "the same external metacell set for all three cell types")

tb <- rbindlist(rows)
print(tb[, .(block, item, Bn, CD4, CD8)], nrows = 100)
print(data.table(ct = CELLS, eq_min = sapply(eq, `[[`, "min"), eq_max = sapply(eq, `[[`, "max"), ca_min = sapply(ca, `[[`, "min"), ca_max = sapply(ca, `[[`, "max"), sig_max = sig_max))

OUT_DIR <- file.path(PROJ, "Results/Figures/Figures_20261027")
fwrite(tb, file.path(OUT_DIR, "Table_4.tsv"), sep = "\t")

# ---- xlsx: 单个 sheet, 按块分段 ----
BLK <- c(B = "Single-cell data used for the pseudobulk correlation (after quality control)", C = "cis-QTL summary statistics (CIMA, OMIX009786)",
         D = "Candidate peak\u2013gene pairs and the Mendelian randomization analysis set", E = "Samples for the marginal correlation methods")
st_title <- createStyle(fontSize = 14, textDecoration = "bold")
st_head  <- createStyle(fontSize = 11, textDecoration = "bold", fontColour = "#FFFFFF", fgFill = "#1F3A5F", halign = "center", valign = "center", border = "TopBottomLeftRight", borderColour = "#BBBBBB")
st_sect  <- createStyle(fontSize = 11, textDecoration = "bold", fgFill = "#DCE6F1", border = "TopBottomLeftRight", borderColour = "#BBBBBB")
st_item  <- createStyle(fontSize = 11, valign = "center", wrapText = TRUE, border = "TopBottomLeftRight", borderColour = "#BBBBBB")
st_val   <- createStyle(fontSize = 11, halign = "center", valign = "center", border = "TopBottomLeftRight", borderColour = "#BBBBBB")
wb <- createWorkbook(); sh <- "Table_4"; addWorksheet(wb, sh, gridLines = FALSE)
writeData(wb, sh, "S4 Table. Data overview for the three cell types analyzed (Bn, CD4, CD8)", startRow = 1, startCol = 1); addStyle(wb, sh, st_title, rows = 1, cols = 1)
writeData(wb, sh, t(c("", "Bn (naive B, TCL1A)", "CD4 (naive CD4 T, CCR7)", "CD8 (naive CD8 T, CCR7)")), startRow = 3, startCol = 1, colNames = FALSE)
addStyle(wb, sh, st_head, rows = 3, cols = 1:4, gridExpand = TRUE)
r <- 4
for (b in names(BLK)) {
  writeData(wb, sh, BLK[[b]], startRow = r, startCol = 1); mergeCells(wb, sh, cols = 1:4, rows = r); addStyle(wb, sh, st_sect, rows = r, cols = 1:4, gridExpand = TRUE); r <- r + 1
  x <- tb[block == b]
  writeData(wb, sh, x[, .(item, Bn, CD4, CD8)], startRow = r, startCol = 1, colNames = FALSE)
  addStyle(wb, sh, st_item, rows = r + seq_len(nrow(x)) - 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sh, st_val,  rows = r + seq_len(nrow(x)) - 1, cols = 2:4, gridExpand = TRUE)
  r <- r + nrow(x)
}
foot <- c("The Chinese Immune Multi-Omics Atlas (CIMA) comprises scRNA-seq and scATAC-seq data from 428 Chinese adults annotated into 73 L4 immune cell types; cis-QTL were mapped within 1 Mb of the gene TSS or peak midpoint (Yin et al.).",
          "Three sample definitions are used: QTL sample sizes (summary statistics, used by the Mendelian randomization methods), individuals in the pseudobulk matrices (passing the RNA and ATAC cell-count filters), and 80 metacells from an external PBMC multiome dataset (WNN correlation). caQTL and eQTL come from the same cohort, so the two QTL data sets overlap in samples.",
          "Candidate pairs: All_GenePeak pairs intersected with the LD-clumped caQTL (P < 1e-8) and eQTL variants, at least 3 shared variants per pair, with the gene and the peak present in both the pseudobulk and WNN data; pairs with a missing p-value in any base method (IVW, MR-Egger, SMR, WNN, pseudobulk) or an IVW residual variance above 10 were removed.")
writeData(wb, sh, foot, startRow = r + 1, startCol = 1)
for (j in seq_along(foot)) { mergeCells(wb, sh, cols = 1:4, rows = r + j); addStyle(wb, sh, createStyle(wrapText = TRUE, valign = "top", fontSize = 9), rows = r + j, cols = 1); setRowHeights(wb, sh, rows = r + j, heights = 52) }
setColWidths(wb, sh, cols = 1:4, widths = c(62, 24, 24, 24)); freezePane(wb, sh, firstActiveRow = 4)
saveWorkbook(wb, file.path(OUT_DIR, "Table_4.xlsx"), overwrite = TRUE)
msg("Saved: ", file.path(OUT_DIR, "Table_4.xlsx"), " + Table_4.tsv")
