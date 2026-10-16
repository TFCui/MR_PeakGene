# Fig6: Motif 富集 —— 合并 3 个细胞类型的结果为 Motif_all.tsv, 并把 paper 手工整理表 Fig6_Motif_table.xlsx 里的 raw p 换成 FDR
# 不画 logo、不做 16 个 motif 的表 (xlsx 是本地 Claude Code 结合外部知识手工整理的; 这里只替换 E-J 列数值和高亮)
# 移植自 20261006_paper/Scripts/Plot6_val_Motif.sh 第 1 节
# 输入: Results/Validation/Valid3_Motif/<CT>_motif_enrichment.tsv.gz
#       (对照) paper 目录 Materials/Motif_details/Motif_all.tsv (旧结果), Materials/Fig6_Motif_table.xlsx (源表, 只读)
# 输出: FIG_DIR/Motif_all.tsv, FIG_DIR/Fig6_Motif_table.xlsx (E-J = Valid3 的 BH FDR; 高亮 = FDR < 0.05 且 FE > 1)
# 运行: Rscript Scripts/Visualization/Plot6_Motif.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")

PAPER_DIR  <- file.path(BASE, "Results/RunMR/Materials_runMR/20261006_paper/Materials")
MOTIF_METH <- c("IVW", "MR_Egger", "SMR", "CCT", "WNN", "Pseudo")   # 与 paper 的 Motif_all.tsv 相同 (Motif 不用统一显示名)

### 1. 合并三个细胞类型
Motif_all <- rbindlist(lapply(CELLS, function(ct) fread(file.path(VAL_DIR, "Valid3_Motif", paste0(ct, "_motif_enrichment.tsv.gz")))))
Motif_all[, `:=`(celltype = factor(celltype, levels = CELLS), method = factor(method, levels = MOTIF_METH))]
stopifnot(!anyNA(Motif_all$celltype), !anyNA(Motif_all$method))
setorder(Motif_all, celltype, method, rank)
print(Motif_all[, .(n_motif = .N, n_sig = sum(p.adjust < 0.05 & fold.enrichment > 1), n_fg = n_fg[1]), by = .(celltype, method)])

out_tsv <- file.path(FIG_DIR, "Motif_all.tsv")
fwrite(Motif_all, out_tsv, sep = "\t", quote = FALSE)
msg("Saved: ", out_tsv)

### 2. 与 paper 旧 Motif_all.tsv 对照 (xlsx 中的数字来自旧结果)
old <- fread(file.path(PAPER_DIR, "Motif_details/Motif_all.tsv"))
stopifnot(identical(names(old), names(Motif_all)))
old[, `:=`(celltype = factor(celltype, levels = CELLS), method = factor(method, levels = MOTIF_METH))]
cmp <- merge(old, Motif_all, by = c("celltype", "method", "motif"), suffixes = c("_old", "_new"), all = TRUE)
cat(sprintf("\n----- vs paper Motif_all.tsv: old %d rows, new %d rows, matched %d\n",
            nrow(old), nrow(Motif_all), cmp[!is.na(pvalue_old) & !is.na(pvalue_new), .N]))
print(cmp[, .(n_fg_old = n_fg_old[1], n_fg_new = n_fg_new[1],
              max_abs_dlog10p = max(abs(log10(pvalue_old) - log10(pvalue_new)), na.rm = TRUE),
              n_sig_old = sum(p.adjust_old < 0.05 & fold.enrichment_old > 1, na.rm = TRUE),
              n_sig_new = sum(p.adjust_new < 0.05 & fold.enrichment_new > 1, na.rm = TRUE),
              n_sig_flip = sum((p.adjust_old < 0.05 & fold.enrichment_old > 1) != (p.adjust_new < 0.05 & fold.enrichment_new > 1), na.rm = TRUE)),
          by = .(celltype, method)])

### 3. 手工整理的 xlsx: E-J 列 (raw p) 换成 Valid3 的 FDR (p.adjust), 高亮 = FDR < 0.05 且 FE > 1; 其余 (logo、解读、表达量) 原样保留
suppressMessages(library(openxlsx))
src_xlsx  <- file.path(PAPER_DIR, "Fig6_Motif_table.xlsx")
out_xlsx  <- file.path(FIG_DIR, "Fig6_Motif_table.xlsx")
XL_ROWS   <- 5:20                                                      # 16 个 motif 行
XL_METH   <- c("SMR", "IVW", "MR_Egger", "CCT", "WNN", "Pseudo")       # E-J 列顺序 (Causal SMR/IVW/Egger/Joint, Marginal Multiomics/Pseudobulk)
tb <- as.data.table(read.xlsx(src_xlsx, sheet = 1, rows = XL_ROWS, cols = 1:10, colNames = FALSE, skipEmptyCols = FALSE))
setnames(tb, c("celltype", "tf", "motif", "logo", XL_METH))

# 不用 openxlsx 回写 (它会把 logo 的 drawing XML 嵌套包一层, Excel 可能报错); 直接解压、改 sheet1.xml 的单元格、重新打包, logo 等文件原样保留
# 样式沿用原表: s="5" = 普通 0.0000, s="6" = 高亮 (黄底粗体 0.0000)
td <- file.path(PROJ, "tmp/fig6_xlsx_src"); unlink(td, recursive = TRUE); dir.create(td, recursive = TRUE)
unzip(src_xlsx, exdir = td)
sheet_f <- file.path(td, "xl/worksheets/sheet1.xml")
xml <- readChar(sheet_f, file.size(sheet_f), useBytes = TRUE); Encoding(xml) <- "UTF-8"
hl_old_xml <- regmatches(xml, gregexpr('<c r="[E-J][0-9]+" s="6"', xml))[[1]]
hl_old_xml <- sub('^<c r="([E-J][0-9]+)".*$', "\\1", hl_old_xml)

tgt <- melt(tb[, c("celltype", "motif", XL_METH), with = FALSE], id.vars = c("celltype", "motif"),
            variable.name = "method", value.name = "p_xlsx")
tgt[, row := rep(seq_along(XL_ROWS), times = length(XL_METH))]            # melt 按列展开: 每个方法 16 行
tgt[, `:=`(cell = paste0(LETTERS[4 + match(as.character(method), XL_METH)], XL_ROWS[row]),
           celltype = factor(celltype, levels = CELLS), method = factor(method, levels = MOTIF_METH))]
tgt <- merge(tgt, cmp[, .(celltype, method, motif, p_old = pvalue_old, fdr_old = p.adjust_old, fe_old = fold.enrichment_old,
                          p_new = pvalue_new, fdr_new = p.adjust_new, fe_new = fold.enrichment_new)],
             by = c("celltype", "method", "motif"), all.x = TRUE)
stopifnot(nrow(tgt) == 96, !anyNA(tgt$fdr_new), !anyNA(tgt$fdr_old))
stopifnot(all(abs(tgt$p_xlsx - signif(tgt$p_old, 3)) <= 1e-9))          # 原表 = 旧结果的 raw p (3 位有效数字)
stopifnot(setequal(hl_old_xml, tgt[fdr_old < 0.05 & fe_old > 1, cell]))   # 原表高亮规则 = 旧 FDR < 0.05 且 FE > 1
tgt[, hl := fdr_new < 0.05 & fe_new > 1]
msg("xlsx highlight: old ", length(hl_old_xml), " -> new ", sum(tgt$hl), "; flipped: ",
    paste(tgt[hl != (fdr_old < 0.05 & fe_old > 1), paste(tb$tf[row], tb$motif[row], method, ifelse(hl, "->sig", "->ns"))], collapse = "; "))

# 写入 FDR (3 位有效数字, 同原表) + 高亮样式
for (i in seq_len(nrow(tgt))) {
  pat <- sprintf('<c r="%s" s="[0-9]+" t="n"><v>[^<]*</v></c>', tgt$cell[i])
  stopifnot(length(regmatches(xml, gregexpr(pat, xml))[[1]]) == 1)
  xml <- sub(pat, sprintf('<c r="%s" s="%d" t="n"><v>%s</v></c>', tgt$cell[i], if (tgt$hl[i]) 6L else 5L,
                          format(signif(tgt$fdr_new[i], 3), scientific = FALSE)), xml)
}
# 标题 / 表头 / 脚注 (XML 里 "<" 要转义)
swap <- function(xml, from, to) { stopifnot(grepl(from, xml, fixed = TRUE)); sub(from, to, xml, fixed = TRUE) }
xml <- swap(xml, "raw p-values", "FDR")
xml <- swap(xml, "Raw p-value by method", "FDR (BH) by method")
n_motif <- Motif_all[celltype == CELLS[1] & method == MOTIF_METH[1], .N]
foot_old <- "Shaded bold p-values: FDR &lt; 0.05 for that method. p-values are from motif enrichment of MR-significant peaks (method-specific foreground) vs background."
foot_new <- paste0("Shaded bold FDR: FDR &lt; 0.05 and fold enrichment &gt; 1 for that method. FDR is Benjamini-Hochberg adjusted across all ", n_motif,
                   " JASPAR motifs within each cell type and method; motif enrichment of MR-significant peaks (method-specific foreground) vs matched background.")
xml <- swap(xml, foot_old, foot_new)
writeChar(xml, sheet_f, eos = NULL, useBytes = TRUE)
# 重新打包 (mirror 模式保留目录结构; [Content_Types].xml 放第一个)
unlink(out_xlsx)
zfiles <- unzip(src_xlsx, list = TRUE)$Name
zfiles <- c("[Content_Types].xml", setdiff(zfiles, "[Content_Types].xml"))
zip::zip(out_xlsx, files = zfiles, root = td, mode = "mirror")
msg("Saved: ", out_xlsx)

### 4. xlsx 中的 16 个目标 motif (来自 paper Plot6_val_Motif.sh): 新旧 p.adjust 对照, 看手工表的结论是否仍成立
target_motifs <- data.table(
  celltype  = c(rep("Bn_TCL1A", 11), rep("CD4_Tn_CCR7", 2), rep("CD8_Tn_CCR7", 3)),
  motif     = c("MA0488.1","MA0492.1","MA1133.1","MA1140.2","MA1136.1",
                "MA0475.2","MA0018.4","MA1483.1","MA1466.1","MA0759.1","MA1632.1",
                "MA0492.1","MA1106.1",
                "MA0018.4","MA1466.1","MA1145.1"),
  tf_target = c("JUN","JUND","JUN;JUNB","JUNB","FOSB;JUNB",
                "FLI1","CREB1","ELF2","ATF6","ELK3","ATF2",
                "JUND","HIF1A",
                "CREB1","ATF6","FOSL2;JUND"))
target_motifs[, `:=`(celltype = factor(celltype, levels = CELLS), ord = .I)]
tg <- cmp[target_motifs, on = .(celltype, motif)]
tg[, cell := fifelse(p.adjust_new < 0.05 & fold.enrichment_new > 1, "S", ".")]
tg[, cell_old := fifelse(p.adjust_old < 0.05 & fold.enrichment_old > 1, "S", ".")]
tg[, cell := paste0(cell_old, ">", cell)]   # 旧 > 新, S = 显著 (p.adjust < 0.05 且 FE > 1)
cat("\n----- 16 target motifs: significance old>new (S = sig, . = not)\n")
print(dcast(tg, ord + celltype + tf_target + motif ~ method, value.var = "cell")[, ord := NULL][])
