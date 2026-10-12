# Fig6: Motif 富集 —— 只合并 3 个细胞类型的结果为 Motif_all.tsv, 并原样复制 paper 的手工整理表 Fig6_Motif_table.xlsx
# 不画 logo、不做 16 个 motif 的表、不生成新的 xlsx (xlsx 是本地 Claude Code 结合外部知识手工整理的)
# 移植自 20261006_paper/Scripts/Plot6_val_Motif.sh 第 1 节
# 输入: Results/Validation/Valid3_Motif/<CT>_motif_enrichment.tsv.gz
#       (对照) paper 目录 Materials/Motif_details/Motif_all.tsv (旧结果), Materials/Fig6_Motif_table.xlsx
# 输出: Results/Figures/Motif_all.tsv, Results/Figures/Fig6_Motif_table.xlsx
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

### 3. 复制手工整理的 xlsx (原样)
out_xlsx <- file.path(FIG_DIR, "Fig6_Motif_table.xlsx")
stopifnot(file.copy(file.path(PAPER_DIR, "Fig6_Motif_table.xlsx"), out_xlsx, overwrite = TRUE))
msg("Copied: ", out_xlsx)

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
