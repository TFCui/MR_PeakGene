# 2026-10-16
# Valid7_GO: 导出 MR-Only gene list, 供 DAVID 网站做 GO / KEGG 富集 (https://davidbioinformatics.nih.gov)
# MR-Only link: IVW / Egger / SMR 任一 FDR < 0.1, 且 WNN、Pseudo 的 FDR 都 >= 0.1
# 每个细胞类型 (Bn_TCL1A / CD4_Tn_CCR7 / CD8_Tn_CCR7) 取 MR-Only link 的 unique gene (symbol), 按字母排序
# 输出: Results/Validation/Valid7_GO/DAVID_input/MROnly_genes_<CT>.txt (纯文本, 一行一个 symbol, 无表头)
# DAVID 上传: Select Identifier = OFFICIAL_GENE_SYMBOL; List Type = Gene List; Species = Homo sapiens; Background 保持默认 (Homo sapiens);
#   Functional Annotation Chart 勾 GOTERM_BP_DIRECT / GOTERM_MF_DIRECT / GOTERM_CC_DIRECT / KEGG_PATHWAY, 导出结果放 Results/Validation/Valid7_GO/DAVID_output/
# 运行 (登录节点, 几秒): Rscript Scripts/Validation/Valid7_GO.R

suppressPackageStartupMessages(library(data.table))

PROJ    <- "/scratch/user/tfcui2025/CIMA/Claude_project"
MR_DIR  <- file.path(PROJ, "Results/MR")
OUT_DIR <- file.path(PROJ, "Results/Validation/Valid7_GO/DAVID_input")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

celltypes <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
FDR_CUT   <- 0.1

for (ct in celltypes) {
  mr  <- fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")))
  mro <- mr[(ivw_fdr < FDR_CUT | egger_fdr < FDR_CUT | smr_fdr < FDR_CUT) & wnn_fdr >= FDR_CUT & pseudo_fdr >= FDR_CUT]
  genes <- sort(unique(mro$gene_id))
  stopifnot(!anyNA(genes), all(nzchar(genes)), !anyDuplicated(genes))
  f <- file.path(OUT_DIR, sprintf("MROnly_genes_%s.txt", ct))
  writeLines(genes, f)
  cat(sprintf("%-12s MR-Only link %5d | gene %5d | %s\n", ct, nrow(mro), length(genes), f))
}
