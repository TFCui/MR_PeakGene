# 2026-10-09
# Valid3_Motif 步骤 2: 合并 3 个 CT 的 motif 结果, 分组, 检查 TF 表达 (流程见 Scripts/Pipeline.md 3.3)
# 移植自 Scripts_reference/Motif1_step2_RunAll.sh 第 1-4 部分 (交互式代码改为脚本; 不做第 5 部分的汇总表)
# 输出: Results/Validation/Valid3_Motif/Motif_sig_3ct.tsv, Motif_expression_motif.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Signac)
  library(Matrix)
})
options(width = 300)

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
RES_DIR  <- file.path(BASE_DIR, "Claude_project/Results/Validation/Valid3_Motif")
RNA_DIR  <- file.path(BASE_DIR, "Data/OMIX009786/cell_rna/QC_filtered")

CELLS        <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
FDR_CUT      <- 0.05
EXPR_PCT_CUT <- 0.05
METHODS      <- c("IVW", "MR_Egger", "SMR", "CCT", "WNN", "Pseudo")
MR_METHODS   <- c("IVW", "MR_Egger", "SMR", "CCT")
MG_METHODS   <- setdiff(METHODS, MR_METHODS)

### 合并三个细胞类型
Motif_all <- rbindlist(lapply(CELLS, function(ct) fread(file.path(RES_DIR, paste0(ct, "_motif_enrichment.tsv.gz")))))
Motif_all[, `:=`(celltype = factor(celltype, levels = CELLS), method = factor(method, levels = METHODS))]
stopifnot(!anyNA(Motif_all$method), !anyNA(Motif_all$celltype))

### 每个 celltype x method 的显著 motif
Motif_sig <- Motif_all[p.adjust < FDR_CUT & fold.enrichment > 1]
Motif_sig[, `:=`(mlog10_padj = -log10(pmax(p.adjust, 1e-300)), method_class = fifelse(method %in% MR_METHODS, "MR", "Marginal"))]
setorder(Motif_sig, celltype, method, p.adjust)
fwrite(Motif_sig, file.path(RES_DIR, "Motif_sig_3ct.tsv"), sep = "\t", quote = FALSE)
cat("显著 motif 数 (p.adjust < 0.05 & fold > 1) / 前景 peak 数:\n")
print(dcast(Motif_all[, .(n_fg = n_fg[1], n_sig = sum(p.adjust < FDR_CUT & fold.enrichment > 1)), by = .(celltype, method)],
            method ~ celltype, value.var = c("n_sig", "n_fg")))

### 1. motif -> TF 基因符号: 拆异二聚体 "::", 去掉 "(var.N)", 统一大写
motif_tf <- unique(Motif_all[, .(motif, motif.name)])
motif_tf <- motif_tf[, .(tf = toupper(trimws(strsplit(sub("\\(var\\.[0-9]+\\)$", "", motif.name), "::", fixed = TRUE)[[1]]))), by = .(motif, motif.name)]

### 2. 每个 celltype x motif: 6 个方法各一列, TRUE = 该方法认为显著; 分组
motif_grp <- dcast(Motif_sig, celltype + motif + motif.name ~ method, value.var = "p.adjust", fun.aggregate = length, fill = 0L, drop = c(TRUE, FALSE))
motif_grp[, (METHODS) := lapply(.SD, function(x) x > 0L), .SDcols = METHODS]
n_mr <- rowSums(motif_grp[, ..MR_METHODS]); n_mg <- rowSums(motif_grp[, ..MG_METHODS])
motif_grp[, group := fcase(n_mr == 4L & n_mg == 2L, "Shared_all6", n_mr >= 1L & n_mg == 0L, "MR_only", n_mr == 0L & n_mg >= 1L, "Marginal_only", default = "Partial")]
print(dcast(motif_grp, celltype ~ group, fun.aggregate = length))

### 3. TF 表达: counts layer 实为 log1p(CP10k), expm1 反变换回 CP10k
tf_expr_stats <- function(obj, genes, ct) {
  feat <- Features(obj[["RNA"]], layer = "counts")
  g    <- intersect(genes, feat)
  m    <- drop0(LayerData(obj, assay = "RNA", layer = "counts", fast = TRUE)[match(g, feat), , drop = FALSE])
  n_cell <- ncol(m); n_expr <- tabulate(m@i + 1L, nbins = nrow(m))
  m@x <- expm1(m@x); s_cp10k <- Matrix::rowSums(m)
  res <- data.table(celltype = ct, tf = g, in_rna = TRUE, pct_expr = n_expr / n_cell, mean_cp10k_expr = s_cp10k / pmax(n_expr, 1L), mean_cp10k_all = s_cp10k / n_cell)
  rbind(res, data.table(celltype = ct, tf = setdiff(genes, feat), in_rna = FALSE), fill = TRUE)
}

TF_expr <- rbindlist(lapply(CELLS, function(ct) {
  cat("读 RNA:", ct, "\n")
  obj <- readRDS(file.path(RNA_DIR, paste0("rna_", ct, "_qc.rds")))
  out <- tf_expr_stats(obj, unique(motif_tf$tf), ct)
  rm(obj); invisible(gc()); out
}))
TF_expr[, celltype := factor(celltype, levels = CELLS)]
cat("RNA 中找不到的 TF:", paste(TF_expr[in_rna == FALSE, unique(tf)], collapse = ", "), "\n")

### 4. motif 层结果: 异二聚体取表达比例最低的亚基 (limiting TF) 作为该 motif 的表达值
BY_COLS <- c("celltype", "motif", "motif.name", METHODS, "group")
Motif_group_expr <- TF_expr[motif_tf[motif_grp, on = .(motif, motif.name), allow.cartesian = TRUE], on = .(celltype, tf)][, {
  j <- if (all(in_rna)) which.min(pct_expr) else which(!in_rna)[1]
  .(tfs = paste(tf, collapse = ";"), all_in_rna = all(in_rna), limiting_tf = tf[j], pct_expr = pct_expr[j],
    mean_cp10k_expr = mean_cp10k_expr[j], mean_cp10k_all = mean_cp10k_all[j], pct_expr_each = paste(sprintf("%.3f", pct_expr), collapse = ";"))
}, by = BY_COLS]
Motif_group_expr[, expressed := all_in_rna & pct_expr >= EXPR_PCT_CUT]
setorder(Motif_group_expr, celltype, group, -pct_expr)
fwrite(Motif_group_expr, file.path(RES_DIR, "Motif_expression_motif.tsv"), sep = "\t", quote = FALSE)

cat("\n表达的 motif 比例 (按 CT x group):\n")
print(Motif_group_expr[, .(n_motif = .N, n_expressed = sum(expressed), frac = round(mean(expressed), 3)), keyby = .(celltype, group)])
cat("\n全部完成\n")
