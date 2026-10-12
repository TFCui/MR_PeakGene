# 2026-10-09
# Valid5_GRN: 与 CIMA SCENIC+ GRN 的重叠 (流程见 Scripts/Pipeline.md 3.5)
# 移植自 Scripts_reference/GRN1_peak_gene.sh 第 0-5、7 部分 (交互式代码改为脚本; 去掉未使用的 Motif 读取; 不画图)
# 两套 GRN: CT = 该 L4 细胞类型 active (AUC > 1) 的 eRegulon; CT_HQ = CT ∩ 对应 lineage 的高质量 eRegulon
# 每个方法: FDR < 0.05 显著; 背景 = 该 CT 全部 MR pair; 2x2 Fisher (在 GRN 中 / 不在 x 显著 / 不显著)
# 一个 job 处理 3 个 CT; 输出 Results/Validation/Valid5_GRN/{GRN_sig_counts_by_method.tsv, GRN_summary.tsv, GRN_pairs_inMR.tsv.gz}

suppressPackageStartupMessages(library(data.table))
options(width = 300)

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
PROJ     <- file.path(BASE_DIR, "Claude_project")
MR_DIR   <- file.path(PROJ, "Results/MR")
GRN_DIR  <- file.path(BASE_DIR, "Results/RunMR/Materials_runMR/20260430_LinkPeakGene_more/MRresult_20260911/Valid6_CIMA_GRN/Data")
OUT_DIR  <- file.path(PROJ, "Results/Validation/Valid5_GRN")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

celltypes <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
FDR_CUT   <- 0.05
methods   <- c("ivw", "egger", "smr", "cct", "wnn", "pseudo")
fdr_cols  <- paste0(methods, "_fdr")
HQ_COLS   <- list(Bn_TCL1A = "selected_eRegulons_B", CD4_Tn_CCR7 = "selected_eRegulons_CD4T", CD8_Tn_CCR7 = "selected_eRegulons_CD8T")

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 1. MR results (三个细胞类型合成一张长表)

MR_all <- rbindlist(lapply(celltypes, function(ct) {
  dt <- fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")))
  dt[, cell_type := ct]
}), use.names = TRUE)
setcolorder(MR_all, "cell_type")
setkey(MR_all, cell_type, gene_id, peak_id)
stopifnot(!anyDuplicated(MR_all, by = c("cell_type", "gene_id", "peak_id")))
print(MR_all[, .N, by = cell_type])

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 2. CIMA GRN 原始数据

meta <- fread(file.path(GRN_DIR, "CIMA_eRegulons_Metadata.csv"))
if (names(meta)[1] %in% c("", "V1")) meta[, (names(meta)[1]) := NULL]
meta[, eRegulon := sub("_\\([0-9]+r\\)$", "", Region_signature_name)]
meta <- unique(meta[, .(TF, eRegulon, Region, Gene, Consensus_name, R2G_importance, R2G_rho)])

# AUC: 每个 L4 细胞类型 active 的 eRegulons (Activators + Repressors)
auc <- fread(file.path(GRN_DIR, "CIMA_eRegulons_AUC_Large_than_1_in_cell_type_L4.csv"))
auc_long <- unique(auc[cell_type_l4 %chin% celltypes,
                       .(eRegulon = trimws(unlist(strsplit(eRegulon_list, ",", fixed = TRUE)))),
                       by = .(cell_type_l4, Repressor_Activator)])

# HQ: 宽表 -> 长表 (各列长度不同, 空值去掉)
hq <- fread(file.path(GRN_DIR, "CIMA_Selected_High_Quality_eRegulons.csv"))
hq_long <- melt(hq, measure.vars = names(hq), variable.name = "hq_set", value.name = "eRegulon", variable.factor = FALSE)
hq_long <- unique(hq_long[!is.na(eRegulon) & nzchar(trimws(eRegulon))][, eRegulon := trimws(eRegulon)])

cat(sprintf("AUC eRegulon: %d, 在 meta 中: %d | HQ eRegulon: %d, 在 meta 中: %d\n",
            uniqueN(auc_long$eRegulon), sum(unique(auc_long$eRegulon) %chin% meta$eRegulon),
            uniqueN(hq_long$eRegulon),  sum(unique(hq_long$eRegulon)  %chin% meta$eRegulon)))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 3. 两套 eRegulon 集合 -> 两套 GRN

reg_sets <- rbindlist(lapply(celltypes, function(ct) {
  ct_regs <- auc_long[cell_type_l4 == ct, unique(eRegulon)]
  hq_regs <- hq_long[hq_set %chin% HQ_COLS[[ct]], unique(eRegulon)]
  rbind(data.table(grn_set = "CT",    cell_type = ct, eRegulon = ct_regs),
        data.table(grn_set = "CT_HQ", cell_type = ct, eRegulon = intersect(ct_regs, hq_regs)))
}))
print(dcast(reg_sets[, .N, by = .(grn_set, cell_type)], cell_type ~ grn_set, value.var = "N"))

GRN <- meta[reg_sets, on = "eRegulon", allow.cartesian = TRUE, nomatch = 0L]
setcolorder(GRN, c("grn_set", "cell_type", "TF", "eRegulon", "Region", "Gene"))
GRN[, in_mr := !is.na(MR_all[.(GRN$cell_type, GRN$Gene, GRN$Region), which = TRUE, mult = "first"])]

grn_summary <- GRN[, .(n_eRegulon  = uniqueN(eRegulon),
                       n_TF        = uniqueN(TF),
                       n_row       = .N,
                       n_pair      = uniqueN(.SD, by = c("Region", "Gene")),
                       n_pair_inMR = uniqueN(.SD[in_mr == TRUE], by = c("Region", "Gene")),
                       n_TF_inMR   = uniqueN(TF[in_mr])),
                   by = .(grn_set, cell_type), .SDcols = c("Region", "Gene")]
grn_summary[MR_all[, .(n_MR = .N), by = cell_type], n_MR := i.n_MR, on = "cell_type"]
grn_summary[, pct_MR_in_GRN := round(100 * n_pair_inMR / n_MR, 2)]
print(grn_summary)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 4. GRN ∩ MR: pair 层

GRN_pair <- GRN[in_mr == TRUE, .(n_TF = uniqueN(TF), TFs = paste(sort(unique(TF)), collapse = ",")),
                by = .(grn_set, cell_type, gene_id = Gene, peak_id = Region)]
GRN_MR   <- MR_all[GRN_pair, on = .(cell_type, gene_id, peak_id)]
setcolorder(GRN_MR, c("grn_set", "cell_type", "gene_id", "peak_id", "n_TF", "TFs"))

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 5. 6 个方法在 GRN pairs 中各 identify 多少 (FDR < FDR_CUT); 背景 = 该 CT 全部 MR pairs
#    fold_vs_All = prop_GRN / prop_All; OR / fisher_p = GRN vs 非 GRN 的 2x2 Fisher

MR_sig <- MR_all[, c("cell_type", "gene_id", "peak_id", fdr_cols), with = FALSE]
MR_sig[, (methods) := lapply(.SD, function(x) x < FDR_CUT), .SDcols = fdr_cols]
MR_sig[, (fdr_cols) := NULL]
MR_sig_long <- melt(MR_sig, id.vars = c("cell_type", "gene_id", "peak_id"), variable.name = "method", value.name = "sig")

count_set <- function(set_name) {
  d <- copy(MR_sig_long)
  d[, in_grn := FALSE]
  d[GRN_pair[grn_set == set_name], in_grn := TRUE, on = .(cell_type, gene_id, peak_id)]
  d[, {
    a <- sum(sig & in_grn); b <- sum(!sig & in_grn); c_ <- sum(sig & !in_grn); d_ <- sum(!sig & !in_grn)
    ft <- fisher.test(matrix(c(a, c_, b, d_), 2))                 # 行: GRN / 非GRN; 列: sig / 非sig
    .(n_pairs_GRN = a + b, n_sig_GRN = a, prop_GRN = a / (a + b),
      n_pairs_All = .N,    n_sig_All = a + c_, prop_All = (a + c_) / .N,
      fold_vs_All = (a / (a + b)) / ((a + c_) / .N),
      OR = unname(ft$estimate), fisher_p = ft$p.value)
  }, by = .(cell_type, method)][, grn_set := set_name][]
}

sig_tab <- rbindlist(lapply(c("CT", "CT_HQ"), count_set))
setcolorder(sig_tab, c("grn_set", "cell_type", "method"))
setorder(sig_tab, grn_set, cell_type, method)
print(sig_tab, digits = 3)
print(dcast(sig_tab, cell_type + method ~ grn_set, value.var = c("n_pairs_GRN", "n_sig_GRN", "prop_GRN", "OR")), digits = 3)

###----------------------------------------------------------------------------------------------------------------------------------------------------------------------
### 6. 保存

fwrite(sig_tab,     file.path(OUT_DIR, "GRN_sig_counts_by_method.tsv"), sep = "\t")
fwrite(grn_summary, file.path(OUT_DIR, "GRN_summary.tsv"),              sep = "\t")
fwrite(GRN_MR,      file.path(OUT_DIR, "GRN_pairs_inMR.tsv.gz"),        sep = "\t")
cat("\n全部完成\n")
