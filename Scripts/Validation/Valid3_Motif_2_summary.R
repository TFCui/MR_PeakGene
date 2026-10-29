# 2026-10-09 (2026-10-20 更新: 删掉 TF 表达检查和 Shared/MR_only/Marginal_only/Partial 分组; 新增 MR 特异 motif 表)
# Valid3_Motif 步骤 2: 合并 3 个 CT 的 motif 结果, 提取 MR 特异 motif (流程见 Scripts/Pipeline.md 3.3)
#   显著: p.adjust (BH) < 0.05 且 fold.enrichment > 1
#   MR 特异: IVW / MR_Egger / SMR 任一显著, 且 WNN、Pseudo 都不显著 (CCT 不参与筛选, 只列 FDR)
# 输出: Results/Validation/Valid3_Motif/Motif_sig_3ct.tsv, Motif_MR_specific.tsv

suppressPackageStartupMessages(library(data.table))
options(width = 300)

BASE_DIR <- "/scratch/user/tfcui2025/CIMA"
RES_DIR  <- file.path(BASE_DIR, "Claude_project/Results/Validation/Valid3_Motif")

CELLS      <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
FDR_CUT    <- 0.05
METHODS    <- c("IVW", "MR_Egger", "SMR", "CCT", "WNN", "Pseudo")
MR_SELECT  <- c("IVW", "MR_Egger", "SMR")                          # 用于筛选 MR 特异 motif 的方法
MG_METHODS <- c("WNN", "Pseudo")

### 合并三个细胞类型
Motif_all <- rbindlist(lapply(CELLS, function(ct) fread(file.path(RES_DIR, paste0(ct, "_motif_enrichment.tsv.gz")))))
Motif_all[, `:=`(celltype = factor(celltype, levels = CELLS), method = factor(method, levels = METHODS))]
stopifnot(!anyNA(Motif_all$method), !anyNA(Motif_all$celltype), Motif_all[, .N, by = .(celltype, motif)][, all(N == length(METHODS))])
Motif_all[, sig := p.adjust < FDR_CUT & fold.enrichment > 1]

### 每个 celltype x method 的显著 motif
Motif_sig <- Motif_all[sig == TRUE][, sig := NULL]
Motif_sig[, `:=`(mlog10_padj = -log10(pmax(p.adjust, 1e-300)), method_class = fifelse(method %in% c("IVW", "MR_Egger", "SMR", "CCT"), "MR", "Marginal"))]
setorder(Motif_sig, celltype, method, p.adjust)
fwrite(Motif_sig, file.path(RES_DIR, "Motif_sig_3ct.tsv"), sep = "\t", quote = FALSE)
cat("显著 motif 数 (p.adjust < 0.05 & fold > 1) / 前景 peak 数:\n")
print(dcast(Motif_all[, .(n_fg = n_fg[1], n_sig = sum(sig)), by = .(celltype, method)], method ~ celltype, value.var = c("n_sig", "n_fg")))

### MR 特异 motif: 宽表, 每个 celltype x motif 一行
wide <- dcast(Motif_all, celltype + motif + motif.name ~ method,
              value.var = c("percent.observed", "percent.background", "fold.enrichment", "p.adjust", "sig", "n_fg"))
old_nm <- c(outer(c("percent.observed", "percent.background", "fold.enrichment", "p.adjust"), METHODS, paste, sep = "_"))
new_nm <- c(outer(c("pct_obs", "pct_bg", "fold", "FDR"), METHODS, function(v, m) paste(m, v, sep = "_")))
setnames(wide, old_nm, new_nm)

sig_mr <- as.matrix(wide[, paste0("sig_", MR_SELECT), with = FALSE])
sig_mg <- as.matrix(wide[, paste0("sig_", MG_METHODS), with = FALSE])
MR_spec <- wide[rowSums(sig_mr) >= 1 & rowSums(sig_mg) == 0]
MR_spec[, selected_by := apply(.SD, 1, function(x) paste(MR_SELECT[x], collapse = ";")), .SDcols = paste0("sig_", MR_SELECT)]
MR_spec[, tfs := vapply(motif.name, function(x) paste(toupper(trimws(strsplit(sub("\\(var\\.[0-9]+\\)$", "", x), "::", fixed = TRUE)[[1]])), collapse = ";"), "")]
MR_spec[, min_MR_FDR := do.call(pmin, .SD), .SDcols = paste0(MR_SELECT, "_FDR")]

keep <- c("celltype", "motif", "motif.name", "tfs", "selected_by",
          c(outer(c("pct_obs", "pct_bg", "fold", "FDR"), MR_SELECT, function(v, m) paste(m, v, sep = "_"))),
          paste0(setdiff(METHODS, MR_SELECT), "_FDR"), paste0("n_fg_", METHODS))
setorder(MR_spec, celltype, min_MR_FDR)
MR_spec <- MR_spec[, ..keep]
fwrite(MR_spec, file.path(RES_DIR, "Motif_MR_specific.tsv"), sep = "\t", quote = FALSE)

cat("\nMR 特异 motif 数 (IVW/MR_Egger/SMR 任一显著, WNN/Pseudo 都不显著):\n")
print(MR_spec[, .N, keyby = celltype])
print(MR_spec[, .(celltype, motif, motif.name, selected_by, IVW_FDR = signif(IVW_FDR, 3), MR_Egger_FDR = signif(MR_Egger_FDR, 3),
                  SMR_FDR = signif(SMR_FDR, 3), WNN_FDR = signif(WNN_FDR, 3), Pseudo_FDR = signif(Pseudo_FDR, 3))], nrows = 500)
cat("\n全部完成\n")
