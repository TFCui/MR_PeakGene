# 2026-10-20
# Valid3_Motif 步骤 3: MR 特异 motif 表 (论文正文 Table 1; 2026-10-29 前为 S3 Table) —— motif logo + 数值 + 文献注释 合并为 xlsx (流程见 Scripts/Pipeline.md 3.3)
#   logo: JASPAR2020 PFM -> PPM -> ggseqlogo (bits), 沿用 20261006_paper/Scripts/Plot6_val_Motif.sh 第 3 节的画法
#   2026-10-23: 论文表只保留 MR 一侧更强的 motif: IVW / MR_Egger / SMR 任一 FDR < MR_FDR_TABLE (0.01) 且 fold > 1, 并去掉没有明确 Brief interpretation 的行 ("No clear ...");
#     Valid3 步骤 2 的 Motif_MR_specific*.tsv (47 行, FDR < 0.05) 不动, 筛选只在建表阶段做; 被去掉的行写进 xlsx 的 Excluded sheet
#   2026-10-29: 改为每个 CT 取 MR 一侧最小 FDR 的前 TOP_N (2) 个 motif (跳过 "No clear ..."), 共 6 行, 进正文;
#     每个 motif 只列一篇文献: 改写后的一句解释 + PubMed 链接取自 Motif_table_main.tsv (人工整理, PMID 已核对), 合并为 "Brief interpretation & link" 列, 不再有 Reference 列
#   高亮: IVW / MR_Egger / SMR / CCT 为 FDR < MR_FDR_TABLE (0.05) 且 fold > 1, WNN / Pseudo 为 FDR < 0.05 且 fold > 1 (黄底粗体); 表内 FDR 写作 p.adj
# 输入: Results/Validation/Valid3_Motif/Motif_MR_specific.tsv, Motif_MR_specific_literature.tsv, Motif_table_main.tsv, <CT>_motif_enrichment.tsv.gz
# 输出: Results/Validation/Valid3_Motif/logos/<ID>.{png,svg}, Results/Figures/Figures_20261102/Table_1.xlsx

suppressPackageStartupMessages({
  library(data.table)
  library(JASPAR2020)
  library(TFBSTools)
  library(ggplot2)
  library(ggseqlogo)
  library(openxlsx)
})

PROJ     <- "/scratch/user/tfcui2025/CIMA/Claude_project"
RES_DIR  <- file.path(PROJ, "Results/Validation/Valid3_Motif")
LOGO_DIR <- file.path(RES_DIR, "logos")
OUT_XLSX <- file.path(PROJ, "Results/Figures/Figures_20261102/Table_1.xlsx")
dir.create(LOGO_DIR, recursive = TRUE, showWarnings = FALSE)

CELLS     <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
FDR_CUT   <- 0.05          # WNN / Pseudo 的显著阈值 (MR 特异 motif 要求它们不显著)
MR_FDR_TABLE <- 0.05       # 论文表里 MR 一侧 (IVW / MR_Egger / SMR, 以及 CCT) 的高亮阈值
TOP_N     <- 2             # 每个 CT 进表的 motif 数
MR3       <- c("IVW", "MR_Egger", "SMR")
MR3_LABEL <- c(IVW = "IVW", MR_Egger = "MR-Egger", SMR = "SMR")
OTHER     <- c("CCT", "WNN", "Pseudo")
OTHER_LAB <- c(CCT = "CCT", WNN = "WNN", Pseudo = "Pseudobulk")

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

### 1. 读表并合并
num <- fread(file.path(RES_DIR, "Motif_MR_specific.tsv"))
lit <- fread(file.path(RES_DIR, "Motif_MR_specific_literature.tsv"), quote = "")
stopifnot(nrow(num) == nrow(lit), !anyDuplicated(num[, .(celltype, motif)]),
          setequal(paste(num$celltype, num$motif, num$motif.name), paste(lit$celltype, lit$motif, lit$motif.name)))
tb <- merge(num, lit[, !"motif.name"], by = c("celltype", "motif"), sort = FALSE)
tb <- tb[num[, .(celltype, motif)], on = .(celltype, motif)]          # 行顺序同 Motif_MR_specific.tsv
msg("Rows: ", nrow(tb), " | unique motifs: ", uniqueN(tb$motif))

# CCT / WNN / Pseudo 的 fold 从原始富集表取 (只用于判断高亮)
raw <- rbindlist(lapply(CELLS, function(ct) fread(file.path(RES_DIR, paste0(ct, "_motif_enrichment.tsv.gz")),
                                                  select = c("celltype", "method", "motif", "fold.enrichment", "p.adjust"))))
fe_other <- dcast(raw[method %in% OTHER], celltype + motif ~ method, value.var = c("fold.enrichment", "p.adjust"))
tb <- fe_other[tb, on = .(celltype, motif)]
for (m in OTHER) stopifnot(all.equal(tb[[paste0("p.adjust_", m)]], tb[[paste0(m, "_FDR")]]) == TRUE)

### 1b. 论文表的筛选: 每个 CT 内按 MR 一侧最小 FDR (IVW / MR_Egger / SMR 中 fold > 1 的) 排序, 跳过没有明确解释的行, 取前 TOP_N 个; 被去掉的行记录原因
tb[, mr_min_FDR := do.call(pmin, lapply(MR3, function(m) fifelse(get(paste0(m, "_fold")) > 1, get(paste0(m, "_FDR")), 1)))]
tb[, has_interp := !grepl("^No clear", brief_interpretation)]
tb[, rank_ct := frank(fifelse(has_interp, mr_min_FDR, Inf), ties.method = "first"), by = celltype]
tb[, keep := has_interp & rank_ct <= TOP_N]
excl <- tb[keep == FALSE][, reason := fifelse(has_interp, sprintf("not in top %d of the cell type", TOP_N), "no clear interpretation")]
tb <- tb[keep == TRUE][order(match(celltype, CELLS), mr_min_FDR)]
stopifnot(nrow(tb) + nrow(excl) == nrow(num), all(tb$mr_min_FDR < FDR_CUT))
# 正文表用的一句解释 + 单篇文献 (人工整理; 行顺序必须与筛选结果一致)
sel <- fread(file.path(RES_DIR, "Motif_table_main.tsv"), colClasses = c(pmid = "character"), quote = "")
stopifnot(identical(paste(tb$celltype, tb$motif, tb$motif.name), paste(sel$celltype, sel$motif, sel$motif.name)),
          all(mapply(grepl, sel$pmid, tb$reference_url, fixed = TRUE)))   # 选中的 PMID 在原来的文献列表中
tb[, `:=`(table_pmid = sel$pmid, table_brief = sel$brief, table_link = paste0("https://pubmed.ncbi.nlm.nih.gov/", sel$pmid, "/"))]
cat(sprintf("\n===== 每个 CT 前 %d 个 (MR min FDR, 跳过 No clear)\n", TOP_N))
print(tb[, .(celltype, motif, motif.name, mr_min_FDR = signif(mr_min_FDR, 2), table_pmid)])
msg("Kept: ", nrow(tb), " | excluded: ", nrow(excl))

### 2. Motif logo (每个 motif ID 一张)
ids <- unique(tb$motif)
for (id in ids) {
  pfm <- getMatrixByID(JASPAR2020, ID = id)
  stopifnot(!is.null(pfm), name(pfm) == tb[motif == id, motif.name][1])
  mat <- Matrix(pfm); ppm <- sweep(mat, 2, colSums(mat), "/")
  p <- ggseqlogo(ppm, method = "bits") + ggtitle(paste0(name(pfm), " (", id, ")")) +
       theme(plot.title = element_text(size = 10, hjust = 0.5))
  ggsave(file.path(LOGO_DIR, paste0(id, ".png")), p, width = 3, height = 1.6, dpi = 300, bg = "white")
  ggsave(file.path(LOGO_DIR, paste0(id, ".svg")), p, width = 3, height = 1.6, bg = "white", device = svglite::svglite)
}
stopifnot(all(file.exists(file.path(LOGO_DIR, paste0(ids, ".png")))), all(file.exists(file.path(LOGO_DIR, paste0(ids, ".svg")))))
msg("Logos: ", length(ids), " PNG + SVG in ", LOGO_DIR)

### 3. Sheet Motif_table
wb <- createWorkbook()
sh <- "Motif_table"; addWorksheet(wb, sh, gridLines = FALSE)
modifyBaseFont(wb, fontSize = 10, fontName = "Arial")

COL <- list(ct = 1, tf = 2, id = 3, logo = 4)
mr_col <- setNames(5 + 4 * (seq_along(MR3) - 1), MR3)              # 每个 MR 方法 4 列: % obs, % bg, Fold, FDR
oth_col <- setNames(17:19, OTHER)
COL$rel <- 20; COL$brief <- 21
NCOL <- 21
R_TITLE <- 1; R_H1 <- 2; R_H2 <- 3; R0 <- 4                          # 数据从第 4 行开始
rows <- R0 + seq_len(nrow(tb)) - 1

st_head  <- createStyle(textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE,
                        fgFill = "#E7E6E6", border = "TopBottomLeftRight", borderColour = "#A6A6A6")
st_cell  <- createStyle(valign = "center", wrapText = TRUE, border = "TopBottomLeftRight", borderColour = "#D9D9D9")
st_ctr   <- createStyle(valign = "center", halign = "center", border = "TopBottomLeftRight", borderColour = "#D9D9D9")
st_pct   <- createStyle(numFmt = "0.0", valign = "center", halign = "center", border = "TopBottomLeftRight", borderColour = "#D9D9D9")
st_fold  <- createStyle(numFmt = "0.00", valign = "center", halign = "center", border = "TopBottomLeftRight", borderColour = "#D9D9D9")
st_fdr   <- createStyle(numFmt = "GENERAL", valign = "center", halign = "center", border = "TopBottomLeftRight", borderColour = "#D9D9D9")
st_hl    <- createStyle(fgFill = "#FFF2CC", textDecoration = "bold")
st_ctcol <- createStyle(textDecoration = "bold", valign = "center", halign = "center", wrapText = TRUE,
                        border = "TopBottomLeftRight", borderColour = "#A6A6A6")

writeData(wb, sh, "Table 1. MR-specific transcription factor motifs in peaks supported by Mendelian randomization (top two per cell type), with motif logos, enrichment statistics and literature evidence",
          startRow = R_TITLE, startCol = 1)
addStyle(wb, sh, createStyle(textDecoration = "bold", fontSize = 12), rows = R_TITLE, cols = 1)

# 两行表头
single <- c(ct = "Cell type", tf = "TF motif", id = "Motif ID", logo = "Motif logo",
            rel = "Immune relevance", brief = "Brief interpretation & link")
for (k in names(single)) {
  writeData(wb, sh, single[[k]], startRow = R_H1, startCol = COL[[k]])
  mergeCells(wb, sh, cols = COL[[k]], rows = R_H1:R_H2)
}
for (m in MR3) {
  writeData(wb, sh, MR3_LABEL[[m]], startRow = R_H1, startCol = mr_col[[m]])
  mergeCells(wb, sh, cols = mr_col[[m]] + 0:3, rows = R_H1)
  writeData(wb, sh, t(c("% obs", "% bg", "Fold", "p.adj")), startRow = R_H2, startCol = mr_col[[m]], colNames = FALSE)
}
writeData(wb, sh, "p.adj (other methods)", startRow = R_H1, startCol = oth_col[[1]])
mergeCells(wb, sh, cols = oth_col, rows = R_H1)
writeData(wb, sh, t(unname(OTHER_LAB[OTHER])), startRow = R_H2, startCol = oth_col[[1]], colNames = FALSE)
addStyle(wb, sh, st_head, rows = R_H1:R_H2, cols = 1:NCOL, gridExpand = TRUE)

# 数据
writeData(wb, sh, tb$celltype,   startRow = R0, startCol = COL$ct)
writeData(wb, sh, tb$motif.name, startRow = R0, startCol = COL$tf)
writeData(wb, sh, tb$motif,      startRow = R0, startCol = COL$id)
n_hl <- 0L
for (m in MR3) {
  c0 <- mr_col[[m]]
  writeData(wb, sh, tb[[paste0(m, "_pct_obs")]],          startRow = R0, startCol = c0)
  writeData(wb, sh, tb[[paste0(m, "_pct_bg")]],           startRow = R0, startCol = c0 + 1)
  writeData(wb, sh, tb[[paste0(m, "_fold")]],             startRow = R0, startCol = c0 + 2)
  writeData(wb, sh, signif(tb[[paste0(m, "_FDR")]], 3),   startRow = R0, startCol = c0 + 3)
  addStyle(wb, sh, st_pct,  rows = rows, cols = c0 + 0:1, gridExpand = TRUE)
  addStyle(wb, sh, st_fold, rows = rows, cols = c0 + 2)
  addStyle(wb, sh, st_fdr,  rows = rows, cols = c0 + 3)
  hl <- which(tb[[paste0(m, "_FDR")]] < MR_FDR_TABLE & tb[[paste0(m, "_fold")]] > 1)
  if (length(hl)) addStyle(wb, sh, st_hl, rows = rows[hl], cols = c0 + 3, stack = TRUE)
  n_hl <- n_hl + length(hl)
}
for (m in OTHER) {
  c0 <- oth_col[[m]]
  writeData(wb, sh, signif(tb[[paste0(m, "_FDR")]], 3), startRow = R0, startCol = c0)
  addStyle(wb, sh, st_fdr, rows = rows, cols = c0)
  thr <- if (m == "CCT") MR_FDR_TABLE else FDR_CUT
  hl <- which(tb[[paste0(m, "_FDR")]] < thr & tb[[paste0("fold.enrichment_", m)]] > 1)
  if (length(hl)) addStyle(wb, sh, st_hl, rows = rows[hl], cols = c0, stack = TRUE)
  n_hl <- n_hl + length(hl)
}
writeData(wb, sh, tb$immune_relevance,     startRow = R0, startCol = COL$rel)
writeData(wb, sh, paste0(tb$table_brief, "\n", tb$table_link), startRow = R0, startCol = COL$brief)
addStyle(wb, sh, st_ctr,  rows = rows, cols = c(COL$id, COL$logo, COL$rel), gridExpand = TRUE)
addStyle(wb, sh, st_cell, rows = rows, cols = c(COL$tf, COL$brief), gridExpand = TRUE)

# Cell type 列: 同一 CT 纵向合并
for (ct in CELLS) {
  rr <- rows[tb$celltype == ct]
  mergeCells(wb, sh, cols = COL$ct, rows = range(rr))
  addStyle(wb, sh, st_ctcol, rows = rr, cols = COL$ct)
}

# logo: 1.5 x 0.8 in, 放在 logo 列每一行
LOGO_W <- 1.8; LOGO_H <- 0.96   # 1.5 x 0.8 -> 1.8 x 0.96 in (原 3 x 1.6 in 的 PNG 按比例缩放)
for (i in seq_len(nrow(tb))) {
  insertImage(wb, sh, file.path(LOGO_DIR, paste0(tb$motif[i], ".png")),
              startRow = rows[i], startCol = COL$logo, width = LOGO_W, height = LOGO_H, units = "in")
}

# 行高 / 列宽
setRowHeights(wb, sh, rows = rows, heights = LOGO_H * 72 + 6)
setRowHeights(wb, sh, rows = R_H1:R_H2, heights = c(18, 30))
setColWidths(wb, sh, cols = 1:NCOL,
             widths = c(13, 17, 10, 27, rep(c(7, 7, 6, 9), 3), 9, 9, 10, 11, 60))
freezePane(wb, sh, firstActiveRow = R0, firstActiveCol = COL$logo + 1)

# 脚注
notes <- c(
  "MR-specific motif: significant (p.adj < 0.05 and fold enrichment > 1) under at least one of IVW, MR-Egger or SMR, and not significant (p.adj >= 0.05) under either marginal method (WNN multiome metacells, pseudobulk). Of the 47 MR-specific motifs, the two with the smallest MR p.adj in each cell type that have a clear literature-based interpretation are listed; the others are in the 'Excluded' sheet. CCT is listed for reference and was not used for selection.",
  "Foreground: peaks of links with FDR < 0.1 for the given method. Background: 50,000 peaks drawn from the candidate-peak universe, matched on GC content and accessibility (Signac MatchRegionStats). Motifs: JASPAR2020 CORE vertebrates (746 motifs); enrichment by Signac FindMotifs (hypergeometric test).",
  "% obs / % bg: percentage of foreground / background peaks containing the motif. Fold = % obs / % bg. p.adj: Benjamini-Hochberg across the 746 motifs within each cell type and method. Shaded bold p.adj: p.adj < 0.05 and fold > 1.",
  "Immune relevance: Strong = direct experimental evidence in the matched cell type (B cells for Bn_TCL1A, CD4+ or CD8+ T cells) or in a related immune disease; Moderate = immune evidence in related cell types or contexts; Weak = indirect evidence; None = no immune evidence found.",
  "Brief interpretation & link: one-sentence interpretation based on a single paper, followed by its PubMed link. For each motif, the literature on the corresponding transcription factor(s) was searched in PubMed for functions in the same or a related immune cell type, or in immune-mediated disease, and the paper was checked against its PubMed record (title, year and abstract).",
  "Motif enrichment does not establish TF binding at these peaks. Full interpretations and all references are in the 'Details' sheet.")
r_note <- max(rows) + 2
writeData(wb, sh, notes, startRow = r_note, startCol = 1)
for (j in seq_along(notes)) mergeCells(wb, sh, cols = 1:NCOL, rows = r_note + j - 1)
addStyle(wb, sh, createStyle(wrapText = TRUE, valign = "top", fontSize = 9), rows = r_note + seq_along(notes) - 1, cols = 1)
setRowHeights(wb, sh, rows = r_note + seq_along(notes) - 1, heights = 30)

### 4. Sheet Details
det <- tb[, c(list(celltype = celltype, motif = motif, motif.name = motif.name, tfs = tfs, selected_by = selected_by,
                   immune_relevance = immune_relevance, evidence_context = evidence_context,
                   table_pmid = table_pmid, table_brief = table_brief, brief_interpretation = brief_interpretation, biological_interpretation = biological_interpretation,
                   reference_url = gsub(" \\| ", "\n", reference_url), reference_note = gsub(" \\| ", "\n", reference_note),
                   consistency_with_old_Fig6 = consistency_with_old),
              .SD), .SDcols = patterns("^n_fg_")]
addWorksheet(wb, "Details")
writeData(wb, "Details", det, headerStyle = st_head)
addStyle(wb, "Details", createStyle(wrapText = TRUE, valign = "top"), rows = 1 + seq_len(nrow(det)), cols = seq_len(ncol(det)), gridExpand = TRUE)
setColWidths(wb, "Details", cols = seq_len(ncol(det)),
             widths = c(13, 10, 17, 12, 14, 11, 30, 10, 40, 40, 80, 42, 60, 40, rep(9, ncol(det) - 14)))
freezePane(wb, "Details", firstActiveRow = 2, firstActiveCol = 4)

### 4b. Sheet Excluded: Valid3 的 47 个 MR 特异 motif (FDR < 0.05) 中没有进入论文表的行及原因 (41 行)
exc <- excl[, .(celltype, motif, motif.name, selected_by, mr_min_FDR = signif(mr_min_FDR, 3), immune_relevance, brief_interpretation, reason)]
setorder(exc, celltype, mr_min_FDR)
addWorksheet(wb, "Excluded")
writeData(wb, "Excluded", exc, headerStyle = st_head)
addStyle(wb, "Excluded", createStyle(wrapText = TRUE, valign = "top"), rows = 1 + seq_len(nrow(exc)), cols = seq_len(ncol(exc)), gridExpand = TRUE)
setColWidths(wb, "Excluded", cols = seq_len(ncol(exc)), widths = c(13, 10, 20, 16, 12, 12, 70, 34))
freezePane(wb, "Excluded", firstActiveRow = 2)

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
msg("Saved: ", OUT_XLSX)

### 5. 检查
exp_hl <- sum(sapply(MR3, function(m) sum(tb[[paste0(m, "_FDR")]] < MR_FDR_TABLE & tb[[paste0(m, "_fold")]] > 1))) +
          sum(sapply(OTHER, function(m) sum(tb[[paste0(m, "_FDR")]] < (if (m == "CCT") MR_FDR_TABLE else FDR_CUT) & tb[[paste0("fold.enrichment_", m)]] > 1)))
stopifnot(n_hl == exp_hl)
back <- read.xlsx(OUT_XLSX, sheet = sh, startRow = R_H2, colNames = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)
back <- back[seq_len(nrow(tb)) + 1, ]                                # 去掉第 2 行表头
stopifnot(identical(back[[COL$id]], tb$motif), isTRUE(all.equal(as.numeric(back[[mr_col[["IVW"]]]]), tb$IVW_pct_obs)))
med <- unzip(OUT_XLSX, list = TRUE)$Name
msg("Embedded images: ", sum(grepl("^xl/media/", med)), " | highlighted FDR cells: ", n_hl,
    " | rows by CT: ", paste(tb[, .N, by = celltype][, paste(celltype, N)], collapse = ", "))
print(tb[, .N, keyby = .(celltype, immune_relevance)])
msg("Done.")
