# Table 5 (论文 S4 Table): S-LDSC (Fig 4A) 和 Fig 4B (位点轨道图) 用到的 7 个 GWAS 性状汇总 (只读汇总, 登录节点, 几秒)
# 口径:
#   - 性状 / 文件: Scripts/Validation/Valid2_SLDSC_2_ldsc_slurm.sh 的 TRAITS / TFILES; 7 个展示性状 = Plot10_Combined.R 的 SL_TRAITS
#     (S-LDSC 实际跑了 9 个, Asthma / PBC 不展示, 不列入本表)
#   - N / SNP 数: 直接读 LDSC 输入 sumstats (selected_sumstats/*.sumstats.gz) 的 N 列和行数; 不引用论文报告的样本量
#     AS / CD 的 cases / controls 来自 munge 日志 GWAS_data/onk1k_GWAS/{AS,CD}.log (--N-cas / --N-con)
#   - h2 (obs) 和 se: Valid2_SLDSC/SLDSC_all.tsv, Combined 的 pool (All peaks) 行
#   - Fig 4B: Valid6_GWAS/GWAS_MRonly_overlap.tsv (strong MR-only link 且 peak 内有该性状 P < 5e-8 的 SNP), 只统计这 7 个性状;
#     图中展示的基因取自 Plot10_Combined.R 的 PAGES (用正则读出, 避免写两份)
#   - 来源 / PMID: 脚本里人工整理 (PMID 已用 PubMed 核对: 题目、期刊、年份)
# 输出: Results/Figures/Figures_20261102/Table_S4.xlsx (单 sheet) + Table_S4.tsv
# 运行: Rscript Scripts/Visualization/Table5_GWASTraits.R

source("/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/00_common.R")
suppressMessages({ library(openxlsx); library(scales) })

GWAS_DIR <- file.path(BASE, "R_package/S_LDSC/GWAS_data/selected_sumstats")
LOG_DIR  <- file.path(BASE, "R_package/S_LDSC/GWAS_data/onk1k_GWAS")
OUT_DIR  <- file.path(PROJ, "Results/Figures/Figures_20261102")
fmt <- function(x) comma(x)

# ---- 1. 性状 (顺序 = Plot10_Combined.R 的 SL_TRAITS) ----
src_plot10 <- readLines(file.path(PROJ, "Scripts/Visualization/Plot10_Combined.R"))
SL_TRAITS  <- eval(parse(text = sub("\\s*#.*$", "", sub("^SL_TRAITS <- ", "", grep("^SL_TRAITS <- ", src_plot10, value = TRUE)))))
stopifnot(identical(SL_TRAITS, c("AS", "CD", "IBD", "RA", "SLE", "MS", "T1D")))

TR <- data.table(
  trait = SL_TRAITS,
  full  = c("Ankylosing spondylitis", "Crohn's disease", "Inflammatory bowel disease", "Rheumatoid arthritis",
            "Systemic lupus erythematosus", "Multiple sclerosis", "Type 1 diabetes"),
  file  = c("CUSTOM.AS", "CUSTOM.CD", "PASS.Inflammatory_Bowel_Disease.deLange2017", "PASS.Rheumatoid_Arthritis.Saevarsdottir2022",
            "PASS.Lupus.Bentham2015", "PASS_Multiple_sclerosis", "PASS.Type_1_Diabetes.Chiou2021"),
  source = c("Pan-UK Biobank", "International IBD Genetics Consortium (IIBDGC)", "de Lange et al. 2017, Nat Genet",
             "Saevarsdottir et al. 2022, Ann Rheum Dis", "Bentham et al. 2015, Nat Genet",
             "International Multiple Sclerosis Genetics Consortium, 2011, Nature (Sawcer et al.)", "Chiou et al. 2021, Nature"),
  ref = c("https://pan.ukbb.broadinstitute.org/", "https://www.ibdgenetics.org/", "PMID 28067908", "PMID 35470158", "PMID 26502338",
          "PMID 21833088; ImmunoBase (protected data)", "PMID 34012112"),
  access = c("Pan-UK Biobank website", "IIBDGC website", "Zenodo 10515792 (sumstats_indep107)", "Zenodo 10515792 (sumstats_indep107)",
             "Zenodo 10515792 (sumstats_indep107)", "ImmunoBase protected download (approval required)", "Zenodo 10515792 (sumstats_indep107)"))
stopifnot(identical(TR$file[1:2], c("CUSTOM.AS", "CUSTOM.CD")))

# ---- 2. LDSC 输入文件: SNP 数 / N ----
# 文件里有 Z / N 缺失的行 (AS、CD 约 3 万行), LDSC 读入时会丢掉; 这里只数完整行, 并与 h2 日志的 "Read summary statistics for X SNPs" 核对
rd_gwas <- function(f) { g <- fread(file.path(GWAS_DIR, paste0(f, ".sumstats.gz")), select = c("SNP", "Z", "N")); g[!is.na(Z) & !is.na(N)] }
log_nsnp <- function(tr) as.numeric(sub(".*for ([0-9]+) SNPs.*", "\\1", grep("Read summary statistics for", readLines(file.path(VAL_DIR, "Valid2_SLDSC/h2/Combined", paste0("Combined_pool_", tr, ".log"))), value = TRUE)[1]))
st <- rbindlist(lapply(seq_len(nrow(TR)), function(i) { g <- rd_gwas(TR$file[i])
  data.table(file = TR$file[i], n_snp = nrow(g), n_min = min(g$N), n_med = median(g$N), n_max = max(g$N)) }))
TR <- merge(TR, st, by = "file", sort = FALSE)
stopifnot(all(TR$n_snp == unname(sapply(TR$trait, log_nsnp))))   # 与 S-LDSC 日志一致
TR[, N_txt := ifelse(n_min == n_max, fmt(n_max), sprintf("%s (median; range %s\u2013%s)", fmt(n_med), fmt(n_min), fmt(n_max)))]
# AS / CD: cases / controls (munge 日志); N = cases + controls 核对
cc <- function(tr) { l <- readLines(file.path(LOG_DIR, paste0(tr, ".log")))
  num <- function(k) as.numeric(sub(".*--N-(con|cas) +([0-9.]+).*", "\\2", grep(k, l, value = TRUE)[1])); c(cas = num("--N-cas"), con = num("--N-con")) }
for (tr in c("AS", "CD")) { x <- cc(tr); stopifnot(!anyNA(x), x[["cas"]] + x[["con"]] == TR[trait == tr, n_max])
  TR[trait == tr, N_txt := sprintf("%s (%s cases / %s controls)", fmt(n_max), fmt(x[["cas"]]), fmt(x[["con"]]))] }

# ---- 3. S-LDSC h2 (Combined, All peaks) ----
SLD <- fread(file.path(VAL_DIR, "Valid2_SLDSC/SLDSC_all.tsv"))
H2  <- SLD[cell_type == "Combined" & set == "pool" & trait %in% SL_TRAITS, .(trait, h2g, h2g_se)]
stopifnot(nrow(H2) == 7)
TR <- merge(TR, H2, by = "trait", sort = FALSE)
TR[, h2_txt := sprintf("%.4f (%.4f)", h2g, h2g_se)]

# ---- 4. Fig 4B: strong MR-only link + peak 内 GWAS P < 5e-8 SNP ----
OV <- fread(file.path(VAL_DIR, "Valid6_GWAS/GWAS_MRonly_overlap.tsv"))
stopifnot(all(OV$gwas_p < 5e-8))
LK <- unique(OV[trait %in% SL_TRAITS, .(cell_type, gene_id, peak_id, trait)])
CT_SH <- c(Bn_TCL1A = "Bn", CD4_Tn_CCR7 = "CD4 T", CD8_Tn_CCR7 = "CD8 T")
pg <- regmatches(paste(src_plot10[76:77], collapse = ""), gregexpr('c\\("[A-Za-z0-9_]+", "[A-Za-z0-9]+"\\)', paste(src_plot10[76:77], collapse = "")))[[1]]
SEL <- rbindlist(lapply(pg, function(s) { v <- strsplit(gsub('c\\(|\\)|"| ', "", s), ",")[[1]]; data.table(cell_type = v[1], gene_id = v[2]) }))
stopifnot(nrow(SEL) == 5, all(SEL$cell_type %in% CELLS))
shown <- LK[SEL, on = .(cell_type, gene_id), nomatch = NULL]   # Fig 4B 的 pair 里属于这 7 个性状的
f4b <- LK[, .(n_links = .N, n_genes = uniqueN(gene_id), n_peaks = uniqueN(peak_id),
              by_ct = paste(sprintf("%s %d", CT_SH[CELLS], sapply(CELLS, function(ct) sum(cell_type == ct))), collapse = " / ")), by = trait]
f4b <- merge(data.table(trait = SL_TRAITS), f4b, by = "trait", all.x = TRUE, sort = FALSE)
f4b[is.na(n_links), `:=`(n_links = 0L, n_genes = 0L, n_peaks = 0L, by_ct = "–")]
sh <- shown[, .(shown = paste(sprintf("%s (%s)", unique(gene_id), CT_SH[unique(cell_type)]), collapse = ", ")), by = trait]
f4b <- merge(f4b, sh, by = "trait", all.x = TRUE, sort = FALSE); f4b[is.na(shown), shown := "–"]
TR <- merge(TR, f4b, by = "trait", sort = FALSE)
# 核对: 9 个性状时的 pair 数 = GWAS_pair_counts.tsv 的 MR_only n_pairs_GWAS
PC <- fread(file.path(VAL_DIR, "Valid6_GWAS/GWAS_pair_counts.tsv"))[category == "MR_only"]
chk <- unique(OV[, .(cell_type, gene_id, peak_id)])[, .N, by = cell_type]
stopifnot(all(PC[match(chk$cell_type, cell_type), n_pairs_GWAS] == chk$N))

TR[, order := match(trait, SL_TRAITS)]; setorder(TR, order)
TR[, `:=`(Trait = sprintf("%s (%s)", full, trait), `In Fig 4A (S-LDSC)` = "Yes")]
out <- TR[, .(Trait, `Study / data source` = source, Reference = ref, `Data access` = access, `S-LDSC input file` = paste0(file, ".sumstats.gz"),
              `N (S-LDSC input file)` = N_txt, `SNPs (S-LDSC input file)` = fmt(n_snp), `h2 (obs) (SE)` = h2_txt,
              `In Fig 4A (S-LDSC)`)]   # 2026-11-05: 去掉 Fig 4B 的 link 列 (n_links / n_genes / n_peaks / by_ct / shown 仍在上面算好并核对, 不输出)
print(out, width = 300)
fwrite(out, file.path(OUT_DIR, "Table_S4.tsv"), sep = "\t")

# ---- xlsx ----
st_title <- createStyle(fontSize = 14, textDecoration = "bold")
st_head  <- createStyle(fontSize = 11, textDecoration = "bold", fontColour = "#FFFFFF", fgFill = "#1F3A5F", halign = "center", valign = "center", wrapText = TRUE, border = "TopBottomLeftRight", borderColour = "#BBBBBB")
st_cell  <- createStyle(fontSize = 11, valign = "center", wrapText = TRUE, border = "TopBottomLeftRight", borderColour = "#BBBBBB")
st_ctr   <- createStyle(fontSize = 11, halign = "center", valign = "center", wrapText = TRUE, border = "TopBottomLeftRight", borderColour = "#BBBBBB")
wb <- createWorkbook(); sh_n <- "Table_S4"; addWorksheet(wb, sh_n, gridLines = FALSE)
writeData(wb, sh_n, "S4 Table. GWAS traits used in the S-LDSC analysis (Fig 4A) and in the GWAS-overlap tracks (Fig 4B)", startRow = 1, startCol = 1); addStyle(wb, sh_n, st_title, rows = 1, cols = 1)
nc <- ncol(out); writeData(wb, sh_n, out, startRow = 3, startCol = 1, headerStyle = st_head)
rr <- 4:(3 + nrow(out)); addStyle(wb, sh_n, st_cell, rows = rr, cols = 1:nc, gridExpand = TRUE)
addStyle(wb, sh_n, st_ctr, rows = rr, cols = 6:9, gridExpand = TRUE); setRowHeights(wb, sh_n, rows = 3, heights = 48); setRowHeights(wb, sh_n, rows = rr, heights = 48)
foot <- c("S-LDSC: each GWAS was analyzed with the 97 baselineLD v2.2 annotations plus the MR-only peak annotation or the all-candidate-peak annotation (HapMap3 SNPs, 1000 Genomes Phase 3 European LD reference); h2 (obs) is the observed-scale SNP heritability estimated by S-LDSC (All peaks model; SE in parentheses). N and the number of SNPs are those of the summary statistics read by S-LDSC (HapMap3 SNPs; SNPs with missing Z or N are dropped). For T1D the sample size varies across SNPs (median and range given).",
          "The multiple sclerosis summary statistics (ImmunoBase) are protected data that require approval from the IMSGC.")
writeData(wb, sh_n, foot, startRow = 4 + nrow(out), startCol = 1)
for (j in seq_along(foot)) { r <- 3 + nrow(out) + j; mergeCells(wb, sh_n, cols = 1:nc, rows = r); addStyle(wb, sh_n, createStyle(wrapText = TRUE, valign = "top", fontSize = 9), rows = r, cols = 1); setRowHeights(wb, sh_n, rows = r, heights = 40) }
setColWidths(wb, sh_n, cols = 1:nc, widths = c(34, 36, 30, 32, 44, 30, 16, 16, 12)); freezePane(wb, sh_n, firstActiveRow = 4, firstActiveCol = 2)
saveWorkbook(wb, file.path(OUT_DIR, "Table_S4.xlsx"), overwrite = TRUE)
msg("Saved: ", file.path(OUT_DIR, "Table_S4.xlsx"), " + Table_S4.tsv")
