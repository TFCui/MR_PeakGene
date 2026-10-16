# 2026-10-09
# Valid2_SLDSC 步骤 3: 合并 S-LDSC 结果 (流程见 Scripts/Pipeline.md 3.2)
# 移植自 Scripts_reference/D3_merge_plot.R 的 Part 1 (不画图)
# 每行 = CT (3 个 + Combined) x 集合 (MROnly + pool; 2026-10-14 起) x GWAS; 第 1 行 (L2_0) = 我们的 annotation
# 输出: Results/Validation/Valid2_SLDSC/SLDSC_all.tsv

suppressPackageStartupMessages(library(data.table))
options(width = 300)

PROJ    <- "/scratch/user/tfcui2025/CIMA/Claude_project"
OUT_DIR <- file.path(PROJ, "Results/Validation/Valid2_SLDSC")
H2_DIR  <- file.path(OUT_DIR, "h2")

M_TOT  <- 5961159                                                 # 1000G EUR MAF >= 5% SNP 数
CELLS  <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7", "Combined")
SETS   <- c("MROnly", "pool")
TRAITS <- c("AS", "CD", "IBD", "RA", "SLE", "MS", "T1D", "Asthma", "PBC")

# log: 总 h2 和 intercept
read_log <- function(f) {
  x  <- readLines(f)
  h2 <- regmatches(x, regexec("Total Observed scale h2: ([-0-9.e]+) \\(([-0-9.e]+)\\)", x))
  h2 <- h2[lengths(h2) == 3][[1]]
  ic <- regmatches(x, regexec("^Intercept: ([-0-9.e]+) \\(([-0-9.e]+)\\)", x))
  ic <- ic[lengths(ic) == 3][[1]]
  list(h2g = as.numeric(h2[2]), h2g_se = as.numeric(h2[3]), intercept = as.numeric(ic[2]), intercept_se = as.numeric(ic[3]))
}

# 一个 .results 的第 1 行 (CT / set / trait 由循环给出, 不从文件名解析: CT 名里有下划线)
read_one <- function(ct, set, trait) {
  f <- file.path(H2_DIR, ct, sprintf("%s_%s_%s.results", ct, set, trait))
  if (!file.exists(f)) return(NULL)
  r  <- fread(f, nrows = 1)
  lg <- read_log(sub("\\.results$", ".log", f))
  data.table(cell_type = ct, set = set, trait = trait, category = r$Category,
             prop_snp = r$Prop._SNPs,
             prop_h2  = r$Prop._h2,     prop_h2_se = r$Prop._h2_std_error,
             enrich   = r$Enrichment,   enrich_se  = r$Enrichment_std_error, enrich_p = r$Enrichment_p,
             tau      = r$Coefficient,  tau_se     = r$Coefficient_std_error, tau_z   = r[["Coefficient_z-score"]],
             h2g = lg$h2g, h2g_se = lg$h2g_se, intercept = lg$intercept, intercept_se = lg$intercept_se)
}

full <- CJ(cell_type = CELLS, set = SETS, trait = TRAITS, sorted = FALSE)
res  <- rbindlist(lapply(seq_len(nrow(full)), function(i) read_one(full$cell_type[i], full$set[i], full$trait[i])))
cat("results:", nrow(res), "(应为 4 CT x 2 集合 x 9 GWAS = 72)\n")
stopifnot(all(res$category == "L2_0"))
res[, category := NULL]

miss <- full[!res, on = .(cell_type, set, trait)]
cat("缺失组合:", nrow(miss), "\n")
if (nrow(miss)) print(miss)

# 派生统计量
res[, tau_p  := 2 * pnorm(-abs(tau_z))]                           # 双侧
res[, tau_p1 := pnorm(-tau_z)]                                    # 单侧, H1: tau > 0
res[, sd_annot    := sqrt(prop_snp * (1 - prop_snp))]
res[, tau_star    := tau    * sd_annot * M_TOT / h2g]
res[, tau_star_se := tau_se * sd_annot * M_TOT / h2g]             # 忽略 h2g 的不确定性
res[, annot := fifelse(set == "pool", "All", "MROnly")]

res[, `:=`(cell_type = factor(cell_type, levels = CELLS), set = factor(set, levels = SETS), trait = factor(trait, levels = TRAITS))]
setcolorder(res, c("cell_type", "set", "annot", "trait"))
setorder(res, trait, cell_type, set)
fwrite(res, file.path(OUT_DIR, "SLDSC_all.tsv"), sep = "\t")
cat("Saved:", file.path(OUT_DIR, "SLDSC_all.tsv"), "\n")
print(res[, .(trait, cell_type, set, prop_snp, enrich, enrich_p, tau_star, tau_z)], digits = 3, nrows = 200)

if (nrow(miss)) quit(status = 1)
