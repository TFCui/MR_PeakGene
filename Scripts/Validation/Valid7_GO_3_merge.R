# 2026-10-15
# Valid7_GO 第 3 步: 合并观测 GSEA 与随机结果, 算经验 P
# 经验 P (每个 gene set x term): (1 + #{随机 NES >= 观测 NES}) / (1 + 该 term 被检验到的随机次数); emp_fdr = 经验 P 在组内做 BH
# 输入: Results/Validation/Valid7_GO/GO_GSEA_MROnly.tsv.gz (第 1 步), tmp/Valid7_GO/perm/<set>_perm_*.rds (第 2 步)
# 输出: 覆盖 GO_GSEA_MROnly.tsv.gz, 加 emp_p, emp_fdr, n_exceed, n_null, n_perm, null_NES_mean
# 环境变量: N_PERM (默认 5000, 用于检查随机次数); OUT_FILE (默认覆盖上面的文件; 冒烟测试时指向 tmp/)

suppressPackageStartupMessages(library(data.table))
options(width = 300)

PROJ     <- "/scratch/user/tfcui2025/CIMA/Claude_project"
OBS_FILE <- file.path(PROJ, "Results/Validation/Valid7_GO/GO_GSEA_MROnly.tsv.gz")
PERM_DIR <- file.path(PROJ, "tmp/Valid7_GO/perm")
OUT_FILE <- Sys.getenv("OUT_FILE", OBS_FILE)
N_PERM   <- as.integer(Sys.getenv("N_PERM", "5000"))
gsets    <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7", "Combined")

obs_all <- fread(OBS_FILE)
obs_all[, intersect(c("emp_p", "emp_fdr", "n_exceed", "n_null", "n_perm", "null_NES_mean"), names(obs_all)) := NULL]   # 重跑时去掉旧的经验 P 列

out <- rbindlist(lapply(gsets, function(s) {
  obs  <- obs_all[gene_set == s]
  fs   <- list.files(PERM_DIR, pattern = sprintf("^%s_perm_[0-9]+-[0-9]+\\.rds$", s), full.names = TRUE)
  if (!length(fs)) { warning(sprintf("%s: 没有随机结果, 跳过", s)); return(NULL) }
  null <- rbindlist(lapply(fs, readRDS))
  n_perm <- uniqueN(null$perm)
  stopifnot(!anyDuplicated(null, by = c("perm", "ID")))
  if (n_perm != N_PERM) warning(sprintf("%s: 随机 %d 次 (期望 %d)", s, n_perm, N_PERM))
  e <- merge(obs[, .(ID, NES_obs = NES)], null[, .(ID, NES)], by = "ID")[
    , .(n_null = .N, n_exceed = sum(NES >= NES_obs), null_NES_mean = mean(NES)), by = ID]
  obs <- merge(obs, e, by = "ID", all.x = TRUE)
  obs[is.na(n_null), `:=`(n_null = 0L, n_exceed = 0L)]
  obs[, emp_p := (1 + n_exceed) / (1 + n_null)]
  obs[, emp_fdr := p.adjust(emp_p, "BH")]
  obs[, n_perm := n_perm]
  cat(sprintf("%-12s term %d | 随机 %d 次 | n_null 中位 %d (最小 %d) | FDR < 0.1: %d | emp_p < 0.05: %d (期望约 %d) | emp_p < 0.01: %d | emp_fdr < 0.1: %d | p 与 emp_p 的 Spearman %.2f | 两者都 < 0.05: %d\n",
              s, nrow(obs), n_perm, as.integer(median(obs$n_null)), min(obs$n_null), sum(obs$p.adjust < 0.1), sum(obs$emp_p < 0.05),
              round(0.05 * nrow(obs)), sum(obs$emp_p < 0.01), sum(obs$emp_fdr < 0.1),
              cor(obs$pvalue, obs$emp_p, method = "spearman"), sum(obs$pvalue < 0.05 & obs$emp_p < 0.05)))
  obs
}), use.names = TRUE)

setcolorder(out, c("gene_set", "n_genes", "ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue",
                   "emp_p", "emp_fdr", "n_exceed", "n_null", "n_perm", "null_NES_mean"))
out <- out[order(match(gene_set, gsets), pvalue)]

cat("\n每组经验 P 最小的 10 个 term\n")
print(out[order(emp_p, pvalue), head(.SD, 10), by = gene_set,
          .SDcols = c("Description", "setSize", "NES", "null_NES_mean", "pvalue", "p.adjust", "emp_p", "emp_fdr")][
  , `:=`(Description = substr(Description, 1, 60), NES = round(NES, 2), null_NES_mean = round(null_NES_mean, 2),
         pvalue = signif(pvalue, 2), p.adjust = signif(p.adjust, 2), emp_p = signif(emp_p, 2), emp_fdr = signif(emp_fdr, 2))][], nrows = 200)

fwrite(out, OUT_FILE, sep = "\t")
cat("\n写出:", OUT_FILE, "|", nrow(out), "行\n")
