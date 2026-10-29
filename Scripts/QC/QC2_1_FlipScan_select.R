# QC2_1_FlipScan_select: 从全基因组 flip 扫描结果中挑 "Egger 不 flip 显著、flip 后不显著" 的候选 pair (登录节点, 轻量)
# 输入: Results/QC/2_FlipPoints/<CT>_Egger_FlipScan.tsv.gz (QC2_1_FlipScan.R 输出)
# 候选: egger_fdr_noflip < 0.05 且 egger_fdr_flip >= 0.1; 按 -log10 p 的落差 (不 flip - flip) 排序
# 输出: Results/QC/2_FlipPoints/FlipScan_candidates.tsv (每个 CT 前 TOP_N 个)
#       Results/QC/2_FlipPoints/FlipScan_candidates_nsnp50.tsv (同样的规则, 另加 n_snp >= MIN_SNP_CAND; Figure S2 "差别大" 的例子从这里挑)
# 运行: Rscript Scripts/QC/QC2_1_FlipScan_select.R

library(data.table)
options(width = 250)
PROJ  <- "/scratch/user/tfcui2025/CIMA/Claude_project"
DIR   <- file.path(PROJ, "Results/QC/2_FlipPoints")
CELLS <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
SIG <- 0.05; NS <- 0.1; TOP_N <- 20; MIN_SNP_CAND <- 50

res <- rbindlist(lapply(CELLS, function(ct) fread(file.path(DIR, paste0(ct, "_Egger_FlipScan.tsv.gz")))[, cell_type := ct]))[status == "ok"]

cat("===== 每个 CT: Egger 显著 (p.adj < 0.05) 不 flip x flip\n")
print(dcast(res[, .N, by = .(cell_type, noflip = fifelse(egger_fdr_noflip < SIG, "noflip_sig", "noflip_ns"),
                             flip = fifelse(egger_fdr_flip < SIG, "flip_sig", "flip_ns"))], cell_type ~ noflip + flip, value.var = "N", fill = 0))

cand <- res[egger_fdr_noflip < SIG & egger_fdr_flip >= NS]
cand[, drop_log10p := -log10(egger_p_noflip) + log10(egger_p_flip)]
setorder(cand, cell_type, -drop_log10p)
cat(sprintf("\n===== 候选 (不 flip p.adj < %g, flip p.adj >= %g): ", SIG, NS)); print(cand[, .N, by = cell_type])

top <- cand[, head(.SD, TOP_N), by = cell_type]
fwrite(top, file.path(DIR, "FlipScan_candidates.tsv"), sep = "\t")
cat("\n===== 每个 CT 前 10 个\n")
print(top[, head(.SD, 10), by = cell_type][, .(cell_type, gene_id, peak_id, n_snp, n_neg_x,
  ivw_fdr = signif(ivw_fdr, 2), b_noflip = signif(egger_beta_noflip, 3), p_noflip = signif(egger_p_noflip, 2), fdr_noflip = signif(egger_fdr_noflip, 2),
  b_flip = signif(egger_beta_flip, 3), p_flip = signif(egger_p_flip, 2), fdr_flip = signif(egger_fdr_flip, 2),
  a_noflip = signif(egger_alpha_noflip, 2), a_flip = signif(egger_alpha_flip, 2))])
cat("\nSaved: ", file.path(DIR, "FlipScan_candidates.tsv"), "\n")

# 另存 SNP 数多的候选 (n_snp >= MIN_SNP_CAND): SNP 太少的 pair 画成散点图说服力不够
cand50 <- cand[n_snp >= MIN_SNP_CAND]
cat(sprintf("\n===== 候选 (n_snp >= %d): ", MIN_SNP_CAND)); print(cand50[, .N, by = cell_type])
top50 <- cand50[, head(.SD, TOP_N), by = cell_type]
fwrite(top50, file.path(DIR, "FlipScan_candidates_nsnp50.tsv"), sep = "\t")
cat("Saved: ", file.path(DIR, "FlipScan_candidates_nsnp50.tsv"), "\n")
