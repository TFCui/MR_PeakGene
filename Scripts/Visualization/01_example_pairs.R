# Example pairs: 每个细胞类型 24 个 (gene, peak), 每 6 个一类, 共 4 类
# 来自 20261006_paper (Plot1 / Plot2_scripts.R / Plot8 三处副本), 这里只保留一份, Plot1 (检查) 和 Plot2 (画图) 共用
# 4 类的筛选条件 (paper Plot1, 按 ivw_p 排序取前 6):
#   1 All methods significant        : 6 个方法 FDR 全 < 0.05
#   2 Not significant in Marginal    : pseudo_p >= 0.05 | wnn_p >= 0.05
#   3 Opposite direction             : ivw_fdr < 0.05 且 (wnn_fdr < 0.05 且 sign(wnn_rho) != sign(ivw_beta) 或 pseudo 同理)
#   4 Significant only in MR         : ivw_fdr < 0.05 & egger_fdr < 0.05 & wnn_fdr >= 0.05 & pseudo_fdr >= 0.05

CATS <- c("All methods significant", "Not significant in Marginal Multiomics / Pseudobulk",
          "Opposite direction (MR vs. marginal)", "Significant only in MR")

EXAMPLE_PAIRS <- list(
  Bn_TCL1A = data.table(
    gene_id = c("SIDT2", "TAGLN", "ERAP2", "NBPF26", "ANKRD36B", "ZNF107",
                "PDK1", "LAMTOR5", "ACCS", "RRP7A", "LAMTOR5", "RRP7A",
                "ERAP2", "NBPF26", "NBPF26", "YBEY", "RPS26", "PDK1",
                "ADD3", "RITA1", "BLK", "CATSPER2", "AFF1", "MCM3AP"),
    peak_id = c("chr11:117185715-117186216", "chr11:117185715-117186216", "chr5:96885894-96886395",
                "chr1:120913924-120914425", "chr2:97589361-97589862", "chr7:64686502-64687003",
                "chr2:172510056-172510557", "chr1:110384972-110385473", "chr11:44065924-44066425",
                "chr22:42537036-42537537", "chr1:110384416-110384917", "chr22:42551722-42552223",
                "chr5:96885894-96886395", "chr1:120913924-120914425", "chr1:120844072-120844573",
                "chr21:46253648-46254149", "chr12:56041553-56042054", "chr2:172513783-172514284",
                "chr10:110153435-110153936", "chr12:113135593-113136094", "chr8:11417989-11418490",
                "chr15:43789350-43789851", "chr4:86998008-86998509", "chr21:46253648-46254149")),
  CD4_Tn_CCR7 = data.table(
    gene_id = c("ACYP2", "LCLAT1", "TMEM8A", "ERAP2", "RNASET2", "ANKRD36B",
                "EMB", "LCLAT1", "COG5", "RNASET2", "RNASET2", "RNASET2",
                "ERAP2", "ERAP2", "NBPF26", "NBPF26", "MFSD14C", "ZNF429",
                "ARPC1B", "ZNF117", "PLEC", "ASL", "MC1R", "SNRPC"),
    peak_id = c("chr2:53994688-53995189", "chr2:30447755-30448256", "chr16:380051-380552",
                "chr5:96885894-96886395", "chr6:166951108-166951609", "chr2:97589361-97589862",
                "chr5:50489402-50489903", "chr2:30446497-30446998", "chr7:107560897-107561398",
                "chr6:166935293-166935794", "chr6:166952979-166953480", "chr6:166942727-166943228",
                "chr5:96885894-96886395", "chr5:96896259-96896760", "chr1:120913924-120914425",
                "chr1:120844072-120844573", "chr9:96777893-96778394", "chr19:21474531-21475032",
                "chr7:99376204-99376705", "chr7:65081104-65081605", "chr8:143966529-143967030",
                "chr7:65770609-65771110", "chr16:89837638-89838139", "chr6:34791707-34792208")),
  CD8_Tn_CCR7 = data.table(
    gene_id = c("ACYP2", "LCLAT1", "DIP2A", "DIP2A", "PTK7", "PRMT2",
                "RNASET2", "COG5", "RNASET2", "CARS", "POLI", "FGL2",
                "ERAP2", "NBPF26", "CORO1B", "ZNF429", "RPS26", "ERAP2",
                "ARPC1B", "MC1R", "S100B", "BANF1", "BUD31", "BUD31"),
    peak_id = c("chr2:53994688-53995189", "chr2:30447755-30448256", "chr21:46550539-46551040",
                "chr21:46561121-46561622", "chr6:43142078-43142579", "chr21:46550539-46551040",
                "chr6:166959199-166959700", "chr7:107560897-107561398", "chr6:166952979-166953480",
                "chr11:3006844-3007345", "chr18:54269292-54269793", "chr7:77199591-77200092",
                "chr5:96885894-96886395", "chr1:120913924-120914425", "chr11:67446256-67446757",
                "chr19:21474531-21475032", "chr12:56041553-56042054", "chr5:96896259-96896760",
                "chr7:99376204-99376705", "chr16:89784313-89784814", "chr21:46633535-46634036",
                "chr11:66001982-66002483", "chr7:99380850-99381351", "chr7:99376204-99376705"))
)

# 每 6 个 = 1 类; 同一 pair 出现在多类时合并类别标签, 只画一次 (Plot2 用)
prep_examples <- function(dt) {
  dt <- copy(dt)[, `:=`(ord = .I, category = rep(CATS, each = 6, length.out = .N))]
  dt[, .(category = paste(unique(category), collapse = " / "), ord = min(ord)), by = .(gene_id, peak_id)][order(ord)]
}

# 每个 pair 在给定 MR 表下是否满足自己所属类别的条件
check_examples <- function(ex, mr) {
  ex <- copy(ex)[, `:=`(ord = .I, cat_id = rep(1:4, each = 6, length.out = .N))]
  d  <- mr[ex, on = .(gene_id, peak_id)]
  fdr6 <- c("ivw_fdr", "egger_fdr", "smr_fdr", "wnn_fdr", "pseudo_fdr", "cct_fdr")
  d[, in_table := !is.na(ivw_p)]
  d[, pass := fcase(
    !in_table, NA,
    cat_id == 1, rowSums(as.matrix(.SD) < 0.05) == 6,
    cat_id == 2, pseudo_p >= 0.05 | wnn_p >= 0.05,
    cat_id == 3, ivw_fdr < 0.05 & ((wnn_fdr < 0.05 & sign(wnn_rho) != sign(ivw_beta)) |
                                    (pseudo_fdr < 0.05 & sign(pseudo_rho) != sign(ivw_beta))),
    cat_id == 4, ivw_fdr < 0.05 & egger_fdr < 0.05 & wnn_fdr >= 0.05 & pseudo_fdr >= 0.05), .SDcols = fdr6]
  d[order(ord), c("ord", "cat_id", "gene_id", "peak_id", "in_table", "pass", fdr6), with = FALSE]
}
