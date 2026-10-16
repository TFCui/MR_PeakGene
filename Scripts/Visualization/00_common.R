# 可视化共用设置: 路径 / 细胞类型 / 方法名与配色 / 主题 / 读 MR 结果
# 用法: 在每个 Plot 脚本开头 source(file.path(dirname(this_file), "00_common.R"))
# 方法名、顺序、配色与 20261006_paper 完全一致 (见 Scripts/Visualization.md)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
options(width = 300)

# ---------- 路径 ----------
PROJ    <- "/scratch/user/tfcui2025/CIMA/Claude_project"
BASE    <- "/scratch/user/tfcui2025/CIMA"
MR_DIR  <- file.path(PROJ, "Results/MR")
VAL_DIR <- file.path(PROJ, "Results/Validation")
FIG_DIR <- file.path(PROJ, "Results/Figures/Figures_20261019")
TMP_DIR <- file.path(PROJ, "tmp")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------- 细胞类型 ----------
CELLS     <- c("Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7")
CELL_DESC <- c(Bn_TCL1A = "Naive B cells", CD4_Tn_CCR7 = "Naive CD4+ T cells", CD8_Tn_CCR7 = "Naive CD8+ T cells")

# ---------- 方法: 显示名 / 数据 key / 配色 ----------
METHODS    <- c("Causal SMR", "Causal IVW", "Causal Egger", "Causal Joint", "Marginal Multiomics", "Marginal Pseudobulk")
METHOD_KEY <- setNames(c("smr", "ivw", "egger", "cct", "wnn", "pseudo"), METHODS)
pal_method <- c(`Causal IVW` = "#E64B35", `Causal Egger` = "#3C5488", `Causal SMR` = "#F0A14A",
                `Causal Joint` = "#8491B4", `Marginal Multiomics` = "#4DBBD5", `Marginal Pseudobulk` = "#00A087")
pal_key    <- setNames(pal_method[METHODS], METHOD_KEY)   # 按小写 key 取色
lab_key    <- setNames(METHODS, METHOD_KEY)               # 小写 key -> 显示名
COL_CONC <- "#D5DDE5"; COL_DISC <- "#7E6148"; COL_INK <- "#2F2F2F"

# ---------- 主题 ----------
theme_paper <- function(base = 10) theme_classic(base_size = base) + theme(
  text = element_text(color = "black"), axis.text = element_text(color = "black", size = base - 1),
  axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
  plot.title = element_text(face = "bold", size = base + 1, hjust = 0), plot.subtitle = element_text(size = base - 1, color = "grey30"),
  legend.title = element_blank(), legend.text = element_text(size = base - 1), legend.key.size = unit(0.35, "cm"),
  plot.tag = element_text(face = "bold", size = base + 3))

# ---------- 读 MR 结果 (已删 NA、ivw_sigma2 <= 10、gene 内 BH; 不再过滤) ----------
read_mr <- function(ct) fread(file.path(MR_DIR, paste0(ct, "_AllPairs_results.tsv.gz")))

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
