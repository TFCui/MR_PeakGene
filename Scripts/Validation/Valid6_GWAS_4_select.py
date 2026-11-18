#!/usr/bin/env python3
# 2026-11-09
# Valid6_GWAS_4_select: Fig 4B 例子的证据汇总与筛选 (登录节点, 几十秒)
#   GeneHancer Elite 元件和 ENCODE cCRE 几乎所有候选 peak 都有, 单列为 Q1_generic, 不算突出证据
#   围绕三个问题, 每个 pair 列出 "突出证据" (highlight); 不要求每项都强, 有突出证据即可写
#     Q1 peak 为什么可能有调控功能 (peak 层)
#        EnhancerAtlas / FANTOM5 金标准支持这个 pair; peak 内 SNP 的 FAVOR aPC-Epigenetics-Active >= 20 (全基因组前 1%) 或 H3K27ac_sum >= 100; CAGE enhancer
#     Q2 其中的 SNP 为什么值得关注 (SNP 层)
#        GWAS P <= 1e-10; RegulomeDB ranking 1a-1d / 2a-2c 或 probability >= 0.6; FAVOR aPC-TF >= 20; CADD >= 10;
#        GWAS SNP 是 MR 用到的 SNP (或 r2 >= 0.8 的代理, 最好有); GWAS SNP 本身是该 peak 的 caQTL (p < 1e-5)
#     Q3 这个区域可能调控哪个基因 (gene 层)
#        MR: >= 3 个 MR 方法显著且 min FDR < 1e-3, 或在 >= 2 个 CT 中是 MR-only;
#        3D 接触: 对应 CT 的 HiChIP / PCHiC, 或 4DGenome; GTEx eQTL 金标准; GeneHancer 把元件连到 target (第 1 名, 或依据含 eQTL / C-HiC);
#        GWAS Catalog 中 target 为同一性状的 reported / mapped gene; 文献 (只覆盖已查的 13 个 gene, 见 GWAS_example_literature.tsv)
#   另列: 同一 CT 同一 peak 上 MR 与 marginal 都显著的其他基因 (co_linked_both), 及其中属于 GWAS Catalog 7 个免疫性状基因的 (co_linked_known)
#   入选条件 (qualify): 不在扩展 MHC (chr6:25.5-34 Mb), Q3 >= 2 条, 且 Q1 >= 1 条或 Q2 >= 2 条; 同一 peak 多个 gene 时按 Q3 条数 -> 总条数取 1 个代表 (peak_best)
#   同一 peak 的代表 (peak_best): 入选 pair 中 causal FDR 最小的 (作者偏向 causal 方法更显著的例子)
#   2026-11-10 起加入 Open Targets QTL 共定位 (Valid6_GWAS_3b_opentargets.py): target 基因 H4 >= 0.8 计入 Q3;
#     最终清单另给两个维度: lit_support (基因-疾病文献强度, 人工判断, 见 LIT_SUPPORT) 和 lead_specific (证据是否落在我们的 lead 信号上)
# 输入 (中间表, tmp/Valid6_GWAS/annot/): GWAS_example_candidates.tsv, _snps.tsv, GWAS_example_external_{peak,snp}.tsv, GWAS_example_opentargets.tsv,
#       GWAS_example_literature.tsv (人工整理的文献表, PMID 已用 PubMed / NCBI E-utilities 核对);
#       Results/MR/<CT>_AllPairs_results.tsv.gz (6 个方法的 FDR); tmp/Valid6_GWAS/gwas_catalog/ (GWAS Catalog 2026-10-29)
# 输出: tmp/Valid6_GWAS/annot/GWAS_example_evidence.tsv (全部 61 个 pair 的证据)
#       Results/Validation/Valid6_GWAS/GWAS_Fig4B_examples.tsv (整理后的推荐 / 备选例子; 分档见 TIERS)
# 运行: python3 Scripts/Validation/Valid6_GWAS_4_select.py

import csv, gzip, os, re

PROJ = "/scratch/user/tfcui2025/CIMA/Claude_project"
VDIR = os.path.join(PROJ, "tmp/Valid6_GWAS/annot")                       # 中间表
FINAL = os.path.join(PROJ, "Results/Validation/Valid6_GWAS/GWAS_Fig4B_examples.tsv")
MR_DIR = os.path.join(PROJ, "Results/MR")
GWASCAT = os.path.join(PROJ, "tmp/Valid6_GWAS/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv")
OUT = os.path.join(VDIR, "GWAS_example_evidence.tsv")
CTS = ["Bn_TCL1A", "CD4_Tn_CCR7", "CD8_Tn_CCR7"]
AB = {"Bn_TCL1A": "Bn", "CD4_Tn_CCR7": "CD4 T", "CD8_Tn_CCR7": "CD8 T"}
CAT_TRAITS = {"ankylosing spondylitis", "Crohn disease", "inflammatory bowel disease", "ulcerative colitis", "rheumatoid arthritis",
              "systemic lupus erythematosus", "multiple sclerosis", "type 1 diabetes mellitus"}

# 文献 (只覆盖 Step 2 查过的 13 个 gene; 依据见 GWAS_example_literature.tsv, PMID 均已核对); 只列支持 target gene 的突出证据
LIT = {
    "GPR65":   "lit: GPR65 up in inflamed IBD mucosa, expressed in CD4+ T (35343079); Th17 regulator (26607794)",
    "ATG16L1": "lit: CD GWAS gene (17200669); T-cell Atg16l1 deletion -> colitis (26910010, mouse)",
    "UBA7":    "lit: colon eQTL names UBA7 as IBD-variant target (25766683); ISGylation drives colitis (26515094, mouse)",
    "TMBIM1":  "lit: rs2271543 (r2=1.00 with lead) TMBIM1 sQTL colocalizes with IBD (30975994); IBD-locus SNP is blood eQTL for TMBIM1 (27005424); TMBIM1 degrades TLR4 (28481357)",
    "C5orf56": "lit: CARINH (C5orf56) KO worsens DSS colitis, fine-mapped IBD variant lowers CARINH / IRF1 (37055591); IRF1-AS1 hub of IBD inflammation module (38329124)",
    "GNAS":    "lit: Gnas in CD4+ T required for Th1/Th17 colitis (22326954, mouse)",
    "ERBB3":   "lit: 12q13 genotype sets ERBB3 in PBMC, APC function (20668683)",
    "THBS3":   "lit: blood eQTL coloc with IBD (PPH4 0.96), replicated in FinnGen (41386873)",
    "DR1":     "lit: colon eQTL coloc in 3 cohorts (41688467); blood eQTL coloc (40208878); CD-PD locus labelled DR1 (36687396)",
    "HCN3":    "lit: IBD candidate gene by MAGMA + DEPICT + SMR (41308126); mQTL SMR across IBD / UC / CD (41386873)",
    "ERAP1":   "lit: 5q15 AS variants change ERAP1 expression (29108111); ERAP1-/- worse DSS colitis (35246034, mouse)",
    "MUC1":    "lit: MUC1 overexpressed in IBD mucosa (19025220); costimulatory on activated CD4+ T (22425740); blood SMR names MUC1 CD-causal (37170220); mQTL SMR (41386873)",
    "TTC33":   "lit: CD SNP rs4631223 (r2=0.75 with lead) is eQTL for TTC33 and PTGER4 (39143063)",
}

# 基因-疾病文献强度 (2026-11-10, 依据见 GWAS_example_literature.tsv; 只看文献把 target 基因和疾病联系起来的证据, 不看位点):
#   strong   = 人类遗传证据 + 功能实验 (人或小鼠)
#   moderate = 人类 QTL 共定位 (重复) 或多方法统计排序, 无功能实验
#   weak     = 只有 eQTL 候选或位点归别的基因;  none = 没有基因特异的疾病文献
LIT_SUPPORT = {"C5orf56": "strong", "GPR65": "strong", "ATG16L1": "strong", "GNAS": "strong", "MUC1": "strong", "UBA7": "strong",
               "ERAP1": "strong", "THBS3": "moderate", "DR1": "moderate", "TMBIM1": "moderate", "HCN3": "moderate", "ERBB3": "moderate",
               "TMPRSS3": "weak", "TTC33": "weak", "ZSWIM8": "none", "PLEKHH2": "none"}
# lead_specific: 外部分子证据落在我们的 lead 信号上 = Open Targets target 共定位 H4 >= 0.8, 或文献中的 QTL 变异与 lead r2 >= 0.8 (下表)
LEAD_QTL_LIT = {"TMBIM1": "TMBIM1 sQTL rs2271543 r2=1.00 (30975994)"}
LIT_DONE = set(LIT) | {"DR1", "NKD1", "TMPRSS3", "POLR2E", "TTC33", "ZSWIM8", "PLEKHH2", "HCN3"}

# 写正文时要交代的问题 (依据见文献表)
CAVEAT = {
    "C5orf56": "Cell Res fine-mapped variant rs2188962 has r2=0.24 with our lead rs2548993; Carinh acts in myeloid cells (T-cell Carinh not required in mouse DSS, 37055591); IRF1 also proposed (CRISPRi, 34314424); our data also link IRF1 (MR + marginal)",
    "ERAP1":   "CD locus usually attributed to ERAP2 (21102463, 25766683); Open Targets: lead colocalizes with ERAP2 (H4 0.96), ERAP1 only 0.65; our data also link ERAP2 (MR + marginal)",
    "ATG16L1": "canonical mechanism is coding T300A (17200669); weaker association in East Asians (26192919); no QTL colocalization with CD at this lead (Open Targets)",
    "GPR65":   "rs8005161 in perfect LD with missense I231L (34135065); GWAS SNP is not a caQTL of this peak; Open Targets GPR65 H4 0.68 vs GALC 0.99",
    "UBA7":    "3p21 signal attributed to MST1 R689C (19079170); Open Targets: no UBA7 colocalization at this lead (H4 0.03; MST1R 0.93)",
    "TMBIM1":  "multi-gene eQTL / sQTL locus: PNKD shares the eQTL (27005424); ARPC2, PNKD also colocalize (Open Targets)",
    "THBS3":   "rs2990245 lies in the GBAP1 promoter (30951202); peak links many genes (HCN3, MUC1 also colocalize)",
    "GNAS":    "rs259956 is a ZNF831 missense variant; ZNF831 and EDN3 T-cell eQTLs colocalize as strongly as GNAS (Open Targets H4 0.995); our data also link ZNF831 (MR + marginal)",
    "TMPRSS3": "21q22.3 T1D locus attributed to UBASH3A (28607106); Open Targets: lead colocalizes only with UBASH3A (H4 0.999); our data also link UBASH3A (MR + marginal)",
    "ERBB3":   "RPS26 / SUOX are stronger eQTLs at 12q13 (22144904); Open Targets: SUOX (0.997), RPS26 (0.994), IKZF4 (0.987) colocalize as strongly as ERBB3 (0.990); our data also link RPS26, SUOX (MR + marginal)",
    "DR1":     "no DR1 functional / perturbation data; peak also linked to CCDC18 (promoter), which also colocalizes (Open Targets)",
    "MUC1":    "rs2990245 in GBAP1 promoter (30951202); 1q22 peak links many genes (MUC1 not the top-FDR link in CD4 T)",
    "TTC33":   "5p13.1 locus attributed to PTGER4 (17447842, 23300802); TTC33 only an eQTL candidate (39143063); lead in no QTL credible set (Open Targets)",
    "ZSWIM8":  "lead in LD with PLAU missense rs2227564 (r2=0.64); PLAU enhancer / eQTL (39314475, preprint); Open Targets ZSWIM8 H4 0.19 vs PLAU 0.97; no ZSWIM8 disease literature",
    "PLEKHH2": "2p21 locus attributed to THADA; no PLEKHH2 immune / disease literature; lead in no QTL credible set (Open Targets)",
    "HCN3":    "HCN3-IBD support is statistical prioritization only (41308126, 41386873; no functional data); same peak and lead as THBS3 / MUC1 (gene-dense 1q22, SCAMP3 also implicated)",
}

# 分档 (作者确认后只改这里再重跑): (cell_type, gene) -> "recommended" / "alternate"; 不在字典里的不进最终清单
TIERS = {   # 2026-11-10 第 2 版 (加入 CodeX 文献审阅 + Open Targets 后重选):
    #   recommended     = lit_support strong / moderate + min causal FDR < 0.01 + lead_specific; 同一 peak 只取 1 个 (1q22 取 THBS3: 外部文献专门提名且重复)
    #   alternate       = 基因-疾病文献 strong / moderate, 但 causal FDR >= 0.01、lead 上没有 target 共定位, 或与推荐例子同一 peak
    #   not_recommended = 文献 weak / none 且 lead 上没有 target 共定位 (位点归别的基因)
    ("CD4_Tn_CCR7", "DR1"): "recommended", ("CD8_Tn_CCR7", "GNAS"): "recommended", ("CD8_Tn_CCR7", "TMBIM1"): "recommended",
    ("CD8_Tn_CCR7", "THBS3"): "recommended", ("CD8_Tn_CCR7", "C5orf56"): "recommended",
    ("CD4_Tn_CCR7", "ERBB3"): "recommended",
    ("CD4_Tn_CCR7", "HCN3"): "alternate", ("CD4_Tn_CCR7", "MUC1"): "alternate", ("Bn_TCL1A", "ERAP1"): "alternate",
    ("Bn_TCL1A", "UBA7"): "alternate", ("CD4_Tn_CCR7", "GPR65"): "alternate",
    ("Bn_TCL1A", "ATG16L1"): "alternate",
    ("CD4_Tn_CCR7", "TMPRSS3"): "not_recommended", ("CD8_Tn_CCR7", "TTC33"): "not_recommended",
    ("CD4_Tn_CCR7", "ZSWIM8"): "not_recommended", ("CD4_Tn_CCR7", "PLEKHH2"): "not_recommended",
}



def rd(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def f(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


cand = rd(os.path.join(VDIR, "GWAS_example_candidates.tsv"))
lit_rows = rd(os.path.join(VDIR, "GWAS_example_literature.tsv"))
key = lambda r: (r["cell_type"], r["gene_id"], r["peak_id"])
ext_peak = {key(r): r for r in rd(os.path.join(VDIR, "GWAS_example_external_peak.tsv"))}
ext_snp = {}
for r in rd(os.path.join(VDIR, "GWAS_example_external_snp.tsv")):
    ext_snp.setdefault(key(r), []).append(r)
snps = {}
for r in rd(os.path.join(VDIR, "GWAS_example_candidates_snps.tsv")):
    snps.setdefault(key(r), []).append(r)
assert set(map(key, cand)) == set(ext_peak) == set(ext_snp), "外部注释与候选表的 pair 不一致, 先重跑 Valid6_GWAS_3_external.py"
ot = {key(r): r for r in rd(os.path.join(VDIR, "GWAS_example_opentargets.tsv"))}
assert set(map(key, cand)) == set(ot), "Open Targets 表与候选表的 pair 不一致, 先重跑 Valid6_GWAS_3b_opentargets.py"
OT_H4 = 0.8

# ---- GWAS Catalog: 7 个免疫性状 (P <= 5e-8) 的 reported / mapped gene ----
cat_genes = set()
with open(GWASCAT, newline="") as fh:
    rdr = csv.reader(fh, delimiter="\t")
    h = next(rdr)
    i_tr, i_ml, i_rep, i_map = h.index("MAPPED_TRAIT"), h.index("PVALUE_MLOG"), h.index("REPORTED GENE(S)"), h.index("MAPPED_GENE")
    for row in rdr:
        if len(row) <= max(i_tr, i_ml, i_rep, i_map):
            continue
        if not (set(row[i_tr].split(", ")) & CAT_TRAITS) or (f(row[i_ml], 0) < 7.30103):
            continue
        for g in re.split(r"[,;]| - | x ", row[i_rep] + "," + row[i_map]):
            if g.strip():
                cat_genes.add(g.strip())

# ---- 同一 CT 同一 peak 上 MR 与 marginal 都显著的基因 ----
peaks = {r["peak_id"] for r in cand}
both, fdr_all = {}, {}
for ct in CTS:
    with gzip.open(os.path.join(MR_DIR, ct + "_AllPairs_results.tsv.gz"), "rt") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            if r["peak_id"] not in peaks:
                continue
            fdr_all[(ct, r["gene_id"], r["peak_id"])] = {m + "_fdr": r[m + "_fdr"] for m in ("ivw", "egger", "smr", "cct", "wnn", "pseudo")}
            mr = any(f(r[m + "_fdr"], 1) < 0.05 for m in ("ivw", "egger", "smr", "cct"))
            mg = any(f(r[m + "_fdr"], 1) < 0.05 for m in ("wnn", "pseudo"))
            if mr and mg:
                both.setdefault((ct, r["peak_id"]), []).append(r["gene_id"])

rows = []
for c in cand:
    k = key(c)
    ep, es = ext_peak[k], ext_snp[k]
    q1, q2, q3 = [], [], []
    # ---- Q1 peak ----
    if c["gs_enhancerAtlas"] == "TRUE":
        q1.append("EnhancerAtlas enhancer for target")
    if c["gs_FANTOM5"] == "TRUE":
        q1.append("FANTOM5 enhancer for target")
    q1g = []   # 普遍存在的证据 (几乎所有候选 peak 都有), 单列, 不计入 n_Q1
    els = ep["genehancer_elements"]
    if "Elite" in els and "Enhancer" in els:
        q1g.append("GeneHancer elite enhancer (" + els.split("(")[0] + ")")
    if re.search(r"\((dELS|pELS|PLS)", ep["encode_ccre"]):
        q1g.append("ENCODE cCRE " + "/".join(sorted(set(re.findall(r"\((dELS|pELS|PLS)", ep["encode_ccre"])))))
    epi = max(f(s["favor_apc_epigenetics_active"], 0) for s in es)
    k27 = max(f(s["favor_encodeh3k27ac_sum"], 0) for s in es)
    if epi >= 20 or k27 >= 100:
        q1.append("active chromatin (aPC-EpiActive %.1f, H3K27ac %.0f)" % (epi, k27))
    if any(s["favor_cage_enhancer"] for s in es):
        q1.append("FANTOM CAGE enhancer")
    # ---- Q2 SNP ----
    gp = f(c["lead_gwas_p"])
    if gp is not None and gp <= 1e-10:
        q2.append("GWAS %s P=%.1e" % (c["lead_trait"], gp))
    for s in es:
        rk, pr = s["regulome_rank"], f(s["regulome_prob"], 0)
        if re.match(r"^(1[a-d]|2[a-c])$", rk) or pr >= 0.6:
            q2.append("RegulomeDB %s (%s, p=%.2f)" % (rk, s["rsid"], pr))
        if f(s["favor_apc_transcription_factor"], 0) >= 20:
            q2.append("aPC-TF %.1f (%s)" % (f(s["favor_apc_transcription_factor"]), s["rsid"]))
        if f(s["favor_cadd_phred"], 0) >= 10:
            q2.append("CADD %.1f (%s)" % (f(s["favor_cadd_phred"]), s["rsid"]))
    if c["gwas_is_iv"] == "TRUE":
        q2.append("GWAS SNP is an MR instrument")
    elif f(c["max_r2_gwas_iv"], 0) >= 0.8:
        q2.append("GWAS SNP tagged by MR instrument (r2=%.2f)" % f(c["max_r2_gwas_iv"]))
    if f(c["lead_caqtl_p"], 1) < 1e-5:
        q2.append("GWAS SNP is caQTL of peak (p=%.1e)" % f(c["lead_caqtl_p"]))
    # ---- Q3 gene ----
    if (int(c["n_mr_sig"]) >= 3 and f(c["min_mr_fdr"]) < 1e-3) or int(c["n_ct_mronly"]) >= 2:
        q3.append("MR %s (FDR %.1e; MR-only in %s CT)" % (c["mr_methods_sig"], f(c["min_mr_fdr"]), c["n_ct_mronly"]))
    if c["gs_HiChIP"] == "TRUE":
        q3.append("HiChIP loop (%s)" % AB[c["cell_type"]])
    if c["gs_PCHiC"] == "TRUE":
        q3.append("PCHiC contact (%s)" % AB[c["cell_type"]])
    if c["gs_4DGenome"] == "TRUE":
        q3.append("4DGenome contact")
    if c["gs_GTEx"] == "TRUE":
        q3.append("GTEx eQTL for target")
    if ep["gh_target_connected"] == "True" and (ep["gh_target_rank"] == "1" or re.search("eQTL|C-HiC", ep["gh_target_methods"])):
        q3.append("GeneHancer links element to target (rank %s/%s, %s)" % (ep["gh_target_rank"], ep["gh_n_connected_genes"], ep["gh_target_methods"]))
    if c["gwascat_trait_gene"] == "TRUE":
        q3.append("GWAS Catalog gene for " + c["traits"])
    o = ot[k]
    if f(o["ot_target_h4"], 0) >= OT_H4:
        q3.append("Open Targets QTL coloc with %s GWAS (H4 %.2f, %s)" % (c["traits"], f(o["ot_target_h4"]), o["ot_target_qtl"]))
    if c["gene_id"] in LIT:
        q3.append(LIT[c["gene_id"]])
    co = sorted(set(both.get((c["cell_type"], c["peak_id"]), [])) - {c["gene_id"]})
    co_known = [g for g in co if g in cat_genes]
    rows.append({
        "cell_type": c["cell_type"], "gene_id": c["gene_id"], "peak_id": c["peak_id"], "traits": c["traits"],
        "lead_rsid": c["lead_rsid"], "lead_gwas_p": c["lead_gwas_p"], "distance": c["distance"], "location": c["location"],
        "mr_methods_sig": c["mr_methods_sig"], "min_mr_fdr": c["min_mr_fdr"], "n_ct_mronly": c["n_ct_mronly"],
        "MHC": c["MHC"], "n_genes_peak_mronly": len(c["multi_gene_peak_genes"].split(",")),
        "n_Q1": len(q1), "n_Q2": len(q2), "n_Q3": len(q3), "n_total": len(q1) + len(q2) + len(q3),
        "Q1_peak": "; ".join(q1), "Q1_generic": "; ".join(q1g), "Q2_snp": "; ".join(q2), "Q3_gene": "; ".join(q3),
        "co_linked_both": ",".join(co), "co_linked_known": ",".join(co_known),
        "literature_searched": c["gene_id"] in LIT_DONE,
        **fdr_all[k],
        "min_causal_fdr": c["min_mr_fdr"],
        "gwas_snps": ";".join(sorted({"%s|%s|%.2e|%s" % (x["rsid"], x["trait"], f(x["gwas_p"]), x["pos_hg38"]) for x in snps[k]})),
        "gwas_is_iv": c["gwas_is_iv"], "max_r2_gwas_iv": c["max_r2_gwas_iv"],
        "ot_target_h4": o["ot_target_h4"], "ot_target_qtl": o["ot_target_qtl"], "ot_coloc_genes": o["ot_coloc_genes"],
    })

def in_xmhc(pid):   # 扩展 MHC (chr6:25.5-34.0 Mb, hg38), LD 复杂
    ch, rest = pid.split(":")
    return ch == "chr6" and 25500000 <= int(rest.split("-")[0]) <= 34000000


for r in rows:
    r["xMHC"] = in_xmhc(r["peak_id"])
    r["qualify"] = (not r["xMHC"]) and r["n_Q3"] >= 2 and (r["n_Q1"] >= 1 or r["n_Q2"] >= 2)
best = {}
for r in rows:
    if not r["qualify"]:
        continue
    pk = (r["cell_type"], r["peak_id"])
    if pk not in best or f(r["min_causal_fdr"]) < f(best[pk]["min_causal_fdr"]):
        best[pk] = r
for r in rows:
    r["peak_best"] = r["qualify"] and best.get((r["cell_type"], r["peak_id"])) is r
rows.sort(key=lambda r: (not r["peak_best"], not r["qualify"], -r["n_Q3"], -r["n_total"], f(r["min_mr_fdr"])))

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("Saved:", OUT, len(rows), "rows | qualify:", sum(r["qualify"] for r in rows), "| peak_best:", sum(r["peak_best"] for r in rows))

# ---- 整理后的例子清单 (Results): TIERS 里的 pair, 推荐在前, 各档内按 causal FDR 排 ----
lit_by_gene = {}
for x in lit_rows:
    lit_by_gene.setdefault(x["gene"], []).append("%s: %s [PMID %s]" % (x["item"], x["evidence"], x["pmid"]))
byk = {(r["cell_type"], r["gene_id"]): r for r in rows}
miss = [k for k in TIERS if k not in byk]
assert not miss, "TIERS 中的 pair 不在候选表: %s" % miss
final_cols = ["tier", "cell_type", "gene_id", "peak_id", "traits", "distance", "location",
              "mr_methods_sig", "min_causal_fdr", "ivw_fdr", "egger_fdr", "smr_fdr", "cct_fdr", "wnn_fdr", "pseudo_fdr", "n_ct_mronly",
              "lead_rsid", "lead_gwas_p", "gwas_snps", "gwas_is_iv", "max_r2_gwas_iv",
              "lit_support", "lead_specific", "ot_target_h4", "ot_target_qtl", "ot_coloc_genes",
              "Q1_peak", "Q1_generic", "Q2_snp", "Q3_gene", "co_linked_known", "co_linked_both", "caveat", "literature"]
TIER_ORDER = {"recommended": 0, "alternate": 1, "not_recommended": 2}
assert set(TIERS.values()) <= set(TIER_ORDER) and all(g in LIT_SUPPORT for _, g in TIERS)
final = []
for (ct, g), tier in TIERS.items():
    r = dict(byk[(ct, g)])
    basis = (["Open Targets target H4 %.2f" % f(r["ot_target_h4"])] if f(r["ot_target_h4"], 0) >= OT_H4 else []) + \
            ([LEAD_QTL_LIT[g]] if g in LEAD_QTL_LIT else [])
    r.update(tier=tier, caveat=CAVEAT.get(g, ""), literature=" || ".join(lit_by_gene.get(g, [])),
             lit_support=LIT_SUPPORT[g], lead_specific=("yes: " + "; ".join(basis)) if basis else "no")
    final.append({c: r[c] for c in final_cols})
final.sort(key=lambda r: (TIER_ORDER[r["tier"]], f(r["min_causal_fdr"])))
with open(FINAL, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=final_cols, delimiter="\t")
    w.writeheader()
    w.writerows(final)
print("Saved:", FINAL, len(final), "rows |", {t: sum(r["tier"] == t for r in final) for t in TIER_ORDER})
