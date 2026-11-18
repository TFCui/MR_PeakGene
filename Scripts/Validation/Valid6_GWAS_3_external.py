#!/usr/bin/env python3
# 2026-11-09
# Valid6_GWAS_3_external: Fig 4B 候选 pair 的外部数据库注释 (网上查询, 只能在登录节点跑; 计算节点无网络)
#   输入: tmp/Valid6_GWAS/annot/GWAS_example_candidates.tsv + GWAS_example_candidates_snps.tsv (Valid6_GWAS_2_annot.R)
#   范围: GWAS_example_candidates.tsv 的全部 pair (2026-11-09 起; 最初只查 13 个 gene 的短名单), 每个 CT x pair 一行
#   每个 peak (hg38):
#     ENCODE cCRE      UCSC API track encodeCcreCombined (与 peak 重叠的 cCRE 及类别 PLS / pELS / dELS / CA-CTCF ...)
#     GeneHancer       UCSC API track geneHancerRegElements (与 peak 重叠的 GH 元件, 类型, Elite)
#                      + geneHancerInteractions (这些元件连到的基因及 score; target gene 是否在内、score、在连到的基因中的名次)
#   每个 peak 内的 GWAS SNP (P < 5e-8):
#     FAVOR            api.genohub.org/v1/rsids/<rsid> (aPC-Epigenetics-Active / aPC-TF / CADD phred / ENCODE H3K27ac、H3K4me1 / CAGE enhancer /
#                      super-enhancer / FAVOR 自带的 GeneHancer 注释)
#     RegulomeDB       regulomedb.org/api/search (ranking, probability, blood 组织特异分数)
#   原始返回存 tmp/Valid6_GWAS/external/ (JSON 缓存, 存在即跳过)
# 输出: tmp/Valid6_GWAS/annot/GWAS_example_external_peak.tsv  (每个 pair 一行: cCRE / GeneHancer)
#       tmp/Valid6_GWAS/annot/GWAS_example_external_snp.tsv   (每个 pair x SNP 一行: FAVOR / RegulomeDB)
# 运行 (登录节点): python3 Scripts/Validation/Valid6_GWAS_3_external.py

import csv, json, os, time, urllib.request, urllib.parse

PROJ = "/scratch/user/tfcui2025/CIMA/Claude_project"
VDIR = os.path.join(PROJ, "tmp/Valid6_GWAS/annot")   # 中间表目录 (输入和输出都在这里)
CACHE = os.path.join(PROJ, "tmp/Valid6_GWAS/external")
os.makedirs(CACHE, exist_ok=True)

FAVOR_KEYS = ["genecode_comprehensive_category", "apc_epigenetics_active", "apc_transcription_factor", "apc_conservation_v2",
              "cadd_phred", "encodeh3k27ac_sum", "encodeh3k4me1_sum", "encodeh3k4me3_sum", "encode_dnase_sum",
              "cage_enhancer", "cage_promoter", "super_enhancer", "remap_overlap_tf", "genehancer"]


def get_json(url, key):
    path = os.path.join(CACHE, key + ".json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "python-urllib"})
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.loads(r.read().decode())
            with open(path, "w") as f:
                json.dump(d, f)
            time.sleep(0.5)
            return d
        except Exception as e:
            print("  retry", attempt + 1, url, e)
            time.sleep(3)
    return None


def ucsc(track, chrom, start, end):
    url = "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=%s;chrom=%s;start=%d;end=%d" % (track, chrom, start, end)
    d = get_json(url, "ucsc_%s_%s_%d_%d" % (track, chrom, start, end))
    return (d or {}).get(track, [])


def parse_peak(pid):
    chrom, rest = pid.split(":")
    s, e = rest.split("-")
    return chrom, int(s), int(e)


cand = {(r["cell_type"], r["gene_id"]): r for r in csv.DictReader(open(os.path.join(VDIR, "GWAS_example_candidates.tsv")), delimiter="\t")}
snps = list(csv.DictReader(open(os.path.join(VDIR, "GWAS_example_candidates_snps.tsv")), delimiter="\t"))
SHORTLIST = sorted(cand.keys(), key=lambda k: int(cand[k]["rank"]))   # 全部候选, 按 rank

peak_rows, snp_rows = [], []
for ct, gene in SHORTLIST:
    r = cand[(ct, gene)]
    pid = r["peak_id"]
    chrom, s, e = parse_peak(pid)
    print("==", ct, gene, pid)

    # ---- ENCODE cCRE ----
    cc = ucsc("encodeCcreCombined", chrom, s, e)
    ccre = "; ".join("%s(%s,z=%.2f)" % (x["name"], x["ccre"], float(x["zScore"])) for x in cc) or "none"

    # ---- GeneHancer 元件 + 连到的基因 ----
    els = ucsc("geneHancerRegElements", chrom, s, e)
    el_ids = {x["name"] for x in els}
    el_str = "; ".join("%s(%s,%s)" % (x["name"], x["elementType"], x["eliteness"]) for x in els) or "none"
    inter = [x for x in ucsc("geneHancerInteractions", chrom, s, e) if x["geneHancerIdentifier"] in el_ids]
    conn = {}
    for x in inter:
        g, sc = x["geneName"], float(x["value"])
        if g not in conn or sc > conn[g][0]:
            conn[g] = (sc, x["geneAssociationMethods"])
    ranked = sorted(conn.items(), key=lambda kv: -kv[1][0])
    tgt = conn.get(gene)
    tgt_rank = [g for g, _ in ranked].index(gene) + 1 if tgt else None
    peak_rows.append({
        "cell_type": ct, "gene_id": gene, "peak_id": pid, "rank_in_candidates": r["rank"],
        "encode_ccre": ccre,
        "genehancer_elements": el_str,
        "gh_target_connected": bool(tgt),
        "gh_target_score": "%.2f" % tgt[0] if tgt else "",
        "gh_target_methods": tgt[1] if tgt else "",
        "gh_target_rank": tgt_rank or "",
        "gh_n_connected_genes": len(conn),
        "gh_top_genes": "; ".join("%s(%.1f)" % (g, v[0]) for g, v in ranked[:6]),
    })

    # ---- 每个 SNP: FAVOR + RegulomeDB ----
    rs_list = []
    for x in snps:
        if x["cell_type"] == ct and x["gene_id"] == gene and x["rsid"] not in [y[0] for y in rs_list]:
            rs_list.append((x["rsid"], x["chr"], x["pos_hg38"]))
    for rsid, schr, spos in rs_list:
        trs = sorted({"%s:%.1e" % (x["trait"], float(x["gwas_p"])) for x in snps
                      if x["cell_type"] == ct and x["gene_id"] == gene and x["rsid"] == rsid})
        fav = get_json("https://api.genohub.org/v1/rsids/" + rsid, "favor_" + rsid) or []
        fav = [v for v in fav if str(v.get("position")) == spos] or fav
        f0 = fav[0] if fav else {}
        reg = get_json("https://regulomedb.org/api/search?" + urllib.parse.urlencode({"regions": rsid, "genome": "GRCh38"}), "regulome_" + rsid) or {}
        v0 = (reg.get("variants") or [{}])[0]
        rs = v0.get("regulome_score", {}) or {}
        row = {"cell_type": ct, "gene_id": gene, "peak_id": pid, "rsid": rsid, "pos_hg38": "%s:%s" % (schr, spos),
               "gwas": ",".join(trs), "is_lead": rsid == r["lead_rsid"],
               "regulome_rank": rs.get("ranking", "") or ("not_scored" if v0 else "not_found"),   # 返回里没有 regulome_score = RegulomeDB 未打分
               "regulome_prob": rs.get("probability", ""),
               "regulome_blood": (rs.get("tissue_specific_scores") or {}).get("blood", ""),
               "regulome_features": ",".join(k for k, val in (v0.get("features") or {}).items() if val is True)}
        for k in FAVOR_KEYS:
            val = f0.get(k)
            row["favor_" + k] = "" if val is None else (("%.2f" % val) if isinstance(val, float) else str(val))
        snp_rows.append(row)


def write(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("Saved:", path, len(rows), "rows")


write(os.path.join(VDIR, "GWAS_example_external_peak.tsv"), peak_rows)
write(os.path.join(VDIR, "GWAS_example_external_snp.tsv"), snp_rows)
