#!/usr/bin/env python3
# 2026-11-10
# Valid6_GWAS_3b_opentargets: Fig 4B 候选 pair 的 lead GWAS SNP 在 Open Targets Platform 中的 QTL 共定位 (登录节点, 需联网; 计算节点无网络)
#   输入: tmp/Valid6_GWAS/annot/GWAS_example_candidates.tsv (Valid6_GWAS_2_annot.R): lead_rsid, traits, gene_id
#   做法 (Open Targets GraphQL API v4, https://api.platform.opentargets.org/api/v4/graphql):
#     1. rsID -> Open Targets variant ID (search)
#     2. 该变异所在的全部 QTL credible set (eQTL / sQTL / tuQTL / pQTL / sc-eQTL ..., 不含 GWAS)
#     3. 每个 QTL credible set 的 colocalisation 中, 对方为该 pair 性状 GWAS (按 traitFromSource 关键词) 的最大 H4
#     4. 每个基因取所有 QTL credible set 中的最大 H4
#   注意: 只看包含 lead SNP 本身的 QTL credible set; 基因名按 Open Targets approvedSymbol (C5orf56 = CARINH)
#   原始返回缓存 tmp/Valid6_GWAS/external/opentargets/ (JSON, 存在即跳过)
# 输出: tmp/Valid6_GWAS/annot/GWAS_example_opentargets.tsv (每个 pair 一行)
#   ot_variant_id, ot_n_qtl_credsets, ot_target_h4 / ot_target_qtl (target 基因的最大 H4 及来源), ot_coloc_genes (H4 >= 0.8 的基因, 按 H4 排)
# 运行 (登录节点): python3 Scripts/Validation/Valid6_GWAS_3b_opentargets.py

import csv, hashlib, json, os, time, urllib.request

PROJ = "/scratch/user/tfcui2025/CIMA/Claude_project"
VDIR = os.path.join(PROJ, "tmp/Valid6_GWAS/annot")
CACHE = os.path.join(PROJ, "tmp/Valid6_GWAS/external/opentargets")
OUT = os.path.join(VDIR, "GWAS_example_opentargets.tsv")
API = "https://api.platform.opentargets.org/api/v4/graphql"
os.makedirs(CACHE, exist_ok=True)

H4_MIN = 0.8
ALIAS = {"C5orf56": "CARINH"}                                    # 项目基因名 -> Open Targets approvedSymbol
TRAIT_KW = {"CD": ["crohn", "inflammatory bowel", "ibd"],
            "IBD": ["crohn", "inflammatory bowel", "ibd", "ulcerative colitis"],
            "T1D": ["type 1 diabetes", "type i diabetes"],
            "RA": ["rheumatoid arthritis"],
            "AS": ["ankylosing spondylitis"], "SLE": ["lupus"], "MS": ["multiple sclerosis"]}


def gql(query):
    path = os.path.join(CACHE, hashlib.md5(query.encode()).hexdigest() + ".json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    body = json.dumps({"query": query}).encode()
    for attempt in range(3):
        try:
            req = urllib.request.Request(API, data=body, headers={"Content-Type": "application/json", "User-Agent": "python-urllib"})
            with urllib.request.urlopen(req, timeout=120) as r:
                d = json.loads(r.read().decode())
            if "errors" in d:
                raise RuntimeError(d["errors"])
            with open(path, "w") as f:
                json.dump(d, f)
            time.sleep(0.2)
            return d
        except Exception as e:
            print("  retry", attempt + 1, e)
            time.sleep(3)
    raise RuntimeError("Open Targets query failed: " + query[:120])


def qtl_coloc(rsid, kws):
    hits = gql('query { search(queryString:"%s", entityNames:["variant"]) { hits { id } } }' % rsid)["data"]["search"]["hits"]
    if not hits:
        return None, 0, {}
    vid = hits[0]["id"]
    cs_rows, page = [], 0
    while True:                                                                  # 分页 (每页 500)
        cs = gql('query { variant(variantId:"%s") { credibleSets(page:{index:%d,size:500}) { count rows { studyLocusId '
                 'study { studyType biosample { biosampleName } target { approvedSymbol } } } } } }' % (vid, page))["data"]["variant"]["credibleSets"]
        cs_rows += cs["rows"]
        page += 1
        if len(cs_rows) >= cs["count"] or not cs["rows"]:
            break
    assert len(cs_rows) == cs["count"], "credible set 数不一致: " + rsid
    best, n_qtl = {}, 0
    for r in cs_rows:
        st = r["study"]
        if st["studyType"] == "gwas" or not st.get("target"):
            continue
        n_qtl += 1
        g = st["target"]["approvedSymbol"]
        src = "%s:%s" % (st["studyType"], (st.get("biosample") or {}).get("biosampleName") or "NA")
        co = gql('query { credibleSet(studyLocusId:"%s") { colocalisation(page:{index:0,size:500}) { rows { h4 '
                 'otherStudyLocus { study { traitFromSource } } } } } }' % r["studyLocusId"])["data"]["credibleSet"]["colocalisation"]["rows"]
        for x in co:
            t = (x["otherStudyLocus"]["study"]["traitFromSource"] or "").lower()
            if x["h4"] is not None and any(k in t for k in kws) and x["h4"] > best.get(g, (0, ""))[0]:
                best[g] = (x["h4"], src)
    return vid, n_qtl, best


cand = list(csv.DictReader(open(os.path.join(VDIR, "GWAS_example_candidates.tsv")), delimiter="\t"))
rows = []
for c in cand:
    kws = sorted({k for tr in c["traits"].split(",") for k in TRAIT_KW[tr]})
    vid, n_qtl, best = qtl_coloc(c["lead_rsid"], kws)
    tg = ALIAS.get(c["gene_id"], c["gene_id"])
    th4, tsrc = best.get(tg, (None, ""))
    top = sorted(((h, g, s) for g, (h, s) in best.items() if h >= H4_MIN), reverse=True)
    rows.append({"cell_type": c["cell_type"], "gene_id": c["gene_id"], "peak_id": c["peak_id"], "lead_rsid": c["lead_rsid"],
                 "traits": c["traits"], "ot_variant_id": vid or "", "ot_n_qtl_credsets": n_qtl,
                 "ot_target_h4": "" if th4 is None else "%.3f" % th4, "ot_target_qtl": tsrc,
                 "ot_coloc_genes": "; ".join("%s %.3f [%s]" % (g, h, s) for h, g, s in top)})
    print(c["cell_type"], c["gene_id"], c["lead_rsid"], "| target H4:", rows[-1]["ot_target_h4"] or "NA", "| n H4>=0.8:", len(top))

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print("Saved:", OUT, len(rows), "rows")
