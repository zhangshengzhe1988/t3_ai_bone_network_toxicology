#!/usr/bin/env python3
"""终审：正文↔图↔图注 数据一致性、叙事顺序、引用顺序、文献对应。"""
import sys, re, glob, os
sys.path.insert(0, "code")
import docx, fitz
from openpyxl import load_workbook

DOCX = "选题3_网络毒理学_AI骨质疏松/manuscript/submission/01_Manuscript_Environmental_Pollution.docx"
d = docx.Document(DOCX)
paras = [p.text.strip() for p in d.paragraphs if p.text.strip()]
full = " ".join(paras)
body = " ".join(paras[paras.index("1. Introduction"):paras.index("Funding")])

# 图文本
fig = {}
for i in range(1, 8):
    dd = fitz.open(f"投稿材料/选题3_正图/Figure{i}.pdf"); fig[i] = dd[0].get_text().replace("\n", " "); dd.close()
si = {}
for f in sorted(glob.glob("选题3_网络毒理学_AI骨质疏松/figures/S*.pdf")):
    dd = fitz.open(f); si[os.path.basename(f)] = dd[0].get_text().replace("\n", " "); dd.close()

def in_fig(fig_id, v):
    t = fig[fig_id]
    for m in re.findall(r"-?\d+(?:\.\d+)?", t.replace(",", "")):
        try: x = float(m)
        except: continue
        if abs(x - v) <= max(abs(v)*0.006, 1e-9): return True
    return False

# (断言, 图号, 值, 图注应包含的关键词, 正文片段关键词)
rows = [
 ("Fig1B ROR letrozole", 1, 2.77, "ROR forest", "2.77"),
 ("Fig1B ROR anastrozole", 1, 4.08, "ROR forest", "4.08"),
 ("Fig1B ROR exemestane", 1, 3.61, "Exemestane FAERS", "3.61"),
 ("Fig1B ROR let frac", 1, 4.93, "ROR forest", "4.93"),
 ("Fig1B ROR ana frac", 1, 2.86, "ROR forest", "2.86"),
 ("Fig1E age 253", 1, 253, "Age and sex stratification", "253"),
 ("Fig1E female 390", 1, 390, "Age and sex stratification", "390"),
 ("Fig2A genes 3713", 2, 3713, "42 drug targets intersect 3,713", "3,713"),
 ("Fig2D GO terms", 2, 256, "GO enrichment", "256"),
 ("Fig3A DEG 1447", 3, 1447, "1,447 DEGs", "1,447"),
 ("Fig3D DEG 28", 3, 28, "28 DEGs", "28 DEGs"),
 ("Fig3F DEG 391", 3, 391, "DEG counts", "391"),
 ("Fig4C M4 r 0.363", 4, 0.363, "M4 and M5", "0.363"),
 ("Fig4C padj 0.0066", 4, 0.0066, "padj = 0.0066", "0.0066"),
 ("Fig4B M4 size 1743", 4, 1743, "Module sizes", "1,743"),
 ("Fig4B M5 size 931", 4, 931, "Module sizes", "931"),
 ("Fig5B lambda 0.060", 5, 0.060, "lambda.min = 0.060", "0.060"),
 ("Fig5B lambda 0.146", 5, 0.146, "lambda.1se = 0.146", "0.146"),
 ("Fig5E FPKM 12.97", 5, 12.97, "Median FPKM", "12.97"),
 ("Fig6A dock 7.22", 6, 7.22, "7.22 kcal/mol", "7.22"),
 ("Fig6E Vina 7.63", 6, 7.63, "six target-ligand combinations", "7.63"),
 ("Fig7A dock 7.20", 7, 7.20, "Docking scores", "7.20"),
 ("Fig7B identity 51.3", 7, 51.3, "9/10 identical", "51.3"),
 ("Fig7D EAR 0.704", 7, 0.704, "exposure-activity ratio", "0.70"),
]
print(f"{'断言':<26}{'图中':<6}{'图注':<6}{'正文':<6}")
bad = []
for lab, fid, v, leg, txt in rows:
    f_ok = in_fig(fid, v)
    l_ok = leg in full
    t_ok = txt in body
    if not (f_ok and l_ok and t_ok): bad.append(lab)
    print(f"{lab:<26}{'Y' if f_ok else 'N':<6}{'Y' if l_ok else 'N':<6}{'Y' if t_ok else 'N':<6}")

# 叙事顺序：各节内图号递增
print("\n=== 各节图表引用顺序 ===")
cur = None
for t in paras:
    if re.match(r"^3\.\d ", t): cur = t; seq = []
    elif cur and re.match(r"^3\.\d ", t) is None and re.match(r"^4\.", t): cur = None
    if cur:
        for m in re.finditer(r"Fig\. (\d)", t):
            if int(m.group(1)) not in seq: seq.append(int(m.group(1)))
    if cur and (re.match(r"^4\.", t) or t == paras[-1]):
        print(f"  {cur[:45]:47s}{seq}")

print("\n=== 图内面板顺序 vs 图注顺序 ===")
for i in range(1, 8):
    key = [t for t in paras if t.startswith(f"Fig. {i}.")]
    if not key: continue
    idx = paras.index(key[0]); leg = paras[idx+1]
    panels = re.findall(r"\(([A-G])\)", leg)
    order_ok = panels == sorted(panels)
    print(f"  Fig {i}: {''.join(panels)} {'OK' if order_ok else 'OUT OF ORDER'}")

print("\n未通过：", bad or "无")
