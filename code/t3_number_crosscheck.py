#!/usr/bin/env python3
"""全文数字 — 图表交叉核对。

对正文中每个数值，检查是否能在 (a) 主图 PDF 文本、(b) 补充图 PDF 文本、
(c) 附表 xlsx 单元格 中找到（支持四舍五入容差），并输出未找到/不一致清单。
"""
import sys, re, glob, os
sys.path.insert(0, "code")
import docx, fitz
from openpyxl import load_workbook

DOCX = sys.argv[1] if len(sys.argv) > 1 else \
    "选题3_网络毒理学_AI骨质疏松/manuscript/04. Manuscript_选题3_Environmental_Pollution版.docx"
FIGD = "投稿材料/选题3_正图"
SIFD = "选题3_网络毒理学_AI骨质疏松/figures"
XLSX = "投稿材料/选题3_附表.xlsx"

d = docx.Document(DOCX)
paras = [p.text.strip() for p in d.paragraphs if p.text.strip()]
start = paras.index("1. Introduction"); end = paras.index("Funding")
body = " ".join(paras[start:end])

# 图/补充图文本
fig_text = {}
for f in sorted(glob.glob(FIGD + "/Figure*.pdf")):
    dd = fitz.open(f); fig_text[os.path.basename(f)] = dd[0].get_text().replace("\n", " "); dd.close()
si_text = {}
for f in sorted(glob.glob(SIFD + "/S*.pdf")):
    dd = fitz.open(f); si_text[os.path.basename(f)] = dd[0].get_text().replace("\n", " "); dd.close()

# 附表数值集合
wb = load_workbook(XLSX, data_only=True)
xl_vals = []
for ws in wb:
    for row in ws.iter_rows(values_only=True):
        for c in row:
            if isinstance(c, (int, float)):
                xl_vals.append(float(c))
            elif isinstance(c, str):
                for m in re.findall(r"-?\d+(?:\.\d+)?", c):
                    try: xl_vals.append(float(m))
                    except: pass

def match(num_str, text):
    if num_str in text: return True
    try: v = float(num_str.replace(",", ""))
    except: return False
    for m in re.findall(r"-?\d+(?:,\d{3})*(?:\.\d+)?", text):
        try: x = float(m.replace(",", ""))
        except: continue
        if v == 0:
            if abs(x) < 1e-12: return True
        elif abs(x - v) <= max(abs(v)*0.006, 1e-9): return True
    return False

def in_xl(num_str):
    if not xl_vals: return False
    try: v = float(num_str.replace(",", ""))
    except: return False
    for x in xl_vals:
        if v == 0:
            if abs(x) < 1e-12: return True
        elif abs(x - v) <= max(abs(v)*0.006, 1e-9): return True
    return False

# 提取正文数字（排除年份、编号、引用）
tokens = []
for m in re.finditer(r"-?\d[\d,]*(?:\.\d+)?(?:×10-?\d+)?", body):
    s = m.group(0)
    ctx = body[max(0, m.start()-45):m.end()+25].replace("\n", " ")
    tokens.append((s, ctx))

seen = set(); rows = []
for s, ctx in tokens:
    key = (s, ctx[:40])
    if key in seen: continue
    seen.add(key)
    base = s.split("×")[0]
    if re.fullmatch(r"(19|20)\d{2}", base.replace(",", "")): continue      # 年份
    if re.fullmatch(r"\d{1,2}", base) and float(base) <= 10 and "Fig" not in ctx and "Table" not in ctx:
        pass
    hit_fig = [k for k, v in fig_text.items() if match(base, v)]
    hit_si  = [k for k, v in si_text.items() if match(base, v)]
    hit_xl  = in_xl(base)
    rows.append((s, ctx, hit_fig, hit_si, hit_xl))

print(f"正文数值 token: {len(rows)}")
print(f"{'数值':>10}  {'图':<18} {'SI':<10} {'表':<4} 上下文")
missing = []
for s, ctx, hf, hs, hx in rows:
    f = (",".join(x.replace("Figure", "F").replace(".pdf", "") for x in hf))[:18] or "-"
    si = (",".join(x[:9] for x in hs))[:10] or "-"
    print(f"{s:>10}  {f:<18} {si:<10} {'Y' if hx else '-':<4} {ctx[:78]}")
    if not (hf or hs or hx): missing.append((s, ctx))

print("\n=== 未在任何图/表/补充图中找到的数值 ===")
for s, ctx in missing:
    print(f"  {s:>10} | {ctx[:100]}")
