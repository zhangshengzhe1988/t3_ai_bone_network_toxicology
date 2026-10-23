#!/usr/bin/env python3
"""重绘 SI 面板 S1/S2/S11 与图形摘要（GA），统一 100 mm 宽、字号 ≥6 pt。"""
import sys, json
sys.path.insert(0, "code")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle, Polygon, FancyArrowPatch, Rectangle
from nature_style import save, C_HIGH as C_HI, C_LOW as C_LO

T3 = "选题3_网络毒理学_AI骨质疏松/"
OUT = T3 + "figures/"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False})

# ---------------- S1 靶点数/药 ----------------
st = json.load(open(T3 + "results/tables/S4_target_sets.json"))
sw = st["st_prob_ge_0.1"]
labels = ["Merged\ntargets", "Letrozole\n(SwissTarget)", "Anastrozole\n(SwissTarget)"]
vals = [len(set(st["drug_targets"])), len(set(sw["letrozole"])), len(set(sw["anastrozole"]))]
fig, ax = plt.subplots(figsize=(100/25.4, 62/25.4))
bars = ax.bar(range(3), vals, color=["#3C5488", C_HI, "#F39B7F"], width=0.6)
for b, v in zip(bars, vals):
    ax.text(b.get_x()+b.get_width()/2, v+0.8, str(v), ha="center", fontsize=7.5)
ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("Target genes (n)", fontsize=7.5); ax.tick_params(axis="y", labelsize=7)
ax.set_ylim(0, max(vals)*1.18)
ax.set_title("Target counts per drug", fontsize=8)
fig.subplots_adjust(left=0.18, right=0.97, bottom=0.22, top=0.88)
save(fig, OUT + "S6_target_counts")

# ---------------- S2 靶点来源构成 ----------------
ch = json.load(open(T3 + "results/tables/S3_chembl_targets.json"))
core = set(pd.read_csv(T3 + "results/tables/S4_intersection.csv")["gene"])
drugs = ["Letrozole", "Anastrozole"]
chembl_n = [len(ch[d.lower()]["targets"]) for d in drugs]
swiss = [set(sw[d.lower()]) for d in drugs]
x = np.arange(2); w = 0.26
fig, ax = plt.subplots(figsize=(100/25.4, 62/25.4))
b1 = ax.bar(x-w, chembl_n, w, color="#8491B4", label="ChEMBL")
b2 = ax.bar(x, [len(s) for s in swiss], w, color=C_LO, label="SwissTarget (≥ 0.1)")
b3 = ax.bar(x+w, [len(s & core) for s in swiss], w, color=C_HI, label="core-gene hits")
for bars in (b1, b2, b3):
    for b in bars:
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.8, f"{int(b.get_height())}",
                ha="center", fontsize=7.0)
ax.set_xticks(x); ax.set_xticklabels(drugs, fontsize=7.5)
ax.set_ylabel("Targets (n)", fontsize=7.5); ax.tick_params(axis="y", labelsize=7)
ax.set_ylim(0, max(chembl_n)*1.30)
ax.legend(fontsize=6.5, frameon=False, loc="upper right", handlelength=1.0, labelspacing=0.35)
ax.set_title("Target sources per drug", fontsize=8)
fig.subplots_adjust(left=0.18, right=0.97, bottom=0.16, top=0.88)
save(fig, OUT + "S7_target_sources")

# ---------------- S11 vKO 扰动基因 ----------------
d = pd.read_csv(T3 + "results/tables/S10_vko_VPS35_diffRegulation.csv")
d = d[d["gene"] != "VPS35"].copy()
d = d.sort_values("distance", ascending=False).head(15)
fig, ax = plt.subplots(figsize=(100/25.4, 78/25.4))
y = np.arange(len(d))[::-1]
ax.hlines(y, 0, d["distance"]*1e6, color="#B0B0B0", lw=1.0)
ax.scatter(d["distance"]*1e6, y, s=18, color=C_HI, zorder=3)
for yi, v in zip(y, d["distance"]*1e6):
    ax.text(v+0.08, yi, f"{v:.1f}", va="center", fontsize=6.2)
ax.set_yticks(y); ax.set_yticklabels(d["gene"], fontsize=6.8)
ax.set_xlabel("Network perturbation (×10$^{-6}$)", fontsize=7.2)
ax.set_xlim(0, (d["distance"]*1e6).max()*1.25)
ax.tick_params(axis="x", labelsize=6.6)
ax.set_title("Top perturbed genes after VPS35 knockout", fontsize=7.8)
fig.subplots_adjust(left=0.28, right=0.97, bottom=0.15, top=0.90)
save(fig, OUT + "S11_vko_top_genes")
print("S1/S2/S11 done")

# ---------------- 图形摘要（130 × 52 mm） ----------------
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle, Polygon, FancyArrowPatch, Rectangle

RED, BLUE, GREEN, PURPLE, PINK, GREY = "#D55E00", "#0072B2", "#009E73", "#7B6AA8", "#B85C8A", "#5c5f66"
W, H = 130, 52
fig = plt.figure(figsize=(W/25.4, H/25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

def arrow(x1, y1, x2, y2, col=GREY, lw=1.2, ms=7):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=ms,
                                 lw=lw, color=col, shrinkA=0, shrinkB=0))

# 三列标题带
for x0, lab, col in [(2, "HUMAN", BLUE), (45, "SHARED TARGET", PURPLE), (88, "FISH", GREEN)]:
    ax.add_patch(FancyBboxPatch((x0, 47.2), 40, 3.8, boxstyle="round,pad=0.2,rounding_size=0.8",
                                fc=col, ec="none"))
    ax.text(x0+20, 48.5, lab, ha="center", va="center", color="white", fontsize=6.6, fontweight="bold")

# ---- 左列：人 ----
ax.add_patch(FancyBboxPatch((8, 40.5), 28, 4.2, boxstyle="round,pad=0.2,rounding_size=2",
                            fc="#F4E9F7", ec=PURPLE, lw=0.7))
ax.text(22, 42.6, "aromatase inhibitor", ha="center", va="center", fontsize=6.2, color="#4a3a6b")
arrow(22, 40.2, 22, 36.6)
ax.add_patch(Ellipse((22, 33), 26, 6.4, fc="#D9E8F5", ec=BLUE, lw=0.8))
ax.text(22, 33, "CYP19A1", ha="center", va="center", fontsize=6.6, fontweight="bold", color="#123c5a")
arrow(22, 29.6, 22, 26.2)
for dx in (-3.2, 0, 3.2):
    ax.add_patch(Circle((22+dx, 24.4), 1.5, fc="#8B5A72", ec="white", lw=0.4))
ax.text(22, 20.6, "oestradiol ↓", ha="center", va="center", fontsize=6.2, color="#6b3a52")
arrow(22, 18.8, 22, 15.6)
# 骨图标
ax.add_patch(Rectangle((17.5, 12.0), 9, 2.6, fc="#F4EBDC", ec="#B9A98F", lw=0.6))
for cx, cy in [(17.5, 12.5), (17.5, 14.1), (26.5, 12.5), (26.5, 14.1)]:
    ax.add_patch(Circle((cx, cy), 1.1, fc="#F4EBDC", ec="#B9A98F", lw=0.6))
ax.text(22, 8.6, "bone loss", ha="center", va="center", fontsize=6.4, fontweight="bold", color="#7a4a1a")
ax.text(22, 5.4, "osteoporosis risk", ha="center", va="center", fontsize=6.0, color=BLUE)

# ---- 中列：共同靶点 ----
ax.add_patch(Ellipse((65, 30), 34, 15, fc="#EDE7F4", ec=PURPLE, lw=0.8))
ax.text(65, 30, "CYP19A1", ha="center", va="center", fontsize=7.0, fontweight="bold", color="#4a3a6b")
ax.add_patch(Circle((49, 33.5), 2.0, fc="#7a4a1a", ec="white", lw=0.5))
ax.add_patch(Circle((81, 33.5), 2.0, fc="#1f5f4a", ec="white", lw=0.5))
arrow(51.4, 33.2, 58, 31.6, col="#4a3a6b", lw=0.9, ms=6)
arrow(78.6, 33.2, 72, 31.6, col="#4a3a6b", lw=0.9, ms=6)
ax.text(49, 37.6, "letrozole", ha="center", fontsize=6.0, color="#7a4a1a")
ax.text(81, 37.6, "anastrozole", ha="center", fontsize=6.0, color="#1f5f4a")
ax.text(65, 20.2, "9/10 active-site residues identical", ha="center", fontsize=6.2, color="#4a3a6b")
ax.text(65, 16.6, "−7.63 / −7.20 kcal/mol (fish)", ha="center", fontsize=6.0, color=RED)
ax.text(65, 13.0, "−7.22 / −7.12 kcal/mol (human)", ha="center", fontsize=6.0, color=GREY)
ax.text(65, 8.2, "conserved initiating event", ha="center", fontsize=6.2, fontweight="bold", color=PURPLE)

# ---- 右列：鱼 ----
ax.text(108, 45.2, "surface water", ha="center", fontsize=6.0, color="#1f5f4a")
for i, y0 in enumerate([41.6, 43.0]):
    xs = np.linspace(99, 117, 60)
    ax.plot(xs, y0 + 0.7*np.sin((xs-99)*1.6 + i*1.2), color="#4DBBD5", lw=0.9)
ax.add_patch(Ellipse((108, 36.2), 13, 5.2, fc="#CFE7F0", ec="#2C7A8C", lw=0.8))
ax.add_patch(Polygon([[114.3, 36.2], [120.5, 39.0], [120.5, 33.4]], closed=True, fc="#CFE7F0", ec="#2C7A8C", lw=0.8))
ax.add_patch(Circle((102.8, 36.8), 0.55, fc="#1f2f33"))
arrow(108, 33.4, 108, 30.6)
ax.add_patch(Ellipse((108, 27.4), 24, 5.6, fc="#E7F3EC", ec=GREEN, lw=0.8))
ax.text(108, 27.4, "cyp19a1a", ha="center", va="center", fontsize=6.6, fontweight="bold", color="#1f5f4a")
arrow(108, 24.6, 108, 21.8)
for dx in (-3.2, 0, 3.2):
    ax.add_patch(Circle((108+dx, 20.0), 1.5, fc="#3f5b66", ec="white", lw=0.4))
ax.text(108, 16.2, "sex steroids ↓", ha="center", fontsize=6.2, color="#2C4A55")
arrow(108, 14.4, 108, 11.6)
ax.add_patch(Ellipse((108, 8.6), 13, 5.0, fc="#F7DCE8", ec=PINK, lw=0.8))
ax.add_patch(Circle((108, 8.6), 1.2, fc="#F2C3D6", ec=PINK, lw=0.5))
ax.text(108, 4.6, "gonadal disruption", ha="center", fontsize=6.2, fontweight="bold", color="#8a3a5e")

ax.text(65, 1.6, "CYP19A1 inhibition → steroidogenic shift → bone loss (human) / gonadal disruption (fish)",
        ha="center", va="center", fontsize=6.2, color="#333333")
save(fig, OUT + "GA_130x52")
print("GA done")
