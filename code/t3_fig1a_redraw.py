#!/usr/bin/env python3
"""重绘 Fig 1A（工作流程图）：170 × 88 mm，Arial 可编辑文本（PDF fonttype 42）。"""
import sys
sys.path.insert(0, "code")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Ellipse, Circle, Polygon, FancyArrowPatch
from nature_style import save

T3 = "选题3_网络毒理学_AI骨质疏松/"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.linewidth": 0.5})

BLUE, PURPLE, ORANGE, GREEN = "#0072B2", "#7B6AA8", "#E69F00", "#009E73"
GREY = "#4D4D4D"

W, H = 170, 88
fig = plt.figure(figsize=(W/25.4, H/25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

def txt(x, y, s, size=6.0, bold=False, color=GREY, ha="left", va="center"):
    ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va,
            fontweight="bold" if bold else "normal", family="Arial")

def arrow(x1, y1, x2, y2, col="#8C8C8C"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=6,
                                 lw=0.8, color=col, shrinkA=0, shrinkB=0))

# ---- 四列标题 ----
cols = [(2, "DATA SOURCES", BLUE), (44, "NETWORK MAPPING", PURPLE),
        (86, "KEY EVIDENCE", ORANGE), (128, "TRANSLATION", GREEN)]
for x0, lab, col in cols:
    ax.add_patch(FancyBboxPatch((x0, 81.2), 38, 5.2, boxstyle="round,pad=0.15,rounding_size=0.7",
                                fc=col, ec="none"))
    txt(x0 + 19, 83.8, lab, size=6.4, bold=True, color="white", ha="center")

# ---- 图标 ----
def cylinder(x, y, col):
    ax.add_patch(Rectangle((x, y - 2.2), 4.6, 4.4, fc=col, ec="none", alpha=0.85))
    ax.add_patch(Ellipse((x + 2.3, y + 2.2), 4.6, 1.5, fc=col, ec="white", lw=0.4))
    ax.add_patch(Ellipse((x + 2.3, y - 2.2), 4.6, 1.5, fc=col, ec="white", lw=0.4))

def cell(x, y, col):
    ax.add_patch(Circle((x + 2.3, y), 2.4, fc="#EAF2F8", ec=col, lw=0.7))
    ax.add_patch(Circle((x + 2.6, y + 0.3), 0.9, fc=col, ec="none"))

def venn(x, y):
    ax.add_patch(Circle((x + 1.9, y), 2.0, fc="#E64B35", ec="none", alpha=0.45))
    ax.add_patch(Circle((x + 3.5, y), 2.0, fc="#4DBBD5", ec="none", alpha=0.45))

def network(x, y):
    pts = [(x + 2.3, y + 2.3), (x + 0.4, y - 0.8), (x + 2.3, y - 2.2), (x + 4.2, y - 0.8)]
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            ax.plot([pts[i][0], pts[j][0]], [pts[i][1], pts[j][1]], color="#B0B0B0", lw=0.5, zorder=1)
    for i, (px, py) in enumerate(pts):
        ax.add_patch(Circle((px, py), 0.75, fc=PURPLE if i == 0 else "#D9D9D9", ec="white", lw=0.3, zorder=2))

def protein(x, y, col):
    ax.add_patch(Ellipse((x + 2.3, y), 4.8, 3.0, fc="#E8F1F8", ec=col, lw=0.7))
    ax.add_patch(Circle((x + 2.3, y), 0.7, fc="#D55E00", ec="none"))

def fish(x, y, col):
    ax.add_patch(Ellipse((x + 2.2, y), 4.0, 2.2, fc="#CFE7F0", ec=col, lw=0.6))
    ax.add_patch(Polygon([[x + 4.0, y], [x + 5.4, y + 1.2], [x + 5.4, y - 1.2]],
                         closed=True, fc="#CFE7F0", ec=col, lw=0.6))
    ax.add_patch(Circle((x + 1.1, y + 0.4), 0.22, fc="#1f2f33"))

def bone(x, y, col):
    ax.add_patch(Rectangle((x + 0.8, y - 0.6), 3.4, 1.2, fc="#F4EBDC", ec=col, lw=0.5))
    for cx in (x + 0.8, x + 4.2):
        ax.add_patch(Circle((cx, y + 0.6), 0.65, fc="#F4EBDC", ec=col, lw=0.5))
        ax.add_patch(Circle((cx, y - 0.6), 0.65, fc="#F4EBDC", ec=col, lw=0.5))

# ---- 行内容 ----
rows = [75.0, 58.0, 41.0, 24.0]
# 数据源
cylinder(4, rows[0], BLUE);  txt(11, rows[0]+1.2, "FAERS", 6.2, True);  txt(11, rows[0]-1.3, "2013–2026; three AIs", 5.8)
cylinder(4, rows[1], BLUE);  txt(11, rows[1]+1.2, "ChEMBL + SwissTarget", 6.2, True); txt(11, rows[1]-1.3, "drug targets (42)", 5.8)
cylinder(4, rows[2], BLUE);  txt(11, rows[2]+1.2, "GeneCards", 6.2, True); txt(11, rows[2]-1.3, "osteoporosis genes (3,713)", 5.8)
cell(4, rows[3], BLUE);      txt(11, rows[3]+1.2, "GEO ×4", 6.2, True);   txt(11, rows[3]-1.3, "MSC, B cell, PBMC cohorts", 5.8)

# 网络映射
venn(46, rows[0]);    txt(53, rows[0]+1.2, "42 ∩ 3,713", 6.2, True); txt(53, rows[0]-1.3, "→ nine core genes", 5.8)
network(46, rows[1]); txt(53, rows[1]+1.2, "STRING PPI", 6.2, True);  txt(53, rows[1]-1.3, "CYP19A1 hub", 5.8)
ax.plot([46, 47.5, 48.5, 50.8], [rows[2]-1.5, rows[2]+1.2, rows[2]-1.0, rows[2]+0.6],
        color=PURPLE, lw=0.9)  # penalty path
txt(53, rows[2]+1.2, "WGCNA + LASSO", 6.2, True); txt(53, rows[2]-1.3, "VPS35 candidate", 5.8)
protein(46, rows[3], PURPLE); txt(53, rows[3]+1.2, "Blind docking", 6.2, True); txt(53, rows[3]-1.3, "Vina + Vinardo", 5.8)

# 证据
txt(88, rows[0]+1.2, "Disproportionality", 6.2, True); txt(88, rows[0]-1.3, "ROR 2.77–4.93", 5.8)
txt(88, rows[1]+1.2, "Two cohorts", 6.2, True);        txt(88, rows[1]-1.3, "no core-gene DE (n = 9, 20)", 5.8)
txt(88, rows[2]+1.2, "Docking", 6.2, True);            txt(88, rows[2]-1.3, "−7.22 / −6.03 kcal/mol", 5.8)
txt(88, rows[3]+1.2, "Virtual knockout", 6.2, True);   txt(88, rows[3]-1.3, "25 genes perturbed", 5.8)

def droplet(x, y, col):
    ax.add_patch(Ellipse((x + 2.3, y - 0.6), 3.2, 3.4, fc="#D6ECF3", ec=col, lw=0.7))
    ax.add_patch(Polygon([[x + 2.3, y + 2.6], [x + 0.9, y + 0.4], [x + 3.7, y + 0.4]],
                         closed=True, fc="#D6ECF3", ec=col, lw=0.7))

def aop(x, y, col):
    for k, px in enumerate((x + 0.9, x + 2.9, x + 4.9)):
        ax.add_patch(Circle((px, y), 0.85, fc="#E7F3EC", ec=col, lw=0.7))
        if k < 2:
            ax.plot([px + 0.9, px + 1.1], [y, y], color=col, lw=0.6)

# 转化
fish(130, rows[0], GREEN); txt(137, rows[0]+1.2, "Zebrafish Cyp19a1a", 6.2, True); txt(137, rows[0]-1.3, "9/10 active-site residues", 5.8)
droplet(130, rows[1], GREEN); txt(137, rows[1]+1.2, "Exposure margin", 6.2, True); txt(137, rows[1]-1.3, "EAR 0.70 at 100 ng/L", 5.8)
bone(130, rows[2], GREEN); txt(137, rows[2]+1.2, "Human bone risk", 6.2, True); txt(137, rows[2]-1.3, "ROR 2.77–4.93", 5.8)
aop(130, rows[3], GREEN); txt(137, rows[3]+1.2, "AOP reading", 6.2, True); txt(137, rows[3]-1.3, "bone loss / gonadal disruption", 5.8)

# 列间箭头
for y in (rows[0], rows[2]):
    arrow(40.5, y, 43.6, y)
    arrow(82.5, y, 85.6, y)
    arrow(124.5, y, 127.6, y)

# 底部 AOP 条
ax.add_patch(FancyBboxPatch((2, 8.0), 166, 6.4, boxstyle="round,pad=0.2,rounding_size=0.8",
                            fc="#F2F5FA", ec="none"))
txt(85, 11.2, "AI exposure → CYP19A1 inhibition → oestrogen ↓ → bone loss (human) / gonadal disruption (fish)",
    size=6.2, ha="center", color="#333333")

save(fig, T3 + "figures/Fig1A_workflow_a4")
fig.savefig(T3 + "figures/Fig1A_workflow_a4.svg", format="svg")
print("Fig1A redrawn")
