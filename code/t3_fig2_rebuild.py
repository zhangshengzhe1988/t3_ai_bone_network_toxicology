#!/usr/bin/env python3
"""重绘 Fig 2 的 B/C/D/E/F 面板：修图例遮挡、标签重叠、字号过小。"""
import sys
sys.path.insert(0, "code")
import numpy as np, pandas as pd, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch
from nature_style import save, C_HIGH, C_LOW

T3 = "选题3_网络毒理学_AI骨质疏松/"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})

# ---------- B 靶点来源 ----------
src = pd.read_csv(T3 + "results/tables/S4b_target_source_sensitivity.csv")
fig, ax = plt.subplots(figsize=(80/25.4, 52/25.4))
x = np.arange(len(src)); w = 0.26
b1 = ax.bar(x - w, src["chembl_n"], w, color="#3C5488", label="ChEMBL")
b2 = ax.bar(x, src["swisstarget_n"], w, color=C_HIGH, label="SwissTarget (p ≥ 0.1)")
b3 = ax.bar(x + w, src["merged_n"], w, color="#9DB4C9", label="merged")
for bars in (b1, b2, b3):
    for b in bars:
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.7, f"{int(b.get_height())}",
                ha="center", fontsize=5.8)
ax.set_xticks(x); ax.set_xticklabels([d.capitalize() for d in src["drug"]], fontsize=6.2)
ax.set_ylabel("Targets (n)", fontsize=6.4); ax.tick_params(axis="y", labelsize=5.8)
ax.set_ylim(0, max(src["merged_n"])*1.55)
ax.legend(fontsize=5.4, frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.20), handlelength=0.9, columnspacing=1.0)

fig.subplots_adjust(left=0.17, right=0.98, bottom=0.16, top=0.72)
save(fig, T3 + "figures/Fig3B_targetdb_venn")

# ---------- C PPI 网络 ----------
e = pd.read_csv(T3 + "results/tables/S5_string_network.tsv", sep="\t")[["preferredName_A", "preferredName_B", "score"]]
nodes = sorted(set(e["preferredName_A"]) | set(e["preferredName_B"]))
nodes = [n for n in nodes if n != "CYP19A1"]
ang = np.linspace(90, 90-360, len(nodes), endpoint=False)
pos = {n: (np.cos(np.deg2rad(a)), np.sin(np.deg2rad(a))) for n, a in zip(nodes, ang)}
deg = pd.concat([e["preferredName_A"], e["preferredName_B"]]).value_counts()
fig, ax = plt.subplots(figsize=(62/25.4, 56/25.4))
for r in e.itertuples():
    a, b = pos[r.preferredName_A], pos[r.preferredName_B]
    ax.plot([a[0], b[0]], [a[1], b[1]], color="#B8B8B8", lw=0.5 + 0.7*r.score, zorder=1)
for n, (px, py) in pos.items():
    r_ = 0.085 + 0.012*deg[n]
    ax.add_patch(Circle((px, py), r_, fc=C_LOW if n != "CYP19A1" else C_HIGH, ec="white", lw=0.6, zorder=3))
    ax.text(px*1.30, py*1.30, n, ha="center", va="center", fontsize=5.6, zorder=4)
# CYP19A1 在 STRING 中无 ≥0.4 的边，单独画出
ax.add_patch(Circle((-1.30, -1.28), 0.085, fc="#F2F2F2", ec="#999999", lw=0.7, ls="--", zorder=3))
ax.text(-1.30, -1.10, "CYP19A1", ha="center", va="center", fontsize=5.6, color="#777777", zorder=4)
ax.text(-1.30, -1.45, "no edge ≥ 0.4", ha="center", va="center", fontsize=4.6, color="#777777", zorder=4)
ax.set_xlim(-1.55, 1.55); ax.set_ylim(-1.55, 1.45); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("STRING PPI network (7 edges, score ≥ 0.4)", fontsize=6.4)
save(fig, T3 + "figures/Fig3C_PPI_network")

# ---------- D 化合物-靶点-通路 ----------
st = json.load(open(T3 + "results/tables/S4_target_sets.json"))
core = ["CYP19A1", "HDAC1", "HDAC6", "ACHE", "ADRB1", "ADRB2", "ADRB3", "HTR2A", "HTR2C"]
go = pd.read_csv(T3 + "results/tables/S6_go_bp.csv").sort_values("p.adjust").head(4)
import textwrap
terms = ["\n".join(textwrap.wrap(t, 26)) for t in go["Description"]]
fig, ax = plt.subplots(figsize=(70/25.4, 56/25.4))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
# 药物
for i, (drug, col) in enumerate([("Letrozole", "#7a4a1a"), ("Anastrozole", "#1f5f4a")]):
    yy = 7.6 - i*1.3
    ax.add_patch(Circle((0.7, yy), 0.28, fc=col, ec="none"))
    ax.text(1.15, yy, drug, va="center", fontsize=5.6)
    for j, g in enumerate(core):
        yy2 = 8.8 - j*0.95
        ax.plot([1.0, 3.2], [yy, yy2], color=col, lw=0.35, alpha=0.55, zorder=1)
# 基因
for j, g in enumerate(core):
    yy2 = 8.8 - j*0.95
    ax.add_patch(Circle((3.35, yy2), 0.13, fc="#3C5488", ec="none"))
    ax.text(3.6, yy2, g, va="center", fontsize=5.2)
# 通路
for k, term in enumerate(terms):
    yk = 8.7 - k*1.9
    ax.add_patch(Circle((7.0, yk), 0.13, fc=C_HIGH, ec="none"))
    ax.text(7.1, yk, term, va="center", fontsize=5.2)
    for j in range(len(core)):
        ax.plot([4.9, 6.8], [8.8 - j*0.95, yk], color="#D0D0D0", lw=0.25, zorder=0)
ax.set_title("Compound–target–pathway network", fontsize=6.6)
save(fig, T3 + "figures/Fig4C_compound_target_pathway")

# ---------- E 证据矩阵 ----------
ch = json.load(open(T3 + "results/tables/S3_chembl_targets.json"))
genes = sorted(core)
M = np.zeros((len(genes), 4))
for i, g in enumerate(genes):
    M[i, 0] = 1 if g in set(st["st_prob_ge_0.1"]["letrozole"]) else 0
    M[i, 1] = 1 if g in set(st["st_prob_ge_0.1"]["anastrozole"]) else 0
    M[i, 2] = 1 if any(g == str(v) for v in ch["letrozole"]["targets"].values()) else 0
    M[i, 3] = deg.get(g, 0)
fig, ax = plt.subplots(figsize=(52/25.4, 56/25.4))
ax.imshow(M[:, :3], cmap=plt.matplotlib.colors.ListedColormap(["#F2F2F2", C_HIGH]), aspect="auto")
for i in range(len(genes)):
    for j in range(3):
        ax.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, fill=False, ec="white", lw=0.8))
    ax.text(3.0, i, f"{int(M[i,3])}", va="center", ha="center", fontsize=5.6)
ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(["Let", "Ana", "ChEMBL", "PPI"], fontsize=5.6)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=5.6)
ax.tick_params(length=0); ax.set_xlim(-0.5, 3.5)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.set_title("Core-gene evidence", fontsize=6.6)
fig.subplots_adjust(left=0.30, right=0.97, bottom=0.14, top=0.88)
save(fig, T3 + "figures/Fig3F_evidence_matrix")

# ---------- F GO 点图 ----------
cats = [("GO-BP", "S6_go_bp.csv", "#E64B35"), ("GO-CC", "S6_go_CC.csv", "#4DBBD5"), ("GO-MF", "S6_go_MF.csv", "#00A087")]
fig, ax = plt.subplots(figsize=(150/25.4, 60/25.4))
y = 0; labels = []; ypos = []
for name, f, col in cats:
    g = pd.read_csv(T3 + "results/tables/" + f).sort_values("p.adjust").head(8)
    ax.axhspan(y-0.5, y+7.5, color=col, alpha=0.07, lw=0)

    for _, r in g.iterrows():
        ax.scatter(-np.log10(r["p.adjust"]), y, s=r["Count"]*7, c=r["Count"],
                   cmap="YlOrRd", vmin=2, vmax=5, edgecolors="grey", linewidths=0.35, zorder=3)
        labels.append(r["Description"][:34] + ("…" if len(r["Description"]) > 34 else ""))
        ypos.append(y); y += 1
ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=5.8)
cat_cols = []
for name, f, col in cats:
    cat_cols += [col]*8
for lab, col in zip(ax.get_yticklabels(), cat_cols):
    lab.set_color(col)
ax.set_ylim(-0.8, y-0.2); ax.set_xlabel("−log10(adj. P)", fontsize=6.6)
ax.tick_params(axis="x", labelsize=5.8); ax.tick_params(axis="y", length=0)
sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=plt.Normalize(2, 5)); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02); cb.set_label("Genes (n)", fontsize=6.4)
cb.ax.tick_params(labelsize=5.8)
handles = [plt.Line2D([], [], marker="s", ls="", color=c, label=n, ms=4) for n, f, c in cats]
ax.legend(handles=handles, fontsize=6.0, frameon=False, ncol=3, loc="lower right",
          bbox_to_anchor=(1.0, 1.02), handletextpad=0.3, columnspacing=1.2)
fig.subplots_adjust(left=0.40, right=0.90, bottom=0.14, top=0.90)
save(fig, T3 + "figures/Fig4A_GO_three")
print("Fig2 B/C/D/E/F rebuilt")
