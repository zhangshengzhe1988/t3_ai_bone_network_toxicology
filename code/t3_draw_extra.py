import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import matplotlib.text as _mtext
def FONT_FLOOR(fig, floor=6.5):
    for _t in fig.findobj(_mtext.Text):
        if _t.get_fontsize() < floor: _t.set_fontsize(floor)
    return fig

"""T3 Fig3E/F、Fig6E/F 扩充面板"""
import sys; sys.path.insert(0, "code")
import pandas as pd, numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from nature_style import save
T3 = "选题3_网络毒理学_AI骨质疏松/"
C_HI=(0.902,0.294,0.208); C_LO=(0.302,0.686,0.835); GREY='#7E8AA2'
plt.rcParams.update({"font.family":"Arial","font.size":6,"pdf.fonttype":42,"axes.linewidth":0.5,
                     "axes.spines.top":False,"axes.spines.right":False})

# ---------- Fig3E 每药靶点来源构成 ----------
st = json.load(open(T3+"results/tables/S4_target_sets.json"))
ch = json.load(open(T3+"results/tables/S3_chembl_targets.json"))
core = set(pd.read_csv(T3+"results/tables/S4_intersection.csv")["gene"])
drugs = ["Letrozole", "Anastrozole"]
chembl_n = [len(ch[d.lower()]["targets"]) for d in drugs]
swiss = [st["st_prob_ge_0.1"][d.lower()] for d in drugs]
swiss_n = [len(s) for s in swiss]
corehit_n = [len(set(s) & core) for s in swiss]
x = np.arange(2); w = 0.26
fig, ax = plt.subplots(figsize=(2.2, 2.0))
b1 = ax.bar(x-w, chembl_n, w, color="#8491B4", label="ChEMBL")
b2 = ax.bar(x, swiss_n, w, color=C_LO, label="SwissTarget\n(≥ 0.1)")
b3 = ax.bar(x+w, corehit_n, w, color=C_HI, label="core-gene hits")
for bars in (b1,b2,b3):
    for b in bars:
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.6, f"{int(b.get_height())}",
                ha='center', fontsize=5.5)
ax.set_xticks(x); ax.set_xticklabels(drugs, fontsize=5.8)
ax.set_ylabel("Targets (n)"); ax.set_ylim(0, max(chembl_n)*1.18)
ax.legend(fontsize=4.8, frameon=False, loc='upper right', handlelength=1.0, labelspacing=0.3)
ax.set_title("Target sources per drug", fontsize=6.5)
fig.tight_layout(); save(fig, T3+"figures/Fig3E_target_sources")

# ---------- Fig3F 核心基因证据矩阵 ----------
ppi = pd.read_csv(T3+"results/tables/S5_string_network.tsv", sep="\t")
deg = pd.concat([ppi["preferredName_A"], ppi["preferredName_B"]]).value_counts()
genes = sorted(core)
M = np.zeros((len(genes), 4))
for i,g in enumerate(genes):
    M[i,0] = 1 if g in set(st["st_prob_ge_0.1"]["letrozole"]) else 0
    M[i,1] = 1 if g in set(st["st_prob_ge_0.1"]["anastrozole"]) else 0
    M[i,2] = 1 if any(g == str(v) for v in ch["letrozole"]["targets"].values()) else 0
    M[i,3] = deg.get(g, 0)
fig, ax = plt.subplots(figsize=(2.1, 2.5))
for i in range(len(genes)):
    for j in range(3):
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=C_HI if M[i,j] else "#F2F2F2", edgecolor="white", lw=0.8))
    ax.text(3.1, i+0.5, f"{int(M[i,3])}", va='center', fontsize=5.2)
ax.add_patch(Rectangle((3, 0), 0.9, len(genes), facecolor="none", edgecolor="#dddddd", lw=0.4))
ax.set_xlim(0, 4.6); ax.set_ylim(0, len(genes))
ax.set_yticks(np.arange(len(genes))+0.5); ax.set_yticklabels(genes, fontsize=5.2)
ax.set_xticks([0.5,1.5,2.5,3.45])
ax.set_xticklabels(["ST\nletro","ST\nana","ChEMBL","PPI\ndeg"], fontsize=4.8)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_title("Core-gene evidence", fontsize=6.5)
fig.tight_layout(); save(fig, T3+"figures/Fig3F_evidence_matrix")

# ---------- Fig6E GSE201774 FPKM 排名 ----------
fp = pd.read_csv(T3+"results/tables/S2b_GSE201774_core_fpkm.csv", index_col=0)
med = fp.median(axis=0).sort_values()
med8 = med.tail(8)
fig, ax = plt.subplots(figsize=(56/25.4, 52/25.4))
cols = [C_HI if g == "VPS35" else "#B8C4CE" for g in med8.index]
ax.barh(range(len(med8)), med8.values, color=cols, height=0.62)
for i, (g, v) in enumerate(med8.items()):
    ax.text(v + med8.max()*0.02, i, f"{v:.1f}", va='center', fontsize=6.0)
ax.set_yticks(range(len(med8))); ax.set_yticklabels(med8.index, fontsize=6.5)
ax.set_xlabel("Median FPKM (n = 16 MSCs)", fontsize=6.5)
ax.set_xlim(0, med8.max()*1.22); ax.tick_params(axis="x", labelsize=6.0)
ax.set_title("GSE201774 expression ranking", fontsize=6.5)
fig.subplots_adjust(left=0.30, bottom=0.24, right=0.97, top=0.90)
save(fig, T3+"figures/Fig6E_GSE201774_fpkm")

# ---------- Fig6F VPS35 共表达模块相关矩阵 ----------
co = pd.read_csv(T3+"results/tables/S6_VPS35_coexpression.csv")
col = co.columns.tolist(); print("co 列:", col)
gene_col = "gene" if "gene" in col else col[0]
top = co.sort_values(co.columns[-1], ascending=False).head(8) if co.shape[1]>1 else co.head(8)
topgenes = ["VPS35"] + [g for g in top[gene_col].tolist() if g != "VPS35"][:7]
mat = pd.read_csv(T3+"results/tables/S2_GSE35958_gene_matrix.csv", index_col=0)
have = [g for g in topgenes if g in mat.columns]
X = mat[have]
C = X.corr(method="spearman").values
fig, ax = plt.subplots(figsize=(2.4, 2.4))
cmap = plt.get_cmap("RdBu_r"); norm = Normalize(-1, 1)
for i in range(len(have)):
    for j in range(len(have)):
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(norm(C[i,j])), edgecolor="none"))
ax.set_xlim(0, len(have)); ax.set_ylim(0, len(have))
ax.set_xticks(np.arange(len(have))+0.5); ax.set_yticks(np.arange(len(have))+0.5)
ax.set_xticklabels(have, rotation=90, fontsize=6.0); ax.set_yticklabels(have, fontsize=6.0)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_title("VPS35 co-expression module (Spearman)", fontsize=6.5)
fig.subplots_adjust(left=0.30, bottom=0.28, right=0.98, top=0.90)
save(fig, T3+"figures/Fig6F_coexpression_matrix")
print("四个新面板完成")
