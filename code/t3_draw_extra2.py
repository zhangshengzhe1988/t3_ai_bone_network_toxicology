import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import matplotlib.text as _mtext
def FONT_FLOOR(fig, floor=6.5):
    for _t in fig.findobj(_mtext.Text):
        if _t.get_fontsize() < floor: _t.set_fontsize(floor)
    return fig

"""T3：Fig2A/D 矢量火山、Fig4C 矢量热图、Fig4E/F 与 Fig5E/F 新增面板"""
import sys; sys.path.insert(0, "code")
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from nature_style import save
T3 = "选题3_网络毒理学_AI骨质疏松/"
C_HI=(0.902,0.294,0.208); C_LO=(0.302,0.686,0.835); GREY='#7E8AA2'
plt.rcParams.update({"font.family":"Arial","font.size":6,"pdf.fonttype":42,"axes.linewidth":0.5,
                     "axes.spines.top":False,"axes.spines.right":False})
core9 = ["CYP19A1","HDAC1","HDAC6","ACHE","ADRB1","ADRB2","ADRB3","HTR2A","HTR2C"]

# ---------- Fig2A / Fig2D 矢量火山（去 rasterized） ----------
def volcano(df, fc_col, p_col, fc_thr, title, out, xlab):
    sig = df[(df[p_col] < 0.05) & (df[fc_col].abs() > fc_thr)]
    fig, ax = plt.subplots(figsize=(2.4, 2.2))
    ax.scatter(df[fc_col], -np.log10(df[p_col]), s=1.2, color='#C9C9C9', linewidths=0)
    ax.scatter(sig[fc_col], -np.log10(sig[p_col]), s=1.2, color='#8491B4', linewidths=0)
    for gname in core9:
        r = df[df['gene'] == gname]
        if len(r):
            ax.scatter(r[fc_col], -np.log10(r[p_col]), s=9, color=C_HI, zorder=3,
                       edgecolor='white', linewidth=0.3)
    ax.axhline(-np.log10(0.05), color='black', lw=0.4, ls='--')
    for v in (fc_thr, -fc_thr): ax.axvline(v, color='black', lw=0.4, ls='--')
    ax.set_xlabel(xlab); ax.set_ylabel("-log10(P)")
    ax.set_title(title, fontsize=6.2)
    fig.tight_layout(); save(fig, out)

g1 = pd.read_csv(T3+"results/tables/S2_GSE35958_limma_all.csv").dropna(subset=['P.Value'])
volcano(g1, 'logFC', 'P.Value', 1.0, "GSE35958: 1,447 DEGs", T3+"figures/Fig2A_GSE35958_volcano", "log2FC (osteoporotic vs control)")
g2 = pd.read_csv(T3+"results/tables/S2b_GSE13850_limma.csv").dropna(subset=['P.Value'])
volcano(g2, 'logFC', 'P.Value', 0.58, "GSE13850: 28 DEGs", T3+"figures/Fig2D_GSE13850_volcano", "log2FC (high vs low BMD)")

# ---------- Fig4C 核心基因热图（全矢量） ----------
mat = pd.read_csv(T3+"results/tables/S2_GSE35958_gene_matrix.csv", index_col=0)
grp = pd.read_csv(T3+"results/tables/S2_GSE35958_groups.csv")
grp = grp[grp.group.isin(["high","low"])].sort_values("group")
genes = [g for g in core9 if g in mat.columns]
X = mat.loc[grp["sample"].tolist(), genes]
Z = ((X - X.mean())/X.std()).T
fig = plt.figure(figsize=(2.0, 2.2))
ax = fig.add_axes([0.20, 0.16, 0.62, 0.76])
cmap = plt.get_cmap("RdBu_r"); norm = Normalize(-2, 2)
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(norm(Z.values[i, j])), edgecolor="none"))
ax.set_xlim(0, Z.shape[1]); ax.set_ylim(0, Z.shape[0])
ax.set_xticks([]); ax.set_yticks(np.arange(len(genes))+0.5); ax.set_yticklabels(genes, fontsize=5.0)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
axb = fig.add_axes([0.20, 0.115, 0.62, 0.032])
for j in range(Z.shape[1]):
    c = "#E64B35" if grp["group"].values[j] == "low" else "#4DBBD5"
    axb.add_patch(Rectangle((j, 0), 1, 1, facecolor=c, edgecolor="none"))
axb.set_xlim(0, Z.shape[1]); axb.set_ylim(0, 1); axb.set_xticks([]); axb.set_yticks([])
axb.set_xlabel("low | high BMD", fontsize=5.0, labelpad=1)
cax = fig.add_axes([0.85, 0.16, 0.035, 0.76]); cax.set_xlim(0,1); cax.set_ylim(-2,2); cax.set_xticks([])
for k in range(40):
    v0 = -2 + 4*k/40
    cax.add_patch(Rectangle((0, v0), 1, 4/40, facecolor=cmap(norm(v0+2/40)), edgecolor="none"))
cax.set_yticks([-2,0,2]); cax.set_yticklabels(["-2","0","2"], fontsize=4.6)
cax.yaxis.set_label_position("right"); cax.yaxis.tick_right()
cax.set_ylabel("z-score", fontsize=5.2, rotation=90, labelpad=4)
for s in ("top","bottom"): cax.spines[s].set_visible(False)
ax.set_title("Core genes (GSE35958)", fontsize=6.2)
FONT_FLOOR(fig)
fig.savefig(T3+"figures/Fig4D_core_heatmap.pdf", format="pdf"); FONT_FLOOR(fig)
fig.savefig(T3+"figures/Fig4D_core_heatmap.png", dpi=300); plt.close(fig)
print("Fig2A/D、Fig4C 矢量重绘完成")

# ---------- Fig4E 核心基因-通路气泡 ----------
go = pd.read_csv(T3+"results/tables/S6_go_bp.csv").head(10)
rows, cols, sizes, colors = [], [], [], []
for i, (_, r) in enumerate(go.iterrows()):
    for g in str(r['geneID']).split('/'):
        if g in core9:
            rows.append(g); cols.append(r['Description'][:38]); sizes.append(-np.log10(r['p.adjust']))
terms = list(dict.fromkeys(cols)); genes_in = [g for g in core9 if g in set(rows)]
Zi = np.zeros((len(genes_in), len(terms)))
for (g, t, s) in zip(rows, cols, sizes):
    Zi[genes_in.index(g), terms.index(t)] = s
fig, ax = plt.subplots(figsize=(2.6, 2.0))
cmap = plt.get_cmap("RdYlBu")
for i, g in enumerate(genes_in):
    for j, t in enumerate(terms):
        if Zi[i, j] > 0:
            ax.scatter(j, i, s=Zi[i, j]*9, color=cmap(Normalize(4, 12)(Zi[i, j])), edgecolor='#333', linewidth=0.3)
short = [t if len(t) <= 26 else t[:24] + "…" for t in terms]
ax.set_xticks(range(len(terms))); ax.set_xticklabels(short, rotation=38, ha='right', fontsize=4.2)
ax.set_yticks(range(len(genes_in))); ax.set_yticklabels(genes_in, fontsize=4.8)
ax.tick_params(length=0)
ax.set_title("Core-gene × GO-BP", fontsize=6.2)
ax.set_xlim(-0.6, len(terms)-0.4); ax.set_ylim(-0.6, len(genes_in)-0.4)
cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02); cb.set_label("−log10(adj. P)", fontsize=6.5)
cb.ax.tick_params(labelsize=6.0)
ax.legend(*sc.legend_elements(num=3), title="genes", loc="lower right", fontsize=6.0,
          title_fontsize=6.5, frameon=False, handletextpad=0.4)
fig.tight_layout(); save(fig, T3+"figures/Fig4E_go_gene_bubble")

# ---------- Fig4F 核心基因表达点图 ----------
fig, ax = plt.subplots(figsize=(2.4, 2.0))
for i, g in enumerate(genes):
    z = Z.values[i]
    ax.scatter(np.arange(len(z)), [i]*len(z), s=np.abs(z)*14, c=z, cmap="RdBu_r", vmin=-2, vmax=2, edgecolor="#555", linewidth=0.25)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=5.5)
ax.set_xticks(range(Z.shape[1])); ax.set_xticklabels([f"S{j+1}" for j in range(Z.shape[1])], fontsize=4.6)
ax.set_xlabel("GSE35958 samples", fontsize=5.0)
ax.set_title("Core-gene expression (z)", fontsize=6.2)
ax.set_xlim(-0.6, Z.shape[1]-0.4)
cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02); cb.set_label("z-score", fontsize=6.5)
cb.ax.tick_params(labelsize=6.0)
fig.tight_layout(); save(fig, T3+"figures/Fig4F_core_dotplot")

# ---------- Fig5E 模块大小与显著性 ----------
mt = pd.read_csv(T3+"results/tables/S7_module_trait.csv").sort_values("size", ascending=True)
fig, ax = plt.subplots(figsize=(2.1, 2.0))
cols = ["#E64B35" if p < 0.001 else ("#F0A07C" if p < 0.05 else "#C9C9C9") for p in mt["p"]]
ax.barh(range(len(mt)), mt["size"], color=cols, height=0.62)
for i, (_, r) in enumerate(mt.iterrows()):
    lab = f"{int(r['size'])}"
    if r["p"] < 0.001: lab += " ***"
    elif r["p"] < 0.05: lab += " *"
    ax.text(r["size"]+20, i, lab, va='center', fontsize=4.6)
ax.set_yticks(range(len(mt))); ax.set_yticklabels([f"ME{int(m)}" for m in mt["module"]], fontsize=4.8)
ax.set_xlabel("Module size (genes)", fontsize=5.2); ax.set_xlim(0, mt["size"].max()*1.22)
ax.set_title("Module size (trait association)", fontsize=6.0)
fig.tight_layout(); save(fig, T3+"figures/Fig5E_module_sizes")

# ---------- Fig5F M4 hub kME ----------
hub = pd.read_csv(T3+"results/tables/S7_module_hub_genes.csv")
h4 = hub[hub.module == 4].sort_values("kME", ascending=True).tail(10)
fig, ax = plt.subplots(figsize=(2.0, 2.0))
ax.barh(range(len(h4)), h4["kME"], color="#8491B4", height=0.6)
for i, (_, r) in enumerate(h4.iterrows()):
    ax.text(r["kME"]+0.008, i, f"{r['kME']:.2f}", va='center', fontsize=4.6)
ax.set_yticks(range(len(h4))); ax.set_yticklabels(h4["gene"], fontsize=5.0)
ax.set_xlabel("kME (M4)", fontsize=5.2); ax.set_xlim(0, 0.95)
ax.set_title("M4 hub genes", fontsize=6.2)
fig.tight_layout(); save(fig, T3+"figures/Fig5F_M4_hubs")
print("Fig4E/F、Fig5E/F 新增完成")
