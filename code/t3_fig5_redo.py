import matplotlib.text as _mtext
def FONT_FLOOR(fig, floor=6.5):
    for _t in fig.findobj(_mtext.Text):
        if _t.get_fontsize() < floor: _t.set_fontsize(floor)
    return fig

"""Fig5 改为 GSE56815（n=80）WGCNA 口径；Fig4D 重叠重算"""
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

sft = pd.read_csv(T3+"results/tables/S7_sft_fit.csv")
tr = pd.read_csv(T3+"results/tables/S7_wgcna_module_trait.csv")
sizes = pd.read_csv(T3+"results/tables/S7b_module_sizes.csv")
hub = pd.read_csv(T3+"results/tables/S7b_hub_genes.csv")
assign = pd.read_csv(T3+"results/tables/S7_wgcna_module_assignment.csv")
expr = pd.read_csv(T3+"results/tables/S2_GSE56815_gene_matrix.csv", index_col=0)
grp = pd.read_csv(T3+"results/tables/S2_GSE56815_groups.csv")
grp = grp[grp.group.isin(["low","high"])].sort_values("group")

# A 软阈值
fig, ax = plt.subplots(figsize=(2.2, 2.0))
ax.plot(sft["Power"], sft["SFT.R.sq"], marker='o', ms=3, lw=0.8, color="#3C5488")
ax.axhline(0.85, color=C_HI, lw=0.6, ls='--')
pw = int(sft.loc[sft["SFT.R.sq"].idxmax(), "Power"])
ax.axvline(pw, color=GREY, lw=0.6, ls='--')
ax.set_xlabel("Soft threshold (power)"); ax.set_ylabel("Scale-free R²")
ax.set_title(f"Scale-free fit (power = {pw})", fontsize=6.2)
fig.tight_layout(); save(fig, T3+"figures/Fig5A_sft_curve")

# B 模块大小
s2 = sizes[sizes.module > 0].sort_values("size", ascending=False)
fig, ax = plt.subplots(figsize=(2.4, 2.0))
cols = [C_HI if m in (4,5) else "#9DB2C8" for m in s2["module"]]
ax.bar(range(len(s2)), s2["size"], color=cols, width=0.62)
ax.set_xticks(range(len(s2))); ax.set_xticklabels([f"M{int(m)}" for m in s2["module"]], fontsize=4.6, rotation=45)
ax.set_ylabel("Genes (n)"); ax.set_title(f"WGCNA modules (n = {len(s2)})", fontsize=6.2)
fig.tight_layout(); save(fig, T3+"figures/Fig5B_module_hist")

# C 模块-性状
t2 = tr.sort_values("r")
fig, ax = plt.subplots(figsize=(2.4, 2.1))
cols = [C_HI if p < 0.01 else ("#F0A07C" if p < 0.05 else "#C9C9C9") for p in t2["padj"]]
ax.barh(range(len(t2)), t2["r"], color=cols, height=0.62)
for i, (_, r) in enumerate(t2.iterrows()):
    lab = "***" if r["padj"] < 0.001 else ("**" if r["padj"] < 0.01 else ("*" if r["padj"] < 0.05 else ""))
    txt = f"{r['r']:.2f}{lab}"
    ax.text(r["r"] + (0.01 if r["r"] >= 0 else -0.01), i, txt, va='center',
            ha='left' if r["r"] >= 0 else 'right', fontsize=4.6,
            fontweight='bold' if r["module"] in (4, 5) else 'normal')
ax.axvline(0, color='#333', lw=0.5)
ax.set_yticks(range(len(t2))); ax.set_yticklabels([f"M{int(m)}" for m in t2["module"]], fontsize=4.8)
ax.set_xlabel("Module–BMD correlation (r)", fontsize=5.2)
ax.set_title("M4/M5: padj = 0.0066", fontsize=6.0)
fig.tight_layout(); save(fig, T3+"figures/Fig5_WGCNA_module_trait")

# D M4/M5 top 基因热图
top = pd.concat([hub[hub.module==m].head(10) for m in (4,5)])["gene"].tolist()
top = [g for g in top if g in expr.columns]
E = np.log2(expr.loc[grp["sample"], top] + 1)
Z = ((E - E.mean())/E.std()).T
fig = plt.figure(figsize=(2.6, 2.4))
ax = fig.add_axes([0.26, 0.16, 0.58, 0.74])
cmap = plt.get_cmap("RdBu_r"); norm = Normalize(-2, 2)
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(norm(Z.values[i, j])), edgecolor="none"))
ax.set_xlim(0, Z.shape[1]); ax.set_ylim(0, Z.shape[0]); ax.set_xticks([])
ax.set_yticks(np.arange(len(top))+0.5); ax.set_yticklabels(top, fontsize=3.6); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
axb = fig.add_axes([0.26, 0.12, 0.58, 0.03])
for j in range(Z.shape[1]):
    c = "#E64B35" if grp["group"].values[j] == "low" else "#4DBBD5"
    axb.add_patch(Rectangle((j, 0), 1, 1, facecolor=c, edgecolor="none"))
axb.set_xlim(0, Z.shape[1]); axb.set_ylim(0, 1); axb.set_xticks([]); axb.set_yticks([])
axb.set_xlabel("low | high BMD (n = 80)", fontsize=4.8, labelpad=1)
cax = fig.add_axes([0.86, 0.16, 0.03, 0.74]); cax.set_xlim(0,1); cax.set_ylim(-2,2); cax.set_xticks([])
for k in range(40):
    v0 = -2 + 4*k/40
    cax.add_patch(Rectangle((0, v0), 1, 4/40, facecolor=cmap(norm(v0+2/40)), edgecolor="none"))
cax.set_yticks([-2,0,2]); cax.set_yticklabels(["-2","0","2"], fontsize=4.4)
for s in ("top","bottom"): cax.spines[s].set_visible(False)
ax.set_title("M4/M5 top genes (GSE56815)", fontsize=6.0)
FONT_FLOOR(fig)
fig.savefig(T3+"figures/Fig5D_module_heatmap.pdf", format="pdf"); FONT_FLOOR(fig)
fig.savefig(T3+"figures/Fig5D_module_heatmap.png", dpi=300); plt.close(fig)

# E 模块大小 + 显著性
mt = tr.merge(sizes, on="module").sort_values("size", ascending=True)
fig, ax = plt.subplots(figsize=(2.1, 2.0))
cols = ["#E64B35" if p < 0.01 else ("#F0A07C" if p < 0.05 else "#C9C9C9") for p in mt["padj"]]
ax.barh(range(len(mt)), mt["size"], color=cols, height=0.62)
for i, (_, r) in enumerate(mt.iterrows()):
    lab = f"{int(r['size'])}"
    if r["padj"] < 0.001: lab += " ***"
    elif r["padj"] < 0.01: lab += " **"
    elif r["padj"] < 0.05: lab += " *"
    ax.text(r["size"]*1.02, i, lab, va='center', fontsize=4.4)
ax.set_yticks(range(len(mt))); ax.set_yticklabels([f"M{int(m)}" for m in mt["module"]], fontsize=4.6)
ax.set_xlabel("Module size (genes)", fontsize=5.0); ax.set_xlim(0, mt["size"].max()*1.25)
ax.set_title("Module size (trait association)", fontsize=6.0)
fig.tight_layout(); save(fig, T3+"figures/Fig5E_module_sizes")

# F M4 hub kME
h4 = hub[hub.module == 4].sort_values("kME").tail(10)
fig, ax = plt.subplots(figsize=(2.0, 2.0))
ax.barh(range(len(h4)), h4["kME"], color="#8491B4", height=0.6)
for i, (_, r) in enumerate(h4.iterrows()):
    ax.text(r["kME"]+0.01, i, f"{r['kME']:.2f}", va='center', fontsize=4.4)
ax.set_yticks(range(len(h4))); ax.set_yticklabels(h4["gene"], fontsize=4.8)
ax.set_xlabel("kME (M4)", fontsize=5.0); ax.set_xlim(0, 1.05)
ax.set_title("M4 hub genes", fontsize=6.2)
fig.tight_layout(); save(fig, T3+"figures/Fig5F_M4_hubs")

# Fig4D M4/M5 与 DEG 重叠（GSE56815 模块 vs GSE35958 DEG）
lim = pd.read_csv(T3+"results/tables/S2_GSE35958_limma_all.csv")
deg = set(lim[(lim['adj.P.Val'] < 0.05) & (lim['logFC'].abs() > 1)]["gene"])
m45 = set(assign[assign.module.isin([4,5])]["gene"])
ov = len(m45 & deg)
fig, ax = plt.subplots(figsize=(2.0, 2.0))
bars = ax.bar([0, 1], [len(m45), len(deg)], color=["#8491B4", "#C9C9C9"], width=0.5)
ax.bar([2], [ov], color=C_HI, width=0.5)
for x, v in zip([0,1,2], [len(m45), len(deg), ov]):
    ax.text(x, v*1.03, f"{v}", ha='center', fontsize=5.0, fontweight='bold')
ax.set_xticks([0,1,2]); ax.set_xticklabels(["M4/M5", "DEGs", "overlap"], fontsize=5.0)
ax.set_ylabel("Genes (n)", fontsize=5.2)
ax.set_title(f"Module–DEG overlap: {ov}", fontsize=6.2)
fig.tight_layout(); save(fig, T3+"figures/Fig4D_module_deg_overlap")
print("Fig5 A-F 已改为 GSE56815 口径；Fig4D 重叠重算 =", ov)
