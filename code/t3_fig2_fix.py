"""T3 Fig2 修正（log2 口径）+ 扩充 F/G 面板"""
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

mat = pd.read_csv(T3+"results/tables/S2_GSE35958_gene_matrix.csv", index_col=0)
LX = np.log2(mat + 1)                      # 正确口径：log2(x+1)
grp = pd.read_csv(T3+"results/tables/S2_GSE35958_groups.csv")
grp = grp[grp.group.isin(["high","low"])].sort_values("group")
lim = pd.read_csv(T3+"results/tables/S2_GSE35958_limma_all.csv")   # 已替换为 log2 版

# ---- Fig2A 火山（log2 口径，1,447 DEGs） ----
sig = lim[(lim['adj.P.Val'] < 0.05) & (lim['logFC'].abs() > 1)]
fig, ax = plt.subplots(figsize=(2.4, 2.2))
ax.scatter(lim['logFC'], -np.log10(lim['P.Value']), s=1.2, color='#C9C9C9', linewidths=0)
ax.scatter(sig['logFC'], -np.log10(sig['P.Value']), s=1.2, color='#8491B4', linewidths=0)
for g in core9:
    r = lim[lim['gene'] == g]
    if len(r): ax.scatter(r['logFC'], -np.log10(r['P.Value']), s=9, color=C_HI, zorder=3, edgecolor='white', linewidth=0.3)
ax.axhline(-np.log10(0.05), color='k', lw=0.4, ls='--')
for v in (1, -1): ax.axvline(v, color='k', lw=0.4, ls='--')
ax.set_xlabel("log2FC (low vs high BMD)"); ax.set_ylabel("-log10(P)")
ax.set_title("GSE35958: 1,447 DEGs", fontsize=6.2)
fig.tight_layout(); save(fig, T3+"figures/Fig2A_GSE35958_volcano")

# ---- Fig2C 核心基因表达（log2 口径 + adj p 标注） ----
fig, axes = plt.subplots(1, 2, figsize=(2.8, 2.2))
for ax, g in zip(axes, ["CYP19A1", "VPS35"]):
    hi = LX.loc[grp[grp.group=="high"]["sample"], g]
    lo = LX.loc[grp[grp.group=="low"]["sample"], g]
    bp = ax.boxplot([hi, lo], positions=[0,1], widths=0.5, patch_artist=True,
                    medianprops=dict(color='#333', lw=0.8))
    for patch, col in zip(bp['boxes'], [C_LO, C_HI]): patch.set_facecolor(col); patch.set_alpha(0.55)
    rng = np.random.default_rng(3)
    for i, v in enumerate([hi, lo]):
        ax.scatter(i + rng.uniform(-0.1, 0.1, len(v)), v, s=7, color='#333', zorder=4, linewidths=0)
    row = lim[lim['gene'] == g].iloc[0]
    ax.set_xticks([0,1]); ax.set_xticklabels(["high", "low"], fontsize=5.4)
    ax.set_ylabel("log2 expression", fontsize=5.4)
    ax.set_title(f"{g}\nadj. p = {row['adj.P.Val']:.3f}", fontsize=6.0)
    ax.set_xlim(-0.5, 1.5)
fig.tight_layout(); save(fig, T3+"figures/Fig2C_core_expression")

# ---- Fig2E 跨队列核心基因 log2FC ----
g2 = pd.read_csv(T3+"results/tables/S2b_GSE13850_limma.csv")
rows = []
for cohort, df in [("GSE35958", lim), ("GSE13850", g2)]:
    for g in core9 + ["VPS35"]:
        r = df[df['gene'] == g]
        if len(r): rows.append((cohort, g, r.iloc[0]['logFC'], r.iloc[0]['adj.P.Val']))
fc = pd.DataFrame(rows, columns=["cohort","gene","log2FC","adjP"])
piv = fc.pivot(index="gene", columns="cohort", values="log2FC").loc[core9 + ["VPS35"]]
padj = fc.pivot(index="gene", columns="cohort", values="adjP").loc[core9 + ["VPS35"]]
fig, ax = plt.subplots(figsize=(2.4, 2.4))
cmap = plt.get_cmap("RdBu_r"); norm = Normalize(-2.5, 2.5)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = piv.values[i, j]
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(norm(v)), edgecolor="white", lw=0.6))
        star = "*" if padj.values[i, j] < 0.05 else ""
        ax.text(j+0.5, i+0.5, f"{v:.2f}{star}", ha='center', va='center', fontsize=4.6,
                color='white' if abs(v) > 1.4 else '#222')
ax.set_xlim(0, piv.shape[1]); ax.set_ylim(0, piv.shape[0])
ax.set_xticks(np.arange(piv.shape[1])+0.5); ax.set_xticklabels(piv.columns, fontsize=5.0)
ax.set_yticks(np.arange(piv.shape[0])+0.5); ax.set_yticklabels(piv.index, fontsize=4.8)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_title("Core-gene log2FC across cohorts", fontsize=6.0)
fig.tight_layout(); save(fig, T3+"figures/Fig2E_core_logFC")

# ---- Fig2F 四队列 DEG 计数 ----
cohorts = ["GSE35958\n(MSC)", "GSE13850\n(B cells)", "GSE56815\n(PBMC)", "GSE56814\n(PBMC)"]
counts = [1447, 28, 391, 0]
rules = ["|log2FC|>1 & adj.p<0.05", "|log2FC|>0.58 & adj.p<0.05", "adj.p<0.05", "adj.p<0.05"]
fig, ax = plt.subplots(figsize=(2.4, 2.0))
cols = [C_HI, "#F0A07C", "#F0A07C", "#C9C9C9"]
bars = ax.bar(range(4), counts, color=cols, width=0.55)
for i, v in enumerate(counts):
    ax.text(i, v + max(counts)*0.03, f"{v}", ha='center', fontsize=5.0, fontweight='bold')
ax.set_xticks(range(4)); ax.set_xticklabels(cohorts, fontsize=4.6)
ax.set_ylabel("DEGs (n)", fontsize=5.2); ax.set_ylim(0, max(counts)*1.2)
ax.set_title("DEGs across four cohorts", fontsize=6.2)
for i, r in enumerate(rules):
    ax.text(i, -max(counts)*0.16, r if i < 2 else "adj.p<0.05", ha='center', fontsize=3.6, color=GREY)
fig.tight_layout(); save(fig, T3+"figures/Fig2F_cohort_degs")

# ---- Fig4C/4F 用 log2 口径重绘 ----
genes = [g for g in core9 if g in LX.columns]
X = LX.loc[grp["sample"].tolist(), genes]
Z = ((X - X.mean())/X.std()).T
fig = plt.figure(figsize=(2.0, 2.2))
ax = fig.add_axes([0.20, 0.16, 0.62, 0.76])
cmap = plt.get_cmap("RdBu_r"); norm = Normalize(-2, 2)
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(norm(Z.values[i, j])), edgecolor="none"))
ax.set_xlim(0, Z.shape[1]); ax.set_ylim(0, Z.shape[0]); ax.set_xticks([])
ax.set_yticks(np.arange(len(genes))+0.5); ax.set_yticklabels(genes, fontsize=5.0); ax.tick_params(length=0)
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
fig.savefig(T3+"figures/Fig4D_core_heatmap.pdf", format="pdf"); fig.savefig(T3+"figures/Fig4D_core_heatmap.png", dpi=300); plt.close(fig)
fig, ax = plt.subplots(figsize=(2.4, 2.0))
for i, g in enumerate(genes):
    z = Z.values[i]
    ax.scatter(np.arange(len(z)), [i]*len(z), s=np.abs(z)*14, c=z, cmap="RdBu_r", vmin=-2, vmax=2, edgecolor="#555", linewidth=0.25)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=5.0)
ax.set_xticks(range(Z.shape[1])); ax.set_xticklabels([f"S{j+1}" for j in range(Z.shape[1])], fontsize=4.6)
ax.set_xlabel("GSE35958 samples", fontsize=5.0); ax.set_title("Core-gene expression (z)", fontsize=6.2)
ax.set_xlim(-0.6, Z.shape[1]-0.4)
fig.tight_layout(); save(fig, T3+"figures/Fig4F_core_dotplot")
print("Fig2 A/C/E 修正 + F 新增；Fig4C/4F log2 口径重绘完成")
