#!/usr/bin/env python3
"""重建 Fig 4（WGCNA, GSE56815）六个面板，统一配色与排版。

数据来源（与正文一致的一套 WGCNA 结果）：
  S7b_module_sizes.csv, S7_wgcna_module_trait.csv, S7_wgcna_hub_genes.csv,
  S2_GSE56815_gene_matrix.csv + groups, /tmp/fitIndices.csv（pickSoftThreshold 复算）
"""
import sys
sys.path.insert(0, "code")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from nature_style import save, C_HIGH, C_LOW

T3 = "选题3_网络毒理学_AI骨质疏松/"
RED, BLUE, GREY = C_HIGH, "#4DBBD5", "#C8C8C8"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})

sizes = pd.read_csv(T3 + "results/tables/S7b_module_sizes.csv")
trait = pd.read_csv(T3 + "results/tables/S7_wgcna_module_trait.csv")
hubs = pd.read_csv(T3 + "results/tables/S7_wgcna_hub_genes.csv")
fit = pd.read_csv("/tmp/fitIndices.csv")

# ---------------- A 软阈值 ----------------
fig, ax = plt.subplots(figsize=(68/25.4, 52/25.4))
ax.plot(fit["Power"], fit["SFT.R.sq"], "o-", color=BLUE, lw=0.9, ms=3.0)
best = fit.loc[fit["SFT.R.sq"].idxmax()]
ax.scatter([16], [float(fit.loc[fit.Power == 16, "SFT.R.sq"].iloc[0])], s=22, color=RED, zorder=4)
ax.axhline(0.8, color="#7E8AA2", lw=0.6, ls="--")
ax.text(1.2, 0.815, "R² = 0.8", fontsize=5.6, color="#555555")
ax.annotate("power 16", xy=(16, float(fit.loc[fit.Power == 16, "SFT.R.sq"].iloc[0])),
            xytext=(7.5, 0.55), fontsize=5.8, color=RED,
            arrowprops=dict(arrowstyle="-", lw=0.5, color=RED))
ax.set_xlabel("Soft-thresholding power", fontsize=6.4)
ax.set_ylabel("Scale-free fit, signed R²", fontsize=6.4)
ax.tick_params(labelsize=5.8)
ax.set_title("Scale-free topology fit", fontsize=7.0)
ax.set_ylim(0, 1.0)
fig.subplots_adjust(left=0.22, right=0.97, bottom=0.19, top=0.88)
save(fig, T3 + "figures/Fig5A_sft_curve")

# ---------------- B 模块大小 ----------------
s = sizes.copy()
s["label"] = ["M" + str(int(m)) for m in s["module"]]
fig, ax = plt.subplots(figsize=(76/25.4, 52/25.4))
cols = [RED if m in (4, 5) else "#9DB4C9" for m in s["module"]]
bars = ax.bar(range(len(s)), s["size"], color=cols, width=0.72)
for i, (b, v, m) in enumerate(zip(bars, s["size"], s["module"])):
    if m in (4, 5):
        ax.text(b.get_x()+b.get_width()/2, v+60, f"{int(v)}", ha="center", fontsize=5.6, color=RED)
ax.set_xticks(range(len(s))); ax.set_xticklabels(s["label"], fontsize=5.6, rotation=45, ha="right")
ax.set_ylabel("Genes (n)", fontsize=6.4); ax.tick_params(axis="y", labelsize=5.8)
ax.set_title("Module sizes (13 modules)", fontsize=7.0)
fig.subplots_adjust(left=0.17, right=0.97, bottom=0.24, top=0.88)
save(fig, T3 + "figures/Fig5B_module_hist")

# ---------------- C 模块-性状相关 ----------------
t = trait.copy().sort_values("r")
t["label"] = ["M" + str(int(m)) for m in t["module"]]
fig, ax = plt.subplots(figsize=(62/25.4, 56/25.4))
cols = [RED if (p < 0.05) else GREY for p in t["padj"]]
bars = ax.barh(range(len(t)), t["r"], color=cols, height=0.68)
for y_, (r, p) in enumerate(zip(t["r"], t["padj"])):
    star = "**" if p < 0.01 else ("*" if p < 0.05 else "")
    xoff = 0.012 if r >= 0 else -0.012
    ax.text(r + xoff, y_, f"{r:.2f}{star}", va="center", ha="left" if r >= 0 else "right", fontsize=5.2)
ax.axvline(0, color="#333333", lw=0.5)
ax.set_yticks(range(len(t))); ax.set_yticklabels(t["label"], fontsize=5.6)
ax.set_xlabel("Module–trait correlation (r; low BMD = 1)", fontsize=6.0)
ax.tick_params(axis="x", labelsize=5.8)
ax.set_xlim(t["r"].min()-0.16, t["r"].max()+0.16)
ax.set_title("Module–trait association", fontsize=7.0)
ax.text(0.02, 0.94, "M4/M5: BH adj. P = 0.0066", transform=ax.transAxes, ha="left", fontsize=5.2, color="#444444")
fig.subplots_adjust(left=0.20, right=0.97, bottom=0.17, top=0.88)
save(fig, T3 + "figures/Fig5_WGCNA_module_trait")

# ---------------- D M4/M5 hub 基因热图 ----------------
mat = pd.read_csv(T3 + "results/tables/S2_GSE56815_gene_matrix.csv", index_col=0)
grp = pd.read_csv(T3 + "results/tables/S2_GSE56815_groups.csv")
grp = grp[grp["group"].isin(["high", "low"])]
order = list(grp.sort_values("group", ascending=False)["sample"])
genes = list(hubs["gene"])
sub = mat.loc[order, genes].astype(float)
z = (sub - sub.mean()) / sub.std()
fig, ax = plt.subplots(figsize=(66/25.4, 58/25.4))
im = ax.imshow(z.T.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
               extent=[0, len(order), len(genes), 0])
ax.set_yticks(np.arange(len(genes)) + 0.5); ax.set_yticklabels(genes, fontsize=4.6)
n_high = int((grp["group"] == "high").sum()); n_low = len(order) - n_high
ax.axvline(n_high, color="k", lw=0.6)
ax.set_xticks([n_high/2, n_high + n_low/2]); ax.set_xticklabels(["high BMD", "low BMD"], fontsize=5.6)
ax.tick_params(axis="x", length=0)
ax.set_title("M4/M5 hub genes (GSE56815, n = 80)", fontsize=6.8)
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02); cb.set_label("z-score", fontsize=5.8)
cb.ax.tick_params(labelsize=5.2)
fig.subplots_adjust(left=0.26, right=0.92, bottom=0.14, top=0.88)
save(fig, T3 + "figures/Fig5D_module_heatmap")

# ---------------- E 模块大小 vs 性状显著性 ----------------
m = sizes.rename(columns={"size": "n_genes"}).merge(trait[["module", "padj"]], on="module", how="left")
m["label"] = ["M" + str(int(x)) for x in m["module"]]
m = m.sort_values("n_genes", ascending=True)
fig, ax = plt.subplots(figsize=(62/25.4, 52/25.4))
cols = [RED if (p < 0.05) else GREY for p in m["padj"].fillna(1)]
ax.barh(range(len(m)), m["n_genes"], color=cols, height=0.68)
for y_, (v, p) in enumerate(zip(m["n_genes"], m["padj"].fillna(1))):
    star = "**" if p < 0.01 else ("*" if p < 0.05 else "")
    ax.text(v + 40, y_, f"{int(v)}{star}", va="center", fontsize=5.2)
ax.set_yticks(range(len(m))); ax.set_yticklabels(m["label"], fontsize=5.6)
ax.set_xlabel("Module size (genes)", fontsize=6.2); ax.tick_params(axis="x", labelsize=5.8)
ax.set_xlim(0, m["n_genes"].max()*1.22)
ax.set_title("Module size vs trait association", fontsize=6.8)
ax.text(0.98, 0.05, "* adj. P < 0.05; ** < 0.01", transform=ax.transAxes, ha="right", fontsize=4.8, color="#444444")
fig.subplots_adjust(left=0.20, right=0.97, bottom=0.17, top=0.88)
save(fig, T3 + "figures/Fig5E_module_sizes")

# ---------------- F M4 hub kME ----------------
h4 = hubs[hubs["module"] == 4].sort_values("kME")
fig, ax = plt.subplots(figsize=(62/25.4, 52/25.4))
ax.barh(range(len(h4)), h4["kME"], color="#8491B4", height=0.66)
for y_, v in enumerate(h4["kME"]):
    ax.text(v + 0.006, y_, f"{v:.2f}", va="center", fontsize=5.4)
ax.set_yticks(range(len(h4))); ax.set_yticklabels(h4["gene"], fontsize=5.6)
ax.set_xlabel("kME (M4)", fontsize=6.2); ax.tick_params(axis="x", labelsize=5.8)
ax.set_xlim(0, 0.92)
ax.set_title("Top M4 hub genes", fontsize=7.0)
fig.subplots_adjust(left=0.24, right=0.97, bottom=0.17, top=0.88)
save(fig, T3 + "figures/Fig5F_M4_hubs")
print("Fig4 panels rebuilt")
