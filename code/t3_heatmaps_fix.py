#!/usr/bin/env python3
"""重绘 Fig 4D（M4/M5 hub 热图）与 Fig 5F（VPS35 共表达 Spearman 矩阵），提高清晰度并加图例。"""
import sys
sys.path.insert(0, "code")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from nature_style import save, C_HIGH

T3 = "选题3_网络毒理学_AI骨质疏松/"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "axes.linewidth": 0.5})

# ---------------- Fig 4D：M4/M5 hub 热图（每模块 6 个，共 12 行） ----------------
hubs = pd.read_csv(T3 + "results/tables/S7_wgcna_hub_genes.csv")
mat = pd.read_csv(T3 + "results/tables/S2_GSE56815_gene_matrix.csv", index_col=0)
grp = pd.read_csv(T3 + "results/tables/S2_GSE56815_groups.csv")
grp = grp[grp["group"].isin(["high", "low"])].copy()

top = (hubs.sort_values(["module", "kME"], ascending=[True, False])
            .groupby("module").head(6).reset_index(drop=True))
genes = list(top["gene"]); modules = list(top["module"])

sub = mat.loc[grp["sample"], genes].astype(float)
z = (sub - sub.mean()) / sub.std()
grp["score"] = z[genes[:6]].mean(axis=1).values
grp["group"] = pd.Categorical(grp["group"], categories=["high", "low"], ordered=True)
grp = grp.sort_values(["group", "score"], ascending=[True, True])
order = list(grp["sample"])
Z = ((mat.loc[order, genes].astype(float) - sub.mean()) / sub.std()).T.values
n_high = int((grp["group"] == "high").sum()); n = len(order)

fig = plt.figure(figsize=(70/25.4, 56/25.4))
ax = fig.add_axes([0.36, 0.13, 0.48, 0.66])
im = ax.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
               extent=[0, n, len(genes), 0], interpolation="nearest")
ax.set_yticks(np.arange(len(genes)) + 0.5); ax.set_yticklabels(genes, fontsize=5.8)
ax.set_xticks([]); ax.tick_params(axis="y", length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
# 按模块给基因名着色
for lab, m in zip(ax.get_yticklabels(), modules):
    lab.set_color("#3C5488" if m == 4 else "#E64B35")

# 顶部分组条（蓝=high, 红=low）
axb = fig.add_axes([0.36, 0.805, 0.48, 0.028])
grp_val = np.array([[0 if g == "high" else 1 for g in grp["group"]]])
axb.imshow(grp_val, aspect="auto", cmap=plt.matplotlib.colors.ListedColormap(["#4DBBD5", C_HIGH]),
           extent=[0, n, 0, 1], vmin=0, vmax=1)
axb.set_xticks([]); axb.set_yticks([])
for s_ in axb.spines.values(): s_.set_visible(False)
axb.text(n_high/2, 1.9, "high BMD", ha="center", fontsize=5.2, color="#2C7A8C")
axb.text(n_high + (n-n_high)/2, 1.9, "low BMD", ha="center", fontsize=5.2, color="#B4432F")

# 色标
cax = fig.add_axes([0.845, 0.20, 0.022, 0.52])
cb = fig.colorbar(im, cax=cax); cb.set_label("z-score", fontsize=5.6, labelpad=1.5)
cb.ax.tick_params(labelsize=5.2); cb.outline.set_linewidth(0.4)

fig.text(0.60, 0.94, "M4/M5 hub genes (GSE56815, n = 80)", ha="center", fontsize=6.8)
fig.text(0.36, 0.045, "M4 labels blue, M5 labels red", fontsize=4.8, color="#444444")
save(fig, T3 + "figures/Fig5D_module_heatmap")

# ---------------- Fig 5F：VPS35 共表达 Spearman 矩阵 ----------------
co = pd.read_csv(T3 + "results/tables/S6_VPS35_coexpression.csv")
gcol = co.columns[0]
partners = [g for g in co.sort_values(co.columns[-1], ascending=False)[gcol] if g != "VPS35"][:7]
genes5 = ["VPS35"] + partners
m5 = pd.read_csv(T3 + "results/tables/S2_GSE35958_gene_matrix.csv", index_col=0)
have = [g for g in genes5 if g in m5.columns]
C = m5[have].corr(method="spearman").values

fig = plt.figure(figsize=(62/25.4, 54/25.4))
ax = fig.add_axes([0.24, 0.25, 0.54, 0.60])
cmap = plt.get_cmap("Reds"); norm = Normalize(0.3, 1.0)
im = ax.imshow(C, cmap=cmap, norm=norm, aspect="equal", interpolation="nearest")
ax.set_xticks(range(len(have))); ax.set_yticks(range(len(have)))
ax.set_xticklabels(have, rotation=45, ha="right", fontsize=5.8)
ax.set_yticklabels(have, fontsize=5.8)
ax.tick_params(length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
# VPS35 行列高亮
ax.add_patch(Rectangle((-0.5, -0.5), len(have), 1, fill=False, edgecolor="#1f3b73", lw=0.9))
ax.add_patch(Rectangle((-0.5, -0.5), 1, len(have), fill=False, edgecolor="#1f3b73", lw=0.9))
ax.set_title("VPS35 co-expression (Spearman r)", fontsize=6.8, pad=4)
cax = fig.add_axes([0.815, 0.30, 0.028, 0.40])
cb = fig.colorbar(im, cax=cax); cb.set_label("Spearman r", fontsize=5.4, labelpad=3)
cb.ax.tick_params(labelsize=5.2); cb.outline.set_linewidth(0.4)

save(fig, T3 + "figures/Fig6F_coexpression_matrix")
print("Fig4D / Fig5F rebuilt")

# ---------------- Fig 5D：VPS35 跨数据集（重绘，去掉杂散虚线） ----------------
import pandas as pd
dat = []
for acc in ["GSE35958", "GSE56814", "GSE56815"]:
    m = pd.read_csv(T3 + f"results/tables/S2_{acc}_gene_matrix.csv", index_col=0)
    v = m["VPS35"].astype(float) if "VPS35" in m.columns else m.loc["VPS35"].astype(float)
    dat.append(pd.DataFrame(dict(dataset=acc, z=(v - v.mean()) / v.std())))
dall = pd.concat(dat)
order = ["GSE35958", "GSE56814", "GSE56815"]
labels = ["GSE35958\n(MSC)", "GSE56814\n(PBMC)", "GSE56815\n(PBMC)"]
fig, ax = plt.subplots(figsize=(70/25.4, 56/25.4))
bp = ax.boxplot([dall[dall.dataset == d]["z"] for d in order], positions=[0, 1, 2], widths=0.5,
                patch_artist=True, showfliers=False, medianprops=dict(color="#333333", lw=0.8))
for patch in bp["boxes"]:
    patch.set_facecolor("#8491B4"); patch.set_alpha(0.85); patch.set_edgecolor("#333333"); patch.set_linewidth(0.5)
rng = np.random.default_rng(1)
for i, d in enumerate(order):
    vals = dall[dall.dataset == d]["z"].values
    ax.scatter(i + rng.uniform(-0.13, 0.13, len(vals)), vals, s=3.0, color="#333333", alpha=0.55, linewidths=0)
    ax.text(i, 2.55, f"n = {len(vals)}", ha="center", fontsize=5.6, color="#444444")
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(labels, fontsize=6.0)
ax.set_ylabel("VPS35 (z-score)", fontsize=6.6); ax.tick_params(axis="y", labelsize=5.8)
ax.set_ylim(-3.6, 2.9)
ax.set_title("VPS35 across datasets", fontsize=7.0)
fig.subplots_adjust(left=0.22, right=0.97, bottom=0.20, top=0.88)
save(fig, T3 + "figures/Fig6D_VPS35_datasets")
print("Fig5D rebuilt")
