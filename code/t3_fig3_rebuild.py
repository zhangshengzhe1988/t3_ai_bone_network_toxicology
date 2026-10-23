#!/usr/bin/env python3
"""重建 Fig 3（转录组）六个面板，统一配色与排版。

判定口径沿用各队列已报告的数字（见面板标题），但统一视觉：显著 #E64B35、非显著 #C8C8C8。
"""
import sys
sys.path.insert(0, "code")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nature_style import save, C_HIGH

T3 = "选题3_网络毒理学_AI骨质疏松/"
SIG, NSIG = C_HIGH, "#C8C8C8"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6})

def mannwhitney_exact(x, y):
    """双侧精确 Mann-Whitney U（枚举所有分组，样本量小可行）。"""
    from itertools import combinations
    n1, n2 = len(x), len(y)
    allv = np.concatenate([x, y])
    ranks = np.empty(len(allv)); order = allv.argsort()
    ranks[order] = np.arange(1, len(allv)+1)
    obs = ranks[:n1].sum()
    # 平局取平均秩
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    if (counts > 1).any():
        for k in np.where(counts > 1)[0]:
            idx = np.where(inv == k)[0]
            ranks[idx] = ranks[idx].mean()
        obs = ranks[:n1].sum()
    n = len(allv)
    def U(idx):
        s = ranks[list(idx)].sum()
        return s - n1*(n1+1)/2
    u_obs = U(range(n1))
    dist = [U(c) for c in combinations(range(n), n1)]
    dist = np.array(dist)
    mu = n1*n2/2
    if u_obs > mu:
        p = np.mean(dist >= u_obs - 1e-9)
    else:
        p = np.mean(dist <= u_obs + 1e-9)
    return min(1.0, 2*p)

def load(f):
    d = pd.read_csv(T3 + "results/tables/" + f)
    for c in ("logFC", "P.Value", "adj.P.Val"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["logFC", "P.Value", "adj.P.Val"])

def volcano(d, title, out, fc_cut, note):
    y = -np.log10(d["P.Value"].clip(lower=1e-300))
    sig = (d["adj.P.Val"] < 0.05) & (d["logFC"].abs() > fc_cut)
    fig, ax = plt.subplots(figsize=(84/25.4, 54/25.4))
    ax.scatter(d.loc[~sig, "logFC"], y[~sig], s=2.4, c=NSIG, linewidths=0)
    ax.scatter(d.loc[sig, "logFC"], y[sig], s=2.6, c=SIG, linewidths=0)
    ax.axhline(-np.log10(0.05), color="#7E8AA2", lw=0.5, ls="--")
    for x in (-fc_cut, fc_cut):
        ax.axvline(x, color="#7E8AA2", lw=0.5, ls="--")
    ax.set_xlabel("log$_2$ fold change", fontsize=6.6)
    ax.set_ylabel("−log$_{10}$(P)", fontsize=6.6)
    ax.tick_params(labelsize=5.8)
    ax.set_title(title, fontsize=7.2)
    ax.text(0.02, 0.96, note, transform=ax.transAxes, fontsize=5.8, va="top", color="#444444")
    lim = max(abs(d["logFC"]).max()*1.05, fc_cut*1.5)
    ax.set_xlim(-lim, lim); ax.set_ylim(0, y.max()*1.06)
    fig.subplots_adjust(left=0.20, right=0.97, bottom=0.19, top=0.88)
    save(fig, T3 + "figures/" + out)

a = load("S2_GSE35958_limma_all.csv")
b = load("S2b_GSE13850_limma.csv")

# A GSE35958
volcano(a, "GSE35958 (MSC, n = 9): 1,447 DEGs", "Fig2A_GSE35958_volcano", 1.0,
        "adj. P < 0.05, |log$_2$FC| > 1")

# B core genes on the same volcano
core = ["CYP19A1", "HDAC1", "HDAC6", "ACHE", "ADRB1", "ADRB2", "ADRB3", "HTR2A", "HTR2C"]
d_core = a[a["gene"].isin(core)].set_index("gene")
y = -np.log10(a["P.Value"].clip(lower=1e-300))
sig = (a["adj.P.Val"] < 0.05) & (a["logFC"].abs() > 1.0)
fig, ax = plt.subplots(figsize=(84/25.4, 54/25.4))
ax.scatter(a.loc[~sig, "logFC"], y[~sig], s=2.2, c=NSIG, linewidths=0)
ax.scatter(a.loc[sig, "logFC"], y[sig], s=2.4, c="#E8AFA8", linewidths=0)
ax.scatter(d_core["logFC"], -np.log10(d_core["P.Value"]), s=16, facecolor="none",
           edgecolor="#1f3b73", linewidths=0.7, zorder=4)
offsets = {"CYP19A1": (6, 6), "HDAC1": (-30, 4), "HDAC6": (6, -8), "ACHE": (6, 4),
           "ADRB1": (-34, -2), "ADRB2": (6, -10), "ADRB3": (-34, -10), "HTR2A": (6, 6), "HTR2C": (-36, 6)}
for g, r in d_core.iterrows():
    dx, dy = offsets.get(g, (5, 4))
    ax.annotate(g, (r["logFC"], -np.log10(max(r["P.Value"], 1e-300))),
                textcoords="offset points", xytext=(dx, dy), fontsize=5.2, color="#1f3b73")
ax.axhline(-np.log10(0.05), color="#7E8AA2", lw=0.5, ls="--")
for x in (-1.0, 1.0):
    ax.axvline(x, color="#7E8AA2", lw=0.5, ls="--")
ax.set_xlabel("log$_2$ fold change", fontsize=6.6); ax.set_ylabel("−log$_{10}$(P)", fontsize=6.6)
ax.tick_params(labelsize=5.8)
ax.set_title("Core genes (open circles)", fontsize=7.2)
ax.set_xlim(-5.2, 5.2); ax.set_ylim(0, y.max()*1.06)
fig.subplots_adjust(left=0.20, right=0.97, bottom=0.19, top=0.88)
save(fig, T3 + "figures/Fig2B_volcano_core_genes")

# C 表达箱线图（线性强度；标注 Mann-Whitney P）
mat = pd.read_csv(T3 + "results/tables/S2_GSE35958_gene_matrix.csv", index_col=0)
grp = pd.read_csv(T3 + "results/tables/S2_GSE35958_groups.csv")
grp = grp[grp["group"].isin(["high", "low"])]
fig, axes = plt.subplots(1, 2, figsize=(84/25.4, 54/25.4))
for ax, g in zip(axes, ["CYP19A1", "VPS35"]):
    hi = mat.loc[grp[grp.group == "high"]["sample"], g].astype(float)
    lo = mat.loc[grp[grp.group == "low"]["sample"], g].astype(float)
    bp = ax.boxplot([hi, lo], widths=0.55, patch_artist=True, showfliers=False,
                    medianprops=dict(color="#333333", lw=0.8))
    for patch, col in zip(bp["boxes"], ["#4DBBD5", C_HIGH]):
        patch.set_facecolor(col); patch.set_alpha(0.85); patch.set_edgecolor("#333333"); patch.set_linewidth(0.5)
    rng = np.random.default_rng(42)
    for i, vals in enumerate([hi, lo], start=1):
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(vals)), vals, s=7, color="#333333", zorder=3, linewidths=0)
    from statistics import mean
    ax.plot([0.8, 1.2], [hi.mean()]*2, color="#333333", lw=0.8)
    ax.plot([1.8, 2.2], [lo.mean()]*2, color="#333333", lw=0.8)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["High\nBMD", "Low\nBMD"], fontsize=5.8)
    ax.set_title(g, fontsize=7.2); ax.tick_params(axis="y", labelsize=5.8)
    if g == "CYP19A1":
        ax.set_ylabel("Normalised expression", fontsize=6.4)
    adjp = {"CYP19A1": 0.0868, "VPS35": 0.3508}[g]   # limma 校正 P（与正文/Table S10 一致）
    ax.text(0.5, 0.90, f"adj. P = {adjp:.3f}",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.4, color="#333333")
fig.subplots_adjust(left=0.14, right=0.97, bottom=0.20, top=0.88, wspace=0.42)
save(fig, T3 + "figures/Fig2C_core_expression")

# D GSE13850
volcano(b, "GSE13850 (B cell, n = 20): 28 DEGs", "Fig2D_GSE13850_volcano", np.log2(1.5),
        "adj. P < 0.05, |log$_2$FC| > log$_2$1.5")

# E 核心基因 log2FC 热图
fc = pd.DataFrame({"GSE13850": [b.set_index("gene")["logFC"].get(g, np.nan) for g in core],
                   "GSE35958": [a.set_index("gene")["logFC"].get(g, np.nan) for g in core]}, index=core)
fig, ax = plt.subplots(figsize=(84/25.4, 54/25.4))
im = ax.imshow(fc.values, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_xticks([0, 1]); ax.set_xticklabels(["GSE13850", "GSE35958"], fontsize=6.0)
ax.set_yticks(range(len(core))); ax.set_yticklabels(core, fontsize=5.8)
for i in range(len(core)):
    for j in range(2):
        v = fc.values[i, j]
        ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=5.2,
                color="white" if abs(v) > 1.6 else "#222222")
ax.set_title("Core-gene log$_2$FC across cohorts", fontsize=7.0)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03); cb.set_label("log$_2$FC", fontsize=6.0)
cb.ax.tick_params(labelsize=5.6)
fig.subplots_adjust(left=0.22, right=0.90, bottom=0.16, top=0.90)
save(fig, T3 + "figures/Fig2E_core_logFC")

# F DEG 计数（横向，避免标签重叠）
labels = ["GSE35958\n(MSC, n=9)", "GSE13850\n(B cell, n=20)", "GSE56815\n(PBMC, n=80)", "GSE56814\n(PBMC, n=73)"]
vals = [1447, 28, 391, 0]
fig, ax = plt.subplots(figsize=(84/25.4, 54/25.4))
bars = ax.barh(range(3, -1, -1), vals, color=[C_HIGH, "#F39B7F", "#9DB4C9", "#C8C8C8"], height=0.62)
for y_, v in zip(range(3, -1, -1), vals):
    ax.text(v + max(vals)*0.02, y_, str(v), va="center", fontsize=6.0)
ax.set_yticks(range(3, -1, -1)); ax.set_yticklabels(labels, fontsize=5.6)
ax.set_xlabel("DEGs (n)", fontsize=6.6); ax.tick_params(axis="x", labelsize=5.8)
ax.set_xlim(0, max(vals)*1.18)
ax.set_title("DEGs across four cohorts", fontsize=7.0)
fig.subplots_adjust(left=0.34, right=0.97, bottom=0.17, top=0.88)
save(fig, T3 + "figures/Fig2F_cohort_degs")
print("Fig3 panels rebuilt")
