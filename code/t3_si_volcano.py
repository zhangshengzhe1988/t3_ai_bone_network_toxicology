#!/usr/bin/env python3
"""统一配色的补充图火山图（S7 GSE56814 / S8 GSE56815）。

风格与主图一致：显著基因 #E64B35，非显著 #C8C8C8，阈值虚线，Arial，字号 ≥6.5 pt。
"""
import sys
sys.path.insert(0, "code")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nature_style import save, C_HIGH

T3 = "选题3_网络毒理学_AI骨质疏松/"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})

SIG, NSIG = C_HIGH, "#C8C8C8"

def volcano(csv, title, out, padj=0.05, lfc=1.0):
    d = pd.read_csv(T3 + "results/tables/" + csv)
    d = d.dropna(subset=["logFC", "P.Value"])
    d["P.Value"] = pd.to_numeric(d["P.Value"], errors="coerce")
    d["adj.P.Val"] = pd.to_numeric(d["adj.P.Val"], errors="coerce")
    d = d.dropna(subset=["P.Value"])
    y = -np.log10(d["P.Value"].clip(lower=1e-300))
    sig = d["adj.P.Val"] < padj   # 与正文一致：仅 BH 校正后 P < 0.05，不设 fold-change 截断
    fig, ax = plt.subplots(figsize=(100/25.4, 74/25.4))
    ax.scatter(d.loc[~sig, "logFC"], y[~sig], s=3.0, c=NSIG, linewidths=0, rasterized=False)
    if sig.any():
        ax.scatter(d.loc[sig, "logFC"], y[sig], s=3.4, c=SIG, linewidths=0)
    ax.axhline(-np.log10(0.05), color="#7E8AA2", lw=0.6, ls="--")
    ax.set_xlabel("log$_2$ fold change", fontsize=7.4)
    ax.set_ylabel("−log$_{10}$(P)", fontsize=7.4)
    ax.tick_params(labelsize=6.8)
    ax.set_title(title, fontsize=7.8)
    n = int(sig.sum())
    minadj = float(np.nanmin(d["adj.P.Val"])) if d["adj.P.Val"].notna().any() else float("nan")
    note = (f"{n} DEGs (BH adj. P < 0.05)" if n
            else f"0 DEGs; min BH adj. P = {minadj:.2f}")
    ax.text(0.02, 0.96, note, transform=ax.transAxes, fontsize=6.6,
            color=SIG if n else "#555555", va="top")
    ax.set_xlim(-abs(d["logFC"]).max()*1.10, abs(d["logFC"]).max()*1.10)
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.16, top=0.90)
    save(fig, T3 + "figures/" + out)

volcano("S2_GSE56814_limma_all.csv", "GSE56814 (PBMC, n = 73): no significant DEGs", "S7_GSE56814_volcano")
volcano("S2_GSE56815_limma_all.csv", "GSE56815 (PBMC, n = 80): 391 significant DEGs", "S8_GSE56815_volcano")
print("volcano panels done")

# ---------------- S6 核心基因表达点图（100×74 mm，重新排版避免标签被裁） ----------------
mat = pd.read_csv(T3 + "results/tables/S2_GSE35958_gene_matrix.csv", index_col=0)
grp = pd.read_csv(T3 + "results/tables/S2_GSE35958_groups.csv")
grp = grp[grp["group"].isin(["high", "low"])]
genes = ["CYP19A1", "HDAC1", "HDAC6", "ACHE", "ADRB1", "ADRB2", "ADRB3", "HTR2A", "HTR2C"]
samples = list(grp["sample"])
Z = mat.loc[samples, genes].T
Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1), axis=0)

fig, ax = plt.subplots(figsize=(100/25.4, 74/25.4))
for i, g in enumerate(genes):
    z = Z.loc[g].values
    sc = ax.scatter(np.arange(len(z)), [i]*len(z), s=np.abs(z)*16, c=z, cmap="RdBu_r",
                    vmin=-2, vmax=2, edgecolor="#555555", linewidth=0.25)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=7.0)
ax.set_xticks(range(Z.shape[1])); ax.set_xticklabels([f"S{j+1}" for j in range(Z.shape[1])], fontsize=6.4)
ax.set_xlabel("GSE35958 samples (high → low BMD)", fontsize=7.0)
ax.set_title("Core-gene expression (z-score)", fontsize=7.8)
ax.set_xlim(-0.6, Z.shape[1]-0.4); ax.set_ylim(-0.6, len(genes)-0.4)
ax.tick_params(length=0)
cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
cb.set_label("z-score", fontsize=7.0); cb.ax.tick_params(labelsize=6.4)
fig.subplots_adjust(left=0.20, right=0.94, bottom=0.17, top=0.90)
save(fig, T3 + "figures/S11_core_dotplot")
print("dotplot done")
