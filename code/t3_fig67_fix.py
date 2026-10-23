#!/usr/bin/env python3
"""修复 Fig 6B/6D/6E/6F 与 Fig 7A/7D 的标签重叠、数值缺失与轴截断。"""
import sys
sys.path.insert(0, "code")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from nature_style import save, C_HIGH, C_LOW

T3 = "选题3_网络毒理学_AI骨质疏松/"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})

# ---------- 6B：vKO top 15（排除 VPS35 自身） ----------
d = pd.read_csv(T3 + "results/tables/S10_vko_VPS35_diffRegulation.csv")
d = d[d["gene"] != "VPS35"].copy()
top = d.reindex(d["Z"].abs().sort_values(ascending=False).index).head(15).iloc[::-1]
fig, ax = plt.subplots(figsize=(78/25.4, 54/25.4))
cols = [C_HIGH if z > 0 else C_LOW for z in top["Z"]]
ax.hlines(range(len(top)), 0, top["Z"], color=cols, lw=1.0)
ax.scatter(top["Z"], range(len(top)), c=cols, s=12, zorder=3)
ax.axvline(0, color="#333333", lw=0.5)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top["gene"], fontsize=5.8)
ax.set_xlabel("Perturbation Z score", fontsize=6.4); ax.tick_params(axis="x", labelsize=5.8)
ax.set_title("VPS35 knockout, top 15 genes", fontsize=6.8)
fig.subplots_adjust(left=0.24, right=0.97, bottom=0.17, top=0.88)
save(fig, T3 + "figures/Fig7C_VPS35_lollipop")

# ---------- 6D：vKO GO ----------
go = pd.read_csv(T3 + "results/tables/S10_vko_GO.csv").sort_values("p.adjust").head(10)
fig, ax = plt.subplots(figsize=(72/25.4, 58/25.4))
ax.scatter(-np.log10(go["p.adjust"]), range(len(go))[::-1], s=go["Count"]*10,
           c=go["Count"], cmap="YlOrRd", linewidths=0.3, edgecolors="grey")
ax.set_yticks(range(len(go))[::-1])
ax.set_yticklabels([t if len(t) <= 30 else t[:28] + "…" for t in go["Description"]], fontsize=5.6)
ax.set_xlabel("−log10(adj. P)", fontsize=6.4); ax.tick_params(axis="x", labelsize=5.8)
ax.set_title("GO enrichment of VPS35 vKO", fontsize=6.8)
fig.subplots_adjust(left=0.56, right=0.97, bottom=0.15, top=0.88)
save(fig, T3 + "figures/Fig7C3_VPS35_vko_GO")

# ---------- 6E：Vina vs Vinardo（横向分组条形） ----------
scores = pd.read_csv(T3 + "results/tables/S14_vinardo_rescore_best.csv")
order = [("human_CYP19A1", "letrozole", "Letrozole\nhuman"), ("human_CYP19A1", "anastrozole", "Anastrozole\nhuman"),
         ("zebrafish_Cyp19a1a", "letrozole", "Letrozole\nfish"), ("zebrafish_Cyp19a1a", "anastrozole", "Anastrozole\nfish"),
         ("human_VPS35", "letrozole", "Letrozole\nVPS35"), ("human_VPS35", "anastrozole", "Anastrozole\nVPS35")]
vina = [scores[(scores.target == t) & (scores.ligand == l)]["vina"].iloc[0] for t, l, _ in order]
vind = [scores[(scores.target == t) & (scores.ligand == l)]["vinardo_total"].iloc[0] for t, l, _ in order]
lab = [o[2] for o in order]
fig, ax = plt.subplots(figsize=(76/25.4, 60/25.4))
y = np.arange(len(lab))[::-1]; h = 0.36
b1 = ax.barh(y + h/2, vina, h, color=C_HIGH, label="Vina")
b2 = ax.barh(y - h/2, vind, h, color="#9DB4C9", edgecolor="#5c6b7a", linewidth=0.4, hatch="///", label="Vinardo")
for bars, vals, col in [(b1, vina, "white"), (b2, vind, "#222222")]:
    for b, v in zip(bars, vals):
        ax.text(v - 0.12, b.get_y()+b.get_height()/2, f"{v:.1f}", va="center", ha="right",
                fontsize=5.2, color=col)
ax.axvline(-7, color="grey", lw=0.6, ls="--")
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=5.2)
ax.set_xlim(-9.6, 0); ax.set_xlabel("Docking score (kcal/mol)", fontsize=6.0)
ax.tick_params(axis="x", labelsize=5.6)
ax.legend(fontsize=5.4, frameon=False, loc="lower left", handlelength=1.2)
ax.set_title("Vina vs Vinardo", fontsize=6.6)
fig.subplots_adjust(left=0.32, right=0.97, bottom=0.20, top=0.90)
save(fig, T3 + "figures/Fig7E_docking_energies2")

# ---------- 6F：VPS35 contacts ----------
c = pd.read_csv(T3 + "results/tables/S9b_VPS35_contacts.csv")
fig, ax = plt.subplots(figsize=(64/25.4, 58/25.4))
y = np.arange(len(c))
colors = [C_HIGH if not np.isnan(v) else "#9DB4C9" for v in c["letrozole_A"]]
ax.barh(y, c["dmin"], color=colors, height=0.62)
for i, r in c.iterrows():
    who = "L" if not np.isnan(r["letrozole_A"]) else "A"
    ax.text(r["dmin"] + 0.03, i, f"{r['dmin']:.2f} ({who})", va="center", fontsize=5.2)
ax.axvline(4.0, color="#333333", lw=0.6, ls="--")
ax.set_yticks(y); ax.set_yticklabels(c["residue"], fontsize=5.6)
ax.set_ylim(len(c)-0.4, -1.6); ax.set_xlim(2.8, 4.75)
ax.set_xlabel("Min. ligand–residue distance (Å)", fontsize=6.0)
ax.tick_params(axis="x", labelsize=5.6)
ax.text(4.04, -0.9, "4 Å cutoff", fontsize=5.0, color="#333333")
ax.set_title("VPS35 contact residues", fontsize=6.8)
fig.subplots_adjust(left=0.26, right=0.97, bottom=0.16, top=0.90)
save(fig, T3 + "figures/Fig7G_VPS35_contacts")

# ---------- 7A：跨物种对接能 ----------
labels = ["Let.\nhuman", "Let.\nzebrafish", "Ana.\nhuman", "Ana.\nzebrafish"]
vals = [-7.22, -7.632, -7.12, -7.201]
fig, ax = plt.subplots(figsize=(70/25.4, 54/25.4))
bars = ax.bar(range(4), vals, color=[C_HIGH, "#F39B7F", "#3C5488", "#9DB4C9"], width=0.6)
for b, v in zip(bars, vals):
    ax.text(b.get_x()+b.get_width()/2, v - 0.08, f"{v:.2f}", ha="center", va="top", fontsize=5.8)
ax.axhline(-7.0, color="grey", lw=0.6, ls="--")
ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=5.6)
ax.set_ylim(-8.3, -6.4); ax.set_ylabel("Docking score (kcal/mol)", fontsize=6.2)
ax.tick_params(axis="y", labelsize=5.6)
ax.set_title("CYP19 orthologues", fontsize=6.8)
fig.subplots_adjust(left=0.22, right=0.97, bottom=0.22, top=0.88)
save(fig, T3 + "figures/Fig8A_cross_species_energy")

# ---------- 7D：EAR ----------
ear = pd.read_csv(T3 + "results/tables/S17_Exposure_EAR.csv")
fig, ax = plt.subplots(figsize=(74/25.4, 56/25.4))
for drug, col in [("Letrozole", "#3C5488"), ("Anastrozole", C_HIGH), ("Exemestane", "#00A087")]:
    sub = ear[ear["drug"] == drug]
    ax.plot(sub["Cwater_ngL"], sub["EAR"], "o-", color=col, lw=0.9, ms=2.6, label=drug)
ax.axhline(1.0, color="k", ls="--", lw=0.6)
ax.axvspan(20, 200, color="#B0B0B0", alpha=0.18, lw=0)
ax.text(63, 0.55, "medaka effect range 20–200 ng/L", fontsize=4.8, ha="center", color="#444444")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"], fontsize=5.6)
ax.set_yticks([0.001, 0.01, 0.1, 1]); ax.set_yticklabels(["0.001", "0.01", "0.1", "1"], fontsize=5.6)
ax.minorticks_off()
ax.set_xlabel("Water concentration (ng/L)", fontsize=6.2)
ax.set_ylabel("Exposure–activity ratio (EAR)", fontsize=6.2)
ax.legend(fontsize=4.9, frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.98), handlelength=1.1)
ax.set_title("Fish plasma-model EAR", fontsize=6.8)
ax.set_ylim(0.001, 3)
fig.subplots_adjust(left=0.20, right=0.97, bottom=0.20, top=0.88)
save(fig, T3 + "figures/S12_exposure_EAR")
print("Fig6/7 panels fixed")
