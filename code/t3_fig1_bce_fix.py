#!/usr/bin/env python3
"""重绘 Fig 1B/1C/1E：修掉轴标签被裁、刻度重叠；统一 Arial 与配色。"""
import sys
sys.path.insert(0, "code")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nature_style import save, C_HIGH, C_LOW

T3 = "选题3_网络毒理学_AI骨质疏松/"
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})

# ---------------- B：ROR 森林图 ----------------
d = pd.read_csv(T3 + "results/tables/S1_FAERS_ROR_20261013.csv", skiprows=2, nrows=5)
labels = [f"{r.drug}\n{r.event.replace('osteoporotic fracture', 'fracture')}" for r in d.itertuples()]
fig, ax = plt.subplots(figsize=(70/25.4, 54/25.4))
y = np.arange(len(d))[::-1]
ax.hlines(y, d["ROR_lo"], d["ROR_hi"], color="#7E8AA2", lw=1.0)
ax.scatter(d["ROR"], y, s=18, color=C_HIGH, zorder=3)
ax.axvline(1, color="#333333", lw=0.6, ls="--")
for yi, r in zip(y, d.itertuples()):
    ax.text(r.ROR_hi + 0.25, yi, f"{r.ROR:.2f} [{r.ROR_lo:.2f}–{r.ROR_hi:.2f}]",
            va="center", fontsize=5.8)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.8)
ax.set_xlabel("ROR (95% CI), FAERS 2013–2026", fontsize=6.6)
ax.set_xlim(0.5, 12.4); ax.tick_params(axis="x", labelsize=5.8)
ax.set_title("Disproportionality signals", fontsize=7.2)
fig.subplots_adjust(left=0.30, right=0.97, bottom=0.17, top=0.88)
save(fig, T3 + "figures/Fig1B_FAERS_ROR_forest")

# ---------------- C：年度趋势 ----------------
trend = pd.read_csv(T3 + "results/tables/S1_FAERS_ROR_20261013.csv", skiprows=9)
trend = trend.dropna(subset=["year"])
fig, ax = plt.subplots(figsize=(70/25.4, 54/25.4))
ax.plot(trend["year"], trend["letrozole_x_osteoporosis_reports"], "o-",
        color=C_HIGH, lw=1.0, ms=2.8)
ax.set_xlabel("Year", fontsize=6.8)
ax.set_ylabel("Reports (letrozole × osteoporosis)", fontsize=6.4)
ax.set_xticks([2013, 2016, 2019, 2022, 2025]); ax.tick_params(labelsize=5.8)
ax.set_xlim(2012.5, 2026.5)
ax.set_title("Annual reports", fontsize=7.2)
fig.subplots_adjust(left=0.26, right=0.97, bottom=0.17, top=0.88)
save(fig, T3 + "figures/Fig1C_FAERS_trend")

# ---------------- E：年龄/性别分层 ----------------
fig, axes = plt.subplots(1, 2, figsize=(70/25.4, 54/25.4))
age_v = [130, 108, 15]; age_l = ["<65", "65–80", ">80"]
axes[0].bar(range(3), age_v, color=[C_LOW, C_HIGH, "#9DB4C9"], width=0.6)
for i, v in enumerate(age_v):
    axes[0].text(i, v + 3, str(v), ha="center", fontsize=5.8)
axes[0].set_xticks(range(3)); axes[0].set_xticklabels(age_l, fontsize=6.0)
axes[0].set_title("Age (n = 253)", fontsize=6.8)
axes[0].tick_params(axis="y", labelsize=5.8); axes[0].set_ylim(0, 155)
sex_v = [4, 390]; sex_l = ["Male", "Female"]
axes[1].bar(range(2), sex_v, color=["#9DB4C9", C_HIGH], width=0.6)
for i, v in enumerate(sex_v):
    axes[1].text(i, v + 8, str(v), ha="center", fontsize=5.8)
axes[1].set_xticks(range(2)); axes[1].set_xticklabels(sex_l, fontsize=6.0)
axes[1].set_title("Sex (n = 394)", fontsize=6.8)
axes[1].tick_params(axis="y", labelsize=5.8); axes[1].set_ylim(0, 440)
fig.text(0.035, 0.5, "Reports", rotation=90, va="center", ha="center", fontsize=6.4)
fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.86, wspace=0.45)
save(fig, T3 + "figures/Fig1E_stratification")
print("Fig1 B/C/E rebuilt")
