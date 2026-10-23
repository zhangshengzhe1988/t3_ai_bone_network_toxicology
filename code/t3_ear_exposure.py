#!/usr/bin/env python3
"""鱼血浆模型暴露-活性比（EAR）情景分析。

方法：read-across / fish plasma model（Rand-Weaver 2013; Margiotta-Casaluci 2014）
  log K(plasma) = 0.73 * logKow - 0.16
  C_fish,plasma  = C_water * K(plasma)          [ng/mL = ng/L * L/kg / 1000]
  EAR            = C_fish,plasma / C_human,plasma
输入均为可引用来源，不做外推假设：
  logKow  : PubChem XLogP（letrozole 2.7; anastrozole 2.1; exemestane 3.1）
  C_human : letrozole 129 nmol/L (36.8 ng/mL, Cmax, SmPC); anastrozole 25.7 ng/mL
            (steady-state trough, label); exemestane 18 ng/mL (Cmax, label)
  C_water : 1–100 ng/L 情景（覆盖已报道 ng/L–µg/L 与 medaka 20–200 ng/L 效应区间）
输出：results/tables/S17_Exposure_EAR.csv, figures/S12_exposure_EAR.pdf
"""
import sys
sys.path.insert(0, "code")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nature_style import save, C_HIGH, C_LOW

OUT = "选题3_网络毒理学_AI骨质疏松/"
drugs = [("Letrozole", 2.7, 36.8, "#3C5488"),
         ("Anastrozole", 2.1, 25.7, "#E64B35"),
         ("Exemestane", 3.1, 18.0, "#00A087")]
waters = [1.0, 10.0, 20.0, 100.0]

rows = []
for name, logk, chum, _ in drugs:
    kpl = 10 ** (0.73 * logk - 0.16)
    for cw in waters:
        cpl = cw * kpl / 1000.0
        rows.append(dict(drug=name, logKow=logk, K_plasma=round(kpl, 1),
                         Cwater_ngL=cw, Cfish_plasma_ngmL=round(cpl, 4),
                         Chuman_ngmL=chum, EAR=round(cpl / chum, 4)))
df = pd.DataFrame(rows)
df.to_csv(OUT + "results/tables/S17_Exposure_EAR.csv", index=False)

fig, ax = plt.subplots(figsize=(90/25.4, 62/25.4))
for name, logk, chum, col in drugs:
    d = df[df.drug == name]
    ax.plot(d.Cwater_ngL, d.EAR, "o-", color=col, lw=1.0, ms=3.2, label=name)
ax.axhline(1.0, color="k", ls="--", lw=0.6)
ax.axvspan(20, 200, color="#B0B0B0", alpha=0.18, lw=0)
ax.text(63, 0.012, "medaka effect range\n20–200 ng/L", fontsize=5.6, ha="center", color="#444444")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Water concentration (ng/L)", fontsize=6.6)
ax.set_ylabel("Exposure–activity ratio (EAR)", fontsize=6.6)
ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"], fontsize=6.4)
ax.set_yticks([0.001, 0.01, 0.1, 1]); ax.set_yticklabels(["0.001", "0.01", "0.1", "1"], fontsize=6.4)
ax.minorticks_off()
ax.set_title("Fish plasma-model EAR", fontsize=7.0)
ax.legend(fontsize=5.8, frameon=False, loc="lower right")
ax.set_ylim(1e-3, 3)
fig.subplots_adjust(left=0.16, right=0.97, bottom=0.17, top=0.90)
save(fig, OUT + "figures/S12_exposure_EAR")
print(df.to_string(index=False))
print("saved S17 csv + S12 pdf")
