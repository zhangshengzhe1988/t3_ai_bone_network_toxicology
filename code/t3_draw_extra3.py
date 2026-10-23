import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import matplotlib.text as _mtext
def FONT_FLOOR(fig, floor=6.5):
    for _t in fig.findobj(_mtext.Text):
        if _t.get_fontsize() < floor: _t.set_fontsize(floor)
    return fig

"""T3 新增：Fig7G VPS35 接触残基、Fig3G STRING 边置信度"""
import sys; sys.path.insert(0, "code")
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nature_style import save
T3 = "选题3_网络毒理学_AI骨质疏松/"
C_HI=(0.902,0.294,0.208); C_LO=(0.302,0.686,0.835)
plt.rcParams.update({"font.family":"Arial","font.size":6,"pdf.fonttype":42,"axes.linewidth":0.5,
                     "axes.spines.top":False,"axes.spines.right":False})

# Fig7G 接触残基
c = pd.read_csv(T3+"results/tables/S9b_VPS35_contacts.csv")
fig, ax = plt.subplots(figsize=(62/25.4, 58/25.4))
y = np.arange(len(c))
ax.barh(y, c["dmin"], color=[C_HI if not np.isnan(v) else "#B0B0B0" for v in c["letrozole_A"]], height=0.62)
ax.axvline(4.0, color='k', lw=0.5, ls='--')
ax.text(3.05, len(c)-0.4, "dashed: 4 Å cutoff", fontsize=6.0, color="#444444")
for i, (_, r) in enumerate(c.iterrows()):
    who = "L" if not np.isnan(r["letrozole_A"]) else "A"
    ax.text(r["dmin"]+0.03, i, f"{r['dmin']:.2f} ({who})", va='center', fontsize=6.0)
ax.set_yticks(y); ax.set_yticklabels(c["residue"], fontsize=6.0)
ax.set_xlabel("Min. distance (Å)", fontsize=6.5)
ax.set_xlim(2.8, 4.6); ax.invert_yaxis()
ax.set_title("VPS35 contact residues", fontsize=6.8)
fig.subplots_adjust(left=0.40, right=0.97, bottom=0.22, top=0.90)
save(fig, T3+"figures/Fig7G_VPS35_contacts")

# Fig3G STRING 边置信度
e = pd.read_csv(T3+"results/tables/S5_string_network.tsv", sep="\t")[["preferredName_A","preferredName_B","score"]]
e = e.sort_values("score")
fig, ax = plt.subplots(figsize=(100/25.4, 48/25.4))
cols = [C_HI if s >= 0.9 else ("#F0A07C" if s >= 0.7 else "#B0B0B0") for s in e["score"]]
ax.barh(range(len(e)), e["score"], color=cols, height=0.6)
for i, (_, r) in enumerate(e.iterrows()):
    ax.text(r["score"]+0.01, i, f"{r['score']:.2f}", va='center', fontsize=6.0)
ax.set_yticks(range(len(e)))
ax.set_yticklabels([f"{a}–{b}" for a, b in zip(e["preferredName_A"], e["preferredName_B"])], fontsize=6.0)
ax.set_xlabel("STRING combined score", fontsize=6.6); ax.set_xlim(0, 1.12)
ax.set_title("PPI edge confidence", fontsize=7.0)
ax.tick_params(axis="x", labelsize=6.0)
fig.subplots_adjust(left=0.34, right=0.98, bottom=0.20, top=0.90)
save(fig, T3+"figures/Fig3G_string_edges")
print("Fig7G + Fig3G 完成")
