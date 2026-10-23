#!/usr/bin/env python3
"""Nature 风格共享绘图模块：Arial、语义色板、PDF 矢量主输出 + PNG 300dpi 预览
用户规范：字体统一 Arial、Nature 系列杂志风格、PDF 矢量。"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 8
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['axes.titlesize'] = 9
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

C_HIGH = "#E64B35"   # 高风险组
C_LOW  = "#4DBBD5"   # 低风险组
C_GREY = "#808080"

def save(fig, path_stem, dpi=None):
    # 按组图槽位设定物理尺寸：使面板以 1:1 放入槽位，避免整体缩小后字号 <5 pt
    try:
        import os, sys
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from panel_slots import SLOT_MM
        key = os.path.basename(path_stem)
        if key in SLOT_MM:
            w, h = SLOT_MM[key]
            fig.set_size_inches(w/25.4, h/25.4)
            _SLOT = True
    except Exception:
        pass
    # 字号下限：所有文本至少 5.5 pt（组图后仍 ≥5 pt）
    try:
        import matplotlib.text as mtext
        for t in fig.findobj(mtext.Text):
            if t.get_fontsize() < 6.5:
                t.set_fontsize(6.5)
    except Exception:
        pass
    _bbox = None if "_SLOT" in dir() else "tight"
    fig.savefig(f"{path_stem}.pdf", format="pdf", bbox_inches=_bbox, dpi=dpi)
    fig.savefig(f"{path_stem}.png", dpi=300, bbox_inches=_bbox)
    plt.close(fig)
    print("saved:", path_stem + ".pdf / .png")

def km_curve(T, E):
    """Kaplan-Meier 步进曲线（numpy 实现，处理删失）"""
    order = np.argsort(T)
    t = np.asarray(T)[order]; e = np.asarray(E)[order]
    n = len(t)
    surv = 1.0
    ts, ss = [0.0], [1.0]
    i = 0
    while i < n:
        j = i
        while j < n and t[j] == t[i]:
            j += 1
        d = int(np.sum(e[i:j]))
        r = n - i
        if r > 0:
            surv *= (1 - d / r)
        ts.append(t[i]); ss.append(surv)
        i = j
    return np.array(ts), np.array(ss)
