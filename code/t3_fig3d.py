"""3D 插画风格矢量绘图工具（matplotlib，纯矢量：分层椭圆/多边形模拟光照）"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon, FancyBboxPatch, Circle, Rectangle, FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches

def rgb(h):
    if isinstance(h, (tuple, list)):
        return tuple(float(v) for v in h[:3])
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))

def mix(c1, c2, t):
    a, b = rgb(c1), rgb(c2)
    return tuple(a[i] + (b[i]-a[i])*t for i in range(3))

def shade(c, amt):
    """amt>0 变亮, <0 变暗"""
    return mix(c, '#ffffff', amt) if amt >= 0 else mix(c, '#000000', -amt)

def sphere(ax, x, y, r, color, z=5, n=14, light=(-0.35, -0.38)):
    """3D 球：分层椭圆模拟径向光照（纯矢量）"""
    for k in range(n, 0, -1):
        t = k/n
        rr = r*t
        col = mix(shade(color, 0.62), color, (1-t)*0.9)
        col = mix(col, shade(color, -0.45), max(0, t-0.55)/0.45*0.85)
        ax.add_patch(Ellipse((x + light[0]*r*(1-t)*0.9, y + light[1]*r*(1-t)*0.9),
                             rr*2, rr*2, facecolor=col, edgecolor='none', zorder=z))
    ax.add_patch(Ellipse((x + light[0]*r*0.55, y + light[1]*r*0.55), r*0.42, r*0.30,
                         facecolor='white', alpha=0.55, edgecolor='none', zorder=z+1))
    ax.add_patch(Ellipse((x, y), r*2, r*2, facecolor='none',
                         edgecolor=shade(color, -0.5), linewidth=0.5, zorder=z+2))

def cylinder(ax, x, y, w, h, color, z=5, n=12):
    """3D 数据库/圆柱体（俯视椭圆 + 带渐变的柱身）"""
    ry = w*0.22
    for k in range(n):
        t0, t1 = k/n, (k+1)/n
        col = mix(shade(color, 0.35), shade(color, -0.30), t0)
        ax.add_patch(Polygon([[x-w/2, y-h], [x+w/2, y-h],
                              [x+w/2, y], [x-w/2, y]], closed=True,
                             facecolor='none', edgecolor='none', zorder=z))
        ax.add_patch(Rectangle((x-w/2, y-h*t1), w, h*(t1-t0) + 0.001,
                               facecolor=col, edgecolor='none', zorder=z))
    for k in range(6, 0, -1):
        ax.add_patch(Rectangle((x-w/2 + (w-k*w/12), y-h), w - 2*(w-k*w/12), h,
                               facecolor='none', edgecolor='none', zorder=z))
    ax.add_patch(Ellipse((x, y), w, ry*2, facecolor=shade(color, 0.30),
                         edgecolor=shade(color, -0.5), linewidth=0.5, zorder=z+2))
    ax.add_patch(Ellipse((x-w*0.15, y+ry*0.25), w*0.55, ry*0.7, facecolor='white',
                         alpha=0.35, edgecolor='none', zorder=z+3))

def cell3d(ax, x, y, r, fill='#EAF3FB', rim='#4A7FA5', z=6, nucleus=True, nuc_r=0.38, seed=3):
    rng = np.random.default_rng(seed)
    ang = np.linspace(0, 2*np.pi, 200)
    rr = r*(1 + 0.045*np.sin(3*ang + seed) + 0.03*np.sin(7*ang + seed*2))
    verts = np.c_[x + rr*np.cos(ang), y + rr*np.sin(ang)]
    ax.add_patch(Polygon(verts, closed=True, facecolor=shade(fill, 0.30),
                         edgecolor=shade(fill, -0.25), linewidth=0.5, zorder=z))
    for k in range(8):
        t = (k+1)/8
        ax.add_patch(Ellipse((x - r*0.12, y - r*0.12), 2*r*(1-0.6*t), 2*r*(1-0.6*t),
                             facecolor=mix(shade(fill, 0.30), shade(fill, -0.12), t),
                             edgecolor='none', zorder=z))
    for _ in range(10):
        a = rng.uniform(0, 2*np.pi); d = rng.uniform(0.35, 0.85)*r
        ax.add_patch(Circle((x + d*np.cos(a), y + d*np.sin(a)), r*0.035,
                            facecolor=shade(fill, -0.15), edgecolor='none', alpha=0.7, zorder=z+1))
    if nucleus:
        sphere(ax, x + r*0.15, y + r*0.10, r*nuc_r, '#8FB3D9', z=z+2)
    ax.add_patch(Ellipse((x, y), 2*r, 2*r, facecolor='none',
                         edgecolor=rim, linewidth=0.55, zorder=z+4))

def dna_helix(ax, x0, y0, length, n=40, amp=1.6, period=6, c1='#0072B2', c2='#9CCAF0', z=6, vertical=False):
    t = np.linspace(0, length, n*4)
    ph = 2*np.pi*t/period
    if vertical:
        ax.plot(x0 + amp*np.sin(ph), y0 + t, color=c1, lw=0.9, zorder=z)
        ax.plot(x0 - amp*np.sin(ph), y0 + t, color=c2, lw=0.9, zorder=z)
        idx = np.arange(0, len(t), max(1, len(t)//16))
        for i in idx:
            ax.plot([x0 + amp*np.sin(ph[i]), x0 - amp*np.sin(ph[i])],
                    [y0 + t[i], y0 + t[i]], color=shade(c1, -0.25), lw=0.45, zorder=z)
    else:
        ax.plot(x0 + t, y0 + amp*np.sin(ph), color=c1, lw=0.9, zorder=z)
        ax.plot(x0 + t, y0 - amp*np.sin(ph), color=c2, lw=0.9, zorder=z)
        idx = np.arange(0, len(t), max(1, len(t)//16))
        for i in idx:
            ax.plot([x0 + t[i], x0 + t[i]], [y0 + amp*np.sin(ph[i]), y0 - amp*np.sin(ph[i])],
                    color=shade(c1, -0.25), lw=0.45, zorder=z)

def protein_blob(ax, x, y, w, h, color='#7fd3cf', edge='#1a6b66', z=6, seed=4, label=None, fs=5.2):
    rng = np.random.default_rng(seed)
    ang = np.linspace(0, 2*np.pi, 160)
    rr = 1 + 0.14*np.sin(5*ang + rng.uniform(0, 6)) + 0.09*np.sin(9*ang + rng.uniform(0, 6))
    verts = np.c_[x + (w/2)*rr*np.cos(ang), y + (h/2)*rr*np.sin(ang)]
    ax.add_patch(Polygon(verts, closed=True, facecolor=shade(color, 0.18), edgecolor=edge,
                         linewidth=0.55, zorder=z))
    for k in range(6):
        t = (k+1)/6
        ax.add_patch(Ellipse((x - w*0.10, y - h*0.10), w*(1-0.5*t)*0.92, h*(1-0.5*t)*0.92,
                             facecolor=mix(shade(color, 0.18), shade(color, -0.15), t),
                             edgecolor='none', zorder=z+1))
    if label:
        ax.text(x, y, label, ha='center', va='center', fontsize=fs, zorder=z+3, color='#1a1a1a')

def bone_slice(ax, x, y, w, h, z=3, color='#EFE6D8', edge='#B9A98F'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=1.2",
                                facecolor=color, edgecolor=edge, linewidth=0.6, zorder=z))
    rng = np.random.default_rng(7)
    for _ in range(26):
        cx = x + rng.uniform(0.12, 0.88)*w; cy = y + rng.uniform(0.12, 0.88)*h
        rx = rng.uniform(0.04, 0.11)*w; ry = rng.uniform(0.06, 0.16)*h
        ax.add_patch(Ellipse((cx, cy), rx*2, ry*2, facecolor='white',
                             edgecolor=shade(edge, -0.05), linewidth=0.35, zorder=z+1))

def bilayer_arc(ax, cx, cy, r, a0, a1, color='#A4A9D9', z=4, half=1.6, n=46):
    th = np.linspace(np.deg2rad(a0), np.deg2rad(a1), n)
    ax.plot(cx + r*np.cos(th), cy + r*np.sin(th), color='white', lw=half*2, zorder=z, solid_capstyle='butt')
    for d in (-half, half):
        ax.plot(cx + (r + d*1.35)*np.cos(th), cy + (r + d*1.35)*np.sin(th),
                color=shade(color, -0.25), lw=0.5, zorder=z+1)
        ax.plot(cx + (r + d)*np.cos(th), cy + (r + d)*np.sin(th), color=color, lw=0.9, zorder=z+1)
        ax.scatter(cx + (r + d)*np.cos(th), cy + (r + d)*np.sin(th), s=1.1, color=color,
                   zorder=z+2, linewidths=0)

def receptor3d(ax, x, y, angle, size, color='#0072B2', z=8):
    a = np.deg2rad(angle)
    u = np.array([np.cos(a), np.sin(a)]); nrm = np.array([-u[1], u[0]])
    base = np.array([x, y])
    ax.plot([*(base)], [*(base + u*size*0.32)], color=shade(color, -0.3), lw=1.6, zorder=z)
    for s in (-1, 1):
        c = base + u*size*0.62 + nrm*s*size*0.16
        ax.add_patch(FancyBboxPatch((c[0]-size*0.20, c[1]-size*0.10), size*0.40, size*0.20,
                                    boxstyle="round,pad=0.01,rounding_size=0.12",
                                    facecolor=shade(color, 0.25 if s > 0 else -0.1),
                                    edgecolor=shade(color, -0.4), lw=0.45, zorder=z+1))
    ax.plot([base[0]+u[0]*size*0.32, base[0]+u[0]*size*0.45],
            [base[1]+u[1]*size*0.32, base[1]+u[1]*size*0.45], color=shade(color, -0.4),
            lw=0.5, zorder=z)

def arrow_flat(ax, a, b, color='#5c5f66', lw=0.8, head=3.0, z=10, style='-|>'):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=head,
                                 color=color, lw=lw, shrinkA=0, shrinkB=0, zorder=z,
                                 joinstyle='miter'))

def banner(ax, x, y, w, h, text, fill, fs=6.2, z=12):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.9",
                                facecolor=fill, edgecolor='none', zorder=z))
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fs,
            color='white', fontweight='bold', zorder=z+1)

def card(ax, x, y, w, h, line1, line2, color, fs1=5.0, fs2=6.0, z=8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.7",
                                facecolor='white', edgecolor=color, linewidth=0.5, zorder=z))
    ax.text(x + 1.2, y + h - 2.6, line1, fontsize=fs1, color='#333', va='center', zorder=z+1)
    ax.text(x + 1.2, y + 1.9, line2, fontsize=fs2, color=color, fontweight='bold', va='center', zorder=z+1)
