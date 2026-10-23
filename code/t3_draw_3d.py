"""选题3 Fig1A 工作流 + Fig7D 机制图：3D 插画风格矢量重绘"""
import sys; sys.path.insert(0, "code")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle, Polygon, Rectangle
import t3_fig3d as T

BLUE='#0072B2'; LBLUE='#56B4E9'; ORANGE='#E69F00'; GREEN='#009E73'; RED='#D55E00'; PURPLE='#7B6AA8'; GREY='#5c5f66'
T3 = "选题3_网络毒理学_AI骨质疏松/"
plt.rcParams.update({"font.family":"Arial","font.size":6,"pdf.fonttype":42})

def canvas(w=183, h=110):
    fig = plt.figure(figsize=(w/25.4, h/25.4))
    ax = fig.add_axes([0,0,1,1]); ax.set_xlim(0,w); ax.set_ylim(h,0); ax.axis('off')
    return fig, ax

def save(fig, stem):
    fig.savefig(stem+".pdf", format="pdf"); fig.savefig(stem+".png", dpi=600); plt.close(fig)
    print("saved:", stem)

# ================= Fig 1A 工作流 =================
def fig1a():
    fig, ax = canvas()
    cols = [(2,42,'DATA SOURCES',BLUE),(47,42,'NETWORK MAPPING',GREY),(92,42,'KEY OUTPUTS',ORANGE),(137,42,'VALIDATION',RED)]
    for x,w,t,c in cols: T.banner(ax,x,3,w,7,t,c)

    # --- col1 数据源
    T.cylinder(ax, 20, 22, 16, 6, '#3C5488'); ax.text(29.5,24,'FAERS (2013–2026)',fontsize=5.0,va='center')
    T.cylinder(ax, 20, 38, 16, 5.5, '#3C5488'); ax.text(29.5,40,'ChEMBL/SwissTarget',fontsize=5.0,va='center')
    T.cylinder(ax, 20, 53, 16, 5.5, '#3C5488'); ax.text(29.5,55,'GeneCards (3,713)',fontsize=5.0,va='center')
    T.cell3d(ax, 18, 67, 5.5, '#EAF3FB', '#4A7FA5'); ax.text(26,67,'GSE35958 (MSC)',fontsize=5.0,va='center')
    T.cell3d(ax, 18, 79, 5.5, '#FDECEC', '#C08080', seed=5); ax.text(26,79,'GSE13850 (B cells)',fontsize=5.0,va='center')
    T.cell3d(ax, 18, 91, 5.0, '#E8F4EC', '#5c9a72', seed=9); ax.text(26,91,'GSE201774 / 56814/15',fontsize=5.0,va='center')
    T.dna_helix(ax, 33, 84, 14, amp=1.3, period=4.5, vertical=True, c1=PURPLE, c2='#C3B1DE')
    ax.text(37,104,'single-cell atlas',fontsize=4.8,ha='center')

    # --- col2 网络
    T.card(ax,48,15,40,9,'Disproportionality','ROR / PRR / E-value',BLUE)
    # 3D Venn
    T.sphere(ax,64,32,6.2,'#8FB3D9'); T.sphere(ax,72.5,32,6.2,'#E7B08C')
    T.sphere(ax,68.3,32,2.2,'#C9A9C9')
    ax.text(68.3,32,'9',fontsize=4.8,ha='center',va='center',fontweight='bold')
    ax.text(68.3,40.5,'42 targets ∩ 3,713 genes',fontsize=4.8,ha='center')
    # PPI 网络（分层节点 + 边）
    nodes = [[56,50],[62,49],[68,50],[74,49],[80,50],[59,56],[65,56],[71,56],[77,56],[68,62]]
    for i,j in [(0,1),(1,2),(2,3),(3,4),(0,5),(1,5),(2,6),(3,7),(4,8),(5,9),(6,9),(7,9),(8,9)]:
        ax.plot([nodes[i][0],nodes[j][0]],[nodes[i][1],nodes[j][1]],color='#b9c3cf',lw=0.4,zorder=5)
    for k,(x,y) in enumerate(nodes):
        T.sphere(ax,x,y,1.5,[BLUE,LBLUE,GREY,ORANGE][k%4],z=6)
    ax.text(68,66,'STRING PPI',fontsize=4.8,ha='center')
    # 火山散点
    rng = np.random.default_rng(11)
    ax.plot([48,86],[76,76],color='#bbb',lw=0.4)
    for i in range(30):
        x = rng.uniform(49,85); y = 76 - abs(rng.normal(0,1))*1.6
        ax.scatter(x, max(y,68.5), s=1.6, color=ORANGE if i%4==0 else '#c9c9c9', zorder=6, linewidths=0)
    ax.text(68,80,'limma (two cohorts)',fontsize=4.8,ha='center')
    # LASSO 路径
    t = np.linspace(0,1,80)
    for k in range(9):
        ax.plot(50+t*32, 92 - (t**1.6)*(9-k)*0.6 - k*0.35, color='#c9c9c9', lw=0.45, zorder=6)
    ax.plot(50+t*32, 92 - (t**1.6)*4.2, color=RED, lw=0.8, zorder=7)
    ax.axvline(72, ymin=0.0001, ymax=0.0001, color='k')
    ax.plot([72,72],[84.5,99.5],color='k',lw=0.4,ls='--',zorder=7)
    ax.text(68,102.5,'WGCNA + LASSO',fontsize=4.8,ha='center')

    # --- col3 关键输出
    T.card(ax,93,15,40,9,'FAERS signals','ROR 2.77 – 4.93',ORANGE)
    T.card(ax,93,26,40,9,'E-value range','4.98 – 9.33',ORANGE)
    # 9 基因网络
    g9 = [[112,40],[105,44],[119,44],[105,51],[119,51],[112,48],[108,57],[116,57],[112,62]]
    for i in range(len(g9)):
        for j in range(i+1,len(g9)):
            if np.hypot(g9[i][0]-g9[j][0], g9[i][1]-g9[j][1]) < 8:
                ax.plot([g9[i][0],g9[j][0]],[g9[i][1],g9[j][1]],color='#c8d2de',lw=0.4,zorder=5)
    for k,(x,y) in enumerate(g9):
        T.sphere(ax,x,y,1.9 if k!=8 else 2.4, RED if k==8 else (BLUE if k<5 else GREY), z=6)
    ax.annotate('CYP19A1', xy=(112,62), xytext=(122,66), fontsize=4.8, color=RED,
                arrowprops=dict(arrowstyle='-', lw=0.4, color=RED))
    ax.text(112,72,'9 core genes',fontsize=4.8,ha='center')
    T.card(ax,93,75,40,9,'Transcriptome','DEG null × 2 cohorts',GREY)
    # VPS35 蛋白 + 对接
    T.protein_blob(ax, 102, 90, 16, 9, '#CFE8DF', '#2E8B57', label='VPS35')
    T.sphere(ax, 113, 89, 2.6, PURPLE)
    ax.text(117,90,'candidate',fontsize=4.6,va='center',color=PURPLE)
    for dx,dy,r in [(0,0,3.4),(-4.5,1.8,2.6),(4.5,-1.4,2.9)]:
        T.sphere(ax, 124+dx, 100+dy, r, '#C9D6E5', z=6)
    ax.text(124,105,'docking (AF2)',fontsize=4.6,ha='center')

    # --- col4 验证
    T.card(ax,138,15,42,9,'Exemestane FAERS','ROR 3.61',RED)
    T.card(ax,138,26,42,9,'Demographics','390/394 female',RED)
    T.card(ax,138,37,42,9,'GSE201774','VPS35 FPKM 12.97',RED)
    T.card(ax,138,48,42,9,'Honest reporting','negatives kept',GREY)
    T.arrow_flat(ax,(141,64),(177,64),color=RED,lw=1.0,head=3.2)
    ax.text(159,61.5,'AI → CYP19A1 → estrogen ↓',fontsize=5.0,ha='center',color=RED,fontweight='bold')
    ax.text(159,67.5,'→ VPS35 perturbation',fontsize=5.0,ha='center',color=RED,fontweight='bold')
    T.cell3d(ax,150,84,6.5,'#F7EEE6','#B08A6E',seed=13)
    ax.text(150,93.5,'osteoblast model',fontsize=4.8,ha='center')
    T.receptor3d(ax, 141, 82, -90, 8, BLUE)
    ax.text(139,72,'FZD6',fontsize=4.6,ha='center',color=BLUE)
    T.protein_blob(ax, 172, 86, 11, 7, '#EBD9F0', PURPLE, seed=8)
    ax.text(172,94,'Wnt',fontsize=4.8,ha='center',color=PURPLE)

    # 列间箭头
    for y in (22,38,52,66,80,94):
        T.arrow_flat(ax,(44.5,y),(46.5,y)); T.arrow_flat(ax,(89.5,y),(91.5,y))
    for y in (22,38,52,66,86):
        T.arrow_flat(ax,(133.5,y),(136.5,y))
    save(fig, T3+"figures/Fig1A_workflow3d")

# ================= Fig7D 机制图 =================
def fig7d():
    fig, ax = canvas()
    ax.text(6,6,'Aromatase inhibitor (letrozole / anastrozole)',fontsize=6.4,fontweight='bold')
    ax.text(96,6,'CYP19A1 inhibition → estrogen deficiency',fontsize=6.4,fontweight='bold')
    ax.text(150,6,'VPS35–retromer / Wnt axis in bone',fontsize=6.4,fontweight='bold')

    # --- 左：药物
    ax.add_patch(FancyBboxPatch((10,20),26,8,boxstyle="round,pad=0.02,rounding_size=2.0",
                                facecolor='white',edgecolor=PURPLE,lw=0.7,zorder=8))
    ax.add_patch(FancyBboxPatch((10,20),13,8,boxstyle="round,pad=0.02,rounding_size=2.0",
                                facecolor=T.shade(PURPLE,0.75),edgecolor='none',zorder=9))
    ax.text(23,24,'aromatase inhibitor',fontsize=5.0,ha='center',va='center',zorder=10)
    T.arrow_flat(ax,(30,34),(30,40),color=GREY)
    # CYP19A1 酶
    T.protein_blob(ax, 30, 47, 26, 12, '#FFD9A8', '#B07A2D', label='CYP19A1', fs=5.4)
    ax.text(30,56,'aromatase (target of docking −7.22)',fontsize=4.6,ha='center',color='#7a5418')
    # 雌激素分子（3 个球）+ 下降箭头
    for i,(x,y) in enumerate([(20,66),(30,66),(40,66)]):
        T.sphere(ax,x,y,2.6,'#E58FB1',z=8)
    ax.text(30,72,'estradiol (E2)',fontsize=5.0,ha='center',color='#a04464')
    T.arrow_flat(ax,(30,76),(30,84),color=RED,lw=1.0,head=3.4)
    ax.text(31.5,80,'↓',fontsize=9,color=RED,va='center')

    # --- 中：骨吸收
    T.bone_slice(ax, 8, 90, 52, 14, color='#EFE6D8', edge='#B9A98F')
    ax.text(34,101,'bone: osteoclast activation (RANKL ↑ / OPG ↓)',fontsize=4.8,ha='center',color='#7a6a52')
    # --- 右：成骨细胞膜 + endosome + Wnt
    T.bilayer_arc(ax, 128, 58, 34, -70, 60, '#A4A9D9', half=1.4)
    T.receptor3d(ax, 150, 34, -60, 12, BLUE)
    ax.text(156,26,'FZD6',fontsize=5.0,color=BLUE)
    T.vesicle_wrap = None
    ax.add_patch(Ellipse((112,54),26,17,facecolor='#F6E3F0',edgecolor='#B0619A',lw=0.6,zorder=7))
    ax.add_patch(Ellipse((112,54),20,13,facecolor='#FDF6FA',edgecolor='#D9A5C6',lw=0.45,zorder=8))
    T.protein_blob(ax, 112, 54, 13, 8, '#D9C6E8', PURPLE, seed=6, label='VPS35', fs=4.8)
    ax.text(112,64,'retromer endosome',fontsize=4.6,ha='center',color='#7b4f72')
    for i,(x,y) in enumerate([(140,58),(146,52),(152,46)]):
        T.sphere(ax,x,y,1.9,'#9FD6C8',z=8)
    ax.text(146,64,'Wnt ligands',fontsize=4.6,ha='center',color='#3d8f7c')
    T.arrow_flat(ax,(132,54),(146,50),color=GREY)
    # 核
    ax.add_patch(Ellipse((160,80),26,16,facecolor='#E7EAF6',edgecolor='#7c86b8',lw=0.6,zorder=7))
    ax.add_patch(Ellipse((160,80),20,11,facecolor='#F4F6FD',edgecolor='#A9B2DB',lw=0.4,zorder=8))
    T.dna_helix(ax, 148, 76, 22, amp=1.1, period=4.0, c1='#7c86b8', c2='#B7BEE0')
    ax.text(160,90,'nucleus: Wnt target genes',fontsize=4.8,ha='center',color='#4a5488')
    T.arrow_flat(ax,(150,46),(154,68),color=GREY)
    ax.text(158,58,'β-catenin →\nbone formation',fontsize=4.6,ha='left',va='center',color=GREY)
    # 结论轴
    ax.text(92,106,'AI → CYP19A1 ↓ → E2 ↓ → VPS35/retromer perturbation → attenuated Wnt → bone loss',
            fontsize=5.2,ha='center',color='#333',fontweight='bold')
    save(fig, T3+"figures/Fig7D_mechanism_3d")

if __name__ == "__main__":
    fig1a(); fig7d()
