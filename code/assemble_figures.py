#!/usr/bin/env python3
"""矢量组图引擎 v2：严格按 A→F 阅读顺序排布（从左到右、从上到下）"""
import os
import fitz, os

ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARIALB = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

def assemble(out, panels, page_w=180, page_h=240, fit_pdf=False, y_from_bottom=True,
             target_width=210.0):
    """target_width=210 → 所有版面按 A4 宽度（210 mm）等比放大。"""
    S = 2.8346
    if target_width and abs(target_width - page_w) > 0.1:
        _k = target_width / page_w
        panels = [(s_, x*_k, y*_k, w*_k, h*_k, lab) for (s_, x, y, w, h, lab) in panels]
        page_h = page_h * _k
        page_w = target_width
    if page_h > 296:                     # A4 高度上限 297 mm
        _k = 292.0 / page_h
        panels = [(s_, x*_k, y*_k, w*_k, h*_k, lab) for (s_, x, y, w, h, lab) in panels]
        page_h = 292
    doc = fitz.open()
    page = doc.new_page(width=page_w*S, height=page_h*S)
    page.insert_font(fontname="Arial", fontfile=ARIAL)
    page.insert_font(fontname="ArialB", fontfile=ARIALB)
    for (src, x, y, w, h, label) in panels:
        if not os.path.exists(src):
            print(f"  跳过缺失面板: {src}")
            continue
        # 面板坐标按 matplotlib 习惯（y 从下往上）书写；PDF 原点在左上，故翻转
        y = page_h - y - h if y_from_bottom else y
        rect = fitz.Rect(x*S, y*S, (x+w)*S, (y+h)*S)
        if src.endswith(".pdf"):
            s = fitz.open(src)
            if fit_pdf:
                sr = s[0].rect
                sc = min(rect.width/sr.width, rect.height/sr.height)
                dw, dh = sr.width*sc, sr.height*sc
                dx = rect.x0 + (rect.width-dw)/2
                dy = rect.y0 + (rect.height-dh)/2
                page.show_pdf_page(fitz.Rect(dx, dy, dx+dw, dy+dh), s, 0)
            else:
                page.show_pdf_page(rect, s, 0)
            s.close()
        else:
            # PNG 保持纵横比，白底 letterbox 居中，不拉伸
            img = fitz.open(src)
            pix = img[0]
            iw, ih = pix.rect.width, pix.rect.height
            img.close()
            scale = min(rect.width/iw, rect.height/ih)
            dw, dh = iw*scale, ih*scale
            dx = rect.x0 + (rect.width-dw)/2
            dy = rect.y0 + (rect.height-dh)/2
            page.draw_rect(fitz.Rect(rect.x0, rect.y0, rect.x0+rect.width, rect.y0+rect.height),
                           color=None, fill=(1,1,1))
            page.insert_image(fitz.Rect(dx, dy, dx+dw, dy+dh), filename=src)
        page.insert_text(fitz.Point(x*S+1.0, y*S+5.4), label, fontsize=12.0,
                         fontname="ArialB", color=(0.1,0.1,0.1))
    doc.save(out)
    print("saved:", out)

def crop_to_content(path, margin=8):
    """按内容 bbox 裁剪页面高度，去掉多余白边"""
    src = fitz.open(path); pg = src[0]
    W, H = pg.rect.width, pg.rect.height
    x0 = y0 = 1e9; x1 = y1 = -1e9
    for b in pg.get_text("blocks"):
        x0 = min(x0, b[0]); y0 = min(y0, b[1]); x1 = max(x1, b[2]); y1 = max(y1, b[3])
    for d in pg.get_drawings():
        r = d["rect"]; x0 = min(x0, r.x0); y0 = min(y0, r.y0); x1 = max(x1, r.x1); y1 = max(y1, r.y1)
    for im in pg.get_images(full=True):
        for r in pg.get_image_rects(im[0]):
            x0 = min(x0, r.x0); y0 = min(y0, r.y0); x1 = max(x1, r.x1); y1 = max(y1, r.y1)
    newH = min(H, max(60.0, y1 + margin))
    if H - newH < 5:
        src.close(); return
    out = fitz.open(); npx = out.new_page(width=W, height=newH)
    npx.show_pdf_page(fitz.Rect(0, 0, W, H), src, 0)
    out.save(path, garbage=3, deflate=True); out.close(); src.close()

T1 = "选题1_胞啃Trogocytosis_肺癌/figures/"
T2 = "选题2_机械刺激_卵巢癌/figures/"
T3 = "选题3_网络毒理学_AI骨质疏松/figures/"
os.makedirs("投稿材料", exist_ok=True)

# ============ 选题1（8 主图，严格阅读顺序）============
P = "投稿材料/选题1_正图/"
os.makedirs(P, exist_ok=True)
assemble(P+"Figure1.pdf", [
 (T1+"Fig1A_flowchart_v6.pdf", 5, 100, 170, 135, "A"),
 (T1+"Fig1B_pancan_lollipop.pdf", 5, 5, 82, 90, "B"),
 (T1+"Fig1C_immune_corr.pdf", 95, 5, 80, 90, "C"),
], page_h=240, fit_pdf=True)
assemble(P+"Figure2.pdf", [
 (T1+"Fig3A_univariate_forest.pdf", 5, 225, 55, 90, "A"),
 (T1+"Fig3B_lasso_path.pdf", 65, 225, 55, 90, "B"),
 (T1+"Fig3_multivariate_forest.pdf", 125, 225, 50, 90, "C"),
 (T1+"Fig3D_training_KM.pdf", 5, 125, 55, 90, "D"),
 (T1+"Fig3E_timeROC.pdf", 65, 125, 55, 90, "E"),
 (T1+"Fig3F_model_translation.pdf", 125, 125, 50, 90, "F"),
 (T1+"Fig3I_risk_plot.pdf", 5, 5, 82, 110, "G"),
 (T1+"Fig3J_calibration.pdf", 92, 5, 83, 110, "H"),
 (T1+"Fig3G_benchmark.pdf", 5, 155, 60, 70, "I"),
 (T1+"Fig3H_head2head.pdf", 70, 155, 45, 70, "J"),
], page_h=320, fit_pdf=True)
assemble(P+"Figure3.pdf", [
 (T1+"Fig4_A_GSE15459_KM.pdf", 5, 115, 50, 120, "A"),
 (T1+"Fig4_B_GSE84437_KM.pdf", 60, 115, 50, 120, "B"),
 (T1+"Fig4_C_GSE26253_RFS.pdf", 115, 115, 60, 120, "C"),
 (T1+"Fig4_D_meta_forest.pdf", 5, 5, 60, 110, "D"),
 (T1+"Fig4D_val_score_dist.pdf", 70, 5, 105, 110, "E"),
], fit_pdf=True)
assemble(P+"Figure7.pdf", [
 (T1+"Fig7A_UMAP_mixed.pdf", 5, 205, 60, 130, "A"),
 (T1+"Fig9_scRNA_TrogoScore_violin.pdf", 70, 205, 105, 130, "B"),
 (T1+"Fig7C_inflamed_desert.pdf", 5, 110, 55, 90, "C"),
 (T1+"Fig7D_featureplots_mixed.pdf", 65, 110, 55, 90, "D"),
 (T1+"Fig7E_celltype_composition.pdf", 125, 110, 50, 90, "E"),
 (T1+"Fig7F_gene_dotplot.pdf", 5, 5, 60, 100, "F"),
 (T1+"Fig7G_cellchat_heatmap.pdf", 70, 5, 105, 100, "G"),
], page_h=340, fit_pdf=True)
assemble(P+"Figure8.pdf", [
 (T1+"Fig8A_cellchat_chord_inflamed.pdf", 5, 150, 50, 85, "A"),
 (T1+"Fig8A2_cellchat_chord_desert.pdf", 60, 150, 50, 85, "B"),
 (T1+"Fig8C_NRP1_density.pdf", 115, 150, 60, 85, "C"),
 (T1+"Fig8D_alt1_circular.pdf", 5, 5, 55, 85, "D"),
 (T1+"Fig8D3_NRP1_vko_GO.pdf", 60, 5, 30, 85, "E"),
 (T1+"docking/Fig8E_pose_labeled.pdf", 90, 5, 50, 85, "F"),
 (T1+"docking/Fig8E_combined_bar.pdf", 145, 5, 30, 85, "G"),
 (T1+"Fig8G_NRP1_box.pdf", 5, 155, 25, 70, "H"),
], fit_pdf=True)

assemble(P+"Figure9.pdf", [
 (T1+"Fig9_immunotherapy.pdf", 5, 5, 170, 115, ""),
], fit_pdf=True)

# ============ 选题2（6 主图，严格阅读顺序）============
P = "投稿材料/选题2_正图/"
os.makedirs(P, exist_ok=True)
assemble(P+"Figure1.pdf", [
 (T2+"Fig1A_flowchart_v3.pdf", 5, 100, 170, 135, "A"),
 (T2+"Fig1C_volcano_mech.pdf", 5, 5, 55, 90, "B"),
 (T2+"Fig1C_score_tumor_normal.pdf", 65, 5, 55, 90, "C"),
 (T2+"Fig1D_geneset_nested.pdf", 125, 5, 50, 90, "D"),
], page_h=240, fit_pdf=True)
assemble(P+"Figure2.pdf", [
 (T2+"Fig2A_cluster_heatmap.pdf", 5, 5, 55, 110, "A"),
 (T2+"Fig2B_subtype_KM.pdf", 65, 5, 55, 110, "B"),
 (T2+"Fig2C_mech_GO.pdf", 125, 5, 50, 110, "C"),
], fit_pdf=True)
assemble(P+"Figure3.pdf", [
 (T2+"Fig3A_univariate_forest.pdf", 5, 155, 55, 80, "A"),
 (T2+"Fig3B_lasso_path.pdf", 65, 155, 55, 80, "B"),
 (T2+"Fig3C_forest_12genes.pdf", 125, 155, 50, 80, "C"),
 (T2+"Fig3_KM_training.pdf", 5, 80, 55, 70, "D"),
 (T2+"Fig3E_timeROC.pdf", 65, 80, 55, 70, "E"),
 (T2+"Fig3F_score_dist.pdf", 125, 80, 50, 70, "F"),
], fit_pdf=True)
assemble(P+"Figure4.pdf", [
 (T2+"Fig4_ICON7_PFS_KM.pdf", 5, 5, 45, 110, "A"),
 (T2+"Fig4B_multi_cohort_forest.pdf", 55, 5, 60, 110, "B"),
 (T2+"Fig4C_PFS_other.pdf", 120, 120, 55, 115, "C"),
 (T2+"Fig4D_meta_LOO.pdf", 120, 5, 55, 110, "D"),
 (T2+"Fig4E_arm_forest.pdf", 5, 120, 55, 115, "E"),
 (T2+"Fig4F_head2head_oxc.pdf", 65, 120, 110, 115, "F"),
], page_h=240, fit_pdf=True)
assemble(P+"Figure5.pdf", [
 (T2+"Fig6A_SCISSOR_UMAP.pdf", 5, 125, 85, 110, "A"),
 (T2+"Fig6B_cellchat_chord_high.pdf", 95, 125, 40, 110, "B"),
 (T2+"Fig6B2_cellchat_chord_low.pdf", 140, 125, 35, 110, "C"),
 (T2+"Fig6B3_scissor_proportions.pdf", 5, 5, 55, 110, "D"),
 (T2+"Fig6C_cellchat_diff.pdf", 65, 5, 55, 110, "E"),
 (T2+"Fig5F_scissor_counts.pdf", 125, 5, 50, 110, "F"),
], fit_pdf=True)
assemble(P+"Figure6.pdf", [
 (T2+"Fig7A_key_gene_corr.pdf", 5, 120, 170, 110, "A"),
 (T2+"Fig7B_PPL_vko.pdf", 5, 5, 60, 108, "B"),
 (T2+"Fig6C_PPL_perturbation_bar.pdf", 70, 5, 70, 108, "C"),
 (T2+"Fig7D_gene_cor.pdf", 145, 5, 30, 108, "D"),
 (T2+"Fig6E_ici_pancancer.pdf", 5, 130, 100, 105, "E"),
], page_h=240, fit_pdf=True)

# ============ 选题3（7 主图，严格阅读顺序）============
P = "投稿材料/选题3_正图/"
os.makedirs(P, exist_ok=True)
assemble(P+"Figure1.pdf", [
 (T3+"Fig1A_workflow_a4.pdf", 5.0, 168.0, 170.0, 88.0, "A"),
 (T3+"Fig1B_FAERS_ROR_forest.pdf", 5.0, 113.0, 56.0, 52.0, "B"),
 (T3+"Fig1C_FAERS_trend.pdf", 64.0, 113.0, 56.0, 52.0, "C"),
 (T3+"Fig1D_consistency.pdf", 5.0, 58.0, 56.0, 52.0, "D"),
 (T3+"Fig1E_stratification.pdf", 64.0, 58.0, 56.0, 52.0, "E"),
 (T3+"Fig1F_exemestane_forest.pdf", 5.0, 3.0, 56.0, 52.0, "F"),
], page_h=259, fit_pdf=True)
assemble(P+"Figure2.pdf", [
 (T3+"Fig3A_intersection_venn.pdf", 5.0, 168.0, 84.0, 50.0, "A"),
 (T3+"Fig3B_targetdb_venn.pdf", 94.0, 168.0, 84.0, 50.0, "B"),
 (T3+"Fig3C_PPI_network.pdf", 5.0, 112.0, 56.0, 52.0, "C"),
 (T3+"Fig4C_compound_target_pathway.pdf", 64.0, 112.0, 62.0, 52.0, "D"),
 (T3+"Fig3F_evidence_matrix.pdf", 129.0, 112.0, 49.0, 52.0, "E"),
 (T3+"Fig4A_GO_three.pdf", 5.0, 56.0, 150.0, 52.0, "F"),
 (T3+"Fig3G_string_edges.pdf", 5.0, 4.0, 84.0, 48.0, "G"),
], page_h=222, fit_pdf=True)
assemble(P+"Figure3.pdf", [
 (T3+"Fig2A_GSE35958_volcano.pdf", 5.0, 118.0, 84.0, 52.0, "A"),
 (T3+"Fig2B_volcano_core_genes.pdf", 94.0, 118.0, 84.0, 52.0, "B"),
 (T3+"Fig2C_core_expression.pdf", 5.0, 60.0, 84.0, 52.0, "C"),
 (T3+"Fig2D_GSE13850_volcano.pdf", 94.0, 60.0, 84.0, 52.0, "D"),
 (T3+"Fig2E_core_logFC.pdf", 5.0, 3.0, 84.0, 52.0, "E"),
 (T3+"Fig2F_cohort_degs.pdf", 94.0, 3.0, 84.0, 52.0, "F"),
], page_h=174, fit_pdf=True)
assemble(P+"Figure4.pdf", [
 (T3+"Fig5A_sft_curve.pdf", 5.0, 115.0, 56.0, 52.0, "A"),
 (T3+"Fig5B_module_hist.pdf", 64.0, 115.0, 62.0, 52.0, "B"),
 (T3+"Fig5_WGCNA_module_trait.pdf", 5.0, 58.0, 56.0, 54.0, "C"),
 (T3+"Fig5D_module_heatmap.pdf", 64.0, 58.0, 66.0, 54.0, "D"),
 (T3+"Fig5E_module_sizes.pdf", 5.0, 3.0, 56.0, 52.0, "E"),
 (T3+"Fig5F_M4_hubs.pdf", 64.0, 3.0, 56.0, 52.0, "F"),
], page_h=170, fit_pdf=True)
assemble(P+"Figure5.pdf", [
 (T3+"Fig6A_lasso_path.pdf", 5.0, 160.0, 70.0, 52.0, "A"),
 (T3+"Fig6B_lasso_cv.pdf", 80.0, 160.0, 62.0, 52.0, "B"),
 (T3+"Fig6C_VPS35_coexpression.pdf", 5.0, 75.0, 95.0, 80.0, "C"),
 (T3+"Fig6D_VPS35_datasets.pdf", 105.0, 75.0, 70.0, 80.0, "D"),
 (T3+"Fig6E_GSE201774_fpkm.pdf", 5.0, 8.0, 62.0, 60.0, "E"),
 (T3+"Fig6F_coexpression_matrix.pdf", 72.0, 8.0, 62.0, 60.0, "F"),
], page_h=220, fit_pdf=True)
assemble(P+"Figure6.pdf", [
 (T3+"Fig7A_docking_merged3.pdf", 5.0, 216.0, 110.0, 68.0, "A"),
 (T3+"Fig7C_VPS35_lollipop.pdf", 120.0, 216.0, 60.0, 52.0, "B"),
 (T3+"Fig7C_mechanism_pathway.pdf", 5.0, 128.0, 130.0, 84.0, "C"),
 (T3+"Fig7C3_VPS35_vko_GO.pdf", 5.0, 52.0, 70.0, 72.0, "D"),
 (T3+"Fig7E_docking_energies2.pdf", 80.0, 52.0, 56.0, 72.0, "E"),
 (T3+"Fig7G_VPS35_contacts.pdf", 5.0, 0.0, 56.0, 48.0, "F"),
], page_h=290, fit_pdf=True)
print("全部正图按阅读顺序重组完成")

assemble(P+"Figure7.pdf", [
 (T3+"Fig8A_cross_species_energy.pdf", 5.0, 232.0, 84.0, 56.0, "A"),
 (T3+"Fig8B_active_site_conservation.pdf", 94.0, 232.0, 82.0, 56.0, "B"),
 (T3+"Fig8C_zebrafish_pose.pdf", 5.0, 166.0, 112.0, 62.0, "C"),
 (T3+"S12_exposure_EAR.pdf", 5.0, 96.0, 88.0, 62.0, "D"),
 (T3+"Fig8D_crossspecies_summary.pdf", 5.0, 2.0, 170.0, 102.0, "E"),
], page_h=302, fit_pdf=True)
# ============ T3 七图统一裁白边 ============
for _i in range(1, 8):
    crop_to_content(f"投稿材料/选题3_正图/Figure{_i}.pdf")
print("T3 七图白边裁剪完成")
