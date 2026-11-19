#!/usr/bin/env python3
"""Stack panel PNGs into journal-style composite figures with A/B/C labels."""
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts")
OUT = ROOT / "paper" / "composites"
OUT.mkdir(parents=True, exist_ok=True)

PAD = 28
LABEL_W = 72
GAP = 18
BG = (255, 255, 255)


def font(size: int):
    for p in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
    ):
        if Path(p).exists():
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


def stack(paths: list[Path], labels: list[str], out: Path, max_w: int = 2400) -> None:
    imgs = []
    for p in paths:
        im = Image.open(p).convert("RGB")
        if im.width > max_w:
            h = int(im.height * max_w / im.width)
            im = im.resize((max_w, h), Image.Resampling.LANCZOS)
        imgs.append(im)
    w = max(im.width for im in imgs)
    scaled = []
    for im in imgs:
        if im.width != w:
            h = int(im.height * w / im.width)
            im = im.resize((w, h), Image.Resampling.LANCZOS)
        scaled.append(im)
    total_h = PAD + sum(im.height for im in scaled) + GAP * (len(scaled) - 1) + PAD
    total_w = LABEL_W + w + PAD
    canvas = Image.new("RGB", (total_w, total_h), BG)
    draw = ImageDraw.Draw(canvas)
    fnt = font(42)
    y = PAD
    for lab, im in zip(labels, scaled):
        canvas.paste(im, (LABEL_W, y))
        draw.text((16, y + 8), lab, fill=(0, 0, 0), font=fnt)
        y += im.height + GAP
    canvas.save(out, "PNG", dpi=(200, 200))
    print(f"wrote {out}  {canvas.size}")


def extract_refined_panel_b() -> Path:
    """Panel B of the refined TOPMed-vs-AoU scatter, as-is, for Figure 2C.

    The source PNG is a two-panel figure; ylabel lives only on panel A, so a
    matching y-axis label is drawn in the left margin of the crop.
    """
    src = (
        ROOT
        / "fig1_overview/panel_c_topmed_vs_aou/topmed_vs_new_haplo_prevalence_scatter_refined.png"
    )
    dest_dir = ROOT / "fig2_constraint" / "panel_c_topmed_vs_aou_loeuf"
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / "fig2_topmed_vs_aou_prevalence_loeuf_panelB.png"
    im = Image.open(src).convert("RGB")
    w, h = im.size
    # Skip the shared title (rows ~21-43). Tick labels of panel B start ~960.
    crop = im.crop((960, 48, w, h))
    label_w = 52
    canvas = Image.new("RGB", (label_w + crop.width, crop.height), BG)
    canvas.paste(crop, (label_w, 0))
    draw = ImageDraw.Draw(canvas)
    fnt = font(22)
    label = "-log10 AoU RNA-outlier prevalence"
    bbox = draw.textbbox((0, 0), label, font=fnt)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    txt = Image.new("RGBA", (tw + 8, th + 8), (255, 255, 255, 0))
    ImageDraw.Draw(txt).text((4, 4), label, fill=(0, 0, 0), font=fnt)
    rot = txt.rotate(90, expand=True)
    x = 6
    y = (canvas.height - rot.height) // 2
    canvas.paste(rot, (x, y), rot)
    # Cover leftover inner letter "B." from the two-panel source title
    ImageDraw.Draw(canvas).rectangle((408, 68, 455, 102), fill=BG)
    canvas.save(dest, "PNG", dpi=(200, 200))
    print(f"wrote {dest}  {canvas.size}")
    return dest


def main() -> int:
    fig2c = extract_refined_panel_b()
    stack(
        [
            ROOT / "fig1_overview/panel_a_study_design/fig1_study_design_topmed_aou.png",
            ROOT / "fig1_overview/panel_b_rank_curves/fig1_underlier_tallies_topmed_aou.png",
            ROOT / "fig1_overview/panel_c_topmed_vs_aou/fig1_topmed_vs_aou_prevalence_uncolored.png",
        ],
        ["A", "B", "C"],
        OUT / "Figure1.png",
    )
    stack(
        [
            ROOT / "fig2_constraint/panel_a_loeuf_calibration/loeuf_calibration_z_-3.png",
            ROOT / "fig2_constraint/panel_b_constrained_vs_tolerant/underlier_tallies_loeuf_z_-3_constrained_tolerant.png",
            fig2c,
        ],
        ["A", "B", "C"],
        OUT / "Figure2.png",
    )
    stack(
        [
            ROOT / "fig3_ptv_vs_rna/panel_a_gene_level/ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4.png",
            ROOT / "fig3_ptv_vs_rna/panel_b_carrier_zscore/ptv_vs_nonptv_zscore.png",
            ROOT / "fig3_ptv_vs_rna/panel_c_clingen_examples/candidate_gene_zscore_dist_clingen_z_-3.png",
        ],
        ["A", "B", "C"],
        OUT / "Figure3.png",
        max_w=2200,
    )
    stack(
        [
            ROOT / "fig4_acmg/panel_a_rbm20/fig4_RBM20_4dot_no_ci.png",
            ROOT / "fig4_acmg/panel_b_myh11/fig4_MYH11_4dot_no_ci.png",
        ],
        ["A", "B"],
        OUT / "Figure4.png",
    )
    stack(
        [
            ROOT / "supplement/pc_enrichment/enrichment_npc_sv_vs_ptv.png",
            ROOT / "supplement/ptv_consequences/ptv_consequence_summary_qc_cohort.png",
        ],
        ["A", "B"],
        OUT / "Supplementary_Figure1.png",
        max_w=1600,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
