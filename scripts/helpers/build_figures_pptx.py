#!/usr/bin/env python3
"""Rebuild Figures_1_to_4_panels.pptx for the active PC_SETTING.

The published deck is one panel per slide, 15 slides, with a hand-authored
title above each picture. Nothing about the slide masters, theme or titles is
PC-dependent, so the deck is rebuilt by taking the published file as a
read-only template and swapping only the embedded PNG for the corresponding
panel from final_paper_draft_pc<PC>/source_panels/.

Panel aspect ratios are preserved rather than assumed: every picture box in the
template is width-fitted, so the new height is width x (h/w) of the new PNG and
the box is re-centred on its old vertical midpoint. At 55 PCs only Figure 1d
changes shape at all, but doing it generally keeps the deck correct if a panel
is ever resized.

  PC_SETTING=2005 python3 scripts/helpers/build_figures_pptx.py
"""
from __future__ import annotations

import io
import os
import re
import sys
import zipfile
from pathlib import Path

from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import FINAL_DRAFT, PC_IS_DEFAULT, ROOT  # noqa: E402

PPTX_NAME = "Figures_1_to_4_panels.pptx"
TEMPLATE = ROOT / "final_paper_draft_pc55" / PPTX_NAME
OUT = FINAL_DRAFT / PPTX_NAME
PANELS = FINAL_DRAFT / "source_panels"

# slide number -> panel PNG in source_panels/. Order is the published deck's.
SLIDE_PANEL = {
    1: "Fig1a_cohort_overview_combined.png",
    2: "Fig1b_gene_ranks.png",
    3: "Fig1c_subject_ranks.png",
    4: "Fig1d_topmed_vs_aou_uncolored.png",
    5: "Fig2a_loeuf_forest.png",
    6: "Fig2b_loeuf_decile.png",
    7: "Fig2c_constrained_vs_tolerant.png",
    8: "Fig2d_topmed_vs_aou_loeuf.png",
    9: "Fig3a_ptv_vs_underlier.png",
    10: "Fig3b_ptv_vs_nonptv_zscore.png",
    11: "Fig3c_clingen_hist.png",
    12: "Fig4a_RBM20_ranks.png",
    13: "Fig4b_RBM20_4dot.png",
    14: "Fig4c_MYH11_ranks.png",
    15: "Fig4d_MYH11_4dot.png",
}

PIC_RE = re.compile(r"<p:pic>.*?</p:pic>", re.S)
EMBED_RE = re.compile(r'r:embed="(rId\d+)"')
XFRM_RE = re.compile(
    r'(<a:off x=")(-?\d+)("\s+y=")(-?\d+)("\s*/>\s*<a:ext cx=")(\d+)("\s+cy=")(\d+)(")'
)


def rels_map(xml: str) -> dict[str, str]:
    out = {}
    for m in re.finditer(r"<Relationship\b[^>]*/>", xml):
        tag = m.group(0)
        rid = re.search(r'Id="([^"]+)"', tag)
        tgt = re.search(r'Target="([^"]+)"', tag)
        if rid and tgt:
            out[rid.group(1)] = tgt.group(1)
    return out


def main() -> int:
    if PC_IS_DEFAULT and os.environ.get("ALLOW_PC55_WRITE") != "1":
        raise SystemExit(
            "refusing to rewrite the published 55-PC deck; set ALLOW_PC55_WRITE=1 "
            "if that is really intended"
        )
    if not TEMPLATE.exists():
        raise SystemExit(f"missing deck template: {TEMPLATE}")
    missing = [p for p in SLIDE_PANEL.values() if not (PANELS / p).exists()]
    if missing:
        raise SystemExit("missing panels: " + ", ".join(missing))

    zin = zipfile.ZipFile(TEMPLATE, "r")
    slide_xml: dict[str, str] = {}
    media_new: dict[str, bytes] = {}

    for slide_no, panel in SLIDE_PANEL.items():
        sname = f"ppt/slides/slide{slide_no}.xml"
        rname = f"ppt/slides/_rels/slide{slide_no}.xml.rels"
        xml = zin.read(sname).decode("utf-8")
        rels = rels_map(zin.read(rname).decode("utf-8"))

        png = PANELS / panel
        with Image.open(png) as im:
            aspect = im.height / im.width

        pics = PIC_RE.findall(xml)
        if len(pics) != 1:
            raise SystemExit(f"slide {slide_no}: expected 1 picture, found {len(pics)}")
        block = pics[0]
        emb = EMBED_RE.search(block)
        if not emb:
            raise SystemExit(f"slide {slide_no}: picture has no r:embed")
        target = rels[emb.group(1)]  # e.g. ../media/image7.png
        media_key = "ppt/" + target.replace("../", "")
        media_new[media_key] = png.read_bytes()

        def fix(m: re.Match) -> str:
            x, y = int(m.group(2)), int(m.group(4))
            cx, cy = int(m.group(6)), int(m.group(8))
            new_cy = int(round(cx * aspect))
            new_y = y + (cy - new_cy) // 2  # keep the old vertical midpoint
            return (
                f"{m.group(1)}{x}{m.group(3)}{new_y}{m.group(5)}"
                f"{cx}{m.group(7)}{new_cy}{m.group(9)}"
            )

        new_block, n = XFRM_RE.subn(fix, block, count=1)
        if n != 1:
            raise SystemExit(f"slide {slide_no}: could not rewrite picture geometry")
        slide_xml[sname] = xml.replace(block, new_block)
        print(f"slide {slide_no:2d}  {media_key:22s} <- {panel}")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            name = info.filename
            if name in slide_xml:
                data = slide_xml[name].encode("utf-8")
            elif name in media_new:
                data = media_new[name]
            else:
                data = zin.read(name)
            zout.writestr(info, data)
    zin.close()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(buf.getvalue())
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
