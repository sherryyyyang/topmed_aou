#!/usr/bin/env python3
"""Swap Nature-style composites into the Word draft and lowercase panel tags.

Caption letters and in-text Fig. refs are updated to match
scripts/assemble_nature_figures.R (1a–1d, 2a–2d, 3a–3c, 4a–4d, S1a–S1c).

PC_SETTING selects which tree is written:

  55    edit final_paper_draft_pc55/…docx in place. This is the published
        manuscript, so it additionally requires ALLOW_PC55_WRITE=1.
  2005  take the published docx as a read-only template, swap in the 2005-PC
        composites, and write the result to submission_package_nov_2026/. The
        body prose is inherited unchanged, so a provenance banner is inserted
        at the top of the document saying so.

Either way the five embedded composites are also written to
source_panels/as_in_docx/, which is what that folder is defined to hold.
"""
from __future__ import annotations

import hashlib
import io
import os
import re
import shutil
import sys
import zipfile
from pathlib import Path

from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import FINAL_DRAFT, PC_IS_DEFAULT, PC_SETTING, ROOT  # noqa: E402

DOCX_NAME = "Outliers_in_TOPMed_paper_draft.docx"
PUBLISHED = ROOT / "final_paper_draft_pc55"
# The 55-PC manuscript is always the text template; only the images change.
TEMPLATE = PUBLISHED / DOCX_NAME
DOCX = FINAL_DRAFT / DOCX_NAME

def _first_existing(*candidates: Path) -> Path:
    """Prefer the reorganized figures/{main,supplementary}/ layout."""
    for p in candidates:
        if p.exists():
            return p
    return candidates[0]


def _main_fig(name: str) -> Path:
    return _first_existing(
        FINAL_DRAFT / "figures" / "main" / name,
        FINAL_DRAFT / name,
    )


def _supp_fig(name: str) -> Path:
    return _first_existing(
        FINAL_DRAFT / "figures" / "supplementary" / name,
        FINAL_DRAFT / "supp_figs" / name,
        FINAL_DRAFT / name,
    )


# Main composites live under figures/main/; S1 under figures/supplementary/
# (with legacy package-root / supp_figs/ fallbacks).
MEDIA = {
    "word/media/image1.png": _main_fig("Figure1.png"),
    "word/media/image2.png": _main_fig("Figure2.png"),
    "word/media/image3.png": _main_fig("Figure3.png"),
    "word/media/image4.png": _main_fig("Figure4.png"),
    "word/media/image5.png": _supp_fig("Supplementary_Figure1.png"),
}

BANNER = (
    "PROVENANCE — {pc}-EXPRESSION-PC REPLICATE. The figures embedded in this "
    "document (Figures 1–4 and Supplementary Figure 1) were regenerated at "
    "{pc} residualized expression PCs by PC_SETTING={pc} "
    "scripts/run_pipeline_pc.sh. The body text, captions and all numbers "
    "quoted in them are inherited verbatim from the published 55-PC "
    "manuscript and have NOT been recomputed — do not read them as {pc}-PC "
    "results. The authoritative {pc}-PC values are in README.md, "
    "supplemental_data/ and supp_figs/tables/ in this same folder."
)


def banner_paragraph(pc: str) -> str:
    """A single bold red paragraph, inserted as the first block of the body."""
    text = BANNER.format(pc=pc)
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return (
        '<w:p><w:pPr><w:pBdr>'
        '<w:top w:val="single" w:sz="12" w:space="1" w:color="C00000"/>'
        '<w:left w:val="single" w:sz="12" w:space="4" w:color="C00000"/>'
        '<w:bottom w:val="single" w:sz="12" w:space="1" w:color="C00000"/>'
        '<w:right w:val="single" w:sz="12" w:space="4" w:color="C00000"/>'
        '</w:pBdr><w:spacing w:after="240"/></w:pPr>'
        '<w:r><w:rPr><w:b/><w:color w:val="C00000"/><w:sz w:val="18"/>'
        '<w:szCs w:val="18"/></w:rPr>'
        f'<w:t xml:space="preserve">{text}</w:t></w:r></w:p>'
    )


def insert_banner(xml: str, pc: str) -> str:
    if "PROVENANCE — " in xml:  # already carries one; replace rather than stack
        xml = re.sub(r"<w:p>(?:(?!</w:p>).)*?PROVENANCE — .*?</w:p>", "", xml, count=1, flags=re.S)
    m = re.search(r"<w:body[^>]*>", xml)
    if not m:
        raise SystemExit("no <w:body> in document.xml")
    return xml[: m.end()] + banner_paragraph(pc) + xml[m.end() :]

# Longest keys first.
BODY = [
    ("Supplementary Fig. 1A", "Supplementary Fig. 1a,b"),
    ("Supplementary Fig. 1B", "Supplementary Fig. 1c"),
    ("Fig. 1A", "Fig. 1a"),
    ("Fig. 1B", "Fig. 1b,c"),
    ("Fig. 1C", "Fig. 1d"),
    ("Fig. 2A", "Fig. 2a,b"),
    ("Fig. 2B", "Fig. 2c"),
    ("Fig. 2C", "Fig. 2d"),
    ("Fig. 3A", "Fig. 3a"),
    ("Fig. 3B", "Fig. 3b"),
    ("Fig. 3C", "Fig. 3c"),
    ("Fig. 4A", "Fig. 4a,b"),
    ("Fig. 4B", "Fig. 4c,d"),
]

CAPTION_FRAGMENTS = [
    ("(A) Study design", "(a) Study design"),
    (
        '<w:t>burden</w:t></w:r><w:r w:rsidRPr="00ED7519"><w:rPr><w:b/><w:sz w:val="16"/><w:szCs w:val="16"/></w:rPr><w:t xml:space="preserve">. Left, genes ranked',
        '<w:t>burden</w:t></w:r><w:r w:rsidRPr="00ED7519"><w:rPr><w:b/><w:sz w:val="16"/><w:szCs w:val="16"/></w:rPr><w:t xml:space="preserve"> across genes and subjects. TOPMed n = 5,489; AoU n = 8,895. Left, genes ranked',
    ),
    (
        '<w:t xml:space="preserve"> across genes and subjects. TOPMed n = 5,489; AoU n = 8,895 (haplo + z-score &lt; −3). Left, genes ranked',
        '<w:t xml:space="preserve"> across genes and subjects. TOPMed n = 5,489; AoU n = 8,895. Left, genes ranked',
    ),
    ("Right, TOPMed subjects", "(c) TOPMed subjects"),
    ("(C) Per-gene", "(d) Per-gene"),
    ("Right, the same gradient", "(b) The same gradient"),
    ("(B) Within-tier", "(c) Within-tier"),
    (
        "versus tolerant (LOEUF ≥ 1.0) genes;",
        "versus tolerant (LOEUF ≥ 1.0) genes (2,939 genes with ≥1 underlier);",
    ),
    (
        "; n = 3,319 genes; Spearman = 0.50), on a −log10 scale",
        "; n = 3,319 genes; Spearman = 0.50; Pearson(−log10) = 0.52), on a −log10 scale",
    ),
    (
        " (n = 3,319 genes; Spearman = 0.50), colored by LOEUF",
        " (n = 3,319 genes; Spearman = 0.50; Pearson(−log10) = 0.52), colored by LOEUF",
    ),
    (
        "ranked by PTV-carrier prevalence, with ",
        "ranked by PTV-carrier prevalence (n = 340 genes), with ",
    ),
    ("(C) TOPMed versus", "(d) TOPMed versus"),
    (" marked (left), and the four frequencies for ", ". (b) Four frequencies for "),
    (
        "(A) Enrichment of underliers for rare SVs within 100 kb (AF ≤ 1%) and for rare PTVs in the gene (AF ≤ 0.5%)",
        "(a) Enrichment of underliers for rare SVs within 100 kb (AF ≤ 1%)",
    ),
    ("(B) VEP", "(c) VEP"),
    # Stripped on-figure titles/subtitles/stats → Word legends (Ferraro-style).
    # Fig 1: inner titles Genes/Subjects; 4,671 TOPMed underlier subjects;
    # ≥2-underlier gene set; replicate / grey-red / −log10 meaning from the
    # unused R plot annotation.
    (". (b) Underlier ", ". (b) Genes. Underlier "),
    (
        ". (c) TOPMed subjects ranked by the number of underlier genes (median 5, with a minority above 1",
        ". (c) Subjects. TOPMed-only: 4,671 subjects with ≥1 underlier, ranked by the number of underlier genes (median 5, with a minority above 1",
    ),
    (
        "; n = 3,319 genes; Spearman = 0.50; Pearson(−log10) = 0.52), on a −log10 scale. One marker per distinct (nTOPMed, nAoU) count pair; radius is the number of genes stacked at that point. The boxed top-right corner is the 624 genes with −log10 prevalence ≥ 3 in both cohorts.",
        "; n = 3,319 genes with ≥2 underliers in both cohorts; Spearman = 0.50; Pearson(−log10) = 0.52), on a −log10 scale (higher −log10 = rarer). RNA-outlier prevalence replicates across TOPMed and AoU. One marker per distinct (nTOPMed, nAoU) count pair; radius is the number of genes stacked at that point. Grey points: all genes with ≥2 underliers in both cohorts (no LOEUF coloring). Red labels: genes with high prevalence in both, plus RBM20-like / ACMG examples. The boxed top-right corner is the 624 genes with −log10 prevalence ≥ 3 in both cohorts.",
    ),
    # Fig 2: inner titles; D1–D10; ribbon CI; LOEUF-tracking claim.
    (
        "Figure 2. Underlier prevalence is calibrated by gene constraint. (a",
        "Figure 2. Underlier prevalence is calibrated by gene constraint: observed RNA-underlier prevalence tracks LOEUF. (a",
    ),
    (") Left, prevalence by constraint. Observed", ") By LOEUF tier. Observed"),
    (
        "genes. (b) The same gradient across LOEUF deciles of approximately equal gene count, confirming",
        "genes. (b) By LOEUF decile. The same gradient across LOEUF deciles of approximately equal gene count (D1 = most constrained … D10 = most tolerant; ribbon = Wilson 95% CI), confirming",
    ),
    (". (c) Within-tier rank curves", ". (c) Constrained vs tolerant. Within-tier rank curves"),
    # Fig 3: full gene-level stats block; fold enrichments; 3c strip/rugs.
    (
        " + log2CPM). Each point is an expression-tested gene; cohort n = 5,489. Spearman ρ = 0.148 over genes with ≥ 1 PTV carrier; Pearson over all genes ≈ 0. Dashed line: y = x. (b) Expression ",
        " + log2CPM). Each point is an expression-tested gene, median TPM &gt; 1 (n = 3,531 with ≥1 PTV carrier or ≥1 underlier; cohort n = 5,489). Spearman ρ = 0.148 among genes with ≥1 PTV carrier (n = 2,078; p = 1.1 × 10⁻¹¹); Pearson over all genes = 0.0076 (~0, no monotonic trend). Underlier-only genes (0 PTV carriers): 1,453; PTV-only: 557; both: 1,521. Genes high on y with near-zero x carry expression-outlier signal that PTV genotypes do not reveal. Dashed line: equal frequency (y = x). (b) Expression ",
    ),
    (
        "s of PTV carriers versus non-carriers of the same gene, at AF &lt; 0.5% (top) and AF &lt; 0.1% (bottom). Carriers are shifted negative and enriched in the low tail, more strongly at the lower allele frequency. (c) Per-subject expression ",
        "s of PTV carriers versus non-carriers of the same gene, at AF &lt; 0.1% (enrichments 5×/14×/37× at z &lt; −2/−3/−4). Shared legend: non-carrier versus PTV-carrier. Carriers are shifted negative and enriched in the low tail. (c) Per-subject expression ",
    ),
    (
        "at AF &lt; 0.5% (top; enrichments 4×/12×/31× at z &lt; −2/−3/−4) and AF &lt; 0.1% (bottom; 5×/14×/37×). Shared legend: non-carrier versus PTV-carrier. Carriers are shifted negative and enriched in the low tail, more strongly at the lower allele frequency.",
        "at AF &lt; 0.1% (enrichments 5×/14×/37× at z &lt; −2/−3/−4). Shared legend: non-carrier versus PTV-carrier. Carriers are shifted negative and enriched in the low tail.",
    ),
    (
        " threshold is marked; underliers and PTV carriers are indicated on each histogram.",
        " threshold is marked. Strip: gene · n PTV · n underliers. Purple rugs = PTV carriers; blue rugs = underliers + the z = −3 line.",
    ),
    # Fig 4: dropped zero-PTV count; four-series / loess / log10 / pink arrow.
    (
        "ranked by PTV-carrier prevalence (n = 340 genes), with ",
        "ranked by PTV-carrier prevalence (n = 340 genes; 284 zero-PTV genes dropped). Four series: TOPMed RNA, AoU RNA, TOPMed PTV, and gnomAD AF. Trend lines: loess for RNA and gnomAD; PTV as a line. y-axis is log10 frequency from ~1 in 1,000,000 to 1 in 100. Pink arrow and gene name mark ",
    ),
    # Supp Fig 1: inner titles Rare SVs / Rare PTVs / PTV consequences.
    (
        " (a) Enrichment of underliers for rare SVs within 100 kb (AF ≤ 1%) as a function of the number of expression PCs residualized, at ",
        " (a) Rare SVs. Enrichment of underliers for rare SVs within 100 kb (AF ≤ 1%) as a function of the number of expression PCs residualized, at ",
    ),
    (
        ". (b) The same PC grid for LoF/PTV in the gene (AF ≤ 0.5%). This sweep is the context for choosing 55 expression PCs when calling the outliers used in Figure 1. (c) VEP consequence composition of the rare PTV call set in the 5,489-subject QC cohort.",
        ". (b) Rare PTVs. The same PC grid for LoF/PTV in the gene (AF ≤ 0.5%). This sweep is the context for choosing 55 expression PCs when calling the outliers used in Figure 1. (c) PTV consequences. VEP consequence composition of the rare PTV call set in the 5,489-subject QC cohort.",
    ),
    (" (A)", " (a)"),
    (" (B)", " (b)"),
    (" (C)", " (c)"),
    ("(A)", "(a)"),
    ("(B)", "(b)"),
    ("(C)", "(c)"),
]


def scale_extents(xml: str, images: dict[str, Path], rels: str) -> str:
    rid_to_media = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="(media/[^"]+)"', rels))
    rid_to_media.update(dict(re.findall(r'Target="(media/[^"]+)"[^>]*Id="(rId\d+)"', rels)))

    def repl(m: re.Match) -> str:
        block = m.group(0)
        emb = re.search(r'r:embed="(rId\d+)"', block)
        ext = re.search(r'(<wp:extent[^>]*cx=")(\d+)("[^>]*cy=")(\d+)(")', block)
        if not emb or not ext:
            return block
        media = rid_to_media.get(emb.group(1), "")
        key = "word/" + media if media.startswith("media/") else media
        png = images.get(key)
        if png is None or not png.exists():
            return block
        im = Image.open(png)
        # Fig 2 is a wide a/b/c | d layout; display it wider than the old 5.1".
        target_in = {
            "word/media/image2.png": 7.2,
        }
        cx = (
            int(round(target_in[key] * 914400))
            if key in target_in
            else int(ext.group(2))
        )
        cy = int(round(cx * im.height / im.width))
        return block[:ext.start()] + f'{ext.group(1)}{cx}{ext.group(3)}{cy}{ext.group(5)}' + block[ext.end():]

    return re.sub(r"<w:drawing>.*?</w:drawing>", repl, xml, flags=re.S)


def main() -> int:
    if PC_IS_DEFAULT and os.environ.get("ALLOW_PC55_WRITE") != "1":
        raise SystemExit(
            "refusing to rewrite the published 55-PC manuscript; set "
            "ALLOW_PC55_WRITE=1 if that is really intended"
        )
    for p in MEDIA.values():
        if not p.exists():
            raise SystemExit(f"missing composite: {p}")
    if not TEMPLATE.exists():
        raise SystemExit(f"missing manuscript template: {TEMPLATE}")

    DOCX.parent.mkdir(parents=True, exist_ok=True)
    # The .bak is the document as it stood before this run swapped figures in.
    # At 2005 PCs that is the published 55-PC manuscript.
    bak = DOCX.with_suffix(".docx.bak_nature_panels")
    if not bak.exists():
        shutil.copy2(TEMPLATE, bak)
        print("backup", bak)

    buf = io.BytesIO()
    with zipfile.ZipFile(TEMPLATE, "r") as zin, zipfile.ZipFile(buf, "w") as zout:
        for info in zin.infolist():
            data = zin.read(info.filename)
            if info.filename == "word/document.xml":
                xml = data.decode("utf-8")
                rels = zin.read("word/_rels/document.xml.rels").decode("utf-8")
                for old, new in BODY + CAPTION_FRAGMENTS:
                    n = xml.count(old)
                    if n:
                        print(f"replace x{n}: {old[:72]!r}")
                    xml = xml.replace(old, new)
                xml = scale_extents(xml, MEDIA, rels)
                if not PC_IS_DEFAULT:
                    xml = insert_banner(xml, PC_SETTING)
                data = xml.encode("utf-8")
            elif info.filename in MEDIA:
                data = MEDIA[info.filename].read_bytes()
            zout.writestr(info, data)
    DOCX.write_bytes(buf.getvalue())
    print("wrote", DOCX)

    # source_panels/as_in_docx/ is by definition the composites inside the docx.
    as_in = FINAL_DRAFT / "source_panels" / "as_in_docx"
    as_in.mkdir(parents=True, exist_ok=True)
    for dest, src in MEDIA.items():
        shutil.copyfile(src, as_in / src.name)
        print(f"  {dest} <- {src.name} md5 {hashlib.md5(src.read_bytes()).hexdigest()[:8]}")
    print("as_in_docx refreshed:", as_in)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
