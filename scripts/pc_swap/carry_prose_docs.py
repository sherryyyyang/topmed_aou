#!/usr/bin/env python3
"""Carry the hand-authored prose documents into the PC_SETTING output folder.

Four documents in the published tree have no generator script -- they were
written by hand:

  FIGURE_SOURCES.md                        figure -> panel -> script map
  MOVED.txt                                folder orientation note
  supp_figs/Supplementary_Figures_methods.md     per-figure methods text
  supp_figs/Supplementary_Figures_provenance.md  per-figure input/DAG text

What they describe -- which script reads which file, what each panel shows,
how each statistic is computed -- is a property of the pipeline and is true at
any PC setting. What they *quote* -- "at 55 PCs, SV OR = 1.29 / 1.83 / 2.88",
"variance in PC1-55 is 67.6%" -- are 55-PC results that cannot be rewritten
mechanically without inventing numbers.

So the documents are carried over with two kinds of change and no others:

  1. Mechanical retargeting that is checkable against the filesystem:
     the "this folder" path, and the S6 table / Excel tab names, which really
     are PC-parameterized (S6_enrichment_or_<PC>pcs.tsv).
  2. A banner at the top recording that every numeric claim below it is a
     55-PC value, and pointing at the regenerated tables that hold the
     PC_SETTING values.

Every substitution actually applied is logged, and the retargeted paths are
asserted to exist, so the carryover cannot silently point at nothing.

Supplementary_Figures_methods.docx is a pandoc rendering of the .md and is
regenerated from the retargeted .md rather than copied.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
from _paths import FINAL_DRAFT, PC_IS_DEFAULT, PC_SETTING, ROOT  # noqa: E402

PUBLISHED = ROOT / "final_paper_draft_pc55"
PANDOC = os.environ.get("PANDOC", "pandoc")

DOCS = [
    "FIGURE_SOURCES.md",
    "MOVED.txt",
    "supp_figs/Supplementary_Figures_methods.md",
    "supp_figs/Supplementary_Figures_provenance.md",
]

MD_BANNER = """> **PROVENANCE — {pc}-EXPRESSION-PC REPLICATE.**
> This document was written for the published 55-expression-PC analysis and is
> carried into `final_paper_draft_pc{pc}/` because its description of the
> pipeline — which script reads which input, what each panel shows, how each
> statistic is computed — holds unchanged at {pc} PCs.
>
> **Every number quoted below is a 55-PC value and has not been recomputed.**
> The corresponding {pc}-PC values are in `README.md`, `supplemental_data/` and
> `supp_figs/tables/` in this same folder, all of which *were* regenerated.
>
> Paths of the form `final_paper_draft/...` have been retargeted to this
> folder, and the PC-parameterized S6 table and Excel tab names have been
> renamed to their {pc}-PC equivalents. No other edit was made.
>
> Regenerate with: `PC_SETTING={pc} ./scripts/run_pipeline_pc.sh`

"""

TXT_BANNER = """PROVENANCE -- {pc}-EXPRESSION-PC REPLICATE
This note was written for the published 55-expression-PC analysis and is
carried into final_paper_draft_pc{pc}/ unchanged apart from path retargeting.
Every number quoted below is a 55-PC value and has not been recomputed; the
{pc}-PC values are in README.md, supplemental_data/ and supp_figs/tables/.
Regenerate with: PC_SETTING={pc} ./scripts/run_pipeline_pc.sh

"""

def substitutions(pc: str) -> list[tuple[str, str, str]]:
    """(regex, replacement, human description)."""
    return [
        (r"final_paper_draft(?!_pc)", f"final_paper_draft_pc{pc}", "this-folder path"),
        (r"S6_enrichment_or_55pcs\.tsv", f"S6_enrichment_or_{pc}pcs.tsv", "S6 table name"),
        (r"S6_enrichment_OR_55PC", f"S6_enrichment_OR_{pc}PC", "S6 Excel tab name"),
    ]


def retarget(text: str, pc: str, label: str) -> str:
    for pattern, repl, what in substitutions(pc):
        text, n = re.subn(pattern, repl, text)
        if n:
            print(f"    {label}: {n:3d} x {what}")
    return text


def main() -> int:
    if PC_IS_DEFAULT:
        print("PC_SETTING=55: prose documents are the published originals, nothing to do")
        return 0

    pc = PC_SETTING
    # The retargeted names must actually exist in the tree we are writing into.
    must_exist = [
        FINAL_DRAFT,
        FINAL_DRAFT / "supp_figs" / "tables" / f"S6_enrichment_or_{pc}pcs.tsv",
    ]
    for p in must_exist:
        if not p.exists():
            raise SystemExit(f"retarget would point at a non-existent path: {p}")

    for rel in DOCS:
        src = PUBLISHED / rel
        dst = FINAL_DRAFT / rel
        if not src.exists():
            raise SystemExit(f"missing published document: {src}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        print(f"  {rel}")
        text = src.read_text(encoding="utf-8")
        text = retarget(text, pc, rel)
        banner = (MD_BANNER if rel.endswith(".md") else TXT_BANNER).format(pc=pc)
        dst.write_text(banner + text, encoding="utf-8")

    # Supplementary_Figures_methods.docx: pandoc rendering of the .md beside it.
    md = FINAL_DRAFT / "supp_figs" / "Supplementary_Figures_methods.md"
    docx = FINAL_DRAFT / "supp_figs" / "Supplementary_Figures_methods.docx"
    if shutil.which(PANDOC) is None:
        raise SystemExit(f"pandoc not found ({PANDOC}); cannot render {docx.name}")
    subprocess.run(
        [PANDOC, "--from", "markdown", "--to", "docx", "--output", str(docx), str(md)],
        check=True,
    )
    print(f"  rendered {docx.relative_to(FINAL_DRAFT)} via pandoc")
    print(f"carried {len(DOCS)} prose documents + 1 pandoc rendering into {FINAL_DRAFT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
