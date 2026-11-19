#!/usr/bin/env python3
"""Re-sync Supplementary Figures 11–16 into supp_figs/ from known sources.

Canonical list: submission_package_nov_2026/supp_figs/CANONICAL_LATEST.md
Main Figures 1–4 stay in submission_package_nov_2026/.
S1–S10 are produced by the usual R assemble / supplement scripts.
"""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(
    "/net/topmed11/working/porchard/for-others/julie-and-esther/"
    "2025-11-17/sheyang/pre_Watershed_scripts"
)
DRAFT = Path("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026")
SUPP = DRAFT / "figures" / "supplementary"
ARCHIVE = DRAFT / "_package_archive" / "supp_figs_archive"
LATEX_FIG = DRAFT / "latex" / "figures"

# Prefer live sources; fall back to the cleanup archive if draft-root copies
# were already moved.
ARCHIVE_DRAFT = ARCHIVE / "pre_canonical_20261118" / "draft_root"

KEEPERS = [
    (
        11,
        "FigS11_hvg_allgenes_composition",
        [
            DRAFT / "Fig_hvg_allgenes_composition.png",
            ARCHIVE_DRAFT / "Fig_hvg_allgenes_composition.png",
        ],
        [
            DRAFT / "Fig_hvg_allgenes_composition.pdf",
            ARCHIVE_DRAFT / "Fig_hvg_allgenes_composition.pdf",
        ],
    ),
    (
        12,
        "FigS12_hvg_allgenes_umap",
        [
            DRAFT / "Fig_hvg_allgenes_umap.png",
            ARCHIVE_DRAFT / "Fig_hvg_allgenes_umap.png",
        ],
        [
            DRAFT / "Fig_hvg_allgenes_umap.pdf",
            ARCHIVE_DRAFT / "Fig_hvg_allgenes_umap.pdf",
        ],
    ),
    (
        13,
        "FigS13_coverage_tracks_clingen_underliers",
        [ROOT / "supplement/exploratory/coverage_tracks_clingen_underliers.png"],
        [ROOT / "supplement/exploratory/coverage_tracks_clingen_underliers.pdf"],
    ),
    (
        14,
        "FigS14_pc55_vs_pc2005_outlier_overlap",
        [
            DRAFT / "Supplementary_Figure_PC55_vs_PC2005_outlier_overlap.png",
            ARCHIVE_DRAFT / "Supplementary_Figure_PC55_vs_PC2005_outlier_overlap.png",
            SUPP / "Supplementary_Figure14.png",
        ],
        [
            DRAFT / "Supplementary_Figure_PC55_vs_PC2005_outlier_overlap.pdf",
            ARCHIVE_DRAFT / "Supplementary_Figure_PC55_vs_PC2005_outlier_overlap.pdf",
            SUPP / "Supplementary_Figure14.pdf",
        ],
    ),
    (
        15,
        "FigS15_topmed_vs_aou_prevalence_clinical_class",
        [
            DRAFT / "fig2_topmed_vs_aou_prevalence_clinical_class.png",
            ARCHIVE_DRAFT / "fig2_topmed_vs_aou_prevalence_clinical_class.png",
            SUPP / "archive/pre_canonical_20261118/draft_root/fig2_topmed_vs_aou_prevalence_clinical_class.png",
            SUPP / "archive/redundant_pass2_20261118/source_panels/fig2_topmed_vs_aou_prevalence_clinical_class.png",
            SUPP / "Supplementary_Figure15.png",
        ],
        [
            DRAFT / "fig2_topmed_vs_aou_prevalence_clinical_class.pdf",
            ARCHIVE_DRAFT / "fig2_topmed_vs_aou_prevalence_clinical_class.pdf",
            SUPP / "archive/pre_canonical_20261118/draft_root/fig2_topmed_vs_aou_prevalence_clinical_class.pdf",
            SUPP / "archive/redundant_pass2_20261118/source_panels/fig2_topmed_vs_aou_prevalence_clinical_class.pdf",
            SUPP / "Supplementary_Figure15.pdf",
        ],
    ),
    (
        16,
        "FigS16_constraint_sanity_allgenes",
        [
            DRAFT / "Fig3a_constraint_sanity_allgenes.png",
            ARCHIVE_DRAFT / "Fig3a_constraint_sanity_allgenes.png",
            SUPP / "archive/pre_canonical_20261118/draft_root/Fig3a_constraint_sanity_allgenes.png",
            SUPP / "archive/redundant_pass2_20261118/source_panels/Fig3a_constraint_sanity_allgenes.png",
            SUPP / "Supplementary_Figure16.png",
        ],
        [
            DRAFT / "Fig3a_constraint_sanity_allgenes.pdf",
            ARCHIVE_DRAFT / "Fig3a_constraint_sanity_allgenes.pdf",
            SUPP / "archive/pre_canonical_20261118/draft_root/Fig3a_constraint_sanity_allgenes.pdf",
            SUPP / "archive/redundant_pass2_20261118/source_panels/Fig3a_constraint_sanity_allgenes.pdf",
            SUPP / "Supplementary_Figure16.pdf",
        ],
    ),
]


def first_existing(paths: list[Path]) -> Path | None:
    for p in paths:
        if p.exists():
            return p
    return None


def main() -> None:
    SUPP.mkdir(parents=True, exist_ok=True)
    LATEX_FIG.mkdir(parents=True, exist_ok=True)
    for num, panel, pngs, pdfs in KEEPERS:
        png = first_existing(pngs)
        pdf = first_existing(pdfs)
        if png is None:
            raise SystemExit(f"missing PNG for S{num}: tried {pngs}")
        sfig = f"Supplementary_Figure{num}"
        # Single-panel figures: only write Supplementary_FigureN (no FigS* twin).
        shutil.copy2(png, SUPP / f"{sfig}.png")
        if pdf is not None:
            shutil.copy2(pdf, SUPP / f"{sfig}.pdf")
        link = LATEX_FIG / f"{sfig}.png"
        if link.exists() or link.is_symlink():
            link.unlink()
        link.symlink_to(Path("../../supp_figs") / f"{sfig}.png")
        print(f"S{num}: {png} -> {sfig}.png")
    print("done; see supp_figs/CANONICAL_LATEST.md")


if __name__ == "__main__":
    main()
