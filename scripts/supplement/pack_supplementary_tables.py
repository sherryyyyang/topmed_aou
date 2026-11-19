#!/usr/bin/env python3
"""Pack supplemental TSVs into one Excel workbook (sheet names <= 31 chars)."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(
    "/net/topmed11/working/porchard/for-others/julie-and-esther/"
    "2025-11-17/sheyang/pre_Watershed_scripts"
)
sys.path.insert(0, str(ROOT / "scripts" / "helpers"))
from _paths import FINAL_DRAFT, PC_SETTING  # noqa: E402

TSV = FINAL_DRAFT / "supp_figs/tables"
XLSX = FINAL_DRAFT / "supp_figs/Supplementary_Tables.xlsx"

# (sheet name <=31, filename under TSV, max rows to embed)
SHEETS = [
    ("S1_enrichment_OR", "Source_Data_S1ab_enrichment_or.tsv", None),
    ("S1_enrichment_bootstrap", "Source_Data_S1ab_enrichment_bootstrap.tsv", None),
    ("S1_PTV_consequences", "Source_Data_S1c_ptv_consequences.tsv", None),
    ("S2_QC_tests", "S2_qc_statistical_tests.tsv", None),
    ("S2_connectivity_Z", "S2_connectivity_z.tsv", None),
    ("S2_extreme_counts", "S2_subject_extreme_outlier_counts.tsv", None),
    ("S2_PCA_scree", "S2_pca_scree_variance.tsv", None),
    ("S3_UMAP_coords", "S3_umap_coords_before_after.tsv", None),
    ("S4_count_correlation", "S4_topmed_vs_aou_count_correlation.tsv", None),
    ("S4_counts_long", "S4_topmed_vs_aou_counts_long.tsv", None),
    ("S_gene_list_cutoffs", "S_gene_list_all_cutoffs.tsv", None),
    ("S5_cutoff_landscape", "S5_cutoff_landscape_summary.tsv", None),
    ("S5_LOEUF_by_cutoff", "S5_loeuf_tier_by_cutoff.tsv", None),
    (f"S6_enrichment_OR_{PC_SETTING}PC", f"S6_enrichment_or_{PC_SETTING}pcs.tsv", None),
    ("S7_overlap_stats", "S7_gene_overlap_and_stats.tsv", None),
    ("S7_NNLS_proportions", "S7_celltype_nnls_proportions.tsv", None),
    ("S7_PC_celltype_r", "S7_pc_celltype_pearson_r.tsv", None),
    ("S7_PC_multivar_R2", "S7_pc_multivar_r2.tsv", None),
    ("S7_cum_composition_R2", "S7_cumulative_composition_r2.tsv", None),
    ("S7_SC_cell_counts", "S7_sc_cell_counts.tsv", None),
    ("S8_gene_z_celltype_r", "S8_gene_z_vs_celltype_pearson.tsv", None),
    ("S8_gene_z_celltype_sum", "S8_gene_z_vs_celltype_summary.tsv", None),
    ("S8_gene_cpm_celltype_r", "S8_gene_log2cpm_vs_celltype_pearson.tsv", None),
    ("S8_z_vs_uncorrected_r", "S8_gene_z_vs_uncorrected_r.tsv", None),
    ("S_participants", "S_participant_characteristics.tsv", None),
    ("S2_Wilson_vs_bootstrap", "underlier_loeuf_wilson_vs_subject_bootstrap.tsv", None),
    ("S_gnomAD_incidence", "cohort_vs_gnomad_incidence_per_gene.tsv", None),
    ("MANIFEST", "MANIFEST.tsv", None),
]

ON_DISK_ONLY = [
    "S_underlier_pairs_topmed_haplo.tsv.gz",
    "S_underlier_pairs_topmed_z-2.tsv.gz",
    "S_underlier_pairs_topmed_z-3.tsv.gz",
    "S_underlier_pairs_topmed_z-4.tsv.gz",
    "S_underlier_pairs_aou_haplo.tsv.gz",
    "S_underlier_pairs_aou_z_-2.tsv.gz",
    "S_underlier_pairs_aou_z_-3.tsv.gz",
    "S_underlier_pairs_aou_z_-4.tsv.gz",
    "S_gene_underlier_counts_by_cutoff.tsv",
]


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False, comment="#")


def main() -> None:
    notes = [
        {
            "item": "Ethnicity / race",
            "status": "NOT AVAILABLE",
            "detail": (
                "Neither data/core/metadata.filtered_with_sv.txt nor "
                "2025-06-11/data/metadata.filtered.txt has race/ethnicity/ancestry. "
                "Do not invent it."
            ),
        },
        {
            "item": "AoU participant table",
            "status": "NOT ON DISK",
            "detail": "AoU is gene-level RNA prevalence only (n=8895). No sex/age/ethnicity file.",
        },
        {
            "item": "Age missingness",
            "status": "WHI ages unrecorded",
            "detail": "All 1,215 WHI subjects in the QC cohort lack AGE.",
        },
        {
            "item": "Pair-level underlier lists",
            "status": "TSV.GZ ON DISK ONLY",
            "detail": (
                "Haplo / z-cutoff subject-gene pair lists are too large for Excel. "
                "See supplemental_tables/S_underlier_pairs_*.tsv.gz. "
                "Gene-level counts are in S_gene_list_cutoffs."
            ),
        },
        {
            "item": "S1 vs S6",
            "status": "OK",
            "detail": (
                "S1 is OR vs n expression PCs. S6 is Fisher OR + 95% CI at the "
                "paper setting (55 PCs; z < -2/-3/-4; SV vs PTV)."
            ),
        },
        {
            "item": "S1 overlay / S7 samples",
            "status": "OK",
            "detail": (
                "Cell-type overlay on S1a/b and S7 use n=5560 samples in "
                "expr_pca.rds ∩ expr_cpm.rds (same universe as the S1 enrichment "
                "sweep). Paper QC kept n=5489 is flagged on S7_NNLS_proportions."
            ),
        },
        {
            "item": "S8 / S1c gene-z vs composition",
            "status": "OK",
            "detail": (
                "If residualization worked, gene-z vs NNLS composition |r| "
                "should be small (near 0) compared with uncorrected log2CPM. "
                "S1 panel c is the residualized-z Pearson r distribution; "
                "S8 adds the top-gene heatmap and uncorrected vs residualized |r|."
            ),
        },
    ]
    present = []
    for fn in ON_DISK_ONLY:
        p = TSV / fn
        present.append(
            {
                "file": fn,
                "present": p.exists(),
                "bytes": p.stat().st_size if p.exists() else 0,
                "note": "Full pair/gene table kept as TSV, not embedded in Excel",
            }
        )

    with pd.ExcelWriter(XLSX, engine="openpyxl") as writer:
        pd.DataFrame(notes).to_excel(writer, sheet_name="NOTES", index=False)
        pd.DataFrame(present).to_excel(writer, sheet_name="ON_DISK_NOT_IN_XLSX", index=False)
        for sheet, fn, max_rows in SHEETS:
            path = TSV / fn
            if not path.exists():
                pd.DataFrame({"missing": [fn]}).to_excel(writer, sheet_name=sheet, index=False)
                continue
            df = load_tsv(path)
            if max_rows is not None and len(df) > max_rows:
                df = df.head(max_rows)
            df.to_excel(writer, sheet_name=sheet[:31], index=False)
    print(f"wrote {XLSX}")


if __name__ == "__main__":
    main()
