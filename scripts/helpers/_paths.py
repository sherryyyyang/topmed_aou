"""Shared paths for Python plot scripts after panel-based reorganization.

PC_SETTING mirrors scripts/helpers/_paths.R: unset or "55" reproduces the
published analysis exactly; any other positive integer (e.g. 505, 2005)
re-points the PC-dependent derived data and every output root at that replicate.

Path overrides (optional; cluster defaults preserved when unset):
  TOPMED_OUTLIERS_ROOT / ANALYSIS_ROOT
  TOPMED_SHEYANG
  TOPMED_SUBMISSION_PACKAGE
  TOPMED_PTV_TSV / TOPMED_ACMG_CSV
  PC_SETTING
"""
import os
from pathlib import Path

_DEFAULT_ROOT = Path(
    "/net/topmed11/working/porchard/for-others/julie-and-esther/"
    "2025-11-17/sheyang/pre_Watershed_scripts"
)

ROOT = Path(
    os.environ.get("TOPMED_OUTLIERS_ROOT")
    or os.environ.get("ANALYSIS_ROOT")
    or _DEFAULT_ROOT
)
SHEYANG = Path(os.environ.get("TOPMED_SHEYANG") or ROOT.parent)

PC_SETTING = os.environ.get("PC_SETTING", "55")
if not PC_SETTING.isdigit() or int(PC_SETTING) < 1:
    raise ValueError(f"PC_SETTING must be a positive integer (got {PC_SETTING!r})")
PC_IS_DEFAULT = PC_SETTING == "55"
PC_WORK = ROOT if PC_IS_DEFAULT else ROOT / f"pc{PC_SETTING}_work"

DATA_CORE = ROOT / "data" / "core"
DATA_EXPR = ROOT / "data" / "expression"
DATA_UNDERLIER = (
    ROOT / "data" / "underliers" if PC_IS_DEFAULT
    else ROOT / "data" / f"underliers_pc{PC_SETTING}"
)
DATA_PIPELINE = DATA_UNDERLIER / "pipeline_out"
DATA_ENRICH = ROOT / "data" / "enrichment"
DATA_REF = ROOT / "data" / "references"
AOU_RNA = ROOT / "call-RunRnaUnderlier" / "rna_outputs"

FIG1A = PC_WORK / "fig1_overview" / "panel_a_study_design"
FIG1B = PC_WORK / "fig1_overview" / "panel_b_rank_curves"
FIG1C = PC_WORK / "fig1_overview" / "panel_c_topmed_vs_aou"
FIG2A = PC_WORK / "fig2_constraint" / "panel_a_loeuf_calibration"
FIG2B = PC_WORK / "fig2_constraint" / "panel_b_constrained_vs_tolerant"
FIG2C = PC_WORK / "fig2_constraint" / "panel_c_topmed_vs_aou_loeuf"
FIG4A = PC_WORK / "fig4_acmg" / "panel_a_rbm20"
FIG4B = PC_WORK / "fig4_acmg" / "panel_b_myh11"
FIG4D = PC_WORK / "fig4_acmg" / "data"
SUP_GN = PC_WORK / "supplement" / "gnomad_incidence"
SUP_PC = PC_WORK / "supplement" / "pc_enrichment"
PAPER = PC_WORK / "paper"
PAPER_COMPOSITES = PAPER / "composites"
FINAL_DRAFT = (
    ROOT / "final_paper_draft" if PC_IS_DEFAULT
    else Path(
        os.environ.get("TOPMED_SUBMISSION_PACKAGE")
        or (SHEYANG / "submission_package_nov_2026")
    )
)

PTV_TSV = Path(os.environ.get("TOPMED_PTV_TSV") or (SHEYANG / "ptv_all.tsv"))
ACMG_CSV = Path(
    os.environ.get("TOPMED_ACMG_CSV")
    or (SHEYANG / "annotation" / "variant_annotations" / "acmg_sf_v3.2.csv")
)
KEPT_SUBJECTS = DATA_CORE / "kept_subjects.rds"
GENCODE_SYM = DATA_CORE / "gencode_v30_gene_symbols.tsv"
