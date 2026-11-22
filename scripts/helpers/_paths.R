# Shared paths. Source after ROOT is set, or let this file set it:
#   source(file.path(Sys.getenv("TOPMED_OUTLIERS_ROOT"), "scripts/helpers/_paths.R"))
#
# Required:
#   TOPMED_OUTLIERS_ROOT or DATA_ROOT  -> analysis directory (data/ and scripts/)
#   PC_SETTING                         -> expression PCs removed before z-scores
# Optional:
#   TOPMED_SHEYANG                     -> parent of ROOT (ptv_all.tsv, annotations)
#   TOPMED_SUBMISSION_PACKAGE          -> where assembled figures are written
#   TOPMED_PTV_TSV / TOPMED_ACMG_CSV   -> those inputs

if (!exists("ROOT") || !nzchar(ROOT)) {
  ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
  if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
  if (!nzchar(ROOT)) ROOT <- Sys.getenv("ANALYSIS_ROOT", unset = "")
  if (!nzchar(ROOT)) {
    stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
  }
}

SHEYANG <- Sys.getenv("TOPMED_SHEYANG", unset = "")
if (!nzchar(SHEYANG)) SHEYANG <- dirname(ROOT)

# --------------------------------------------------------------------------
# Expression-PC setting. Required. The selection rule is not in this repository.
# PC_SETTING 55 keeps the unsuffixed data/underliers directory. Any other count
# writes data/underliers_pc<N> and pc<N>_work/, so a second count does not
# overwrite the first.
# --------------------------------------------------------------------------
PC_SETTING <- Sys.getenv("PC_SETTING", unset = "")
if (!grepl("^[0-9]+$", PC_SETTING) || as.integer(PC_SETTING) < 1) {
  stop("Set PC_SETTING to a positive integer (expression PCs removed).")
}
PC_IS_DEFAULT <- identical(PC_SETTING, "55")
# Suffix appended to the directories that hold PC-dependent derived data and
# rendered output. Empty when PC_SETTING is 55.
PC_TAG <- if (PC_IS_DEFAULT) "" else paste0("_pc", PC_SETTING)
# Non-default PC runs write under pc<N>_work/.
PC_WORK <- if (PC_IS_DEFAULT) ROOT else file.path(ROOT, paste0("pc", PC_SETTING, "_work"))

DATA_CORE      <- file.path(ROOT, "data", "core")
DATA_EXPR      <- file.path(ROOT, "data", "expression")
# PC-dependent: underlier calls follow the z-score matrix.
DATA_UNDERLIER <- if (PC_IS_DEFAULT) {
  file.path(ROOT, "data", "underliers")
} else {
  file.path(ROOT, "data", paste0("underliers_pc", PC_SETTING))
}
DATA_PIPELINE  <- file.path(DATA_UNDERLIER, "pipeline_out")
DATA_ENRICH    <- file.path(ROOT, "data", "enrichment")
DATA_REF       <- file.path(ROOT, "data", "references")
AOU_RNA        <- file.path(ROOT, "call-RunRnaUnderlier", "rna_outputs")

FIG1A <- file.path(PC_WORK, "fig1_overview", "panel_a_study_design")
FIG1B <- file.path(PC_WORK, "fig1_overview", "panel_b_rank_curves")
FIG1C <- file.path(PC_WORK, "fig1_overview", "panel_c_topmed_vs_aou")
FIG2A <- file.path(PC_WORK, "fig2_constraint", "panel_a_loeuf_calibration")
FIG2B <- file.path(PC_WORK, "fig2_constraint", "panel_b_constrained_vs_tolerant")
FIG2C <- file.path(PC_WORK, "fig2_constraint", "panel_c_topmed_vs_aou_loeuf")
FIG2X <- file.path(PC_WORK, "fig2_constraint", "extra")
FIG3A <- file.path(PC_WORK, "fig3_ptv_vs_rna", "panel_a_gene_level")
FIG3B <- file.path(PC_WORK, "fig3_ptv_vs_rna", "panel_b_carrier_zscore")
FIG3C <- file.path(PC_WORK, "fig3_ptv_vs_rna", "panel_c_clingen_examples")
FIG3X <- file.path(PC_WORK, "fig3_ptv_vs_rna", "extra")
FIG4A <- file.path(PC_WORK, "fig4_acmg", "panel_a_rbm20")
FIG4B <- file.path(PC_WORK, "fig4_acmg", "panel_b_myh11")
FIG4L <- file.path(PC_WORK, "fig4_acmg", "rbm20_like")
FIG4D <- file.path(PC_WORK, "fig4_acmg", "data")
SUP_PC  <- file.path(PC_WORK, "supplement", "pc_enrichment")
SUP_PTV <- file.path(PC_WORK, "supplement", "ptv_consequences")
SUP_GN  <- file.path(PC_WORK, "supplement", "gnomad_incidence")
SUP_EX  <- file.path(PC_WORK, "supplement", "exploratory")
PAPER   <- file.path(PC_WORK, "paper")
PAPER_COMPOSITES <- file.path(PAPER, "composites")
FINAL_DRAFT <- if (PC_IS_DEFAULT) {
  file.path(ROOT, "final_paper_draft")
} else {
  draft <- Sys.getenv("TOPMED_SUBMISSION_PACKAGE", unset = "")
  if (!nzchar(draft)) draft <- ROOT
  draft
}

PTV_TSV <- Sys.getenv("TOPMED_PTV_TSV", unset = "")
if (!nzchar(PTV_TSV)) PTV_TSV <- file.path(SHEYANG, "ptv_all.tsv")
ACMG_CSV <- Sys.getenv("TOPMED_ACMG_CSV", unset = "")
if (!nzchar(ACMG_CSV)) {
  ACMG_CSV <- file.path(SHEYANG, "annotation", "variant_annotations", "acmg_sf_v3.2.csv")
}

KEPT_SUBJECTS <- file.path(DATA_CORE, "kept_subjects.rds")
EXPR_Z_JOIN   <- if (PC_IS_DEFAULT) {
  file.path(DATA_EXPR, "expr_z_join.rds")
} else {
  file.path(DATA_UNDERLIER, "expr_z_join.rds")
}
EXPR_CPM_ADJ  <- file.path(DATA_EXPR, "expr_cpm_adj.rds")
LOEUF_LOOKUP  <- file.path(DATA_CORE, "loeuf_lookup.rds")
GENCODE_SYM   <- file.path(DATA_CORE, "gencode_v30_gene_symbols.tsv")
METADATA      <- file.path(DATA_CORE, "metadata.filtered_with_sv.txt")
CLINGEN_HI    <- file.path(DATA_REF, "ClinGen_HI_Sufficient_evidence.csv")

# Long-form residualized z-scores (gene_id, subjectid, z), gene-major with
# subject varying fastest. Both files share byte-identical key columns.
Z_LONGFORM <- file.path(
  DATA_EXPR,
  sprintf("expr_scaled_genotype_pc_5_expression_pc_%s_longer_form.tsv", PC_SETTING)
)
