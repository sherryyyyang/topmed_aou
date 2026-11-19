# Assemble final_paper_draft*/supplemental_data/ for the active PC_SETTING.
#
# Each Nature-style Source_Data_*.tsv is a rename of one panel table produced
# upstream in the pipeline. The 55-PC folder was assembled by hand and its
# producer was never committed; the mapping below was recovered by matching the
# published files against the panel tables by checksum (all 19 matched exactly),
# so running this at PC_SETTING=55 reproduces the published folder.
#
# `pc_dependent` records whether a table follows the z-score matrix. FALSE rows
# are identical at both PC settings by construction: they either come from
# inputs that carry no TOPMed z-scores at all (AoU release, gnomAD incidence,
# PTV consequence counts) or from the expression-PC sweep, which spans the whole
# PC grid and so is the same table whichever point on it the paper quotes.

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

suppressPackageStartupMessages(library(tibble))

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

MAP <- tribble(
  ~out,                                                       ~dir,    ~src,                                                              ~pc_dependent,
  "Source_Data_Fig1d_topmed_vs_aou_prevalence_z-3.tsv",       "FIG1C", "topmed_vs_new_z_-3_prevalence_per_gene.tsv",                       TRUE,
  "Source_Data_Fig1d_topmed_vs_aou_prevalence_haplo.tsv",     "FIG1C", "topmed_vs_new_haplo_prevalence_per_gene.tsv",                      TRUE,
  "Source_Data_Fig2d_topmed_vs_aou_loeuf.tsv",                "FIG1C", "topmed_vs_new_z_-3_prevalence_per_gene.tsv",                       TRUE,
  "Source_Data_Fig2a_loeuf_tier_wilson.tsv",                  "FIG2A", "underlier_loeuf_calibration_by_tier_z_-3.tsv",                     TRUE,
  "Source_Data_Fig2b_loeuf_decile_wilson.tsv",                "FIG2A", "underlier_loeuf_calibration_by_decile_z_-3.tsv",                   TRUE,
  "Source_Data_robustness_Fig2_tier_haplo.tsv",               "FIG2A", "underlier_loeuf_calibration_by_tier_haplo.tsv",                    TRUE,
  "Source_Data_robustness_Fig2_tier_z-2.tsv",                 "FIG2A", "underlier_loeuf_calibration_by_tier_z_-2.tsv",                     TRUE,
  # The published folder took this from the unsuffixed alias that
  # plot_loeuf_calibration.R leaves behind at its default cutoff. Reading the
  # explicit _z_-4 file instead is byte-identical at 55 PCs and cannot go stale.
  "Source_Data_robustness_Fig2_tier_z-4.tsv",                 "FIG2A", "underlier_loeuf_calibration_by_tier_z_-4.tsv",                     TRUE,
  "Source_Data_Fig2c_constrained_tolerant_genes.tsv",         "FIG2B", "underlier_tallies_loeuf_z_-3_constrained_tolerant_genes.tsv",      TRUE,
  "Source_Data_robustness_Fig2_wilson_vs_subject_bootstrap.tsv", "FIG2X", "underlier_loeuf_wilson_vs_subject_bootstrap.tsv",              TRUE,
  "Source_Data_Fig3a_ptv_vs_underlier_gene.tsv",              "FIG3X", "ptv_vs_underlier_gene_level_z_-3.tsv",                            TRUE,
  "Source_Data_Fig3b_ptv_carrier_zscore_summary.tsv",         "FIG3B", "ptv_vs_nonptv_zscore_summary.tsv",                                TRUE,
  "Source_Data_Fig3c_clingen_zscore_summary.tsv",             "FIG3C", "candidate_gene_zscore_summary_clingen_z_-3.tsv",                  TRUE,
  "Source_Data_Fig4_acmg_4dot.tsv",                           "FIG4D", "shared_rare_outlier_genes_4dot.tsv",                              TRUE,
  "Source_Data_S1ab_enrichment_or.tsv",                       "SUP_PC", "enrichment_npc_sv_vs_ptv_underlier.tsv",                         FALSE,
  "Source_Data_S1ab_enrichment_bootstrap.tsv",                "SUP_PC", "enrichment_peak_vs_gd_bootstrap.tsv",                            FALSE,
  "Source_Data_S1c_ptv_consequences.tsv",                     "SUP_PTV", "ptv_consequence_counts_qc_cohort.tsv",                          FALSE,
  "Source_Data_extra_cohort_vs_gnomad_incidence.tsv",         "SUP_GN", "cohort_vs_gnomad_incidence_per_gene.tsv",                         FALSE,
  "Source_Data_AoU_rna_outlier_prevalence_z-3.tsv",          "AOU_RNA", "rna_outlier_prevalence_per_gene_haplo_z_-3.tsv",                  FALSE
)

# The PC sweep and the AoU/gnomAD/PTV-consequence inputs are not written by the
# 2005-PC run, so resolve those from the published tree; everything else must
# come from this run's own work tree or the table would be stale.
resolve_dir <- function(key) {
  d <- get(key)
  if (key %in% c("SUP_PC", "SUP_PTV", "SUP_GN") && !PC_IS_DEFAULT && !dir.exists(d)) {
    d <- sub(paste0("^", PC_WORK), ROOT, d)
  }
  d
}

PUBLISHED <- file.path(ROOT, "final_paper_draft_pc55", "supplemental_data")
# SOURCEDATA_OUT lets the 55-PC round-trip check write to scratch. Resolving the
# default through FINAL_DRAFT means PC_SETTING=55 points at the published folder,
# so refuse that outright unless it is asked for by name.
OUT <- Sys.getenv("SOURCEDATA_OUT", unset = file.path(FINAL_DRAFT, "supplemental_data"))
if (normalizePath(OUT, mustWork = FALSE) == normalizePath(PUBLISHED, mustWork = FALSE) &&
    !identical(Sys.getenv("ALLOW_PC55_WRITE"), "1")) {
  stop("refusing to write into the published 55-PC supplemental_data ", PUBLISHED,
       "; set SOURCEDATA_OUT to a scratch path")
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Haplo-definition tables are PC-dependent in principle but come out identical
# at both settings in practice: the haplo call uses only the log2CPM arm, and
# expr_cpm_adj.rds is held fixed across the two runs (see README). They are
# therefore exempt from the "must differ from 55 PCs" check below.
MAP$expect_differs <- MAP$pc_dependent & !grepl("haplo", MAP$out)

MAP$src_path <- file.path(vapply(MAP$dir, resolve_dir, character(1)), MAP$src)
missing <- !file.exists(MAP$src_path)
if (any(missing)) {
  stop("missing upstream panel tables:\n  ", paste(MAP$src_path[missing], collapse = "\n  "))
}

# Staleness rail. The 2005-PC work tree is seeded by copying the 55-PC panel
# directories, so a table that no script rewrote still sits there looking
# plausible while holding 55-PC numbers. Detect that by content rather than by
# timestamp: a table built from a z-score cutoff cannot come out byte-identical
# to its 55-PC counterpart, so if it does, nothing rewrote it. Timestamps are
# not usable here because re-deriving the underliers moves their mtime without
# changing what is in them.
if (!PC_IS_DEFAULT) {
  same_as_pc55 <- function(p) {
    ref <- sub(paste0("^", PC_WORK), ROOT, p)
    file.exists(ref) && identical(tools::md5sum(p)[[1]], tools::md5sum(ref)[[1]])
  }
  check <- MAP$expect_differs
  stale <- check & vapply(MAP$src_path, same_as_pc55, logical(1))
  if (any(stale)) {
    stop("these z-dependent panel tables are byte-identical to their 55-PC ",
         "counterparts, so they are still seeded copies:\n  ",
         paste(MAP$src_path[stale], collapse = "\n  "),
         "\nrun the owning panel script standalone so it rewrites them")
  }
  say("staleness check OK: all ", sum(check),
      " z-dependent tables differ from their 55-PC counterparts")
}

for (i in seq_len(nrow(MAP))) {
  ok <- file.copy(MAP$src_path[i], file.path(OUT, MAP$out[i]), overwrite = TRUE)
  if (!ok) stop("copy failed for ", MAP$out[i])
}
say("wrote ", nrow(MAP), " source-data tables into ", OUT)

# Upstream data lineage: describes the raw BCF/BAM/parquet layers feeding the
# project and is identical whatever residualization is applied downstream.
prov_src <- file.path(PUBLISHED, "PROVENANCE.tsv")
if (file.exists(prov_src) &&
    normalizePath(OUT, mustWork = FALSE) != normalizePath(PUBLISHED, mustWork = FALSE)) {
  file.copy(prov_src, file.path(OUT, "PROVENANCE.tsv"), overwrite = TRUE)
  say("copied PROVENANCE.tsv (upstream lineage, PC-independent)")
}

manifest <- data.frame(
  file        = MAP$out,
  n_rows      = vapply(MAP$src_path, function(p) length(readLines(p, warn = FALSE)) - 1L, integer(1)),
  pc_setting  = ifelse(MAP$pc_dependent, PC_SETTING, "any"),
  provenance  = ifelse(MAP$pc_dependent, "RECOMPUTED", "PC_INDEPENDENT"),
  source      = sub(paste0("^", ROOT, "/"), "", MAP$src_path),
  row.names   = NULL
)
write.table(manifest, file.path(OUT, "MANIFEST.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
say("wrote MANIFEST.tsv (", sum(MAP$pc_dependent), " recomputed, ",
    sum(!MAP$pc_dependent), " PC-independent)")
