ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
# Shared definition-resolution helper for the LOEUF plot scripts.
#
# Supports:
#   strict (legacy z<-4, unsuffixed outputs),
#   z_-4 / z_-3 / z_-2 (strict variants with explicit suffixes),
#   haplo (logCPM-only definition).

resolve_underlier_def <- function() {
    label <- Sys.getenv("DEFINITION_LABEL", unset = "strict")
    label <- tolower(label)

    valid_labels <- c("strict", "z_-4", "z_-3", "z_-2", "haplo")
    if (!label %in% valid_labels) {
        stop("DEFINITION_LABEL must be one of ",
             paste(valid_labels, collapse = ", "),
             " (got '", label, "')")
    }

    if (label %in% c("strict", "z_-4", "z_-3", "z_-2")) {
        # Keep strict legacy outputs unsuffixed when DEFINITION_LABEL=strict.
        suffix <- ifelse(label == "strict", "", paste0("_", label))
        z_txt <- ifelse(label == "strict", "-4", sub("^z_", "", label))

        list(
            label         = label,
            short         = paste0("strict z<", z_txt),
            long          = paste0("z-score < ", z_txt, "  AND  log2CPM_adj < gene_mean - 1"),
            tag           = paste0("[strict z<", z_txt, " + logCPM]"),
            z_cutoff      = as.numeric(z_txt),
            suffix        = suffix,
            rds_underlier = file.path(DATA_UNDERLIER, ifelse(label == "strict", "z_score_underlier_clean.rds",
                                   paste0("z_score_underlier_", label, "_clean.rds"))),
            rds_tally_g   = file.path(DATA_UNDERLIER, ifelse(label == "strict", "z_score_tally_clean.rds",
                                   paste0("z_score_tally_", label, "_clean.rds"))),
            rds_tally_s   = file.path(DATA_UNDERLIER, ifelse(label == "strict", "z_score_tally_subj_clean.rds",
                                   paste0("z_score_tally_subj_", label, "_clean.rds"))),
            tier_tsv      = file.path(FIG2A, ifelse(label == "strict", "underlier_loeuf_calibration_by_tier.tsv",
                                   paste0("underlier_loeuf_calibration_by_tier", suffix, ".tsv"))),
            decile_tsv    = file.path(FIG2A, ifelse(label == "strict", "underlier_loeuf_calibration_by_decile.tsv",
                                   paste0("underlier_loeuf_calibration_by_decile", suffix, ".tsv"))),
            callouts_tsv  = file.path(FIG2X, ifelse(label == "strict", "underlier_loeuf_callouts.tsv",
                                   paste0("underlier_loeuf_callouts", suffix, ".tsv")))
        )
    } else {
        list(
            label         = "haplo",
            short         = "haplo (per-gene z = -1/SD)",
            long          = "log2CPM_adj < gene_mean - 1   (per-gene z < -1/SD)",
            tag           = "[haplo, per-gene z = -1/SD]",
            suffix        = "_haplo",
            rds_underlier = file.path(DATA_UNDERLIER, "z_score_underlier_haplo_clean.rds"),
            rds_tally_g   = file.path(DATA_UNDERLIER, "z_score_tally_haplo_clean.rds"),
            rds_tally_s   = file.path(DATA_UNDERLIER, "z_score_tally_subj_haplo_clean.rds"),
            tier_tsv      = file.path(FIG2A, "underlier_loeuf_calibration_by_tier_haplo.tsv"),
            decile_tsv    = file.path(FIG2A, "underlier_loeuf_calibration_by_decile_haplo.tsv"),
            callouts_tsv  = file.path(FIG2X, "underlier_loeuf_callouts_haplo.tsv")
        )
    }
}
