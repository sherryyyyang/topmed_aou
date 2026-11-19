#!/usr/bin/env Rscript
# Export constrained vs tolerant underlier tables (documentation / downstream use).
#
# Cite this script in methods or supplementary notes (e.g. point e).
#
# Usage:
#   cd .../pre_Watershed_scripts
#   DEFINITION_LABEL=z_-3 Rscript export_underlier_constrained_tolerant_tables.R
#
# Environment:
#   DEFINITION_LABEL  strict | z_-4 | z_-3 | z_-2 | haplo  (default: strict)
#   UNDERLIER_WD      directory with *_clean.rds, kept_subjects.rds, loeuf_lookup.rds
#
# Outputs (suffix follows DEFINITION_LABEL, e.g. _z_-3):
#   underlier_tallies_loeuf{SUFFIX}_constrained_tolerant_genes.tsv
#       Gene-level: LOEUF tier, rank, % subjects, min/mean z and logCPM across underliers.
#   underlier_tallies_loeuf{SUFFIX}_constrained_tolerant_subject_gene.tsv
#       One row per subject-gene underlier with individual expression_zscore, expression_logcpm.
#
# Inputs: same as plot_underlier_tallies_loeuf_constrained_tolerant.R
# Related figure script: plot_underlier_tallies_loeuf_constrained_tolerant.R

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(tidyr)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG2B

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_FOCUS <- c("Constrained", "Tolerant")

pick_wd <- function() {
    env <- Sys.getenv("UNDERLIER_WD", unset = "")
    cands <- unique(c(
        if (nzchar(env)) env else NA_character_,
        "/net/wonderland/home/sheyang/topmed11/pre_Watershed_scripts",
        "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
    ))
    cands <- cands[!is.na(cands)]
    for (d in cands) {
        if (file.exists(file.path(d, "_underlier_def.R")) &&
            file.exists(file.path(d, "kept_subjects.rds"))) {
            return(d)
        }
    }
    stop("Could not find data directory; set UNDERLIER_WD to folder with _underlier_def.R and kept_subjects.rds")
}

build_constrained_tolerant_tables <- function(wd = pick_wd()) {
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(wd)
    source("_underlier_def.R", local = TRUE)
    DEF <- resolve_underlier_def()

    cat("Underlier definition:", DEF$short, DEF$tag, "\n")
    cat("Working directory:", getwd(), "\n")

    cat("Loading kept_subjects.rds and expr_cpm_adj.rds ...\n")
    t0 <- Sys.time()
    kept_subjects <- read_rds(KEPT_SUBJECTS)
    n_subjects <- length(kept_subjects)
    expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
    all_eligible_genes <- rownames(expr_cpm_adj)
    rm(expr_cpm_adj)
    gc(verbose = FALSE)
    cat(
        "  kept subjects =", n_subjects, "; eligible genes =", length(all_eligible_genes),
        " (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s)\n"
    )

    cat("Loading", DEF$rds_underlier, "...\n")
    z_score_underlier <- read_rds(DEF$rds_underlier)
    loeuf_lookup <- read_rds(LOEUF_LOOKUP)

    under_per_gene <- z_score_underlier |>
        group_by(gene_id) |>
        tally() |>
        ungroup()

    gene_df <- tibble(gene_id = all_eligible_genes) |>
        mutate(gene_id_no_version = str_remove(gene_id, "\\.\\d+$")) |>
        left_join(under_per_gene, by = "gene_id") |>
        mutate(
            n             = replace_na(n, 0L),
            pct_underlier = n / n_subjects
        )

    loeuf_full <- loeuf_lookup |>
        distinct(gene_id, .keep_all = TRUE) |>
        select(gene_id, LOEUF_full = LOEUF, tier_full = LOEUF_tier, symbol_full = symbol)
    loeuf_nv <- loeuf_lookup |>
        distinct(gene_id_no_version, .keep_all = TRUE) |>
        select(gene_id_no_version, LOEUF_nv = LOEUF, tier_nv = LOEUF_tier, symbol_nv = symbol)

    gene_df <- gene_df |>
        left_join(loeuf_full, by = "gene_id") |>
        left_join(loeuf_nv, by = "gene_id_no_version") |>
        mutate(
            LOEUF      = coalesce(LOEUF_full, LOEUF_nv),
            symbol     = coalesce(symbol_full, symbol_nv),
            LOEUF_tier = factor(
                ifelse(is.na(LOEUF), "No LOEUF", as.character(coalesce(tier_full, tier_nv))),
                levels = TIER_LEVELS
            )
        ) |>
        select(gene_id, gene_id_no_version, symbol, LOEUF, LOEUF_tier, n, pct_underlier)

    prev_tbl <- gene_df |>
        group_by(LOEUF_tier) |>
        summarise(n_genes_in_tier = n(), .groups = "drop")

    tier_df <- gene_df |>
        filter(n > 0, LOEUF_tier %in% TIER_FOCUS) |>
        group_by(LOEUF_tier) |>
        arrange(desc(pct_underlier), .by_group = TRUE) |>
        mutate(rank_in_tier = row_number()) |>
        ungroup() |>
        left_join(prev_tbl, by = "LOEUF_tier") |>
        mutate(frac_tier_rank = rank_in_tier / n_genes_in_tier)

    tier_df$LOEUF_tier <- factor(tier_df$LOEUF_tier, levels = TIER_FOCUS)

    gene_z_summary <- z_score_underlier |>
        group_by(gene_id) |>
        summarise(
            n_underlier_pairs      = n(),
            expression_zscore_min  = min(expression_zscore, na.rm = TRUE),
            expression_zscore_mean = mean(expression_zscore, na.rm = TRUE),
            expression_logcpm_min  = min(expression_logcpm, na.rm = TRUE),
            expression_logcpm_mean = mean(expression_logcpm, na.rm = TRUE),
            .groups = "drop"
        )

    export_genes <- tier_df |>
        left_join(gene_z_summary, by = "gene_id") |>
        mutate(
            pct_subjects_with_underlier = 100 * pct_underlier,
            gene_label = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id_no_version),
            constrained_vs_tolerant = tolower(as.character(LOEUF_tier))
        ) |>
        select(
            gene_label,
            gene_id,
            gene_id_no_version,
            symbol,
            constrained_vs_tolerant,
            LOEUF_tier,
            LOEUF,
            n_subjects_with_underlier = n,
            pct_subjects_with_underlier,
            rank_in_tier,
            frac_tier_rank,
            n_underlier_pairs,
            expression_zscore_min,
            expression_zscore_mean,
            expression_logcpm_min,
            expression_logcpm_mean
        ) |>
        arrange(LOEUF_tier, rank_in_tier)

    gene_meta <- export_genes |>
        select(
            gene_label, gene_id, gene_id_no_version, symbol,
            constrained_vs_tolerant, LOEUF_tier, LOEUF,
            n_subjects_with_underlier, pct_subjects_with_underlier,
            rank_in_tier, frac_tier_rank
        )

    export_subject_gene <- z_score_underlier |>
        semi_join(tier_df |> select(gene_id), by = "gene_id") |>
        left_join(gene_meta, by = "gene_id") |>
        select(
            gene_label, gene_id, gene_id_no_version, symbol,
            constrained_vs_tolerant, LOEUF_tier, LOEUF,
            subjectid,
            expression_zscore, expression_logcpm,
            n_subjects_with_underlier, pct_subjects_with_underlier,
            rank_in_tier, frac_tier_rank
        ) |>
        arrange(constrained_vs_tolerant, rank_in_tier, gene_id, subjectid)

    list(
        DEF = DEF,
        wd = wd,
        n_subjects = n_subjects,
        gene_df = gene_df,
        tier_df = tier_df,
        export_genes = export_genes,
        export_subject_gene = export_subject_gene
    )
}

write_constrained_tolerant_tables <- function(tbl, wd = tbl$wd) {
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(wd)
    suffix <- tbl$DEF$suffix

    out_genes <- paste0("underlier_tallies_loeuf", suffix, "_constrained_tolerant_genes.tsv")
    out_subj  <- paste0("underlier_tallies_loeuf", suffix, "_constrained_tolerant_subject_gene.tsv")

    write_tsv(tbl$export_genes, out_genes)
    write_tsv(tbl$export_subject_gene, out_subj)

    cat("Wrote:", normalizePath(out_genes), "(", nrow(tbl$export_genes), "genes)\n")
    cat("Wrote:", normalizePath(out_subj), "(", nrow(tbl$export_subject_gene), "subject-gene rows)\n")

    invisible(list(genes = out_genes, subject_gene = out_subj))
}

if (sys.nframe() == 0L) {
    tbl <- build_constrained_tolerant_tables()
    write_constrained_tolerant_tables(tbl)
}
