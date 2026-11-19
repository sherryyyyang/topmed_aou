#!/usr/bin/env Rscript
# Rebuild the per-gene TOPMed-vs-AoU prevalence tables that Fig. 1d, Fig. 2d and
# the Fig. 4 four-dot data all read out of FIG1C.
#
#   topmed_vs_new_z_-3_prevalence_per_gene.tsv
#   topmed_vs_new_haplo_prevalence_per_gene.tsv
#
# The script that originally produced these is gone from the tree (SOURCE.txt in
# panel_c_topmed_vs_aou records "original script missing (.pyc only)"), so this
# is a reimplementation. It is held honest by VALIDATE_AGAINST: under
# PC_SETTING=55 the rebuilt table is diffed against the published one and any
# mismatch is reported loudly.
#
# TOPMed side follows PC_SETTING; the AoU side is a fixed external release
# (n = 8,895) and is identical in both runs.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

N_TOPMED <- 5489L
N_AOU    <- 8895L

# Published tables live here and are read-only; the rebuilt copies go to OUT_DIR.
PUBLISHED_FIG1C <- file.path(ROOT, "fig1_overview", "panel_c_topmed_vs_aou")
OUT_DIR <- Sys.getenv("FIG1C_OUT", unset = FIG1C)
if (normalizePath(OUT_DIR, mustWork = FALSE) ==
    normalizePath(PUBLISHED_FIG1C, mustWork = FALSE)) {
  stop("refusing to write into the published 55-PC panel dir ", PUBLISHED_FIG1C,
       "; set FIG1C_OUT to a scratch path")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

loeuf <- read_rds(LOEUF_LOOKUP) |>
  distinct(gene_id_no_version, .keep_all = TRUE) |>
  transmute(gene_nv = gene_id_no_version, LOEUF, LOEUF_tier, symbol_loeuf = symbol)

build_one <- function(tag, topmed_rds, aou_tsv) {
  tm <- read_rds(file.path(DATA_UNDERLIER, topmed_rds)) |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
    group_by(gene_nv) |>
    summarise(n_topmed = sum(n), .groups = "drop")

  aou <- read_tsv(file.path(AOU_RNA, aou_tsv), show_col_types = FALSE) |>
    group_by(gene_nv) |>
    summarise(symbol_aou = dplyr::first(symbol),
              n_new = sum(n_underlier_subjects), .groups = "drop")

  out <- full_join(tm, aou, by = "gene_nv") |>
    mutate(n_topmed = replace_na(n_topmed, 0L), n_new = replace_na(n_new, 0L)) |>
    left_join(loeuf, by = "gene_nv") |>
    transmute(
      gene_nv,
      symbol = coalesce(na_if(symbol_aou, ""), na_if(symbol_loeuf, ""), gene_nv),
      LOEUF, LOEUF_tier,
      n_topmed, N_topmed = N_TOPMED, prev_topmed = n_topmed / N_TOPMED,
      n_new, N_new = N_AOU, prev_new = n_new / N_AOU,
      ratio_new_over_topmed = ifelse(prev_topmed > 0, prev_new / prev_topmed, NA_real_)) |>
    arrange(gene_nv)

  path <- file.path(OUT_DIR, sprintf("topmed_vs_new_%s_prevalence_per_gene.tsv", tag))
  write_tsv(out, path)
  msg("wrote ", path, "  (", nrow(out), " genes; ",
      sum(out$n_topmed > 0), " TOPMed>=1, ", sum(out$n_new > 0), " AoU>=1)")
  out
}

msg("PC_SETTING = ", PC_SETTING, "; TOPMed tallies from ", DATA_UNDERLIER)
z3 <- build_one("z_-3", "z_score_tally_z_-3_clean.rds", "rna_outlier_prevalence_per_gene_haplo_z_-3.tsv")
hp <- build_one("haplo", "z_score_tally_haplo_clean.rds", "rna_outlier_prevalence_per_gene_haplo.tsv")

# ------------------------------------------------------------------ validate
# Only meaningful for the published setting, where a reference table exists.
if (PC_IS_DEFAULT) {
  ref <- read_tsv(file.path(PUBLISHED_FIG1C, "topmed_vs_new_z_-3_prevalence_per_gene.tsv"),
                  show_col_types = FALSE) |> arrange(gene_nv)
  msg("validation vs published table: ", nrow(ref), " reference rows, ", nrow(z3), " rebuilt rows")
  extra_pub <- nrow(anti_join(ref, z3, by = "gene_nv"))
  extra_new <- nrow(anti_join(z3, ref, by = "gene_nv"))
  both <- inner_join(dplyr::select(ref, gene_nv, ref_tm = n_topmed, ref_new = n_new),
                     dplyr::select(z3, gene_nv, new_tm = n_topmed, new_new = n_new),
                     by = "gene_nv")
  bad_tm <- sum(both$ref_tm != both$new_tm); bad_new <- sum(both$ref_new != both$new_new)
  msg("  genes shared: ", nrow(both), "  n_topmed mismatches: ", bad_tm,
      "  n_new mismatches: ", bad_new)
  msg("  only in published: ", extra_pub, "  only in rebuilt: ", extra_new)
  if (bad_tm || bad_new || extra_pub || extra_new) {
    stop("rebuilt TOPMed-vs-AoU table does not match the published one")
  }
  msg("  ROUND-TRIP OK: reproduces the published prevalence table exactly")
}
msg("DONE")
