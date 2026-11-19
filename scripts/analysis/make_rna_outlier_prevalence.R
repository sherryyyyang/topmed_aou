#!/usr/bin/env Rscript
# Build the RNA-outlier prevalence-per-gene table (haplo + strict z<-3) from the
# already-computed expr_z_join.rds, restricted to the 5,489 kept_subjects cohort
# (same denominator as the released figures). Mirrors run_underlier_pipeline.R.
suppressPackageStartupMessages({ library(tidyverse) })

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- DATA_PIPELINE

OUT <- DATA_PIPELINE; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOGCPM_DROP <- 1

zj   <- read_rds(EXPR_Z_JOIN)
kept <- read_rds(KEPT_SUBJECTS)
sym  <- read_tsv(GENCODE_SYM, col_names = c("gene_nv", "symbol"),
                 show_col_types = FALSE)

# per-gene mean log2CPM is computed over the FULL expression matrix (all tested
# samples), then underliers are counted within the kept_subjects cohort -- this
# reproduces the released per-gene numbers exactly.
zj <- zj |> group_by(gene_id) |>
  mutate(gene_mean_logcpm = mean(expression_logcpm)) |> ungroup() |>
  mutate(pass_logcpm = expression_logcpm < gene_mean_logcpm - LOGCPM_DROP) |>
  filter(subjectid %in% kept)
n_tested <- n_distinct(zj$subjectid)
cat("kept subjects (denominator):", n_tested, "\n")

call_prev <- function(df, tag) {
  prev <- df |>
    distinct(gene_id, subjectid) |>
    count(gene_id, name = "n_underlier_subjects") |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$"),
           n_subjects_tested = n_tested,
           rna_outlier_prevalence = n_underlier_subjects / n_tested) |>
    left_join(sym, by = "gene_nv") |>
    dplyr::select(gene_id, gene_nv, symbol, n_underlier_subjects,
                  n_subjects_tested, rna_outlier_prevalence) |>
    arrange(desc(rna_outlier_prevalence))
  f <- file.path(OUT, paste0("rna_outlier_prevalence_per_gene_", tag, ".tsv"))
  write_tsv(prev, f)
  cat(sprintf("  %-6s -> %s  (%d genes)\n", tag, f, nrow(prev)))
  prev
}

haplo <- call_prev(filter(zj, pass_logcpm), "haplo")
z3    <- call_prev(filter(zj, pass_logcpm, expression_zscore < -3), "z_-3")

# sanity check vs released per-gene numbers
cat("\n[validate] RPS26 haplo n (expect 1455): ",
    filter(haplo, symbol == "RPS26")$n_underlier_subjects, "\n")
cat("[validate] PTPRM z<-3 n (expect 30): ",
    filter(z3, symbol == "PTPRM")$n_underlier_subjects, "\n")
