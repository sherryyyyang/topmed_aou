#!/usr/bin/env Rscript
# Per-gene kept-subject counts at z < -2 ONLY (no log2CPM drop), for the
# de Goede comparison where GTEx Table S6 is also Z-only.
suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(tibble)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

Z_CACHE <- file.path(ROOT, "allgenes_pipeline_out", "expr_zscore_pc2005_allgenes.rds")
OUT_TSV <- file.path(
  ROOT, "allgenes_pipeline_out",
  "underlier_gene_counts_z_-2_zonly_allgenes_pc2005_kept.tsv"
)
Z_CUT <- -2

if (file.exists(OUT_TSV)) {
  cat("exists ", OUT_TSV, "\n", sep = "")
  quit(save = "no", status = 0)
}

kept <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept)
stopifnot(n_kept == 5489L)
nv <- function(x) str_remove(x, "\\..*$")

cat("Recounting kept-subject z < -2 (no CPM drop) ...\n")
z_mat <- read_rds(Z_CACHE)
kept_in_z <- intersect(kept, colnames(z_mat))
if (length(kept_in_z) != n_kept) {
  stop("kept_subjects missing from z grid: ", n_kept - length(kept_in_z))
}
n_hit <- rowSums(z_mat[, kept_in_z, drop = FALSE] < Z_CUT)
out <- tibble(
  gene_id = names(n_hit),
  gene_nv = nv(names(n_hit)),
  n_underlier_subjects = as.integer(n_hit)
)
write_tsv(out, OUT_TSV)
cat("wrote ", OUT_TSV, "  genes=", nrow(out),
    "  with >=1=", sum(out$n_underlier_subjects > 0), "\n", sep = "")
