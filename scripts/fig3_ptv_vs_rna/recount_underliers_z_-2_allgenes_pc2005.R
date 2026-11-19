#!/usr/bin/env Rscript
# Per-gene kept-subject underlier counts at z < -2, same 2005-PC all-gene
# matrices as the z < -3 Fig 3a calls.
suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(tibble)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

Z_CACHE <- file.path(ROOT, "allgenes_pipeline_out", "expr_zscore_pc2005_allgenes.rds")
CPM_RDS <- file.path(ROOT, "allgenes_pipeline_out", "expr_cpm_adj.rds")
OUT_TSV <- file.path(
  ROOT, "allgenes_pipeline_out",
  "underlier_gene_counts_z_-2_allgenes_pc2005_kept.tsv"
)
Z_CUT <- -2
LOGCPM_DROP <- 1

if (file.exists(OUT_TSV)) {
  cat("exists ", OUT_TSV, "\n", sep = "")
  quit(save = "no", status = 0)
}

kept <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept)
stopifnot(n_kept == 5489L)
nv <- function(x) str_remove(x, "\\..*$")

cat("Recounting kept-subject underliers at z < -2 ...\n")
z_mat <- read_rds(Z_CACHE)
cpm <- as.matrix(read_rds(CPM_RDS))
if (!setequal(rownames(z_mat), rownames(cpm))) {
  stop("z and expr_cpm_adj gene IDs do not match")
}
if (!setequal(colnames(z_mat), colnames(cpm))) {
  stop("z and expr_cpm_adj sample IDs do not match")
}
cpm <- cpm[rownames(z_mat), colnames(z_mat), drop = FALSE]
gene_mean <- rowMeans(cpm)
kept_in_z <- intersect(kept, colnames(z_mat))
if (length(kept_in_z) != n_kept) {
  stop("kept_subjects missing from z grid: ", n_kept - length(kept_in_z))
}
n_hit <- rowSums(
  (z_mat[, kept_in_z, drop = FALSE] < Z_CUT) &
    (cpm[, kept_in_z, drop = FALSE] < (gene_mean - LOGCPM_DROP))
)
out <- tibble(
  gene_id = names(n_hit),
  gene_nv = nv(names(n_hit)),
  n_underlier_subjects = as.integer(n_hit)
)
write_tsv(out, OUT_TSV)
cat("wrote ", OUT_TSV, "  genes=", nrow(out),
    "  with >=1=", sum(out$n_underlier_subjects > 0), "\n", sep = "")
