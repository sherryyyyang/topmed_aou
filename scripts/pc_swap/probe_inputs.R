#!/usr/bin/env Rscript
# One-off reconnaissance for the 2005-PC replicate build.
#
# Two questions decide how deep the PC swap can go:
#   1. Does expr_eigenvectors_df.rds actually carry >= 2005 expression PCs?
#   2. Is expr_cpm_adj.rds reproducible as removeBatchEffect(expr_cpm, 55 PCs)?
# If (2) holds, the log2CPM arm of the underlier definition can be re-derived at
# 2005 PCs as well, instead of being frozen at its 55-PC value.
suppressPackageStartupMessages({ library(readr); library(limma) })

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

msg("reading expr_eigenvectors_df.rds ...")
ev <- read_rds(file.path(DATA_EXPR, "expr_eigenvectors_df.rds"))
msg("  dim=", paste(dim(ev), collapse = " x "))
msg("  ID head: ", paste(head(ev$ID, 3), collapse = " | "))
msg("  rownames head: ", paste(head(rownames(ev), 3), collapse = " | "))

msg("reading expr_cpm.rds ...")
cpm_raw <- read_rds(file.path(DATA_EXPR, "expr_cpm.rds"))
msg("  colnames head: ", paste(head(colnames(cpm_raw), 3), collapse = " | "))

msg("reading expr_cpm_adj.rds ...")
cpm_adj <- read_rds(EXPR_CPM_ADJ)
msg("  rownames identical cpm/cpm_adj: ", identical(rownames(cpm_raw), rownames(cpm_adj)))
msg("  colnames identical cpm/cpm_adj: ", identical(colnames(cpm_raw), colnames(cpm_adj)))
msg("  setequal(ev$ID, colnames(cpm)): ", setequal(ev$ID, colnames(cpm_raw)))
msg("  identical order:                ", identical(as.character(ev$ID), colnames(cpm_raw)))

msg("reading expr_genotype_eigenvectors_df.rds ...")
gv <- read_rds(file.path(DATA_EXPR, "expr_genotype_eigenvectors_df.rds"))
msg("  setequal(gv$ID, colnames(cpm)): ", setequal(gv$ID, colnames(cpm_raw)))

# Reproduce the 55-PC adjustment exactly as run_underlier_pipeline.R step 8 did.
E <- as.matrix(ev[, paste0("PC", 1:55), drop = FALSE])
rownames(E) <- as.character(ev$ID)
if (!all(colnames(cpm_adj) %in% rownames(E))) {
  msg("!! sample ids do not line up; cannot test the log2CPM round-trip")
} else {
  msg("running removeBatchEffect with 55 PCs ...")
  rec <- removeBatchEffect(cpm_raw[, colnames(cpm_adj), drop = FALSE],
                           covariates = E[colnames(cpm_adj), , drop = FALSE])
  d <- max(abs(rec - cpm_adj))
  msg("  max |recomputed - saved| = ", format(d, scientific = TRUE))
  msg("  reproduces expr_cpm_adj: ", d < 1e-6)
}
msg("DONE")
