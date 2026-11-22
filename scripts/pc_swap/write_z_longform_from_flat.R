#!/usr/bin/env Rscript
# Write expr_scaled_genotype_pc_5_expression_pc_<N>_longer_form.tsv from
# expr_scaled_flat_list.rds (or residualize if that NPC is absent from the flat
# list). Row order matches expr_cpm_adj so build_underlier_artifacts.R accepts it.
#
#   PC_SETTING=505 Rscript scripts/pc_swap/write_z_longform_from_flat.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
  library(data.table); library(tibble)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/_enrichment_helpers.R"))

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
N_GENO <- 5L
npc <- as.integer(PC_SETTING)
stopifnot(npc >= 1L)
out <- Z_LONGFORM
if (file.exists(out) && !identical(Sys.getenv("FORCE", "0"), "1")) {
  msg("already exists (set FORCE=1 to overwrite): ", out)
  quit(save = "no", status = 0)
}

msg("loading expr_cpm_adj.rds for axis order ...")
cpm <- read_rds(EXPR_CPM_ADJ)
genes <- rownames(cpm); subjects <- colnames(cpm)
msg("  ", length(genes), " genes x ", length(subjects), " subjects")

flat_path <- file.path(DATA_EXPR, "expr_scaled_flat_list.rds")
msg("loading ", flat_path, " (~28 GB) ...")
flat <- read_rds(flat_path)
nm <- names(flat)
npcs_available <- as.integer(str_remove(nm, "^.*_expression_pc_"))
msg("  flat list has ", length(flat), " matrices; npc range ",
    min(npcs_available), "..", max(npcs_available))

i <- match(npc, npcs_available)
if (!is.na(i)) {
  msg("using flat matrix for npc=", npc, " (index ", i, ", name=", nm[i], ")")
  mat <- flat[[i]]
} else {
  msg("npc=", npc, " not in flat list; residualizing from eigenvectors ...")
  METADATA_PATH <- Sys.getenv("TOPMED_METADATA", unset = "")
  if (!nzchar(METADATA_PATH)) METADATA_PATH <- file.path(DATA_CORE, "metadata.filtered.txt")
  expr_qc_normalized <- read_rds(file.path(DATA_EXPR, "expr_qc_normalized.rds"))
  expr_eigenvectors_df <- rekey_expression_pcs(
    read_rds(file.path(DATA_EXPR, "expr_eigenvectors_df.rds")),
    expr_qc_normalized, METADATA_PATH)
  expr_genotype_eigenvectors_df <- read_rds(file.path(DATA_EXPR, "expr_genotype_eigenvectors_df.rds"))
  mat <- residualize_phenotypes(npc, N_GENO, expr_qc_normalized,
                                expression_pcs = expr_eigenvectors_df,
                                genotype_pcs = expr_genotype_eigenvectors_df)
}
rm(flat); gc(verbose = FALSE)

if (is.null(rownames(mat)) || is.null(colnames(mat))) {
  stop("scaled matrix is missing dimnames")
}
strip_ver <- function(x) sub("\\..*$", "", x)
rg <- rownames(mat); cg <- colnames(mat)
if (setequal(rg, genes)) {
  mat <- mat[genes, , drop = FALSE]
} else if (setequal(strip_ver(rg), strip_ver(genes))) {
  map <- match(strip_ver(genes), strip_ver(rg))
  stopifnot(!anyNA(map))
  mat <- mat[map, , drop = FALSE]
  rownames(mat) <- genes
} else {
  stop("gene dimnames of scaled matrix do not match expr_cpm_adj")
}
if (!setequal(cg, subjects)) stop("subject dimnames do not match expr_cpm_adj")
mat <- mat[, subjects, drop = FALSE]
stopifnot(identical(rownames(mat), genes), identical(colnames(mat), subjects))

msg("writing long-form TSV (gene-major, subject-fastest) ...")
t0 <- Sys.time()
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(out, ".tmp")
if (file.exists(tmp)) file.remove(tmp)
writeLines("gene_id\tsubjectid\texpression_zscore", tmp)
block <- 100L
n_blocks <- ceiling(length(genes) / block)
for (b in seq_len(n_blocks)) {
  i1 <- (b - 1L) * block + 1L
  i2 <- min(b * block, length(genes))
  g <- genes[i1:i2]
  m <- mat[g, , drop = FALSE]
  dt <- data.table(
    gene_id = rep(g, each = length(subjects)),
    subjectid = rep(subjects, times = length(g)),
    expression_zscore = as.vector(t(m))
  )
  fwrite(dt, tmp, sep = "\t", col.names = FALSE, quote = FALSE, append = TRUE)
  if (b %% 20L == 0L || b == n_blocks) {
    msg(sprintf("  genes %d/%d (%.1f min)", i2, length(genes),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
file.rename(tmp, out)
msg("wrote ", out)
msg(sprintf("elapsed %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

ref <- file.path(DATA_EXPR, "expr_scaled_genotype_pc_5_expression_pc_55_longer_form.tsv")
if (file.exists(ref)) {
  msg("verifying key columns match 55-PC longform ...")
  a <- fread(cmd = sprintf("head -n %d %s | cut -f1,2 | tail -n +2", length(subjects) + 1, out),
             header = FALSE)
  b <- fread(cmd = sprintf("head -n %d %s | cut -f1,2 | tail -n +2", length(subjects) + 1, ref),
             header = FALSE)
  stopifnot(identical(a[[1]], b[[1]]), identical(a[[2]], b[[2]]))
  msg("  first gene block keys match 55-PC longform")
}
msg("done")
