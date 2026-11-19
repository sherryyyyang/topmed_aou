#!/usr/bin/env Rscript
# All-gene HVG underliers at 2005 expression PCs + 5 genotype PCs, then the
# Fig 3a gene-level PTV vs underlier table.
#
# z is residualize_expression on the FULL qc_normalized HVG matrix
# (~18,062 genes x 6,024 samples), NOT the 8,119 shared-coding subset in
# expr_zscore_pc2005.rds. Log2CPM is the pipeline 55-PC-adjusted expr_cpm_adj.
#
# Underlier: z_2005 < -3 AND log2CPM < gene_mean - 1.
# Gene-mean grid = all samples used for residualize / scale (prefer the
# 6,024-sample z grid). Calls are then restricted to kept_subjects (5,489)
# for the Fig 3a table, matching plot_ptv_vs_underlier_gene_level.R.
#
# Universe = all-gene HVG (rownames of qc_normalized / expr_cpm_adj), then
# median TPM > 1, then n_ptv > 0 OR n_underlier > 0. Does not use the
# published coding-only HVG universe.
#
# Writes:
#   allgenes_pipeline_out/expr_zscore_pc2005_allgenes.rds
#   allgenes_pipeline_out/ptv_vs_underlier_gene_level_z_-3_allgenes_pc2005.tsv
# and a clearly named copy under submission_package_nov_2026/supplemental_data/.
# Does NOT overwrite Source_Data_Fig3a_ptv_vs_underlier_gene.tsv.
# Does NOT write into data/underliers, data/expression, or final_paper_draft_pc55.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Matrix)
  library(limma)
  library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

OUT <- file.path(ROOT, "allgenes_pipeline_out")
DRAFT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
SUPP <- file.path(DRAFT, "supplemental_data")
QC_NORM <- file.path(OUT, "qc_normalized.rds")
EXPR_PC <- file.path(OUT, "expr_pc_df.rds")
GENO_PC <- file.path(OUT, "geno_pc_df.rds")
NEW_CPM <- file.path(OUT, "expr_cpm_adj.rds")
Z_CACHE <- file.path(OUT, "expr_zscore_pc2005_allgenes.rds")
GENE_TSV <- file.path(OUT, "ptv_vs_underlier_gene_level_z_-3_allgenes_pc2005.tsv")
SUPP_TSV <- file.path(SUPP, "Source_Data_Fig3a_ptv_vs_underlier_gene_allgenes_hvg_pc2005.tsv")
PUB_FIG3A_TSV <- file.path(SUPP, "Source_Data_Fig3a_ptv_vs_underlier_gene.tsv")

N_EXPR_PCS <- 2005L
N_GENO_PCS <- 5L
Z_CUT <- -3
LOGCPM_DROP <- 1
TPM_MIN <- 1
TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")

# Refuse to clobber published coding-HVG Fig 3a source data or 55-PC draft.
.protected <- c(
  PUB_FIG3A_TSV,
  file.path(ROOT, "final_paper_draft_pc55", "supplemental_data",
            "Source_Data_Fig3a_ptv_vs_underlier_gene.tsv")
)
.stop_if_protected <- function(path) {
  if (normalizePath(path, mustWork = FALSE) %in%
      vapply(.protected, normalizePath, character(1), mustWork = FALSE)) {
    stop("Refusing to overwrite protected published table: ", path)
  }
}

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

nv <- function(x) str_remove(x, "\\..*$")

# Match run_underlier_pipeline_allgenes.R::residualize_expression exactly.
residualize_expression <- function(normalized_expression, expr_pc_df, geno_pc_df,
                                   n_expr_pcs, n_geno_pcs) {
  epc <- dplyr::select(expr_pc_df, ID, dplyr::one_of(str_c("PC", 1:n_expr_pcs)))
  gpc <- dplyr::select(geno_pc_df, ID, dplyr::one_of(str_c("Genotype_PC", 1:n_geno_pcs)))
  pheno <- left_join(epc, gpc, by = "ID") |> drop_na() |> column_to_rownames("ID") |> as.matrix()
  mat <- as.matrix(normalized_expression[, rownames(pheno)])
  resid <- removeBatchEffect(mat, covariates = pheno)
  t(scale(t(resid)))
}

# --------------------------------------------------------------------------
# Full all-gene HVG z at 2005 expr PCs + 5 geno PCs.
# --------------------------------------------------------------------------
if (file.exists(Z_CACHE)) {
  msg("Loading cached all-gene 2005-PC z: ", Z_CACHE)
  z_mat <- read_rds(Z_CACHE)
  msg("  z matrix: ", comma(nrow(z_mat)), " genes x ",
      comma(ncol(z_mat)), " samples")
} else {
  msg("Computing full-matrix 2005-PC z (removeBatchEffect + per-gene scale) ...")
  msg("  loading qc_normalized / expr_pc_df / geno_pc_df ...")
  qn <- read_rds(QC_NORM)
  expr_pc_df <- read_rds(EXPR_PC)
  geno_pc_df <- read_rds(GENO_PC)
  pc_need <- str_c("PC", seq_len(N_EXPR_PCS))
  gpc_need <- str_c("Genotype_PC", seq_len(N_GENO_PCS))
  if (!all(pc_need %in% names(expr_pc_df))) {
    stop("expr_pc_df is missing PC1..PC", N_EXPR_PCS)
  }
  if (!all(gpc_need %in% names(geno_pc_df))) {
    stop("geno_pc_df is missing Genotype_PC1..", N_GENO_PCS)
  }

  # dgeMatrix: materialize so rownames() is a plain character vector.
  qn <- as.matrix(qn)
  if (is.null(rownames(qn)) || !length(rownames(qn))) {
    stop("qc_normalized has no rownames after as.matrix()")
  }
  msg("  FT matrix: ", comma(nrow(qn)), " genes x ", comma(ncol(qn)), " samples")
  msg("  residualize_expression n_expr_pcs=", N_EXPR_PCS,
      " n_geno_pcs=", N_GENO_PCS, " ...")
  z_mat <- residualize_expression(
    qn, expr_pc_df, geno_pc_df, N_EXPR_PCS, N_GENO_PCS
  )
  rm(qn, expr_pc_df, geno_pc_df)
  gc(verbose = FALSE)
  write_rds(z_mat, Z_CACHE, compress = "gz")
  msg("  wrote ", Z_CACHE, "  ", nrow(z_mat), " x ", ncol(z_mat))
}

# --------------------------------------------------------------------------
# Underlier calls, then Fig 3a gene table.
# --------------------------------------------------------------------------
msg("Loading 55-PC-adjusted log2CPM (allgenes expr_cpm_adj) ...")
cpm <- as.matrix(read_rds(NEW_CPM))
kept <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept)
msg("  kept_subjects: ", comma(n_kept))

if (!setequal(rownames(z_mat), rownames(cpm))) {
  stop("z and expr_cpm_adj gene IDs do not match")
}
eligible_ids <- rownames(z_mat)
eligible_nv <- unique(nv(eligible_ids))
z_subj <- colnames(z_mat)
cpm_subj <- colnames(cpm)
if (!setequal(z_subj, cpm_subj)) {
  stop("z and expr_cpm_adj sample IDs do not match")
}
# Align columns to the z / residualize grid.
cpm <- cpm[rownames(z_mat), z_subj, drop = FALSE]
z_align <- z_mat[, z_subj, drop = FALSE]

# Gene mean on the same samples used for residualize / scale (document this).
gene_mean_grid <- "allgene_qc_z_samples"
n_gene_mean <- ncol(z_align)
gene_mean <- rowMeans(cpm)
msg(sprintf(
  "  gene_mean grid: %s (n=%s); underlier calls then restricted to kept_subjects",
  gene_mean_grid, comma(n_gene_mean)
))

under_mat <- (z_align < Z_CUT) & (cpm < (gene_mean - LOGCPM_DROP))
storage.mode(under_mat) <- "logical"
kept_in_z <- intersect(kept, z_subj)
msg(sprintf(
  "  kept in z grid: %s / %s; kept missing from z: %s",
  comma(length(kept_in_z)), comma(n_kept),
  comma(length(setdiff(kept, z_subj)))
))
under_kept <- under_mat[, kept_in_z, drop = FALSE]
hit <- which(under_kept, arr.ind = TRUE)
under <- tibble(
  gene_id = rownames(under_kept)[hit[, 1]],
  subjectid = colnames(under_kept)[hit[, 2]]
) |>
  mutate(gene_nv = nv(gene_id)) |>
  distinct(subjectid, gene_nv)
msg("  underlier pairs in kept_subjects: ", comma(nrow(under)),
    "  genes: ", comma(n_distinct(under$gene_nv)))
rm(cpm, z_align, under_mat, under_kept, hit, z_mat)
gc(verbose = FALSE)
msg("  all-gene HVG universe: ", comma(length(eligible_nv)))

tpm <- read_tsv(file.path(DATA_CORE, "gene_median_tpm_kept.tsv"), show_col_types = FALSE) |>
  mutate(gene_nv = nv(gene_nv)) |>
  distinct(gene_nv, .keep_all = TRUE)
tpm_ok_nv <- tpm |> filter(median_tpm > TPM_MIN) |> pull(gene_nv)
msg("  genes with median TPM > ", TPM_MIN, ": ", comma(length(tpm_ok_nv)))

msg("Loading PTV table (kept_subjects only) ...")
ptv <- read_tsv(PTV_TSV, show_col_types = FALSE) |>
  filter(sample %in% kept) |>
  mutate(gene_nv = nv(gene)) |>
  distinct(subjectid = sample, gene_nv) |>
  count(gene_nv, name = "n_ptv_carriers")

under_gene <- under |>
  count(gene_nv, name = "n_underlier_subjects")

loeuf_nv <- read_rds(LOEUF_LOOKUP) |>
  distinct(gene_id_no_version, .keep_all = TRUE) |>
  select(gene_nv = gene_id_no_version, symbol, LOEUF, LOEUF_tier)
gencode_nv <- read_tsv(
  GENCODE_SYM, col_names = c("gene_nv", "gencode_symbol"), show_col_types = FALSE
) |>
  mutate(gene_nv = nv(gene_nv)) |>
  distinct(gene_nv, .keep_all = TRUE)

genes <- full_join(ptv, under_gene, by = "gene_nv") |>
  mutate(
    n_ptv_carriers = replace_na(n_ptv_carriers, 0L),
    n_underlier_subjects = replace_na(n_underlier_subjects, 0L)
  ) |>
  filter(gene_nv %in% eligible_nv) |>
  filter(gene_nv %in% tpm_ok_nv) |>
  filter(n_ptv_carriers > 0 | n_underlier_subjects > 0) |>
  left_join(loeuf_nv, by = "gene_nv") |>
  left_join(gencode_nv, by = "gene_nv") |>
  mutate(
    symbol = coalesce(na_if(gencode_symbol, ""), na_if(symbol, ""), gene_nv),
    gene_label = symbol,
    pct_ptv = 100 * n_ptv_carriers / n_kept,
    pct_under = 100 * n_underlier_subjects / n_kept,
    LOEUF_tier = factor(
      ifelse(is.na(LOEUF_tier), "No LOEUF", as.character(LOEUF_tier)),
      levels = TIER_LEVELS
    )
  )

out_tbl <- genes |>
  select(
    gene_nv, gene_label, symbol, LOEUF, LOEUF_tier,
    n_ptv_carriers, n_underlier_subjects, pct_ptv, pct_under
  ) |>
  arrange(desc(pct_under))

.stop_if_protected(GENE_TSV)
write_tsv(out_tbl, GENE_TSV)
msg("Fig 3a all-gene HVG table: ", comma(nrow(out_tbl)), " genes  -> ", GENE_TSV)

dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)
.stop_if_protected(SUPP_TSV)
file.copy(GENE_TSV, SUPP_TSV, overwrite = TRUE)
msg("copied source TSV: ", SUPP_TSV)
if (!file.exists(PUB_FIG3A_TSV)) {
  stop("Published Source_Data_Fig3a TSV is missing (should have been left alone)")
}
msg("published Source_Data_Fig3a left alone: ", PUB_FIG3A_TSV)
msg(sprintf(
  "PROVENANCE gene_mean_grid=%s n=%s; call_subjects=kept_subjects n=%s; z=2005expr+5geno; log2CPM=allgenes_55PC",
  gene_mean_grid, n_gene_mean, n_kept
))
msg("DONE")
