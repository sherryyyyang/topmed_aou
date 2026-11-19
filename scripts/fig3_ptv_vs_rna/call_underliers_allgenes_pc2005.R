#!/usr/bin/env Rscript
# Z-score underliers for every gene in gene_median_tpm_kept.tsv (2005 PCs).
#
# The published matrix is GENCODE v30 protein-coding ∩ HVG (8,124 genes).
# This run keeps those official calls and adds z < -3 + log2CPM-drop calls for
# every other quantified gene (lncRNA, pseudogene, non-HVG coding, etc.).
#
# Extra-gene assay (not in expr_cpm_adj):
#   log2CPM          = edgeR::cpm(log=TRUE) on gene_counts.csv, same 5,560
#                    QC samples as expr_cpm_adj
#   expression_logcpm = log2CPM residualized on 55 expression PCs
#                    (same log2CPM arm as the published definition)
#   expression_zscore = per-gene scale of log2CPM residualized on
#                    2005 expression PCs + 5 genotype PCs
#   underlier         = z < -3 AND log2CPM < gene_mean - 1, then kept_subjects
#
# Writes under pc2005_work/fig3_ptv_vs_rna/allgenes/ (does not touch published
# 55-PC or 2005-PC underlier folders).
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(edgeR)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
Sys.setenv(PC_SETTING = "2005")
source(file.path(ROOT, "scripts/helpers/_paths.R"))

OUT <- file.path(PC_WORK, "fig3_ptv_vs_rna", "allgenes")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
COUNTS_CSV <- file.path(SHEYANG, "counts", "gene_counts.csv")
N_EXPR_Z <- as.integer(PC_SETTING)
N_EXPR_LOG <- 55L
N_GENO <- 5L
Z_CUT <- -3
LOGCPM_DROP <- 1
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}
cache_rds <- function(path, expr) {
  if (file.exists(path)) {
    msg("cache load ", path)
    return(read_rds(path))
  }
  val <- force(expr)
  write_rds(val, path)
  msg("wrote ", path)
  val
}

proj_out <- function(Y, Q) {
  Y - tcrossprod(Y %*% Q, Q)
}
cov_Q <- function(C) {
  C <- scale(C, center = TRUE, scale = FALSE)
  qr.Q(qr(C))
}

msg("loading axes, PCs, TPM, official underliers ...")
cpm_pub <- read_rds(EXPR_CPM_ADJ)
samples <- colnames(cpm_pub)
pub_genes <- rownames(cpm_pub)
pub_nv <- str_remove(pub_genes, "\\..*$")
kept <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept)

meta <- read_tsv(METADATA, show_col_types = FALSE)
ev <- read_rds(file.path(DATA_EXPR, "expr_eigenvectors_df.rds"))
gv <- read_rds(file.path(DATA_EXPR, "expr_genotype_eigenvectors_df.rds"))
ev$NWD <- meta$WGS[match(ev$ID, meta$RNA)]
stopifnot(setequal(ev$NWD, samples), setequal(gv$ID, samples))

epc_all <- as.matrix(ev[match(samples, ev$NWD), paste0("PC", seq_len(N_EXPR_Z))])
rownames(epc_all) <- samples
gpc <- as.matrix(gv[match(samples, gv$ID), paste0("Genotype_PC", seq_len(N_GENO))])
rownames(gpc) <- samples
Q_log <- cov_Q(epc_all[, seq_len(N_EXPR_LOG), drop = FALSE])
Q_z <- cov_Q(cbind(epc_all, gpc))
rm(ev, gv, epc_all, gpc)
gc(verbose = FALSE)

tpm <- read_tsv(file.path(DATA_CORE, "gene_median_tpm_kept.tsv"), show_col_types = FALSE) |>
  mutate(gene_nv = str_remove(gene_nv, "\\..*$")) |>
  distinct(gene_nv, .keep_all = TRUE)
tpm_nv <- tpm$gene_nv

orig_under <- read_rds(file.path(DATA_UNDERLIER, "z_score_underlier_z_-3_clean.rds")) |>
  filter(subjectid %in% kept)

EXTRA_UNDER <- file.path(OUT, "z_score_underlier_z_-3_extra_genes.rds")
extra_under <- cache_rds(EXTRA_UNDER, {
  msg("reading gene_counts.csv ...")
  dt <- fread(COUNTS_CSV, sep = ",", header = TRUE, data.table = FALSE)
  gene_ids <- dt[[1]]
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_ids
  rm(dt)
  gc(verbose = FALSE)

  nwd <- meta$WGS[match(colnames(mat), meta$RNA)]
  keep_col <- !is.na(nwd) & nwd %in% samples & !duplicated(nwd)
  mat <- mat[, keep_col, drop = FALSE]
  colnames(mat) <- nwd[keep_col]
  stopifnot(setequal(colnames(mat), samples))
  mat <- mat[, samples, drop = FALSE]
  msg("  counts aligned: ", nrow(mat), " genes x ", ncol(mat), " samples")

  extra_idx <- which(!str_remove(rownames(mat), "\\..*$") %in% pub_nv &
                       str_remove(rownames(mat), "\\..*$") %in% tpm_nv)
  extra <- mat[extra_idx, , drop = FALSE]
  rm(mat)
  gc(verbose = FALSE)
  msg("  extra genes in TPM table, not in expr_cpm_adj: ", nrow(extra))

  dge <- DGEList(extra)
  logcpm <- cpm(dge, log = TRUE)
  rm(dge, extra)
  gc(verbose = FALSE)

  msg("  residualizing log2CPM (55 PCs for drop, 2005+5 for z) ...")
  log_adj <- proj_out(logcpm, Q_log)
  z_resid <- proj_out(logcpm, Q_z)
  rm(logcpm)
  gc(verbose = FALSE)
  z <- t(scale(t(z_resid)))
  rm(z_resid)
  gc(verbose = FALSE)

  gene_mean <- rowMeans(log_adj)
  ok <- is.finite(rowMeans(z)) & is.finite(gene_mean)
  log_adj <- log_adj[ok, , drop = FALSE]
  z <- z[ok, , drop = FALSE]
  gene_mean <- gene_mean[ok]
  msg("  finite extra genes: ", nrow(z))

  pass <- (z < Z_CUT) & (log_adj < (gene_mean - LOGCPM_DROP))
  kept_col <- colnames(z) %in% kept
  pass[, !kept_col] <- FALSE
  hit <- which(pass, arr.ind = TRUE)
  tibble(
    gene_id = rownames(z)[hit[, 1]],
    subjectid = colnames(z)[hit[, 2]],
    expression_zscore = z[hit],
    expression_logcpm = log_adj[hit]
  )
})

under <- bind_rows(
  orig_under |> transmute(gene_id, subjectid, expression_zscore, expression_logcpm),
  extra_under
) |>
  mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
  distinct(subjectid, gene_nv, .keep_all = TRUE)
msg("combined underlier pairs: ", nrow(under),
    "  genes: ", n_distinct(under$gene_nv))

ptv <- read_tsv(PTV_TSV, show_col_types = FALSE) |>
  filter(sample %in% kept) |>
  mutate(gene_nv = str_remove(gene, "\\..*$")) |>
  distinct(subjectid = sample, gene_nv) |>
  count(gene_nv, name = "n_ptv_carriers")

under_gene <- under |>
  distinct(subjectid, gene_nv) |>
  count(gene_nv, name = "n_underlier_subjects")

loeuf_nv <- read_rds(LOEUF_LOOKUP) |>
  distinct(gene_id_no_version, .keep_all = TRUE) |>
  select(gene_nv = gene_id_no_version, symbol, LOEUF, LOEUF_tier)
gencode_nv <- read_tsv(
  GENCODE_SYM, col_names = c("gene_nv", "gencode_symbol"), show_col_types = FALSE
) |>
  mutate(gene_nv = str_remove(gene_nv, "\\..*$")) |>
  distinct(gene_nv, .keep_all = TRUE)

tpm_ok_nv <- tpm |> filter(median_tpm > 1) |> pull(gene_nv)
TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")

genes <- full_join(ptv, under_gene, by = "gene_nv") |>
  mutate(
    n_ptv_carriers = replace_na(n_ptv_carriers, 0L),
    n_underlier_subjects = replace_na(n_underlier_subjects, 0L)
  ) |>
  filter(gene_nv %in% tpm_nv) |>
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

out_tsv <- file.path(OUT, "ptv_vs_underlier_gene_level_z_-3_allgenes.tsv")
write_tsv(
  genes |>
    select(
      gene_nv, gene_label, symbol, LOEUF, LOEUF_tier,
      n_ptv_carriers, n_underlier_subjects, pct_ptv, pct_under
    ) |>
    arrange(desc(pct_under)),
  out_tsv
)
msg("Fig 3a all-gene table: ", nrow(genes), " genes  -> ", out_tsv)
write_rds(under, file.path(OUT, "z_score_underlier_z_-3_allgenes.rds"))
msg("DONE")
