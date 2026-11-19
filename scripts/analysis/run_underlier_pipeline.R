#!/usr/bin/env Rscript
# =============================================================================
# run_underlier_pipeline.R
#
# End-to-end, reproducible expression-outlier ("underlier") pipeline distilled
# from outlier_processing.R. Produces per-(gene,subject) residualized expression
# z-scores and underlier calls, starting from raw gene counts + PLINK genotypes.
#
# RUN IN ONE LINE:
#     Rscript run_underlier_pipeline.R                 # uses the CONFIG below
#     Rscript run_underlier_pipeline.R my_config.R     # optional: override CONFIG
#
# To adapt to a new dataset you should ONLY need to edit the CONFIG block below
# (input/output paths + a few parameters). Every intermediate is cached to
# CFG$out_dir as .rds, so the run is restartable and each step is inspectable.
#
# PIPELINE (faithful to the released z-scores: 55 expression PCs + 5 genotype PCs)
#   1. Load counts + metadata; rename RNA sample ids -> WGS ids (id-harmonization)
#   2. Restrict to protein-coding genes (GENCODE)
#   3. filter_expression(): highly-variable-gene selection (Poisson & CV2)
#   4. Freeman-Tukey normalization (corral)
#   5. Sample QC: network connectivity outlier removal (WGCNA), Z < CONN_Z dropped
#   6. Re-filter + re-normalize the QC-passing samples
#   7. Expression PCA (PCAtools) -> take first N_EXPR_PCS PCs
#   8. log2-CPM (edgeR) and PC-adjust it   -> expression_logcpm  (2-fold criterion)
#   9. Genotype PCA (PLINK BED -> GDS -> snpgdsPCA) -> first N_GENO_PCS PCs
#  10. Residualize normalized expression on [expr PCs + genotype PCs], scale
#      per gene -> expression_zscore
#  11. Join z-score + log2CPM; call underliers per z-cutoff
#      (expression_zscore < z  AND  expression_logcpm < gene_mean - LOGCPM_DROP)
# =============================================================================

# ------------------------------- CONFIG (EDIT) -------------------------------
CFG <- list(
  # --- inputs ---
  counts_csv   = "../counts/gene_counts.csv",        # genes x samples; a gene-id column (CFG$gene_id_col) + one column per sample
  gene_id_col  = "Name",                             # name of the gene-id column in counts_csv
  # OPTIONAL sample-id harmonization. Set metadata_tsv = NULL to skip renaming
  # (use this when your counts columns ALREADY match the PLINK/genotype sample ids).
  # If set, counts columns (meta_rna_col) are renamed to genotype ids (meta_wgs_col).
  metadata_tsv = NULL,                               # e.g. ".../metadata.filtered.txt", or NULL to skip
  meta_rna_col = "RNA",                              # metadata column matching counts sample names
  meta_wgs_col = "WGS",                              # metadata column matching PLINK/genotype sample ids
  gencode_gff  = "../annotation/gencode.v30.annotation.gff3",
  plink_dir    = "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/genotypes/plink",  # holds chr{1..22,X}.{bed,bim,fam}
  plink_chroms = c(1:22, "X"),

  # --- output ---
  out_dir      = "./underlier_pipeline_out",

  # --- parameters ---
  # NOTE: the number of PCs to regress out should be OPTIMIZED for each dataset
  # (e.g. via Gavish-Donoho / an outlier-enrichment sweep). My original TOPMed
  # run used 55 expression PCs and 5 genotype PCs.
  n_expr_pcs   = 55,                                 # expression PCs regressed out (my original: 55); optimize per dataset
  n_geno_pcs   = 5,                                  # genotype PCs regressed out (my original: 5); optimize per dataset
  conn_z       = -3,                                 # drop samples with network-connectivity Z below this
  # Emit BOTH underlier definitions as SEPARATE, labelled outputs (not either/or):
  #   haplo      -> haploinsufficient (log2CPM-only) underliers
  #   strict z<-3 -> z-score + log2CPM underliers
  do_haplo     = TRUE,                               # haplo: log2CPM < gene_mean - logcpm_drop (no z threshold)
  z_cutoffs    = c(-3),                              # strict: expression_zscore < z AND the log2CPM drop (add more if needed)
  logcpm_drop  = 1,                                  # log2CPM < gene_mean - logcpm_drop (>=2-fold drop); used by both
  ld_maf       = 0.05, ld_missing = 0.10, ld_thresh = sqrt(0.1), ld_slide_bp = 10e6,  # genotype LD-pruning
  n_threads    = 16L
)
# -----------------------------------------------------------------------------

# optional: Rscript run_underlier_pipeline.R my_config.R  (sourced to override CFG)
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1 && file.exists(.args[1])) { source(.args[1], local = TRUE); message("Applied config overrides from ", .args[1]) }

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(corral); library(PCAtools)
  library(scran); library(scuttle); library(WGCNA); library(SNPRelate); library(gdsfmt)
  library(rtracklayer); library(vroom); library(tidyverse); library(magrittr)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
setwd_out <- function(f) file.path(CFG$out_dir, f)

# cache helper: compute-once, restartable. Each step writes <name>.rds to out_dir.
cache <- function(name, expr) {
  path <- setwd_out(paste0(name, ".rds"))
  if (file.exists(path)) { message("[cache] load  ", name); return(read_rds(path)) }
  message("[run]   ", name); val <- force(expr); write_rds(val, path); val
}

# ------------------------------- FUNCTIONS -----------------------------------
# HVG selection: intersect Poisson mean-variance HVGs with CV2 HVGs (scran).
filter_expression <- function(counts) {
  gv_pois <- modelGeneVarByPoisson(counts)
  top_pois <- getTopHVGs(gv_pois)
  gv_cv2 <- modelGeneCV2(counts)
  top_cv2 <- getTopHVGs(gv_cv2, var.field = "ratio", var.threshold = 1.0)
  counts[intersect(top_pois, top_cv2), ]
}

# Per-gene z from residualizing normalized expression on expression + genotype PCs.
residualize_expression <- function(normalized_expression, expr_pc_df, geno_pc_df, n_expr_pcs, n_geno_pcs) {
  epc <- dplyr::select(expr_pc_df, ID, dplyr::one_of(str_c("PC", 1:n_expr_pcs)))
  gpc <- dplyr::select(geno_pc_df, ID, dplyr::one_of(str_c("Genotype_PC", 1:n_geno_pcs)))
  pheno <- left_join(epc, gpc, by = "ID") |> drop_na() |> column_to_rownames("ID") |> as.matrix()
  mat <- as.matrix(normalized_expression[, rownames(pheno)])
  resid <- removeBatchEffect(mat, covariates = pheno)   # regress out the PCs
  t(scale(t(resid)))                                     # standardize each gene across samples -> z
}

# Underlier caller. Both definitions require the >=2-fold log2CPM drop vs the
# gene's mean; the strict definition additionally requires expression_zscore < z.
#   z_cutoff = NULL  -> "haplo" (log2CPM-only)
#   z_cutoff = -3    -> "strict z<-3"
call_underliers <- function(expr_z_join, z_cutoff, logcpm_drop) {
  out <- expr_z_join |>
    group_by(gene_id) |>
    mutate(gene_mean_logcpm = mean(expression_logcpm)) |>
    ungroup() |>
    filter(expression_logcpm < gene_mean_logcpm - logcpm_drop)
  if (!is.null(z_cutoff)) out <- filter(out, expression_zscore < z_cutoff)
  out
}

# ------------------------------- PIPELINE ------------------------------------
## 1. counts (+ OPTIONAL sample-id harmonization RNA -> WGS via metadata)
use_meta <- !is.null(CFG$metadata_tsv) && nzchar(CFG$metadata_tsv) && file.exists(CFG$metadata_tsv)
rna2wgs <- NULL
if (use_meta) {
  meta <- read_tsv(CFG$metadata_tsv, show_col_types = FALSE)
  rna2wgs <- setNames(meta[[CFG$meta_wgs_col]], meta[[CFG$meta_rna_col]])
  message("  metadata provided: renaming counts sample ids ", CFG$meta_rna_col, " -> ", CFG$meta_wgs_col)
} else {
  message("  no metadata: assuming counts sample ids already match the genotype/PLINK ids")
}

raw_counts <- cache("raw_counts", {
  d <- vroom(CFG$counts_csv, show_col_types = FALSE)
  m <- column_to_rownames(as.data.frame(d), CFG$gene_id_col) |> as.matrix()
  if (use_meta) {                              # rename RNA ids to genotype (WGS) ids
    keep <- colnames(m)[colnames(m) %in% names(rna2wgs) & !is.na(rna2wgs[colnames(m)])]
    m <- m[, keep, drop = FALSE]
    colnames(m) <- rna2wgs[colnames(m)]
  }
  m[, !is.na(colnames(m)) & !duplicated(colnames(m)), drop = FALSE]
})

## 2. restrict to protein-coding genes (GENCODE); keep a gene_id -> symbol map
gene_meta <- cache("gene_meta", {
  g <- as_tibble(readGFF(CFG$gencode_gff)) |>
    filter(type == "gene", gene_type == "protein_coding") |>
    transmute(gene_id, gene_nv = str_remove(gene_id, "\\..*$"), symbol = gene_name) |>
    distinct(gene_id, .keep_all = TRUE)
  g
})
coding_ids <- unique(gene_meta$gene_id)
coding_counts <- raw_counts[rownames(raw_counts) %in% coding_ids, , drop = FALSE]

## 3-4. HVG selection + Freeman-Tukey normalization (first pass, for QC)
norm1 <- cache("normalized_pass1", corral_preproc(filter_expression(coding_counts), rtype = "freemantukey"))
colnames(norm1) <- colnames(coding_counts)

## 5. connectivity-based sample QC (drop network-outlier samples)
keep_samples <- cache("qc_kept_samples", {
  adj <- 0.5 + 0.5 * cor(norm1)
  conn <- fundamentalNetworkConcepts(adj)$Connectivity
  names(conn) <- colnames(norm1)
  z <- (conn - mean(conn)) / sd(conn)
  kept <- names(z)[z > CFG$conn_z]
  message(sprintf("  connectivity QC: removed %d of %d samples", ncol(norm1) - length(kept), ncol(norm1)))
  kept
})

## 6. re-filter + re-normalize on QC-passing samples
qc_counts <- coding_counts[, keep_samples, drop = FALSE]
qc_normalized <- cache("qc_normalized", {
  n <- corral_preproc(filter_expression(qc_counts), rtype = "freemantukey")
  colnames(n) <- colnames(qc_counts); n
})

## 7. expression PCA -> first N_EXPR_PCS PCs  (this IS expr_eigenvectors_subset)
expr_pc_df <- cache("expr_pc_df", {
  pca_res <- pca(as.matrix(qc_normalized), scale = TRUE)
  ev <- pca_res$rotated
  rownames(ev) <- colnames(qc_normalized)
  as.data.frame(ev) |> rownames_to_column("ID")
})
expr_eigenvectors_subset <- expr_pc_df |>
  dplyr::select(ID, dplyr::one_of(str_c("PC", 1:CFG$n_expr_pcs))) |>
  column_to_rownames("ID") |> as.matrix()

## 8. log2-CPM (edgeR) on raw counts, then PC-adjust  -> expression_logcpm
expr_cpm_adj <- cache("expr_cpm_adj", {
  counts_aligned <- coding_counts[rownames(qc_normalized), rownames(expr_eigenvectors_subset), drop = FALSE]
  cpm_log <- edgeR::cpm(DGEList(counts_aligned), log = TRUE)
  removeBatchEffect(cpm_log, covariates = expr_eigenvectors_subset)
})
expr_cpm_long <- as_tibble(expr_cpm_adj, rownames = "gene_id") |>
  pivot_longer(-gene_id, names_to = "subjectid", values_to = "expression_logcpm")

## 9. genotype PCA:  PLINK BED -> per-chr GDS -> combined GDS -> snpgdsPCA
geno_pc_df <- cache("geno_pc_df", {
  gds_dir <- setwd_out("gds"); dir.create(gds_dir, showWarnings = FALSE)
  for (chr in CFG$plink_chroms) {
    stem <- file.path(CFG$plink_dir, paste0("chr", chr))
    out_gds <- file.path(gds_dir, paste0("chr", chr, ".gds"))
    if (!file.exists(out_gds))
      snpgdsBED2GDS(paste0(stem, ".bed"), paste0(stem, ".fam"), paste0(stem, ".bim"), out_gds)
  }
  combined <- setwd_out("snv_genotypes.gds")
  if (!file.exists(combined))
    snpgdsCombineGeno(list.files(gds_dir, pattern = "\\.gds$", full.names = TRUE), combined,
                      method = "exact", compress.annotation = "", compress.geno = "")
  gds <- snpgdsOpen(combined)
  on.exit(snpgdsClose(gds), add = TRUE)
  shared <- intersect(read.gdsn(index.gdsn(gds, "sample.id")), keep_samples)
  pruned <- snpgdsLDpruning(gds, method = "corr", slide.max.bp = CFG$ld_slide_bp, maf = CFG$ld_maf,
                            missing.rate = CFG$ld_missing, ld.threshold = CFG$ld_thresh,
                            num.thread = CFG$n_threads, sample.id = shared, verbose = FALSE)
  gpca <- snpgdsPCA(gds, snp.id = unlist(pruned), sample.id = shared,
                    algorithm = "exact", num.thread = CFG$n_threads, eigen.cnt = -1L, verbose = FALSE)
  ev <- gpca$eigenvect
  rownames(ev) <- gpca$sample.id
  colnames(ev) <- str_c("Genotype_PC", 1:ncol(ev))
  as.data.frame(ev) |> rownames_to_column("ID")
})

## 10. residualize normalized expression on [expr PCs + genotype PCs] -> z-scores
expr_z_long <- cache("expr_zscore_long", {
  z <- residualize_expression(qc_normalized, expr_pc_df, geno_pc_df, CFG$n_expr_pcs, CFG$n_geno_pcs)
  as.data.frame(z) |> rownames_to_column("gene_id") |>
    pivot_longer(-gene_id, names_to = "subjectid", values_to = "expression_zscore")
})

## 11. join z + log2CPM, then call underliers per z-cutoff
expr_z_join <- cache("expr_z_join", left_join(expr_z_long, expr_cpm_long, by = c("gene_id", "subjectid")))
write_tsv(expr_z_join, setwd_out("expr_z_join.tsv.gz"))

# definitions to emit: haplo (log2CPM-only) + each strict z-cutoff
defs <- list()
if (isTRUE(CFG$do_haplo)) defs <- c(defs, list(list(tag = "haplo", z = NULL)))
for (zc in CFG$z_cutoffs) defs <- c(defs, list(list(tag = paste0("z_", zc), z = zc)))

n_tested_subjects <- n_distinct(expr_z_join$subjectid)   # denominator for prevalence

for (d in defs) {
  u <- call_underliers(expr_z_join, d$z, CFG$logcpm_drop)
  write_rds(u, setwd_out(paste0("underliers_", d$tag, ".rds")))
  write_tsv(u, setwd_out(paste0("underliers_", d$tag, ".tsv.gz")))

  # END-PRODUCT #2: RNA outlier prevalence per gene (one row per gene), keyed by
  # symbol so it merges with the per-gene q^2 incidence table.
  prev <- u |>
    distinct(gene_id, subjectid) |>
    count(gene_id, name = "n_underlier_subjects") |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$"),
           n_subjects_tested = n_tested_subjects,
           rna_outlier_prevalence = n_underlier_subjects / n_tested_subjects) |>
    left_join(dplyr::select(gene_meta, gene_nv, symbol), by = "gene_nv") |>
    dplyr::select(gene_id, gene_nv, symbol, n_underlier_subjects,
                  n_subjects_tested, rna_outlier_prevalence) |>
    arrange(desc(rna_outlier_prevalence))
  write_tsv(prev, setwd_out(paste0("rna_outlier_prevalence_per_gene_", d$tag, ".tsv")))

  message(sprintf("  %s: %d (subject,gene) events across %d genes -> rna_outlier_prevalence_per_gene_%s.tsv",
                  d$tag, nrow(u), n_distinct(u$gene_id), d$tag))
}
message("DONE. Outputs in ", normalizePath(CFG$out_dir))
