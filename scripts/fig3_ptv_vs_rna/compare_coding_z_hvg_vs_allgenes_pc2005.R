#!/usr/bin/env Rscript
# Compare residualized z-scores for the SAME GENCODE v30 protein-coding genes:
#   published 2005-PC coding-only HVG pipeline
#     vs
#   all-gene HVG pipeline residualized at n_expr_pcs=2005, n_geno_pcs=5
#
# Fair PC match: published 2005-PC z, NOT the 55-PC comparison.
# Does NOT overwrite:
#   Fig_coding_zscore_hvg_coding_vs_allgenes.{png,pdf}
#   coding_z_hvg_vs_allgenes_{summary,per_gene,cache}.*
#   published underlier folders
#
# All-gene 2005-PC z is computed here (not previously cached):
#   residualize_expression(qc_normalized, expr_pc_df, geno_pc_df, 2005, 5)
#   = removeBatchEffect(FT HVG, covariates=PC1:2005 + GenoPC1:5) then t(scale(t(resid)))
# Z / scale grid = all samples with both PC sets (all-gene QC, n=6024).
# Comparison subjects = intersect with published 5,560.
#
# Log2CPM arm (both sides): 55-PC-adjusted expr_cpm_adj.
# Underlier = z_2005 < -3 AND log2CPM_adj55 < gene_mean - 1.
# Gene mean is computed on the same 6024-sample z grid used for residualize /
# scale, then calls are intersected to the 5,560 shared subjects for the 2x2.
#
# Writes (new names only):
#   allgenes_pipeline_out/expr_zscore_pc2005.rds
#   allgenes_pipeline_out/Fig_coding_zscore_hvg_coding_vs_allgenes_pc2005.{png,pdf}
#   allgenes_pipeline_out/coding_z_hvg_vs_allgenes_pc2005_summary.tsv
#   allgenes_pipeline_out/coding_z_hvg_vs_allgenes_pc2005_per_gene.tsv
# and copies the figure (new name) into submission_package_nov_2026/.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Matrix)
  library(limma)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

OUT <- file.path(ROOT, "allgenes_pipeline_out")
DRAFT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
STEM <- "Fig_coding_zscore_hvg_coding_vs_allgenes_pc2005"
PUB_Z <- file.path(ROOT, "data/underliers_pc2005/expr_z_join.rds")
PUB_CPM <- file.path(ROOT, "data/expression/expr_cpm_adj.rds")
PUB_UNDER <- file.path(ROOT, "data/underliers_pc2005/z_score_underlier_z_-3_clean.rds")
NEW_CPM <- file.path(OUT, "expr_cpm_adj.rds")
GENE_META <- file.path(OUT, "gene_meta.rds")
QC_NORM <- file.path(OUT, "qc_normalized.rds")
EXPR_PC <- file.path(OUT, "expr_pc_df.rds")
GENO_PC <- file.path(OUT, "geno_pc_df.rds")
Z_CACHE <- file.path(OUT, "expr_zscore_pc2005.rds")
CACHE <- file.path(OUT, "coding_z_hvg_vs_allgenes_pc2005_cache.rds")
N_EXPR_PCS <- 2005L
N_GENO_PCS <- 5L
Z_CUT <- -3
LOGCPM_DROP <- 1
N_SCATTER <- 2e6L
DELTA_LOTS <- c(5L, 10L, 20L)
set.seed(1)

# Refuse to clobber the published 55-PC comparison artifacts.
.protected <- c(
  file.path(OUT, "Fig_coding_zscore_hvg_coding_vs_allgenes.png"),
  file.path(OUT, "Fig_coding_zscore_hvg_coding_vs_allgenes.pdf"),
  file.path(OUT, "coding_z_hvg_vs_allgenes_summary.tsv"),
  file.path(OUT, "coding_z_hvg_vs_allgenes_per_gene.tsv"),
  file.path(OUT, "coding_z_hvg_vs_allgenes_cache.rds"),
  file.path(DRAFT, "Fig_coding_zscore_hvg_coding_vs_allgenes.png"),
  file.path(DRAFT, "Fig_coding_zscore_hvg_coding_vs_allgenes.pdf")
)
.stop_if_protected <- function(path) {
  if (normalizePath(path, mustWork = FALSE) %in%
      vapply(.protected, normalizePath, character(1), mustWork = FALSE)) {
    stop("Refusing to overwrite 55-PC artifact: ", path)
  }
}

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

nv <- function(x) str_remove(x, "\\..*$")

# Match run_underlier_pipeline.R::residualize_expression exactly.
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
# All-gene 2005-PC z (shared coding genes only).
# --------------------------------------------------------------------------
if (file.exists(Z_CACHE)) {
  msg("Loading cached all-gene 2005-PC z ...")
  z_new_mat <- read_rds(Z_CACHE)
  msg("  z matrix: ", comma(nrow(z_new_mat)), " genes x ",
      comma(ncol(z_new_mat)), " samples")
} else {
  msg("Computing all-gene 2005-PC z (removeBatchEffect + per-gene scale) ...")
  msg("  axes from CPM / GENCODE (not the z tables) ...")
  pub_cpm0 <- read_rds(PUB_CPM)
  new_cpm0 <- read_rds(NEW_CPM)
  gm0 <- as.data.table(read_rds(GENE_META))
  pub_ids0 <- rownames(pub_cpm0)
  new_ids0 <- rownames(new_cpm0)
  coding_nv0 <- unique(gm0$gene_nv[gm0$gene_type == "protein_coding"])
  pub_nv0 <- nv(pub_ids0)
  new_nv0 <- nv(new_ids0)
  stopifnot(!anyDuplicated(pub_nv0), !anyDuplicated(new_nv0))
  shared_nv0 <- intersect(pub_nv0, intersect(new_nv0, coding_nv0))
  shared_ids0 <- pub_ids0[pub_nv0 %in% shared_nv0]
  new_shared0 <- new_ids0[new_nv0 %in% shared_nv0]
  if (!setequal(shared_ids0, new_shared0)) {
    stop("Shared coding genes have versioned ID mismatches between CPM matrices")
  }
  rm(pub_cpm0, new_cpm0, gm0)
  gc(verbose = FALSE)
  msg("  shared coding genes to residualize: ", comma(length(shared_ids0)))

  msg("  loading qc_normalized / expr_pc_df / geno_pc_df ...")
  qn <- read_rds(QC_NORM)
  expr_pc_df <- read_rds(EXPR_PC)
  geno_pc_df <- read_rds(GENO_PC)
  pc_need <- str_c("PC", seq_len(N_EXPR_PCS))
  gpc_need <- str_c("Genotype_PC", seq_len(N_GENO_PCS))
  if (!all(pc_need %in% names(expr_pc_df))) {
    miss_n <- sum(!pc_need %in% names(expr_pc_df))
    msg("  expr_pc_df missing ", miss_n, " of PC1..PC", N_EXPR_PCS,
        " — recomputing PCA (PCAtools, scale=TRUE) ...")
    if (!requireNamespace("PCAtools", quietly = TRUE)) {
      stop("PCAtools is required to recompute expression PCs")
    }
    qn_full <- as.matrix(qn)
    pca_res <- PCAtools::pca(qn_full, scale = TRUE)
    ev <- pca_res$rotated
    rownames(ev) <- colnames(qn_full)
    expr_pc_df <- as.data.frame(ev) |> rownames_to_column("ID")
    pca_cache <- file.path(OUT, "expr_pc_df_recomputed_for_pc2005.rds")
    .stop_if_protected(pca_cache)
    write_rds(expr_pc_df, pca_cache)
    msg("  wrote recomputed PCs: ", pca_cache)
    rm(qn_full, pca_res, ev)
    gc(verbose = FALSE)
  } else {
    msg("  expr_pc_df has PC1..PC", N_EXPR_PCS, " (n_PC_cols=",
        sum(grepl("^PC[0-9]+$", names(expr_pc_df))), ")")
  }
  if (!all(gpc_need %in% names(geno_pc_df))) {
    stop("geno_pc_df is missing Genotype_PC1..", N_GENO_PCS)
  }

  # dgeMatrix: materialize before ID checks so rownames() is a plain character.
  qn <- as.matrix(qn)
  qn_ids <- rownames(qn)
  if (is.null(qn_ids) || !length(qn_ids)) {
    stop("qc_normalized has no rownames after as.matrix()")
  }
  miss_g <- setdiff(shared_ids0, qn_ids)
  if (length(miss_g)) {
    stop("Shared coding genes missing from qc_normalized: ",
         paste(head(miss_g, 5), collapse = ", "))
  }
  qn_shared <- qn[shared_ids0, , drop = FALSE]
  rm(qn)
  gc(verbose = FALSE)
  msg("  FT matrix subset: ", nrow(qn_shared), " genes x ",
      ncol(qn_shared), " samples")

  msg("  residualize_expression n_expr_pcs=", N_EXPR_PCS,
      " n_geno_pcs=", N_GENO_PCS, " ...")
  z_new_mat <- residualize_expression(
    qn_shared, expr_pc_df, geno_pc_df, N_EXPR_PCS, N_GENO_PCS
  )
  rm(qn_shared, expr_pc_df, geno_pc_df)
  gc(verbose = FALSE)
  .stop_if_protected(Z_CACHE)
  write_rds(z_new_mat, Z_CACHE, compress = "gz")
  msg("  wrote ", Z_CACHE, "  ", nrow(z_new_mat), " x ", ncol(z_new_mat))
}

# --------------------------------------------------------------------------
if (file.exists(CACHE)) {
  msg("Loading cached pairwise comparison ...")
  .cache <- read_rds(CACHE)
  merged <- .cache$merged
  per_gene <- .cache$per_gene
  gm <- .cache$gm
  env_stats <- .cache$env_stats
  list2env(env_stats, envir = environment())
  rm(.cache)
  gc(verbose = FALSE)
  msg("  cache: ", comma(nrow(merged)), " pairs, ", comma(nrow(per_gene)), " genes")
} else {
  msg("Gene / subject axes (matrices + GENCODE types; not the z tables) ...")
  pub_cpm <- read_rds(PUB_CPM)
  new_cpm <- read_rds(NEW_CPM)
  gm <- as.data.table(read_rds(GENE_META))
  pub_ids <- rownames(pub_cpm)
  new_ids <- rownames(new_cpm)
  pub_subj <- colnames(pub_cpm)
  new_subj <- colnames(new_cpm)

  coding_nv <- unique(gm$gene_nv[gm$gene_type == "protein_coding"])
  pub_nv <- nv(pub_ids)
  new_nv <- nv(new_ids)
  stopifnot(!anyDuplicated(pub_nv), !anyDuplicated(new_nv))

  new_coding_nv <- intersect(new_nv, coding_nv)
  shared_nv <- intersect(pub_nv, new_coding_nv)
  pub_only_nv <- setdiff(pub_nv, new_nv)
  new_only_coding_nv <- setdiff(new_coding_nv, pub_nv)
  shared_ids <- pub_ids[pub_nv %in% shared_nv]
  shared_subj <- intersect(pub_subj, new_subj)

  new_shared_ids <- new_ids[new_nv %in% shared_nv]
  if (!setequal(shared_ids, new_shared_ids)) {
    stop("Shared coding genes have versioned ID mismatches between matrices")
  }
  if (!setequal(shared_ids, rownames(z_new_mat))) {
    stop("Cached 2005-PC z rownames do not match shared coding gene IDs")
  }

  n_shared_genes <- length(shared_nv)
  n_shared_subj <- length(shared_subj)
  n_pairs <- as.numeric(n_shared_genes) * n_shared_subj
  z_n_subj <- ncol(z_new_mat)
  z_subj <- colnames(z_new_mat)
  msg(sprintf(
    "  coding with z: shared=%s  published-only=%s  allgene-only=%s",
    comma(n_shared_genes), comma(length(pub_only_nv)), comma(length(new_only_coding_nv))
  ))
  msg(sprintf(
    "  subjects: shared=%s  published-only=%s  allgene-only=%s  pairs=%s",
    comma(n_shared_subj),
    comma(length(setdiff(pub_subj, new_subj))),
    comma(length(setdiff(new_subj, pub_subj))),
    comma(n_pairs)
  ))
  msg(sprintf(
    "  all-gene z/scale grid: %s samples; gene_mean grid for underliers: same %s",
    comma(z_n_subj), comma(z_n_subj)
  ))
  if (length(pub_only_nv)) {
    pub_only_lab <- gm[gene_nv %in% pub_only_nv, paste(symbol, collapse = ", ")]
    msg("  published-only coding HVG: ", pub_only_lab)
  }

  # Published 2005-PC z: 8,124 coding genes x 5,560 samples. Keep shared genes.
  msg("Loading published 2005-PC expr_z_join (then subset) ...")
  old <- as.data.table(read_rds(PUB_Z))
  old <- old[gene_id %chin% shared_ids,
             .(gene_id, subjectid, z_old = expression_zscore)]
  gc(verbose = FALSE)
  msg("  published rows kept: ", comma(nrow(old)))

  msg("Melting all-gene 2005-PC z on shared subjects ...")
  miss_s <- setdiff(shared_subj, z_subj)
  if (length(miss_s)) {
    stop("Shared subjects missing from all-gene 2005-PC z: ", length(miss_s))
  }
  z_dt <- as.data.table(z_new_mat[, shared_subj, drop = FALSE],
                        keep.rownames = "gene_id")
  new <- melt(z_dt, id.vars = "gene_id", variable.name = "subjectid",
              value.name = "z_new", variable.factor = FALSE)
  rm(z_dt)
  gc(verbose = FALSE)
  msg("  all-gene rows kept: ", comma(nrow(new)))

  msg("Inner-joining shared gene x subject pairs ...")
  setkey(old, gene_id, subjectid)
  setkey(new, gene_id, subjectid)
  merged <- old[new, nomatch = 0]
  rm(old, new)
  gc(verbose = FALSE)
  merged <- merged[is.finite(z_old) & is.finite(z_new)]
  msg("  finite shared pairs: ", comma(nrow(merged)))
  if (nrow(merged) == 0L) stop("No shared finite gene x subject pairs")

  msg("Pairwise and per-gene correlations ...")
  pearson_all <- cor(merged$z_old, merged$z_new, method = "pearson")
  spearman_all <- cor(merged$z_old, merged$z_new, method = "spearman")
  abs_dz <- abs(merged$z_old - merged$z_new)
  mae <- mean(abs_dz)
  med_abs <- median(abs_dz)
  rmse <- sqrt(mean((merged$z_old - merged$z_new)^2))
  frac_gt_0.5 <- mean(abs_dz > 0.5)
  frac_gt_1 <- mean(abs_dz > 1)
  frac_gt_2 <- mean(abs_dz > 2)
  rm(abs_dz)
  gc(verbose = FALSE)

  per_gene <- merged[, {
    n <- .N
    if (n < 10L) {
      .(n = n, pearson = NA_real_, spearman = NA_real_)
    } else {
      .(n = n,
        pearson = cor(z_old, z_new, method = "pearson"),
        spearman = cor(z_old, z_new, method = "spearman"))
    }
  }, by = gene_id]
  per_gene[, gene_nv := nv(gene_id)]
  r_ok <- per_gene$pearson[is.finite(per_gene$pearson)]
  r_q <- quantile(r_ok, probs = c(0.25, 0.5, 0.75), names = FALSE)
  msg(sprintf(
    "  overall Pearson=%.4f  Spearman=%.4f  MAE=%.4f  median|dz|=%.4f",
    pearson_all, spearman_all, mae, med_abs
  ))
  msg(sprintf(
    "  per-gene Pearson median=%.4f  IQR=[%.4f, %.4f]  n_genes=%s",
    r_q[2], r_q[1], r_q[3], comma(length(r_ok))
  ))

  msg("Underlier concordance on shared pairs (z_2005<-3 AND log2CPM_adj55 drop) ...")
  msg("  gene_mean computed on all-gene z grid (n=", comma(z_n_subj),
      "), not the 5560 shared subset")
  cpm_ids <- intersect(shared_ids, rownames(new_cpm))
  cpm_cols <- intersect(z_subj, colnames(new_cpm))
  if (!setequal(cpm_ids, shared_ids) || !setequal(cpm_cols, z_subj)) {
    stop("All-gene expr_cpm_adj does not cover the z gene x sample grid")
  }
  cpm_mat <- as.matrix(new_cpm[shared_ids, z_subj, drop = FALSE])
  z_align <- z_new_mat[shared_ids, z_subj, drop = FALSE]
  gene_mean <- rowMeans(cpm_mat)
  under_mat <- (z_align < Z_CUT) & (cpm_mat < (gene_mean - LOGCPM_DROP))
  storage.mode(under_mat) <- "logical"
  # Restrict the 2x2 to shared subjects (published underliers live on 5560).
  under_shared <- under_mat[, shared_subj, drop = FALSE]
  hit <- which(under_shared, arr.ind = TRUE)
  u_new <- data.table(
    gene_id = rownames(under_shared)[hit[, 1]],
    subjectid = colnames(under_shared)[hit[, 2]]
  )
  u_new <- unique(u_new)
  rm(cpm_mat, z_align, under_mat, under_shared, hit, new_cpm, pub_cpm)
  gc(verbose = FALSE)

  u_old <- as.data.table(read_rds(PUB_UNDER))
  u_old <- unique(u_old[gene_id %chin% shared_ids & subjectid %chin% shared_subj,
                        .(gene_id, subjectid)])
  u_old[, old_under := TRUE]
  u_new[, new_under := TRUE]
  setkey(u_old, gene_id, subjectid)
  setkey(u_new, gene_id, subjectid)
  u_join <- merge(u_old, u_new, all = TRUE)
  u_join[is.na(old_under), old_under := FALSE]
  u_join[is.na(new_under), new_under := FALSE]
  n_both <- as.integer(u_join[old_under == TRUE & new_under == TRUE, .N])
  n_old_only <- as.integer(u_join[old_under == TRUE & new_under == FALSE, .N])
  n_new_only <- as.integer(u_join[old_under == FALSE & new_under == TRUE, .N])
  n_old_u <- n_both + n_old_only
  n_new_u <- n_both + n_new_only
  n_neither <- nrow(merged) - n_both - n_old_only - n_new_only
  jaccard <- if ((n_both + n_old_only + n_new_only) > 0L) {
    n_both / (n_both + n_old_only + n_new_only)
  } else {
    NA_real_
  }

  n_old_g <- u_old[, .(n_under_old = .N), by = gene_id]
  n_new_g <- u_new[, .(n_under_new = .N), by = gene_id]
  per_gene <- merge(per_gene, n_old_g, by = "gene_id", all.x = TRUE)
  per_gene <- merge(per_gene, n_new_g, by = "gene_id", all.x = TRUE)
  per_gene[is.na(n_under_old), n_under_old := 0L]
  per_gene[is.na(n_under_new), n_under_new := 0L]
  per_gene[, delta_n := n_under_new - n_under_old]
  n_cor <- cor(per_gene$n_under_old, per_gene$n_under_new, method = "pearson")
  n_cor_s <- cor(per_gene$n_under_old, per_gene$n_under_new, method = "spearman")
  n_chg <- sapply(DELTA_LOTS, function(k) as.integer(sum(abs(per_gene$delta_n) >= k)))
  names(n_chg) <- paste0("n_genes_abs_delta_ge_", DELTA_LOTS)
  n_merged <- nrow(merged)

  msg(sprintf(
    "  2x2  both=%s  old-only=%s  new-only=%s  neither=%s  Jaccard=%.4f",
    comma(n_both), comma(n_old_only), comma(n_new_only), comma(n_neither), jaccard
  ))
  msg(sprintf(
    "  n_underlier per gene: Pearson=%.4f  Spearman=%.4f  |d|>=5: %s  |d|>=10: %s  |d|>=20: %s",
    n_cor, n_cor_s, comma(n_chg[1]), comma(n_chg[2]), comma(n_chg[3])
  ))

  env_stats <- list(
    pub_nv = pub_nv, new_coding_nv = new_coding_nv,
    n_shared_genes = n_shared_genes, pub_only_nv = pub_only_nv,
    new_only_coding_nv = new_only_coding_nv, pub_subj = pub_subj,
    new_subj = new_subj, n_shared_subj = n_shared_subj,
    pearson_all = pearson_all, spearman_all = spearman_all,
    mae = mae, med_abs = med_abs, rmse = rmse,
    frac_gt_0.5 = frac_gt_0.5, frac_gt_1 = frac_gt_1, frac_gt_2 = frac_gt_2,
    r_ok = r_ok, r_q = r_q, n_both = n_both, n_old_only = n_old_only,
    n_new_only = n_new_only, n_neither = n_neither, jaccard = jaccard,
    n_old_u = n_old_u, n_new_u = n_new_u, n_cor = n_cor, n_cor_s = n_cor_s,
    n_chg = n_chg, n_merged = nrow(merged),
    z_n_subj = z_n_subj, z_subj = z_subj
  )
  sc_keep <- merged[sample.int(nrow(merged), min(N_SCATTER, nrow(merged)))]
  .stop_if_protected(CACHE)
  write_rds(
    list(merged = sc_keep, per_gene = per_gene, gm = gm, env_stats = env_stats,
         n_merged = nrow(merged)),
    CACHE, compress = "gz"
  )
  msg("  wrote cache (scatter subsample + per-gene + stats)")
}

# --------------------------------------------------------------------------
msg("Writing summary / per-gene tables ...")
pub_only_sym <- gm[gene_nv %in% pub_only_nv, paste(sort(symbol), collapse = ";")]
s_pear <- quantile(per_gene$spearman, probs = c(0.25, 0.5, 0.75), na.rm = TRUE,
                   names = FALSE)
d_q <- quantile(per_gene$delta_n, probs = c(0.25, 0.5, 0.75), names = FALSE)
summary_dt <- data.table(
  metric = c(
    "n_published_coding_hvg",
    "n_allgene_coding_hvg",
    "n_shared_coding_genes",
    "n_published_only_coding",
    "n_allgene_only_coding",
    "published_only_symbols",
    "n_published_subjects",
    "n_allgene_subjects",
    "n_shared_subjects",
    "n_allgene_only_subjects",
    "n_allgene_z_scale_subjects",
    "underlier_gene_mean_n_subjects",
    "underlier_gene_mean_grid",
    "n_shared_finite_pairs",
    "pearson_all_pairs",
    "spearman_all_pairs",
    "mae_abs_dz",
    "median_abs_dz",
    "rmse_dz",
    "frac_pairs_abs_dz_gt_0.5",
    "frac_pairs_abs_dz_gt_1",
    "frac_pairs_abs_dz_gt_2",
    "per_gene_pearson_median",
    "per_gene_pearson_q25",
    "per_gene_pearson_q75",
    "per_gene_spearman_median",
    "per_gene_spearman_q25",
    "per_gene_spearman_q75",
    "n_genes_per_gene_cor",
    "underlier_both",
    "underlier_old_only",
    "underlier_new_only",
    "underlier_neither",
    "underlier_jaccard",
    "n_underlier_old_shared",
    "n_underlier_new_shared",
    "per_gene_n_underlier_pearson",
    "per_gene_n_underlier_spearman",
    "median_delta_n_underlier",
    "iqr_delta_n_underlier_q25",
    "iqr_delta_n_underlier_q75",
    names(n_chg),
    "pc_match",
    "underlier_definition"
  ),
  value = c(
    as.character(length(pub_nv)),
    as.character(length(new_coding_nv)),
    as.character(n_shared_genes),
    as.character(length(pub_only_nv)),
    as.character(length(new_only_coding_nv)),
    pub_only_sym,
    as.character(length(pub_subj)),
    as.character(length(new_subj)),
    as.character(n_shared_subj),
    as.character(length(setdiff(new_subj, pub_subj))),
    as.character(z_n_subj),
    as.character(z_n_subj),
    "allgene_qc_6024 (same samples as residualize_expression / scale); 2x2 then intersected to shared 5560",
    as.character(n_merged),
    sprintf("%.6f", pearson_all),
    sprintf("%.6f", spearman_all),
    sprintf("%.6f", mae),
    sprintf("%.6f", med_abs),
    sprintf("%.6f", rmse),
    sprintf("%.6f", frac_gt_0.5),
    sprintf("%.6f", frac_gt_1),
    sprintf("%.6f", frac_gt_2),
    sprintf("%.6f", r_q[2]),
    sprintf("%.6f", r_q[1]),
    sprintf("%.6f", r_q[3]),
    sprintf("%.6f", s_pear[2]),
    sprintf("%.6f", s_pear[1]),
    sprintf("%.6f", s_pear[3]),
    as.character(length(r_ok)),
    as.character(n_both),
    as.character(n_old_only),
    as.character(n_new_only),
    as.character(n_neither),
    sprintf("%.6f", jaccard),
    as.character(n_old_u),
    as.character(n_new_u),
    sprintf("%.6f", n_cor),
    sprintf("%.6f", n_cor_s),
    sprintf("%.6f", d_q[2]),
    sprintf("%.6f", d_q[1]),
    sprintf("%.6f", d_q[3]),
    as.character(n_chg),
    "published_2005PC_vs_allgene_2005PC",
    "z_2005 < -3 AND log2CPM_adj55 < gene_mean - 1"
  )
)
sum_path <- file.path(OUT, "coding_z_hvg_vs_allgenes_pc2005_summary.tsv")
gene_path <- file.path(OUT, "coding_z_hvg_vs_allgenes_pc2005_per_gene.tsv")
.stop_if_protected(sum_path)
.stop_if_protected(gene_path)
fwrite(summary_dt, sum_path, sep = "\t")
if (!"symbol" %in% names(per_gene)) {
  per_gene <- merge(per_gene, gm[, .(gene_nv, symbol, gene_type)],
                    by = "gene_nv", all.x = TRUE)
}
per_gene[, abs_delta := abs(delta_n)]
setorder(per_gene, -abs_delta, gene_nv)
per_gene[, abs_delta := NULL]
fwrite(
  per_gene[, .(gene_id, gene_nv, symbol, n, pearson, spearman,
               n_under_old, n_under_new, delta_n)],
  gene_path,
  sep = "\t"
)

# --------------------------------------------------------------------------
msg("Plotting sanity figure ...")
n_plot <- min(N_SCATTER, nrow(merged))
sc <- merged[sample.int(nrow(merged), n_plot)]
lim <- c(-8, 8)
ann <- sprintf(
  "Pearson r = %.3f\nSpearman ρ = %.3f\n%s pairs (plot: %s)",
  pearson_all, spearman_all, comma(n_merged), comma(n_plot)
)

p_sc <- ggplot(sc, aes(z_old, z_new)) +
  geom_bin2d(bins = 80, boundary = 0) +
  geom_abline(slope = 1, intercept = 0, colour = "white", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, colour = "grey20", linewidth = 0.35,
              linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10", option = "C", name = "pairs") +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  labs(
    x = "Published z (coding HVG, 2005 PCs)",
    y = "All-gene z (2005 PCs)"
  ) +
  annotate(
    "text", x = -7.6, y = 7.4, label = ann, hjust = 0, vjust = 1,
    size = 2.6, family = PAPER_FONT, lineheight = 0.95
  ) +
  theme_paper() +
  theme(legend.position = "right")

p_hist <- ggplot(per_gene[is.finite(pearson)], aes(x = pearson)) +
  geom_histogram(bins = 40, fill = unname(OI["blue"]), colour = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = r_q[2], colour = unname(OI["vermillion"]),
             linewidth = 0.6) +
  scale_x_continuous(limits = c(min(0.5, min(r_ok) - 0.02), 1.001),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Per-gene Pearson r (old z vs new z)",
    y = "Coding genes"
  ) +
  annotate(
    "text",
    x = r_q[2], y = Inf,
    label = sprintf("median = %.3f\nIQR [%s, %s]", r_q[2],
                    format(round(r_q[1], 3), nsmall = 3),
                    format(round(r_q[3], 3), nsmall = 3)),
    hjust = 1.05, vjust = 1.2, size = 2.6, family = PAPER_FONT,
    colour = unname(OI["vermillion"]), lineheight = 0.95
  ) +
  theme_paper()

p_n <- ggplot(per_gene, aes(n_under_old, n_under_new)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_point(alpha = 0.25, size = 0.7, colour = unname(OI["blue"]),
             stroke = 0) +
  labs(
    x = "n underliers, published 2005-PC (shared subjects)",
    y = "n underliers, all-gene 2005-PC"
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf(
      "r = %.3f   Jaccard (pairs) = %.3f\n|Δn|≥10: %s / %s genes",
      n_cor, jaccard, comma(n_chg[2]), comma(nrow(per_gene))
    ),
    hjust = -0.05, vjust = 1.2, size = 2.6, family = PAPER_FONT, lineheight = 0.95
  ) +
  theme_paper()

fig <- (p_sc | p_hist) / p_n +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = PAPER_FONT, face = "bold", size = 14, colour = "black"
    )
  )

fig_stem <- file.path(OUT, STEM)
.stop_if_protected(paste0(fig_stem, ".png"))
.stop_if_protected(paste0(fig_stem, ".pdf"))
paper_ggsave(fig, fig_stem, width = NATURE_DOUBLE_IN, height = 6.4)
if (dir.exists(DRAFT)) {
  draft_png <- file.path(DRAFT, paste0(STEM, ".png"))
  draft_pdf <- file.path(DRAFT, paste0(STEM, ".pdf"))
  .stop_if_protected(draft_png)
  .stop_if_protected(draft_pdf)
  file.copy(paste0(fig_stem, ".png"), draft_png, overwrite = TRUE)
  file.copy(paste0(fig_stem, ".pdf"), draft_pdf, overwrite = TRUE)
  src <- file.path(DRAFT, "source_panels")
  if (dir.exists(src)) {
    src_png <- file.path(src, paste0(STEM, ".png"))
    src_pdf <- file.path(src, paste0(STEM, ".pdf"))
    .stop_if_protected(src_png)
    .stop_if_protected(src_pdf)
    file.copy(paste0(fig_stem, ".png"), src_png, overwrite = TRUE)
    file.copy(paste0(fig_stem, ".pdf"), src_pdf, overwrite = TRUE)
  }
  draft_sum <- file.path(DRAFT, "coding_z_hvg_vs_allgenes_pc2005_summary.tsv")
  .stop_if_protected(draft_sum)
  file.copy(sum_path, draft_sum, overwrite = TRUE)
}

verdict <- {
  r <- pearson_all
  j <- jaccard
  if (is.finite(r) && r >= 0.95 && is.finite(j) && j >= 0.7) {
    paste0(
      "At 2005 expression PCs the coding-gene z-scores from the all-gene HVG ",
      "pipeline still track the published coding-only HVG z-scores very closely ",
      sprintf("(Pearson r=%.3f; per-gene median r=%.3f). ", r, r_q[2]),
      sprintf("Underlier calls (z<-3 and 55-PC log2CPM drop) overlap with Jaccard=%.3f. ", j),
      "Including noncoding HVGs when building PCs does not materially rewrite ",
      "the published coding underlier list at this PC setting."
    )
  } else if (is.finite(r) && r >= 0.85) {
    paste0(
      "At 2005 expression PCs the coding-gene z-scores remain highly correlated ",
      sprintf("(Pearson r=%.3f; per-gene median r=%.3f), ", r, r_q[2]),
      "but underlier membership is only partly shared ",
      sprintf("(Jaccard=%.3f; both=%s, published-only=%s, all-gene-only=%s). ",
              j, comma(n_both), comma(n_old_only), comma(n_new_only)),
      "The continuous scores agree; the z<-3 tail is more sensitive to whether ",
      "PCs were built from coding-only vs all-gene HVGs."
    )
  } else {
    paste0(
      "At 2005 expression PCs the all-gene HVG residualization diverges from ",
      "the published coding-only HVG z-scores ",
      sprintf("(Pearson r=%.3f; Jaccard=%.3f). ", r, j),
      "Noncoding HVGs in the PC basis change coding z-scores enough that ",
      "underlier calls should not be treated as interchangeable."
    )
  }
}

cat("\n==== SUMMARY ====\n")
print(summary_dt, nrows = 90)
cat("\nTop |Δn_underlier| genes:\n")
print(head(per_gene[, .(symbol, gene_nv, n_under_old, n_under_new, delta_n, pearson)], 15))
cat("\n==== VERDICT ====\n")
cat(verdict, "\n")
msg("DONE  figure ", paste0(fig_stem, ".png"))
