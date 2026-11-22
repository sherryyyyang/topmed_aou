# Shared helpers for the RV/underlier enrichment PC sweeps.
#
# residualize_phenotypes() is copied verbatim from outlier_processing.R.
#
# The enrichment functions accumulate the 2x2 table from a logical matrix aligned
# to the expression matrix rather than pivoting 45M rows to long format. That is
# algebraically the same Fisher test that outlier_enrichment() runs; it is
# verified against the published rare-SV odds ratios wherever the sweeps overlap.
suppressPackageStartupMessages({ library(tidyverse); library(limma) })

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))

residualize_phenotypes <- function(n_expression_pcs, n_genotype_pcs, normalized_expression, original_metadata = NULL, expression_pcs, genotype_pcs) {
    expression_pc_names <- str_c("PC", 1:n_expression_pcs)
    expression_pc_subset <- dplyr::select(expression_pcs, ID, one_of(expression_pc_names))
    genotype_pc_names <- str_c("Genotype_PC", 1:n_genotype_pcs)
    genotype_pc_subset <- dplyr::select(genotype_pcs, ID, one_of(genotype_pc_names))
    combined_phenodata <- left_join(expression_pc_subset, genotype_pc_subset)
    combined_phenotype_matrix <- column_to_rownames(combined_phenodata, "ID") |> as.matrix()
    normalized_expression_ordered <- as.matrix(normalized_expression[,rownames(combined_phenotype_matrix)])
    expr_residuals <- removeBatchEffect(normalized_expression_ordered, covariates = combined_phenotype_matrix)
    expr_scaled <- t(expr_residuals) |> scale() |> t()
    expr_scaled
}

# Re-key the expression PC table from RNA IDs (TOR*) to WGS IDs (NWD*) if needed;
# expr_eigenvectors_df.rds is saved with RNA IDs while the expression matrix and
# the variant tables use WGS IDs.
rekey_expression_pcs <- function(expression_pcs, normalized_expression, metadata_tsv) {
    if (any(expression_pcs$ID %in% colnames(normalized_expression))) return(expression_pcs)
    mp <- read_tsv(metadata_tsv, show_col_types = FALSE) |> (\(d) setNames(d$WGS, d$RNA))()
    expression_pcs$ID <- unname(mp[expression_pcs$ID])
    stopifnot(all(expression_pcs$ID %in% colnames(normalized_expression)))
    message("re-keyed expression PC table from RNA to WGS IDs")
    expression_pcs
}

# variant_tab: gene_key + sample_key, one row per carrier gene/sample pair.
build_rv_matrix <- function(scaled_matrix, variant_tab, strip_gene_version) {
    samples <- str_remove_all(colnames(scaled_matrix), "_.*$")
    genes   <- if (strip_gene_version) str_remove(rownames(scaled_matrix), "\\..*$") else rownames(scaled_matrix)
    stopifnot(!any(duplicated(genes)), !any(duplicated(samples)))
    keep_col <- samples %in% variant_tab$sample_key
    m <- scaled_matrix[, keep_col, drop = FALSE]
    samples <- samples[keep_col]
    ri <- match(variant_tab$gene_key, genes)
    ci <- match(variant_tab$sample_key, samples)
    ok <- !is.na(ri) & !is.na(ci)
    rv <- matrix(FALSE, nrow(m), ncol(m))
    rv[cbind(ri[ok], ci[ok])] <- TRUE
    list(m = m, rv = rv, samples = samples,
         n_pairs_variant = sum(ok), n_samples = ncol(m))
}

enrich_from_matrix <- function(m, rv, z_scores) {
    stopifnot(!anyNA(m))
    rvv <- as.vector(rv); mv <- as.vector(m)
    n_rv <- sum(rvv); n_tot <- length(rvv)
    map_dfr(z_scores, function(z) {
        one <- function(hit, lab) {
            n_hit <- sum(hit); n_hit_rv <- sum(hit & rvv)
            tab <- matrix(c(n_tot - n_hit - (n_rv - n_hit_rv), n_hit - n_hit_rv,
                            n_rv - n_hit_rv, n_hit_rv), nrow = 2)
            r <- fisher.test(tab)
            tibble(odds_ratio = unname(r$estimate), p_value = r$p.value,
                   ci_95_lo = r$conf.int[1], ci_95_hi = r$conf.int[2],
                   direction = lab, n_hit = n_hit, n_hit_variant = n_hit_rv)
        }
        bind_rows(one(mv > z, "overlier"), one(mv < -z, "underlier")) |> mutate(cutoff = z)
    })
}

# Per-subject 2x2 contributions, so a bootstrap over subjects is just resampling
# and summing these vectors. Columns of m and rv must be subjects.
per_sample_counts <- function(m, rv, z) {
    hit <- m < -z
    tibble(n_hit_rv = colSums(hit & rv), n_hit = colSums(hit),
           n_rv = colSums(rv), n_tot = nrow(m))
}

or_from_counts <- function(n_hit_rv, n_hit, n_rv, n_tot) {
    a <- n_hit_rv; b <- n_hit - a; c <- n_rv - a
    d <- n_tot - n_hit - n_rv + a
    (a * d) / (b * c)
}
