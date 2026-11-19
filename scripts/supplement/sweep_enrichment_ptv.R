#!/usr/bin/env Rscript
# Same PC sweep as the rare-SV enrichment, but the variant side is LoF/PTVs
# (LOFTEE high-confidence, AF <= 0.5%) from ptv_all.tsv instead of rare SVs
# within 100 kb of the gene.
#
# Unlike the SV table, a PTV is inside its gene, so there is no flanking window:
# a gene x sample pair is "variant positive" if that sample carries a PTV in that
# gene.
#
# Speed note: the odds ratio is computed from 2x2 counts accumulated on a logical
# matrix aligned to the expression matrix, instead of pivoting 45M rows to long
# format and joining. That is algebraically the same test; it is verified against
# the published rare-SV odds ratios in the SV validation step below before any
# PTV numbers are produced.
suppressPackageStartupMessages({
  library(tidyverse); library(limma); library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

CUTOFFS   <- c(2, 3, 3.5, 4, 5)
N_GENO    <- 5
NPCS_LOW  <- seq(5, 50, 5)          # not in expr_scaled_flat_list.rds; recomputed
OUT       <- "expr_enrichments_ptv.rds"
PTV_TSV   <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/ptv_all.tsv"
METADATA  <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/metadata.filtered.txt"

# --- verbatim from outlier_processing.R ---
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

# variant_tab: gene_key + sample_key, one row per carrier gene/sample pair.
# Returns a logical matrix with the same shape as the (column-subset) expression
# matrix, TRUE where that sample carries a qualifying variant in that gene.
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
    list(m = m, rv = rv, n_pairs_variant = sum(ok), n_samples = ncol(m))
}

enrich_from_matrix <- function(m, rv, z_scores) {
    stopifnot(!anyNA(m))
    rvv <- as.vector(rv)
    mv <- as.vector(m)
    n_rv <- sum(rvv)
    n_tot <- length(rvv)
    map_dfr(z_scores, function(z) {
        one <- function(hit, lab) {
            n_hit <- sum(hit); n_hit_rv <- sum(hit & rvv)
            # rows = outlier FALSE/TRUE, cols = variant FALSE/TRUE
            tab <- matrix(c(n_tot - n_hit - (n_rv - n_hit_rv), n_hit - n_hit_rv,
                            n_rv - n_hit_rv, n_hit_rv), nrow = 2)
            r <- fisher.test(tab)
            tibble(odds_ratio = unname(r$estimate), p_value = r$p.value,
                   ci_95_lo = r$conf.int[1], ci_95_hi = r$conf.int[2],
                   direction = lab, n_hit = n_hit, n_hit_variant = n_hit_rv)
        }
        bind_rows(one(mv > z, "overlier"),
                  one(mv < -z, "underlier")) |> mutate(cutoff = z)
    })
}

message("loading expression inputs ...")
expr_qc_normalized <- read_rds("expr_qc_normalized.rds")
expr_eigenvectors_df <- read_rds("expr_eigenvectors_df.rds")
expr_genotype_eigenvectors_df <- read_rds("expr_genotype_eigenvectors_df.rds")
rna_to_wgs <- read_tsv(METADATA, show_col_types = FALSE) |> (\(d) setNames(d$WGS, d$RNA))()
if (!any(expr_eigenvectors_df$ID %in% colnames(expr_qc_normalized))) {
    expr_eigenvectors_df$ID <- unname(rna_to_wgs[expr_eigenvectors_df$ID])
    message("re-keyed expression PC table from RNA to WGS IDs")
}

# --- variant tables ---
sv_tab <- read_rds("subset_genotypes_100000_0.01_select.rds") |>
    transmute(gene_key = gene_id, sample_key = sample) |> distinct()
ptv_tab <- read_tsv(PTV_TSV, show_col_types = FALSE) |>
    transmute(gene_key = str_remove(gene, "\\..*$"), sample_key = sample) |> distinct()
message(sprintf("SV table: %s carrier gene/sample pairs; PTV table: %s pairs (%s genes, %s samples)",
                comma(nrow(sv_tab)), comma(nrow(ptv_tab)),
                comma(n_distinct(ptv_tab$gene_key)), comma(n_distinct(ptv_tab$sample_key))))

message("loading expr_scaled_flat_list.rds (~28 GB, holds npc >= 55) ...")
flat <- read_rds("expr_scaled_flat_list.rds")
npcs_available <- as.numeric(str_remove(names(flat), "^.*_expression_pc_"))

get_matrix <- function(npc) {
    i <- match(npc, npcs_available)
    if (!is.na(i)) return(flat[[i]])
    residualize_phenotypes(npc, N_GENO, expr_qc_normalized,
                           expression_pcs = expr_eigenvectors_df,
                           genotype_pcs = expr_genotype_eigenvectors_df)
}

# --- step 1: verify the fast 2x2 path reproduces the published SV odds ratios ---
saved_sv <- bind_rows(read_rds(file.path(DATA_ENRICH, "expr_enrichments_all_cutoffs.rds")),
                      read_rds("expr_enrichments_all_cutoffs_2_5.rds")) |>
    transmute(npcs_format = as.numeric(npcs_format), cutoff = as.numeric(cutoff),
              direction, or_saved = odds_ratio) |> distinct()
val_npcs <- c(55, 2005)
message("validating fast path against saved SV sweep at npc = ", paste(val_npcs, collapse = ", "))
val <- map_dfr(val_npcs, function(npc) {
    b <- build_rv_matrix(get_matrix(npc), sv_tab, strip_gene_version = FALSE)
    enrich_from_matrix(b$m, b$rv, CUTOFFS) |> mutate(npcs_format = npc)
})
val_chk <- inner_join(val |> dplyr::select(npcs_format, cutoff, direction, or_new = odds_ratio),
                      saved_sv, by = c("npcs_format", "cutoff", "direction")) |>
    mutate(rel_diff = abs(or_new - or_saved) / or_saved)
cat("\nSV validation (fast 2x2 path vs saved sweep):\n")
print(as.data.frame(arrange(val_chk, npcs_format, cutoff, direction)))
stopifnot(nrow(val_chk) > 0, max(val_chk$rel_diff) < 1e-6)
message("fast path reproduces the saved SV odds ratios")

# --- step 2: PTV sweep over the same PC grid as the SV sweep ---
npc_grid <- sort(unique(c(NPCS_LOW, npcs_available)))
message("PTV sweep over ", length(npc_grid), " PC values: ", paste(range(npc_grid), collapse = " to "))
done <- if (file.exists(OUT)) read_rds(OUT) else NULL
for (npc in npc_grid) {
    if (!is.null(done) && npc %in% done$npcs_format) next
    t0 <- Sys.time()
    b <- build_rv_matrix(get_matrix(npc), ptv_tab, strip_gene_version = TRUE)
    res <- enrich_from_matrix(b$m, b$rv, CUTOFFS) |>
        mutate(npcs = sprintf("genotype_pc_%d_expression_pc_%d", N_GENO, npc),
               npcs_format = npc, variant_class = "PTV",
               n_samples = b$n_samples, n_pairs_variant = b$n_pairs_variant)
    done <- bind_rows(done, res)
    write_rds(done, OUT)
    u <- filter(res, direction == "underlier")
    message(sprintf("[%s] npc=%-4d %.1f min  underlier OR  %s", format(Sys.time(), "%H:%M:%S"),
                    npc, as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    paste(sprintf("z<-%s: %.3f", u$cutoff, u$odds_ratio), collapse = "  ")))
}
write_rds(val, "expr_enrichments_sv_fastpath_validation.rds")
message("wrote ", OUT)
