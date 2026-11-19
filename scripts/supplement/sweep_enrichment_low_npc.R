#!/usr/bin/env Rscript
# Extend the RV/underlier enrichment sweep down to low expression-PC counts
# (5..55), which the saved sweeps never covered -- they start at 55.
#
# expr_scaled_flat_list.rds only holds matrices for npc >= 55, so the
# residualized matrices are recomputed here with residualize_phenotypes() copied
# verbatim from outlier_processing.R. npc = 55 is included as an anchor: its odds
# ratios must match expr_enrichments_all_cutoffs*.rds.
suppressPackageStartupMessages({
  library(tidyverse); library(magrittr); library(limma)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

# usage: Rscript sweep_enrichment_low_npc.R [--npcs 5,10,55] [--out FILE]
# Already-computed PC counts present in --out are skipped, so the sweep is
# resumable. Do not run two copies against the same --out at once.
args    <- commandArgs(trailingOnly = TRUE)
argval  <- function(flag, default) {
  i <- match(flag, args); if (is.na(i)) default else args[i + 1]
}
NPCS    <- as.numeric(strsplit(argval("--npcs", "5,10,15,20,25,30,35,40,45,50,55"), ",")[[1]])
OUT     <- argval("--out", "expr_enrichments_low_npc.rds")
CUTOFFS <- c(2, 3, 3.5, 4, 5)
N_GENO  <- 5
# expr_eigenvectors_df.rds is keyed by RNA ID (TOR*) while expr_qc_normalized.rds
# and the SV table are keyed by WGS ID (NWD*), so the expression PC table has to
# be re-keyed before it can be lined up with the expression matrix.
METADATA <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/metadata.filtered.txt"

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

# Same logic as outlier_enrichment(), but the long-format join is built once and
# reused across cutoffs (the join does not depend on the cutoff). The long table
# is materialized directly instead of via pivot_longer purely for speed; row
# order does not affect fisher.test.
outlier_enrichment_multi <- function(scaled_matrix, sv_subset, z_scores) {
    long <- tibble(gene_id = rep(rownames(scaled_matrix), times = ncol(scaled_matrix)),
                   sample  = rep(colnames(scaled_matrix), each = nrow(scaled_matrix)),
                   expression_zscore = as.vector(scaled_matrix))
    long$sample <- str_remove_all(long$sample, "_.*$")
    longer_filter <- filter(long, is_in(sample, unique(sv_subset$sample)))
    rm(long); gc(verbose = FALSE)

    sv_join <- left_join(longer_filter, sv_subset, by = c("gene_id", "sample"))
    rm(longer_filter); gc(verbose = FALSE)
    rv <- !is.na(sv_join$variant_id)

    map_dfr(z_scores, function(z) {
        ft <- function(hit, lab) {
            r <- fisher.test(hit, rv)
            tibble(odds_ratio = r$estimate, p_value = r$p.value,
                   ci_95_lo = r$conf.int[1], ci_95_hi = r$conf.int[2], direction = lab)
        }
        bind_rows(ft(sv_join$expression_zscore > z, "overlier"),
                  ft(sv_join$expression_zscore < -z, "underlier")) |>
            mutate(cutoff = z)
    })
}

message("loading inputs ...")
expr_qc_normalized <- read_rds("expr_qc_normalized.rds")
expr_eigenvectors_df <- read_rds("expr_eigenvectors_df.rds")
expr_genotype_eigenvectors_df <- read_rds("expr_genotype_eigenvectors_df.rds")
sv_subset <- read_rds("subset_genotypes_100000_0.01_select.rds")

rna_to_wgs <- read_tsv(METADATA, show_col_types = FALSE) |> (\(d) setNames(d$WGS, d$RNA))()
if (!any(expr_eigenvectors_df$ID %in% colnames(expr_qc_normalized))) {
    expr_eigenvectors_df$ID <- unname(rna_to_wgs[expr_eigenvectors_df$ID])
    message("re-keyed expression PC table from RNA to WGS IDs")
}
stopifnot(all(expr_eigenvectors_df$ID %in% colnames(expr_qc_normalized)),
          all(expr_genotype_eigenvectors_df$ID %in% colnames(expr_qc_normalized)))
message(sprintf("expression matrix: %d genes x %d samples",
                nrow(expr_qc_normalized), ncol(expr_qc_normalized)))

done <- if (file.exists(OUT)) read_rds(OUT) else NULL
for (npc in NPCS) {
    if (!is.null(done) && npc %in% done$npcs_format) { message("skip npc=", npc); next }
    t0 <- Sys.time()
    message("npc=", npc, ": residualizing ...")
    m <- residualize_phenotypes(npc, N_GENO, expr_qc_normalized,
                                expression_pcs = expr_eigenvectors_df,
                                genotype_pcs = expr_genotype_eigenvectors_df)
    message("npc=", npc, ": enrichment ...")
    res <- outlier_enrichment_multi(m, sv_subset, CUTOFFS) |>
        mutate(npcs = sprintf("genotype_pc_%d_expression_pc_%d", N_GENO, npc),
               npcs_format = npc)
    rm(m); gc(verbose = FALSE)
    done <- bind_rows(done, res)
    write_rds(done, OUT)
    message(sprintf("npc=%d done in %.1f min", npc,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    print(as.data.frame(filter(res, direction == "underlier") |>
                        dplyr::select(cutoff, odds_ratio, ci_95_lo, ci_95_hi)))
}

# --- anchor check: any PC count also present in a previously saved sweep must match ---
ref <- map_dfr(c("expr_enrichments_all_cutoffs_2_5.rds", "expr_enrichments_all_cutoffs.rds",
                 "expr_enrichments_low_npc.rds"),
               function(f) if (file.exists(f) && normalizePath(f) != normalizePath(OUT))
                   read_rds(f) else NULL) |>
    filter(direction == "underlier") |>
    transmute(npcs_format = as.numeric(npcs_format), cutoff = as.numeric(cutoff),
              or_saved = odds_ratio) |> distinct()
chk <- done |> filter(direction == "underlier") |>
    dplyr::select(npcs_format, cutoff, or_new = odds_ratio) |>
    inner_join(ref, by = c("npcs_format", "cutoff")) |>
    mutate(abs_diff = abs(or_new - or_saved))
cat(sprintf("\nanchor check: %d points shared with previously saved sweeps, max |diff| = %.3g\n",
            nrow(chk), if (nrow(chk)) max(chk$abs_diff) else NA))
if (nrow(chk)) print(as.data.frame(arrange(chk, npcs_format, cutoff)))
message("wrote ", OUT)
