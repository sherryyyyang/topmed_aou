#!/usr/bin/env Rscript
# Complete the rare-SV sweep: the saved sweeps only ran z<-2 and z<-3 out to 1055
# expression PCs, so those two cutoffs cannot be compared at the Gavish-Donoho
# choice of 2005. This recomputes all five cutoffs over the full PC grid
# (5..3955) with the fast 2x2 path, and requires an exact match against the saved
# sweeps everywhere they overlap.
suppressPackageStartupMessages({ library(tidyverse); library(scales) })

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/_enrichment_helpers.R"))

CUTOFFS  <- c(2, 3, 3.5, 4, 5)
N_GENO   <- 5
NPCS_LOW <- seq(5, 50, 5)
OUT      <- "expr_enrichments_sv_full.rds"
METADATA <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/metadata.filtered.txt"

message("loading inputs ...")
expr_qc_normalized <- read_rds("expr_qc_normalized.rds")
expr_eigenvectors_df <- rekey_expression_pcs(read_rds("expr_eigenvectors_df.rds"),
                                             expr_qc_normalized, METADATA)
expr_genotype_eigenvectors_df <- read_rds("expr_genotype_eigenvectors_df.rds")
sv_tab <- read_rds("subset_genotypes_100000_0.01_select.rds") |>
    transmute(gene_key = gene_id, sample_key = sample) |> distinct()

message("loading expr_scaled_flat_list.rds (~28 GB) ...")
flat <- read_rds("expr_scaled_flat_list.rds")
npcs_available <- as.numeric(str_remove(names(flat), "^.*_expression_pc_"))
get_matrix <- function(npc) {
    i <- match(npc, npcs_available)
    if (!is.na(i)) return(flat[[i]])
    residualize_phenotypes(npc, N_GENO, expr_qc_normalized,
                           expression_pcs = expr_eigenvectors_df,
                           genotype_pcs = expr_genotype_eigenvectors_df)
}

npc_grid <- sort(unique(c(NPCS_LOW, npcs_available)))
done <- if (file.exists(OUT)) read_rds(OUT) else NULL
for (npc in npc_grid) {
    if (!is.null(done) && npc %in% done$npcs_format) next
    t0 <- Sys.time()
    b <- build_rv_matrix(get_matrix(npc), sv_tab, strip_gene_version = FALSE)
    res <- enrich_from_matrix(b$m, b$rv, CUTOFFS) |>
        mutate(npcs = sprintf("genotype_pc_%d_expression_pc_%d", N_GENO, npc),
               npcs_format = npc, variant_class = "SV",
               n_samples = b$n_samples, n_pairs_variant = b$n_pairs_variant)
    done <- bind_rows(done, res)
    write_rds(done, OUT)
    u <- filter(res, direction == "underlier")
    message(sprintf("[%s] npc=%-4d %.1f min  underlier OR  %s", format(Sys.time(), "%H:%M:%S"),
                    npc, as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    paste(sprintf("z<-%s: %.3f", u$cutoff, u$odds_ratio), collapse = "  ")))
}

saved <- bind_rows(read_rds(file.path(DATA_ENRICH, "expr_enrichments_all_cutoffs.rds")),
                   read_rds("expr_enrichments_all_cutoffs_2_5.rds"),
                   read_rds(file.path(DATA_ENRICH, "expr_enrichments_low_npc.rds"))) |>
    transmute(npcs_format = as.numeric(npcs_format), cutoff = as.numeric(cutoff),
              direction, or_saved = odds_ratio) |> distinct()
chk <- done |> dplyr::select(npcs_format, cutoff, direction, or_new = odds_ratio) |>
    inner_join(saved, by = c("npcs_format", "cutoff", "direction")) |>
    mutate(abs_diff = abs(or_new - or_saved))
cat(sprintf("\nagreement with saved SV sweeps: %d shared points, max |diff| = %.3g\n",
            nrow(chk), max(chk$abs_diff)))
stopifnot(nrow(chk) > 400, max(chk$abs_diff) < 1e-9)
cat(sprintf("new points added beyond the saved sweeps: %d\n",
            nrow(done) - nrow(chk)))
message("wrote ", OUT)
