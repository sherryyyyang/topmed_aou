#!/usr/bin/env Rscript
# Are the odds ratios at different PC counts actually statistically different?
#
# For each variant class and z-cutoff this compares three PC settings pairwise:
#   peak  = the PC count with the largest odds ratio
#   2005  = the Gavish-Donoho rank estimate
#   55    = the value used in the released analysis
#
# Two answers are reported:
#  1. Whether the two Fisher 95% CIs overlap. This is what "do the CIs overlap"
#     usually means, but it is NOT a test of the difference here: both odds ratios
#     are computed on the same subjects and the same gene x sample pairs, so the
#     two estimates are strongly positively correlated and comparing marginal CIs
#     is far too conservative.
#  2. A paired nonparametric bootstrap over subjects (subjects are the independent
#     unit; gene x sample pairs within a subject are not). Both PC settings are
#     recomputed on each resampled subject set, so the correlation between them is
#     preserved, and the CI is placed on the ratio of the two odds ratios.
suppressPackageStartupMessages({ library(tidyverse); library(scales); library(patchwork) })

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- SUP_PC
source(file.path(ROOT, "scripts/helpers/_enrichment_helpers.R"))

B_BOOT   <- 2000
GD_NPC   <- 2005
USED_NPC <- 55
N_GENO   <- 5
METADATA <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/metadata.filtered.txt"
set.seed(1)

sweeps <- list(
  "Rare SV within 100 kb (AF <= 1%)" = list(
      rds = "expr_enrichments_sv_full.rds",
      tab = read_rds("subset_genotypes_100000_0.01_select.rds") |>
            transmute(gene_key = gene_id, sample_key = sample) |> distinct(),
      strip = FALSE),
  "LoF / PTV in gene (AF <= 0.5%)" = list(
      rds = "expr_enrichments_ptv.rds",
      tab = read_tsv("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/ptv_all.tsv",
                     show_col_types = FALSE) |>
            transmute(gene_key = str_remove(gene, "\\..*$"), sample_key = sample) |> distinct(),
      strip = TRUE))

message("loading expression inputs ...")
expr_qc_normalized <- read_rds("expr_qc_normalized.rds")
expr_eigenvectors_df <- rekey_expression_pcs(read_rds("expr_eigenvectors_df.rds"),
                                             expr_qc_normalized, METADATA)
expr_genotype_eigenvectors_df <- read_rds("expr_genotype_eigenvectors_df.rds")
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

results <- list()
for (vc in names(sweeps)) {
    sw <- sweeps[[vc]]
    fisher_tab <- read_rds(sw$rds) |> filter(direction == "underlier") |>
        transmute(cutoff = as.numeric(cutoff), npcs_format = as.numeric(npcs_format),
                  odds_ratio, ci_95_lo, ci_95_hi)
    cutoffs <- sort(unique(fisher_tab$cutoff))
    peaks <- fisher_tab |> group_by(cutoff) |> slice_max(odds_ratio, n = 1) |> ungroup()
    need_npc <- sort(unique(c(peaks$npcs_format, GD_NPC, USED_NPC)))
    message(sprintf("\n== %s ==\nPC settings needed: %s", vc, paste(need_npc, collapse = ", ")))

    # per-subject 2x2 contributions at each needed PC setting, for every cutoff
    counts <- list()
    n_sub <- NULL
    for (npc in need_npc) {
        b <- build_rv_matrix(get_matrix(npc), sw$tab, strip_gene_version = sw$strip)
        if (is.null(n_sub)) { n_sub <- b$n_samples; sub_ids <- b$samples }
        stopifnot(identical(b$samples, sub_ids))   # paired bootstrap needs identical subject sets
        for (z in cutoffs) counts[[paste(npc, z)]] <- as.matrix(per_sample_counts(b$m, b$rv, z))
        message(sprintf("  npc=%-5d per-subject counts done (%d subjects)", npc, b$n_samples))
    }

    # one shared set of bootstrap replicates for every comparison in this class
    W <- matrix(rmultinom(B_BOOT, n_sub, rep(1 / n_sub, n_sub)), nrow = n_sub)
    boot_or <- function(key) {
        A <- crossprod(W, counts[[key]])
        or_from_counts(A[, "n_hit_rv"], A[, "n_hit"], A[, "n_rv"], A[, "n_tot"])
    }
    point_or <- function(key) {
        s <- colSums(counts[[key]])
        or_from_counts(s["n_hit_rv"], s["n_hit"], s["n_rv"], s["n_tot"])
    }

    for (z in cutoffs) {
        pk <- filter(peaks, cutoff == z)$npcs_format
        cmp <- tribble(~comparison,        ~npc_a,   ~npc_b,
                       "peak vs 2005",    pk,       GD_NPC,
                       "55 vs 2005",      USED_NPC, GD_NPC,
                       "55 vs peak",      USED_NPC, pk)
        for (i in seq_len(nrow(cmp))) {
            ka <- paste(cmp$npc_a[i], z); kb <- paste(cmp$npc_b[i], z)
            if (cmp$npc_a[i] == cmp$npc_b[i]) next
            d <- log(boot_or(ka)) - log(boot_or(kb))
            fa <- filter(fisher_tab, cutoff == z, npcs_format == cmp$npc_a[i])
            fb <- filter(fisher_tab, cutoff == z, npcs_format == cmp$npc_b[i])
            results[[length(results) + 1]] <- tibble(
                variant_class = vc, cutoff = z, comparison = cmp$comparison[i],
                npc_a = cmp$npc_a[i], npc_b = cmp$npc_b[i],
                or_a = fa$odds_ratio, a_lo = fa$ci_95_lo, a_hi = fa$ci_95_hi,
                or_b = fb$odds_ratio, b_lo = fb$ci_95_lo, b_hi = fb$ci_95_hi,
                fisher_ci_overlap = fa$ci_95_lo <= fb$ci_95_hi & fb$ci_95_lo <= fa$ci_95_hi,
                or_ratio = point_or(ka) / point_or(kb),
                boot_lo = exp(quantile(d, 0.025, names = FALSE)),
                boot_hi = exp(quantile(d, 0.975, names = FALSE)),
                boot_p = max(2 * min(mean(d <= 0), mean(d >= 0)), 1 / B_BOOT))
        }
    }
}

res <- bind_rows(results) |>
    mutate(boot_signif = boot_lo > 1 | boot_hi < 1,
           cutoff_lab = sprintf("z < -%s", cutoff))
write_tsv(res, "enrichment_peak_vs_gd_bootstrap.tsv")

cat("\n=== Fisher CI overlap vs paired bootstrap ===\n")
print(as.data.frame(res |> transmute(variant_class = str_sub(variant_class, 1, 12), cutoff_lab,
    comparison, npc_a, npc_b, or_a = round(or_a, 2), or_b = round(or_b, 2),
    fisher_ci_overlap, or_ratio = round(or_ratio, 3),
    boot_ci = sprintf("[%.3f, %.3f]", boot_lo, boot_hi),
    boot_p = signif(boot_p, 3), boot_signif)))

pdat <- res |> mutate(cutoff_lab = factor(cutoff_lab, levels = rev(sprintf("z < -%s", sort(unique(res$cutoff))))),
                      comparison = factor(comparison, levels = c("peak vs 2005", "55 vs 2005", "55 vs peak")))
p <- ggplot(pdat, aes(or_ratio, cutoff_lab, color = boot_signif)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = boot_lo, xmax = boot_hi), height = 0.22, linewidth = 0.5) +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("%d vs %d", npc_a, npc_b)), vjust = -1, size = 2.5,
            color = "grey30", show.legend = FALSE) +
  facet_grid(comparison ~ variant_class, scales = "free_x") +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey45"),
                     labels = c(`TRUE` = "differ (bootstrap CI excludes 1)",
                                `FALSE` = "not distinguishable"), name = NULL) +
  scale_x_log10() +
  labs(x = "ratio of odds ratios (log scale), with paired subject-bootstrap 95% CI", y = NULL,
       title = "Is the enrichment optimum distinguishable from 55 or 2005 expression PCs?",
       subtitle = str_wrap(paste0(
         "Ratio of the two underlier odds ratios, with a 95% CI from ", comma(B_BOOT),
         " bootstrap resamples of subjects. Both PC settings are recomputed on the same ",
         "resampled subjects, so this accounts for the fact that the two estimates share data. ",
         "A CI excluding 1 means the two PC settings give genuinely different enrichment."), width = 118)) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, color = "grey30"),
        legend.position = "bottom", strip.text = element_text(size = 8.5))
ggsave(file.path(OUT_DIR, "enrichment_peak_vs_gd_bootstrap.pdf"), p, width = 9, height = 7, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "enrichment_peak_vs_gd_bootstrap.png"), p, width = 9, height = 7, dpi = 190)
message("\nwrote enrichment_peak_vs_gd_bootstrap.{tsv,pdf,png}")
