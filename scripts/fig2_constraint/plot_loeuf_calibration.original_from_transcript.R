#!/usr/bin/env Rscript
# Dedicated calibration figure: does the empirical RNA-underlier prevalence
# track LOEUF the way it should?
#
# Two panels:
#   A. Per-tier forest plot
#        rows = LOEUF tier (Constrained, Mid-low, Mid-high, Tolerant, No LOEUF)
#        x    = underlier prevalence per 10,000 (gene, subject) pairs
#        bars = 95% Wilson CI on the underlying proportion
#               (binomial proportion is an approximation -- pairs are not
#                independent, but this is the standard summary stat)
#   B. Per-decile trend
#        x    = LOEUF decile (D1 = most constrained ... D10 = most tolerant)
#        y    = prevalence per 10,000 (gene, subject) pairs, with 95% Wilson CI
#        secondary text: median LOEUF per decile
#
# Underlier definition matches outlier_processing.R lines 318-322:
#   per (gene, subject):  expression_zscore < -4
#                    AND  expression_logcpm < gene_mean - 1
#
# Inputs:
#   z_score_underlier.rds   pre-filtered (gene_id, subjectid, ...)
#   expr_cpm_adj.rds        only used for ncol() / rownames()
#   loeuf_lookup.rds        from build_loeuf_lookup.R
# Outputs:
#   loeuf_calibration.png / .pdf
#   underlier_loeuf_calibration_by_tier.tsv
#   underlier_loeuf_calibration_by_decile.tsv

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(ggplot2)
    library(scales)
    library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))


TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_COLORS <- c(
    "Constrained" = "#B2182B",
    "Mid-low"     = "#EF8A62",
    "Mid-high"    = "#67A9CF",
    "Tolerant"    = "#2166AC",
    "No LOEUF"    = "grey60"
)

# Wilson 95% CI on a binomial proportion. Returns p, lo, hi.
wilson_ci <- function(k, n, conf = 0.95) {
    z <- qnorm(1 - (1 - conf) / 2)
    p <- ifelse(n > 0, k / n, NA_real_)
    denom <- 1 + z^2 / n
    centre <- (p + z^2 / (2 * n)) / denom
    half   <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
    tibble(p = p, lo = pmax(centre - half, 0), hi = pmin(centre + half, 1))
}

cat("Loading expr_cpm_adj.rds for total subject/gene counts ...\n")
expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
n_subjects <- ncol(expr_cpm_adj)
all_eligible_genes <- rownames(expr_cpm_adj)
n_genes <- length(all_eligible_genes)
rm(expr_cpm_adj); gc(verbose = FALSE)
cat("  total subjects =", n_subjects, "; total eligible genes =", n_genes, "\n")

cat("Loading z_score_underlier.rds ...\n")
z_score_underlier <- read_rds("z_score_underlier.rds")
cat("  ", nrow(z_score_underlier), "underlier (gene, subject) rows\n")

cat("Loading loeuf_lookup.rds ...\n")
loeuf_lookup <- read_rds(LOEUF_LOOKUP)

# Per-gene underlier counts over the FULL eligible-gene set (genes with zero
# underliers contribute zero to the numerator, full n_subjects to the
# denominator).
under_per_gene <- z_score_underlier |>
    group_by(gene_id) |>
    tally(name = "n") |>
    ungroup()

gene_df <- tibble(gene_id = all_eligible_genes) |>
    mutate(gene_id_no_version = str_remove(gene_id, "\\.\\d+$")) |>
    left_join(under_per_gene, by = "gene_id") |>
    mutate(n = replace_na(n, 0L))

loeuf_full <- loeuf_lookup |>
    distinct(gene_id, .keep_all = TRUE) |>
    select(gene_id, LOEUF_full = LOEUF, tier_full = LOEUF_tier)
loeuf_nv <- loeuf_lookup |>
    distinct(gene_id_no_version, .keep_all = TRUE) |>
    select(gene_id_no_version, LOEUF_nv = LOEUF, tier_nv = LOEUF_tier)

gene_df <- gene_df |>
    left_join(loeuf_full, by = "gene_id") |>
    left_join(loeuf_nv,   by = "gene_id_no_version") |>
    mutate(
        LOEUF = coalesce(LOEUF_full, LOEUF_nv),
        LOEUF_tier = factor(
            ifelse(is.na(LOEUF), "No LOEUF",
                   as.character(coalesce(tier_full, tier_nv))),
            levels = TIER_LEVELS
        )
    ) |>
    select(gene_id, LOEUF, LOEUF_tier, n)

# ---------- Per-tier prevalence + Wilson CI ----------
tier_tbl <- gene_df |>
    group_by(LOEUF_tier) |>
    summarise(
        n_genes_in_tier   = n(),
        n_underlier_pairs = sum(n),
        n_total_pairs     = n() * n_subjects,
        .groups = "drop"
    ) |>
    arrange(LOEUF_tier)

ci_tier <- wilson_ci(tier_tbl$n_underlier_pairs, tier_tbl$n_total_pairs)
tier_tbl <- bind_cols(tier_tbl, ci_tier) |>
    mutate(
        prev_per_10k    = 10000 * p,
        prev_per_10k_lo = 10000 * lo,
        prev_per_10k_hi = 10000 * hi
    )

cat("\nPer-tier prevalence (Wilson 95% CI):\n")
print(tier_tbl |> select(LOEUF_tier, n_genes_in_tier, n_underlier_pairs,
                         n_total_pairs, prev_per_10k,
                         prev_per_10k_lo, prev_per_10k_hi))
write_tsv(tier_tbl, "underlier_loeuf_calibration_by_tier.tsv")
cat("Wrote: underlier_loeuf_calibration_by_tier.tsv\n")

# ---------- Per-decile prevalence + Wilson CI ----------
# Use only genes with LOEUF; deciles via qcut so each decile has ~the same
# number of genes.
gene_with_loeuf <- gene_df |> filter(!is.na(LOEUF)) |> arrange(LOEUF)
decile_breaks <- quantile(gene_with_loeuf$LOEUF,
                          probs = seq(0, 1, by = 0.1), na.rm = TRUE)
gene_with_loeuf <- gene_with_loeuf |>
    mutate(LOEUF_decile = cut(LOEUF, breaks = decile_breaks,
                              include.lowest = TRUE,
                              labels = paste0("D", 1:10)))

decile_tbl <- gene_with_loeuf |>
    group_by(LOEUF_decile) |>
    summarise(
        n_genes_in_decile  = n(),
        loeuf_min          = min(LOEUF),
        loeuf_max          = max(LOEUF),
        loeuf_median       = median(LOEUF),
        n_underlier_pairs  = sum(n),
        n_total_pairs      = n() * n_subjects,
        .groups = "drop"
    ) |>
    arrange(LOEUF_decile)

ci_decile <- wilson_ci(decile_tbl$n_underlier_pairs, decile_tbl$n_total_pairs)
decile_tbl <- bind_cols(decile_tbl, ci_decile) |>
    mutate(
        prev_per_10k    = 10000 * p,
        prev_per_10k_lo = 10000 * lo,
        prev_per_10k_hi = 10000 * hi
    )

cat("\nPer-decile prevalence (Wilson 95% CI):\n")
print(decile_tbl |> select(LOEUF_decile, n_genes_in_decile,
                           loeuf_min, loeuf_max, loeuf_median,
                           prev_per_10k, prev_per_10k_lo, prev_per_10k_hi))
write_tsv(decile_tbl, "underlier_loeuf_calibration_by_decile.tsv")
cat("Wrote: underlier_loeuf_calibration_by_decile.tsv\n")

# ---------- Panel A: per-tier forest ----------
tier_for_plot <- tier_tbl |>
    mutate(LOEUF_tier = factor(LOEUF_tier, levels = rev(TIER_LEVELS)),
           label = sprintf("%.2f (%.2f-%.2f)",
                           prev_per_10k, prev_per_10k_lo, prev_per_10k_hi))

p_tier <- ggplot(tier_for_plot,
                 aes(y = LOEUF_tier, x = prev_per_10k,
                     color = LOEUF_tier)) +
    geom_errorbarh(aes(xmin = prev_per_10k_lo, xmax = prev_per_10k_hi),
                   height = 0.18, linewidth = 0.5) +
    geom_point(size = 3) +
    geom_text(aes(label = label),
              nudge_y = 0.32, size = 2.9, color = "grey25",
              show.legend = FALSE) +
    scale_color_manual(values = TIER_COLORS, drop = FALSE) +
    scale_x_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.10))) +
    labs(title    = "A. Per-LOEUF-tier underlier prevalence",
         subtitle = "Observed RNA underliers per 10,000 (gene, subject) pairs (Wilson 95% CI)",
         x        = "Prevalence per 10,000 (gene, subject) pairs",
         y        = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title.position = "plot",
          axis.text.y = element_text(face = "bold"))

# ---------- Panel B: per-decile trend ----------
p_decile <- ggplot(decile_tbl,
                   aes(x = LOEUF_decile, y = prev_per_10k, group = 1)) +
    geom_ribbon(aes(ymin = prev_per_10k_lo, ymax = prev_per_10k_hi),
                fill = "grey75", alpha = 0.45) +
    geom_line(color = "#2166AC", linewidth = 0.7) +
    geom_point(color = "#2166AC", size = 2.4) +
    geom_text(aes(label = sprintf("LOEUF\nmed=%.2f", loeuf_median)),
              vjust = 1.7, size = 2.4, color = "grey30") +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.10))) +
    labs(title    = "B. Per-LOEUF-decile underlier prevalence",
         subtitle = "D1 = most constrained ... D10 = most tolerant; ribbon = Wilson 95% CI",
         x        = "LOEUF decile (gene-level)",
         y        = "Prevalence per 10,000 (gene, subject) pairs") +
    theme_classic(base_size = 11) +
    theme(plot.title.position = "plot")

# ---------- compose + save ----------
n_with_loeuf <- nrow(gene_with_loeuf)
caption_txt <- sprintf(
    "%s underlier (gene, subject) pairs across %s eligible genes (%s with LOEUF) and %s subjects.\nUnderlier per (gene, subject): expression_zscore < -4 AND expression_logcpm < gene_mean - 1.",
    comma(nrow(z_score_underlier)),
    comma(n_genes),
    comma(n_with_loeuf),
    comma(n_subjects)
)

combined <- (p_tier | p_decile) +
    plot_layout(widths = c(1, 1.15)) +
    plot_annotation(
        title    = "Calibration of RNA-underlier prevalence vs LOEUF",
        subtitle = "Question: does observed RNA-underlier prevalence track LOEUF the way it should?",
        caption  = caption_txt,
        theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                         plot.subtitle = element_text(size = 10, color = "grey25"),
                         plot.caption  = element_text(size = 9,  color = "grey35"))
    )

out_png <- "loeuf_calibration.png"
out_pdf <- "loeuf_calibration.pdf"
ggsave(file.path(OUT_DIR, out_png), combined, width = 12, height = 5.0, dpi = 200)
ggsave(file.path(OUT_DIR, out_pdf), combined, width = 12, height = 5.0, device = cairo_pdf)
cat("\nWrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")
