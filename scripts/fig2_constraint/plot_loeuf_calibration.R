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
#        x    = LOEUF decile (D1..D10; lower decile = more constrained)
#        y    = prevalence per 10,000 (gene, subject) pairs, with 95% Wilson CI
#        secondary text: median LOEUF per decile
#
# Underlier definition matches outlier_processing.R lines 318-322:
#   per (gene, subject):  expression_zscore < -4
#                    AND  expression_logcpm < gene_mean - 1
#
# Inputs:
#   z_score_underlier_clean.rds   pre-filtered (gene_id, subjectid, ...) AFTER
#                                 the apply_extreme_subject_filter.R QC.
#   kept_subjects.rds             vector of kept subjectids
#   expr_cpm_adj.rds              only used for rownames() = gene universe
#   loeuf_lookup.rds              from build_loeuf_lookup.R
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
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG2A

source(file.path(ROOT, "scripts/helpers/_underlier_def.R"))
DEF <- resolve_underlier_def()
cat("Underlier definition: ", DEF$short, " ", DEF$tag, "\n", sep = "")

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_LEGEND_LABELS <- c(
    "Constrained" = "Constrained  (LOEUF < 0.35)",
    "Mid-low"     = "Mid-low  (0.35 \u2264 LOEUF < 0.6)",
    "Mid-high"    = "Mid-high  (0.6 \u2264 LOEUF < 1.0)",
    "Tolerant"    = "Tolerant  (LOEUF \u2265 1.0)",
    "No LOEUF"    = "No LOEUF"
)
LEGEND_ORDER <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")

# Wilson 95% CI on a binomial proportion. Returns p, lo, hi.
wilson_ci <- function(k, n, conf = 0.95) {
    z <- qnorm(1 - (1 - conf) / 2)
    p <- ifelse(n > 0, k / n, NA_real_)
    denom <- 1 + z^2 / n
    centre <- (p + z^2 / (2 * n)) / denom
    half   <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
    tibble(p = p, lo = pmax(centre - half, 0), hi = pmin(centre + half, 1))
}

cat("Loading kept_subjects.rds and expr_cpm_adj.rds ...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_subjects <- length(kept_subjects)
expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
all_eligible_genes <- rownames(expr_cpm_adj)
n_genes <- length(all_eligible_genes)
rm(expr_cpm_adj); gc(verbose = FALSE)
cat("  kept subjects =", n_subjects, "; total eligible genes =", n_genes, "\n")

cat("Loading ", DEF$rds_underlier, " ...\n", sep = "")
z_score_underlier <- read_rds(DEF$rds_underlier)
cat("  ", nrow(z_score_underlier), "underlier (gene, subject) rows after QC filter (",
    DEF$short, ")\n", sep = "")

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
write_tsv(tier_tbl, DEF$tier_tsv)
cat("Wrote:", DEF$tier_tsv, "\n")

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
write_tsv(decile_tbl, DEF$decile_tsv)
cat("Wrote:", DEF$decile_tsv, "\n")

# ---------- Panel A: per-tier forest ----------
tier_for_plot <- tier_tbl |>
    mutate(LOEUF_tier = factor(LOEUF_tier, levels = rev(TIER_LEVELS)),
           label = sprintf("%.2f (%.2f-%.2f)",
                           prev_per_10k, prev_per_10k_lo, prev_per_10k_hi))

# The published 55-PC z<-3 panel uses a hand-tuned window starting at 5 per
# 10k. scale_x_continuous(limits=) DROPS rows outside the window rather than
# zooming, so that window is only safe when the data actually lie inside it:
# at 2005 PCs every tier falls below 5 and three of five rows vanished, and at
# z<-4 all five vanished. Keep the hand-tuned axis only where it was tuned to
# fit; otherwise derive it from the data.
TUNED_LIMITS <- c(5, NA)
TUNED_BREAKS <- c(5, 10, 15)
data_lo <- min(tier_for_plot$prev_per_10k_lo, na.rm = TRUE)
data_hi <- max(tier_for_plot$prev_per_10k_hi, na.rm = TRUE)
if (data_lo >= TUNED_LIMITS[1] && data_hi <= max(TUNED_BREAKS)) {
    X_TIER_LIMITS <- TUNED_LIMITS
    X_TIER_BREAKS <- TUNED_BREAKS
} else {
    axis_lo <- floor(data_lo)
    axis_hi <- ceiling(data_hi)
    X_TIER_LIMITS <- c(axis_lo, NA)
    X_TIER_BREAKS <- scales::breaks_pretty(4)(c(axis_lo, axis_hi))
    cat("Tier axis: data (", signif(data_lo, 3), "-", signif(data_hi, 3),
        ") outside the hand-tuned 5-15 window; using a data-derived axis from ",
        axis_lo, "\n", sep = "")
}
# Guard against a future setting silently losing rows again.
dropped <- tier_for_plot |>
    filter(prev_per_10k_lo < X_TIER_LIMITS[1] |
           (!is.na(X_TIER_LIMITS[2]) & prev_per_10k_hi > X_TIER_LIMITS[2]))
if (nrow(dropped) > 0) {
    stop("Tier forest x-limits would drop ", nrow(dropped), " tier(s): ",
         paste(as.character(dropped$LOEUF_tier), collapse = ", "))
}

p_tier <- ggplot(tier_for_plot,
                 aes(y = LOEUF_tier, x = prev_per_10k,
                     color = LOEUF_tier)) +
    geom_errorbar(aes(xmin = prev_per_10k_lo, xmax = prev_per_10k_hi),
                  width = 0.18, linewidth = 0.5,
                  orientation = "y") +
    geom_point(size = 3) +
    geom_text(aes(label = label),
              nudge_y = 0.28, size = 2.6, color = "grey25",
              family = PAPER_FONT, show.legend = FALSE) +
    scale_color_manual(
        values = TIER_COLORS,
        breaks = LEGEND_ORDER,
        labels = TIER_LEGEND_LABELS[LEGEND_ORDER],
        drop = FALSE,
        name = "LOEUF tier"
    ) +
    scale_x_continuous(
        limits = X_TIER_LIMITS,
        breaks = X_TIER_BREAKS,
        expand = expansion(mult = c(0.04, 0.36))
    ) +
    labs(x = "Prevalence per 10,000 pairs", y = NULL,
         title = "By LOEUF tier") +
    coord_cartesian(clip = "off") +
    theme_paper() +
    theme_paper_inner_title() +
    theme(
        legend.position = "none",
        axis.text.y     = element_text(face = "bold"),
        plot.margin     = margin(6, 10, 8, 4)
    )

# ---------- Panel B: per-decile trend ----------
p_decile <- ggplot(decile_tbl,
                   aes(x = LOEUF_decile, y = prev_per_10k, group = 1)) +
    geom_ribbon(aes(ymin = prev_per_10k_lo, ymax = prev_per_10k_hi),
                fill = "grey75", alpha = 0.45) +
    geom_line(color = COL_TOPMED, linewidth = 0.7) +
    geom_point(color = COL_TOPMED, size = 2.4) +
    geom_text(aes(label = sprintf("%.2f", loeuf_median)),
              vjust = -0.7, size = 2.3, color = "grey30", family = PAPER_FONT) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.16))) +
    labs(x = "LOEUF decile (gene-level)",
         y = "Prevalence per 10,000 pairs",
         title = "By LOEUF decile") +
    theme_paper() +
    theme_paper_inner_title() +
    theme(plot.margin = margin(6, 6, 8, 8))

.paper_plots <- list(a = p_tier, b = p_decile)
combined <- (p_tier | p_decile) + plot_layout(widths = c(1.18, 1))
out_stem <- file.path(OUT_DIR, paste0("loeuf_calibration", DEF$suffix))
paper_ggsave(combined, out_stem, width = NATURE_DOUBLE_IN, height = 3.4)
if (!is_paper_assemble()) cat("\nWrote:", out_stem, "\n")
