#!/usr/bin/env Rscript
# Per-gene comparison: observed RNA-underlier % vs the LOEUF-derived baseline.
#
# Framing
# -------
# LOEUF is a gene-level tolerance metric (gnomAD UCL of obs/exp LoF count); it
# does not directly emit a per-subject "expected % underlier" number.
# The defensible per-gene baseline is the empirical mean RNA-underlier
# prevalence among genes that share the gene's LOEUF decile (and, secondarily,
# its LOEUF tier). The decile baseline is from
# underlier_loeuf_calibration_by_decile.tsv produced by plot_loeuf_calibration.R.
#
# For each call-out gene this script reports:
#   observed_pct      = n_underlier_subjects / n_subjects
#   expected_pct_dec  = decile baseline % per (gene, subject) pair
#                       = decile_prev_per_10k / 10000
#   expected_pct_tier = tier baseline %       (tier_prev_per_10k / 10000)
#   fold_decile       = observed_pct / expected_pct_dec
#   fold_tier         = observed_pct / expected_pct_tier
#   p_binom           = one-sided binomial p-value, H0: subject is underlier
#                       w.p. expected_pct_dec, observed = n in n_subjects trials
#
# Definition switch (env var DEFINITION_LABEL via _underlier_def.R):
#   strict  = z<-4 AND logCPM<gene_mean-1            (default)
#   haplo   = logCPM<gene_mean-1 (gene-specific z = -1/SD)
#
# Inputs:
#   loeuf_lookup.rds
#   expr_cpm_adj.rds
#   kept_subjects.rds
#   z_score_underlier_clean.rds  OR  z_score_underlier_haplo_clean.rds
#   underlier_loeuf_callouts.tsv     (or _haplo)        gene list to compare
#   underlier_loeuf_calibration_by_tier.tsv   (or _haplo) tier baselines
#   underlier_loeuf_calibration_by_decile.tsv (or _haplo) decile baselines
#
# Outputs:
#   underlier_loeuf_per_gene_compare(_haplo).tsv  full per-gene table
#   loeuf_per_gene_compare(_haplo).png / .pdf     lollipop figure

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
OUT_DIR <- FIG2X

# Log10 ticks covering [lo, hi] (positive proportions). Returns numeric breaks.
log_pct_breaks <- function(lo, hi, n_target = 8) {
    if (!is.finite(lo) || !is.finite(hi) || lo <= 0 || hi <= 0 || hi < lo) {
        return(c(lo, hi))
    }
    lo10 <- floor(log10(lo))
    hi10 <- ceiling(log10(hi))
    cand <- 10^seq(lo10, hi10)
    mult <- c(1, 2, 5)
    br <- sort(unique(as.vector(outer(cand, mult, `*`))))
    br <- br[br >= lo * 0.5 & br <= hi * 2]
    if (length(br) > n_target) {
        idx <- round(seq(1, length(br), length.out = n_target))
        br <- unique(br[idx])
    }
    br
}

source(file.path(ROOT, "scripts/helpers/_underlier_def.R"))
DEF <- resolve_underlier_def()
cat("Underlier definition: ", DEF$short, " ", DEF$tag, "\n", sep = "")

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_COLORS <- c(
    "Constrained" = "#B2182B",
    "Mid-low"     = "#EF8A62",
    "Mid-high"    = "#67A9CF",
    "Tolerant"    = "#2166AC",
    "No LOEUF"    = "grey60"
)

cat("Loading inputs ...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_subjects <- length(kept_subjects)
expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
all_eligible_genes <- rownames(expr_cpm_adj)
rm(expr_cpm_adj); gc(verbose = FALSE)

loeuf      <- read_rds(LOEUF_LOOKUP)
under      <- read_rds(DEF$rds_underlier)
callouts   <- read_tsv(DEF$callouts_tsv,  show_col_types = FALSE)
tier_tbl   <- read_tsv(DEF$tier_tsv,      show_col_types = FALSE)
decile_tbl <- read_tsv(DEF$decile_tsv,    show_col_types = FALSE)

# Per-gene n (over the FULL eligible-gene set, so genes with zero underliers
# pick up n=0 -- here we only need it for the call-out subset, but we keep
# the same logic for consistency).
under_per_gene <- under |>
    group_by(gene_id) |>
    tally(name = "n") |>
    ungroup()

# Recompute deciles consistently with plot_loeuf_calibration.R so each gene
# gets its own decile assignment (callouts CSV does not store this).
gene_df <- tibble(gene_id = all_eligible_genes) |>
    mutate(gene_id_no_version = str_remove(gene_id, "\\.\\d+$")) |>
    left_join(under_per_gene, by = "gene_id") |>
    mutate(n = replace_na(n, 0L))

loeuf_full <- loeuf |>
    distinct(gene_id, .keep_all = TRUE) |>
    select(gene_id, LOEUF_full = LOEUF, tier_full = LOEUF_tier,
           symbol_full = symbol)
loeuf_nv <- loeuf |>
    distinct(gene_id_no_version, .keep_all = TRUE) |>
    select(gene_id_no_version, LOEUF_nv = LOEUF, tier_nv = LOEUF_tier,
           symbol_nv = symbol)
gene_df <- gene_df |>
    left_join(loeuf_full, by = "gene_id") |>
    left_join(loeuf_nv,   by = "gene_id_no_version") |>
    mutate(
        LOEUF      = coalesce(LOEUF_full, LOEUF_nv),
        symbol     = coalesce(symbol_full, symbol_nv),
        LOEUF_tier = factor(
            ifelse(is.na(LOEUF), "No LOEUF",
                   as.character(coalesce(tier_full, tier_nv))),
            levels = TIER_LEVELS
        )
    ) |>
    select(gene_id, gene_id_no_version, symbol, LOEUF, LOEUF_tier, n)

# Decile assignment matching plot_loeuf_calibration.R.
gene_with_loeuf <- gene_df |> filter(!is.na(LOEUF))
decile_breaks <- quantile(gene_with_loeuf$LOEUF,
                          probs = seq(0, 1, by = 0.1), na.rm = TRUE)
gene_df <- gene_df |>
    mutate(LOEUF_decile = cut(LOEUF, breaks = decile_breaks,
                              include.lowest = TRUE,
                              labels = paste0("D", 1:10)))

# Restrict to call-out genes (top 20 Constrained + top 20 Tolerant).
compare_df <- callouts |>
    select(gene_id, symbol, LOEUF, LOEUF_tier) |>
    left_join(gene_df |> select(gene_id, n, LOEUF_decile),
              by = "gene_id") |>
    left_join(tier_tbl |> select(LOEUF_tier, tier_prev_per_10k = prev_per_10k),
              by = "LOEUF_tier") |>
    left_join(decile_tbl |> select(LOEUF_decile,
                                   decile_prev_per_10k = prev_per_10k,
                                   loeuf_median),
              by = "LOEUF_decile") |>
    mutate(
        observed_pct      = n / n_subjects,
        expected_pct_dec  = decile_prev_per_10k / 10000,
        expected_pct_tier = tier_prev_per_10k   / 10000,
        fold_decile       = observed_pct / expected_pct_dec,
        fold_tier         = observed_pct / expected_pct_tier,
        # One-sided binomial: P(X >= n | n_subjects, expected_pct_dec).
        p_binom = pbinom(q = n - 1L, size = n_subjects,
                         prob = expected_pct_dec, lower.tail = FALSE)
    ) |>
    arrange(LOEUF_tier, desc(observed_pct))

compare_path <- paste0("underlier_loeuf_per_gene_compare", DEF$suffix, ".tsv")
write_tsv(compare_df, compare_path)
cat("Wrote:", compare_path, " (", nrow(compare_df), "rows )\n")

cat("\nTop 5 Constrained call-outs (observed vs LOEUF-decile expected):\n")
print(compare_df |>
        filter(LOEUF_tier == "Constrained") |>
        head(5) |>
        select(symbol, LOEUF, LOEUF_decile, n, observed_pct,
               expected_pct_dec, fold_decile, p_binom))
cat("\nTop 5 Tolerant call-outs (observed vs LOEUF-decile expected):\n")
print(compare_df |>
        filter(LOEUF_tier == "Tolerant") |>
        head(5) |>
        select(symbol, LOEUF, LOEUF_decile, n, observed_pct,
               expected_pct_dec, fold_decile, p_binom))

# ---------- Lollipop figure ----------
# Y order: within each tier panel, sort by observed_pct desc.
plot_df <- compare_df |>
    mutate(
        gene_label = ifelse(!is.na(symbol) & symbol != "", symbol,
                            str_remove(gene_id, "\\.\\d+$")),
        tier_label = paste0(LOEUF_tier, " (LOEUF ",
                            ifelse(LOEUF_tier == "Constrained", "< 0.35",
                                   "\u2265 1.0"), ")")
    ) |>
    arrange(LOEUF_tier, observed_pct) |>
    group_by(LOEUF_tier) |>
    mutate(y = factor(gene_label, levels = gene_label)) |>
    ungroup()

facet_levels <- c(
    "Constrained" = "Constrained  (LOEUF < 0.35)",
    "Tolerant"    = "Tolerant  (LOEUF \u2265 1.0)"
)
plot_df$facet <- factor(facet_levels[as.character(plot_df$LOEUF_tier)],
                        levels = facet_levels)

# Long format for two points per gene (expected and observed).
long <- plot_df |>
    select(gene_id, gene_label, y, LOEUF_tier, facet,
           Observed = observed_pct,
           `Expected (LOEUF decile)` = expected_pct_dec,
           fold_decile, p_binom, n, LOEUF, LOEUF_decile) |>
    pivot_longer(c(Observed, `Expected (LOEUF decile)`),
                 names_to = "what", values_to = "pct") |>
    mutate(what = factor(what,
                         levels = c("Expected (LOEUF decile)", "Observed")))

# Floor for log axis: must accommodate both the smallest expected (decile-mean
# baseline, can be < 1/n_subjects) and the smallest observed.
all_pct <- c(plot_df$expected_pct_dec, plot_df$observed_pct)
all_pct <- all_pct[is.finite(all_pct) & all_pct > 0]
y_floor <- min(all_pct) * 0.75
y_ceiling <- max(all_pct) * 1.15
log_breaks <- log_pct_breaks(y_floor * 0.7, y_ceiling)
log_labels <- percent_format(accuracy = 0.01)(log_breaks)

p_main <- ggplot(long, aes(x = pct, y = y)) +
    geom_segment(data = plot_df,
                 aes(x = expected_pct_dec, xend = observed_pct,
                     y = y, yend = y),
                 color = "grey60", linewidth = 0.4) +
    geom_point(aes(shape = what, color = LOEUF_tier),
               size = 2.6, stroke = 0.9) +
    geom_text(data = plot_df,
              aes(x = observed_pct, y = y,
                  label = sprintf("  %.1fx (n=%d)", fold_decile, n)),
              hjust = 0, size = 2.6, color = "grey25") +
    scale_x_log10(breaks = log_breaks, labels = log_labels,
                  limits = c(y_floor * 0.7, y_ceiling),
                  expand = expansion(mult = c(0.02, 0.08))) +
    scale_shape_manual(values = c("Expected (LOEUF decile)" = 1,
                                  "Observed" = 19),
                       name = NULL) +
    scale_color_manual(values = TIER_COLORS, drop = FALSE, guide = "none") +
    guides(shape = guide_legend(override.aes = list(color = "grey25"))) +
    facet_wrap(~ facet, scales = "free_y", ncol = 2) +
    labs(title    = sprintf("Per-gene observed RNA-underlier %% vs LOEUF-decile expected baseline  %s",
                             DEF$tag),
         subtitle = sprintf("Underlier definition: %s.  Hollow circle = LOEUF-decile mean prevalence; filled = observed pct_underlier; label = fold-enrichment + n.",
                             DEF$long),
         x        = "% of subjects (log scale)",
         y        = NULL,
         caption  = sprintf("n_subjects (kept) = %s; expected = decile_prev_per_10k / 10,000 from %s.",
                            comma(n_subjects), DEF$decile_tsv)) +
    theme_classic(base_size = 11) +
    theme(legend.position    = "top",
          legend.box.spacing = unit(0.1, "lines"),
          plot.title.position = "plot",
          strip.text         = element_text(face = "bold", size = 11),
          strip.background   = element_blank(),
          axis.text.y        = element_text(face = "bold", size = 9),
          panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2))

out_png <- paste0("loeuf_per_gene_compare", DEF$suffix, ".png")
out_pdf <- paste0("loeuf_per_gene_compare", DEF$suffix, ".pdf")
ggsave(file.path(OUT_DIR, out_png), p_main, width = 13, height = 7.5, dpi = 200)
ggsave(file.path(OUT_DIR, out_pdf), p_main, width = 13, height = 7.5, device = cairo_pdf)
cat("\nWrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")
