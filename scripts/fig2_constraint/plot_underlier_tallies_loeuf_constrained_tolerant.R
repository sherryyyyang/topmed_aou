#!/usr/bin/env Rscript
# Two-tier variant of plot_underlier_tallies_loeuf.R:
#   Same logic, but curves + inset boxplot only for **Constrained** and **Tolerant**
#   LOEUF tiers (extremes of constraint).
#
# Pick underlier definition via env: DEFINITION_LABEL=strict|z_-4|z_-3|z_-2|haplo
# Optional: UNDERLIER_WD=/path/to/dir with *_clean.rds, kept_subjects.rds, loeuf_lookup.rds
#
# Outputs:
#   underlier_tallies_loeuf{SUFFIX}_constrained_tolerant.png / .pdf
#   underlier_tallies_loeuf{SUFFIX}_constrained_tolerant_genes.tsv
#   underlier_tallies_loeuf{SUFFIX}_constrained_tolerant_subject_gene.tsv
#
# X-axis semantics (methods text): within each LOEUF tier, only genes with
# ≥1 RNA underlier are ordered by descending underlier rate; x = rank/N_genes_in_tier
# using the full tier gene count N_genes_in_tier (eligible genes assigned to tier).

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(ggplot2)
    library(scales)
    library(patchwork)
    library(ggrepel)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG2B
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
setwd(ROOT)
source(file.path(ROOT, "scripts/helpers/_underlier_def.R"))
DEF <- resolve_underlier_def()
cat("Underlier definition: ", DEF$short, " ", DEF$tag, "\n", sep = "")
cat("Working directory:", getwd(), "\n")

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_FOCUS <- c("Constrained", "Tolerant")
TIER_LABELS <- c(
    "Constrained" = "Constrained",
    "Mid-low"     = "Mid-low",
    "Mid-high"    = "Mid-high",
    "Tolerant"    = "Tolerant",
    "No LOEUF"    = "No LOEUF"
)

cat("Loading kept_subjects.rds and expr_cpm_adj.rds ...\n")
t0 <- Sys.time()
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_subjects <- length(kept_subjects)
expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
all_eligible_genes <- rownames(expr_cpm_adj)
n_genes <- length(all_eligible_genes)
rm(expr_cpm_adj); gc(verbose = FALSE)
cat(
    "  kept subjects =", n_subjects, "; total eligible genes =", n_genes,
    "  (",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    "s)\n"
)

cat("Loading ", DEF$rds_underlier, " ...\n", sep = "")
z_score_underlier <- read_rds(DEF$rds_underlier)

cat("Loading loeuf_lookup.rds ...\n")
loeuf_lookup <- read_rds(LOEUF_LOOKUP)

under_per_gene <- z_score_underlier |>
    group_by(gene_id) |>
    tally() |>
    ungroup()

gene_df <- tibble(gene_id = all_eligible_genes) |>
    mutate(gene_id_no_version = str_remove(gene_id, "\\.\\d+$")) |>
    left_join(under_per_gene, by = "gene_id") |>
    mutate(
        n             = replace_na(n, 0L),
        pct_underlier = n / n_subjects
    )

loeuf_full <- loeuf_lookup |>
    distinct(gene_id, .keep_all = TRUE) |>
    select(gene_id, LOEUF_full = LOEUF, tier_full = LOEUF_tier, symbol_full = symbol)
loeuf_nv <- loeuf_lookup |>
    distinct(gene_id_no_version, .keep_all = TRUE) |>
    select(gene_id_no_version, LOEUF_nv = LOEUF, tier_nv = LOEUF_tier, symbol_nv = symbol)

gene_df <- gene_df |>
    left_join(loeuf_full, by = "gene_id") |>
    left_join(loeuf_nv, by = "gene_id_no_version") |>
    mutate(
        LOEUF      = coalesce(LOEUF_full, LOEUF_nv),
        symbol     = coalesce(symbol_full, symbol_nv),
        LOEUF_tier = factor(
            ifelse(is.na(LOEUF), "No LOEUF", as.character(coalesce(tier_full, tier_nv))),
            levels = TIER_LEVELS
        )
    ) |>
    select(gene_id, gene_id_no_version, symbol, LOEUF, LOEUF_tier, n, pct_underlier)

prev_tbl <- gene_df |>
    group_by(LOEUF_tier) |>
    summarise(
        n_genes_in_tier        = n(),
        n_genes_with_underlier = sum(n > 0),
        n_underlier_pairs      = sum(n),
        n_total_pairs          = n() * n_subjects,
        prevalence_per_10k     = 10000 * n_underlier_pairs / n_total_pairs,
        median_pct_underlier   = median(pct_underlier),
        mean_pct_underlier     = mean(pct_underlier),
        max_pct_underlier      = max(pct_underlier),
        .groups = "drop"
    ) |>
    arrange(LOEUF_tier)

callouts <- gene_df |>
    filter(LOEUF_tier %in% TIER_FOCUS, n > 0) |>
    group_by(LOEUF_tier) |>
    arrange(desc(pct_underlier), .by_group = TRUE) |>
    slice_head(n = 6) |>
    ungroup()
# Do NOT pin extra immune genes (e.g. HLA-DQB2): under PC2005 HLA-DQB2 falls
# to ~rank 1186 / ~35% within-tier, which painted a lone callout mid-ladder.
# IL32 is already inside the Tolerant top-6 under PC2005.

tier_df <- gene_df |>
    filter(n > 0, LOEUF_tier %in% TIER_FOCUS) |>
    group_by(LOEUF_tier) |>
    arrange(desc(pct_underlier), .by_group = TRUE) |>
    mutate(rank_in_tier = row_number()) |>
    ungroup() |>
    left_join(
        prev_tbl |> select(LOEUF_tier, n_genes_in_tier),
        by = "LOEUF_tier"
    ) |>
    mutate(frac_tier_rank = rank_in_tier / n_genes_in_tier)

tier_df$LOEUF_tier <- factor(tier_df$LOEUF_tier, levels = TIER_FOCUS)
tier_label_vec <- TIER_LABELS[TIER_FOCUS]
tier_colors_focus <- TIER_COLORS[TIER_FOCUS]

# Gene-level table for downstream searches (all genes on the rank curves).
gene_z_summary <- z_score_underlier |>
    group_by(gene_id) |>
    summarise(
        n_underlier_pairs      = n(),
        expression_zscore_min  = min(expression_zscore, na.rm = TRUE),
        expression_zscore_mean = mean(expression_zscore, na.rm = TRUE),
        expression_logcpm_min  = min(expression_logcpm, na.rm = TRUE),
        expression_logcpm_mean = mean(expression_logcpm, na.rm = TRUE),
        .groups = "drop"
    )

export_genes <- tier_df |>
    left_join(gene_z_summary, by = "gene_id") |>
    mutate(
        pct_subjects_with_underlier = 100 * pct_underlier,
        gene_label = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id_no_version),
        constrained_vs_tolerant = tolower(as.character(LOEUF_tier))
    ) |>
    select(
        gene_label,
        gene_id,
        gene_id_no_version,
        symbol,
        constrained_vs_tolerant,
        LOEUF_tier,
        LOEUF,
        n_subjects_with_underlier = n,
        pct_subjects_with_underlier,
        rank_in_tier,
        frac_tier_rank,
        n_underlier_pairs,
        expression_zscore_min,
        expression_zscore_mean,
        expression_logcpm_min,
        expression_logcpm_mean
    ) |>
    arrange(LOEUF_tier, rank_in_tier)

out_tsv <- file.path(OUT_DIR, paste0("underlier_tallies_loeuf", DEF$suffix, "_constrained_tolerant_genes.tsv"))
write_tsv(export_genes, out_tsv)
cat("Wrote:", normalizePath(out_tsv), " (", nrow(export_genes), " genes)\n")

gene_meta <- export_genes |>
    select(
        gene_label, gene_id, gene_id_no_version, symbol,
        constrained_vs_tolerant, LOEUF_tier, LOEUF,
        n_subjects_with_underlier, pct_subjects_with_underlier,
        rank_in_tier, frac_tier_rank
    )

export_subject_gene <- z_score_underlier |>
    semi_join(tier_df |> select(gene_id), by = "gene_id") |>
    left_join(gene_meta, by = "gene_id") |>
    select(
        gene_label, gene_id, gene_id_no_version, symbol,
        constrained_vs_tolerant, LOEUF_tier, LOEUF,
        subjectid,
        expression_zscore, expression_logcpm,
        n_subjects_with_underlier, pct_subjects_with_underlier,
        rank_in_tier, frac_tier_rank
    ) |>
    arrange(constrained_vs_tolerant, rank_in_tier, gene_id, subjectid)

out_tsv_subj <- file.path(OUT_DIR, paste0("underlier_tallies_loeuf", DEF$suffix, "_constrained_tolerant_subject_gene.tsv"))
write_tsv(export_subject_gene, out_tsv_subj)
cat("Wrote:", normalizePath(out_tsv_subj), " (", nrow(export_subject_gene), " subject-gene underlier rows)\n")


# Tolerant (blue): straight fan into a name column right of the top-left
# points. Constrained (orange): L-elbow down into the white pocket under
# the orange stairs (vertical drop, then longer horizontal run) so leader
# lines clear the staircase instead of cutting across it.
callouts_xy <- callouts |>
    mutate(
        rank_in_tier = tier_df$rank_in_tier[match(gene_id, tier_df$gene_id)],
        frac_tier_rank = tier_df$frac_tier_rank[match(gene_id, tier_df$gene_id)],
        lbl = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id_no_version)
    )

tol_xy <- callouts_xy |>
    filter(as.character(LOEUF_tier) == "Tolerant") |>
    arrange(desc(pct_underlier)) |>
    mutate(
        x_txt = 0.155,
        y_txt = 10^seq(log10(0.0220), log10(0.0066), length.out = n()),
        x_elbow = frac_tier_rank,
        y_elbow = pct_underlier,
        elbow = FALSE
    )
con_xy <- callouts_xy |>
    filter(as.character(LOEUF_tier) == "Constrained") |>
    arrange(desc(pct_underlier), frac_tier_rank) |>
    mutate(
        # Drop below the orange stairs, then a short horizontal run into the
        # open pocket (0.30 overshot into mid stairs; keep left of the steps).
        x_txt = 0.18,
        y_txt = 10^seq(log10(0.00042), log10(0.00019), length.out = n()),
        x_elbow = frac_tier_rank,
        y_elbow = y_txt,
        elbow = TRUE
    )
callouts_xy <- bind_rows(tol_xy, con_xy)
con_elbow <- callouts_xy |> filter(elbow)

y_floor <- 1 / n_subjects
log_breaks <- c(y_floor, 0.001, 0.002, 0.005, 0.01, 0.02)
log_labels <- percent_format(accuracy = 0.01)(log_breaks)

n_under_genes_focus <- gene_df |> filter(LOEUF_tier %in% TIER_FOCUS, n > 0) |> nrow()

p_gene <- ggplot(tier_df, aes(x = frac_tier_rank, y = pct_underlier, color = LOEUF_tier)) +
    geom_line(linewidth = 0.55, alpha = 0.95) +
    # Straight leaders for Tolerant; L-elbow (vertical then horizontal) for
    # Constrained so orange lines stay clear of the orange staircase.
    geom_segment(
        data = callouts_xy |> filter(!elbow),
        aes(x = frac_tier_rank, y = pct_underlier, xend = x_txt, yend = y_txt, color = LOEUF_tier),
        linewidth = 0.25, alpha = 0.7, show.legend = FALSE
    ) +
    geom_segment(
        data = con_elbow,
        aes(x = frac_tier_rank, y = pct_underlier, xend = x_elbow, yend = y_elbow, color = LOEUF_tier),
        linewidth = 0.25, alpha = 0.7, show.legend = FALSE
    ) +
    geom_segment(
        data = con_elbow,
        aes(x = x_elbow, y = y_elbow, xend = x_txt, yend = y_txt, color = LOEUF_tier),
        linewidth = 0.25, alpha = 0.7, show.legend = FALSE
    ) +
    geom_point(
        data = callouts_xy,
        aes(x = frac_tier_rank, y = pct_underlier, color = LOEUF_tier),
        size = 1.6, show.legend = FALSE
    ) +
    geom_label(
        data = callouts_xy,
        aes(x = x_txt, y = y_txt, label = lbl, color = LOEUF_tier),
        inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 2.4, fontface = "italic",
        family = PAPER_FONT, linewidth = 0, fill = "white",
        label.padding = unit(0.08, "lines"), show.legend = FALSE
    ) +
    scale_color_manual(
        values = tier_colors_focus, labels = tier_label_vec, drop = FALSE,
        name = "LOEUF tier"
    ) +
    scale_x_continuous(labels = label_percent(accuracy = 0.1), expand = expansion(mult = c(0.02, 0.06))) +
    scale_y_log10(breaks = log_breaks, labels = log_labels, limits = c(y_floor * 0.9, NA)) +
    labs(
        x = "Within-tier rank (fraction)",
        y = "% of subjects with underliers (log scale)",
        title = "Constrained vs tolerant"
    ) +
    coord_cartesian(clip = "off") +
    theme_paper() +
    theme_paper_inner_title() +
    theme(
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(6, 8, 4, 4)
    )

p_box <- ggplot(
    tier_df,
    aes(y = LOEUF_tier, x = pct_underlier, fill = LOEUF_tier, color = LOEUF_tier)
) +
    geom_boxplot(width = 0.48, outlier.size = 0.45, outlier.alpha = 0.5,
                 alpha = 0.40, linewidth = 0.4) +
    scale_fill_manual(values = tier_colors_focus, drop = FALSE) +
    scale_color_manual(values = tier_colors_focus, drop = FALSE) +
    scale_x_log10(
        breaks = c(0.002, 0.02),
        labels = c("0.2%", "2%"),
        limits = c(y_floor * 0.9, NA),
        # Extra right pad so the 2% tick / whiskers are not clipped in the inset.
        expand = expansion(mult = c(0.06, 0.22))
    ) +
    scale_y_discrete(
        limits = rev(TIER_FOCUS),
        drop = FALSE,
        # Tiny labels so orange/blue boxes dominate the inset grob.
        labels = c(Constrained = "Constrained", Tolerant = "Tolerant")
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 6, base_family = PAPER_FONT) +
    theme(
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 4.0, color = "black",
                                   margin = margin(t = 1)),
        axis.text.y = element_text(size = 4.8, face = "bold", color = "black",
                                   margin = margin(r = 0.8)),
        axis.ticks = element_line(linewidth = 0.25, colour = "black"),
        axis.ticks.length = unit(1.5, "pt"),
        axis.line = element_line(linewidth = 0.30, colour = "black"),
        plot.margin = margin(3, 6, 3, 1),
        plot.background = element_rect(fill = alpha("white", 0.97),
                                       color = "grey65", linewidth = 0.35),
        panel.background = element_rect(fill = "white", color = NA)
    )

# Keep 2c a single ggplot so patchwork can line up its y-axis with 2a.
# annotation_custom scales the grob: absolute axis widths would crush the
# panel, so convert axis-l / panel widths to null units (small labels, big boxes).
g_box <- ggplotGrob(p_box)
axis_col <- unique(g_box$layout$l[g_box$layout$name == "axis-l"])
panel_col <- unique(g_box$layout$l[g_box$layout$name == "panel"])
if (length(axis_col) == 1L && length(panel_col) == 1L) {
    g_box$widths[[axis_col]] <- unit(0.40, "null")
    g_box$widths[[panel_col]] <- unit(1.00, "null")
}
# Inset kept fully inside the panel (was xmax=0.92 → right edge clipped in the 2x2).
p_gene_combo <- p_gene +
    annotation_custom(
        grob = g_box,
        xmin = 0.34, xmax = 0.78,
        ymin = 0.0026, ymax = 0.022
    )

.paper_plots <- list(c = p_gene_combo)
out_stem <- file.path(OUT_DIR, paste0("underlier_tallies_loeuf", DEF$suffix, "_constrained_tolerant"))
paper_ggsave(p_gene_combo, out_stem, width = 7.2, height = 4.4)
if (!is_paper_assemble()) cat("\nWrote:", out_stem, "\n")
