#!/usr/bin/env Rscript
# LOEUF-stratified gene-rank figure (complements underlier_tallies.png):
#   Per-tier gene curves + inset boxplot -- genes ranked by per-gene RNA-underlier
#   rate within each LOEUF tier.
# Subject-ranked panel is omitted here because it duplicates Panel B of
# underlier_tallies.png (same tally RDS; LOEUF is per-gene only).
#
# Underlier definition (per outlier_processing.R lines 318-322), per
# (gene, subject):
#       expression_zscore  < -4
#   AND expression_logcpm  < gene_mean(expression_logcpm) - 1
#
# Inputs (in working dir):
#   z_score_underlier_clean.rds    pre-filtered (gene_id, subjectid, ...) AFTER
#                                  the apply_extreme_subject_filter.R QC
#                                  (drops subjects with |z|>5 count >= 60).
#   kept_subjects.rds              vector of kept subjectids (for n_subjects)
#   expr_cpm_adj.rds               only used for nrow() = total #genes
#   loeuf_lookup.rds               from build_loeuf_lookup.R
# Outputs:
#   underlier_tallies_loeuf.png / .pdf
#   underlier_loeuf_prevalence_by_tier.tsv
#   underlier_loeuf_callouts.tsv

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

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG2X

source(file.path(ROOT, "scripts/helpers/_underlier_def.R"))
DEF <- resolve_underlier_def()
cat("Underlier definition: ", DEF$short, " ", DEF$tag, "\n", sep = "")

logcpm_off <- 1L

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_COLORS <- c(
    "Constrained" = "#B2182B",
    "Mid-low"     = "#EF8A62",
    "Mid-high"    = "#67A9CF",
    "Tolerant"    = "#2166AC",
    "No LOEUF"    = "grey60"
)
TIER_LABELS <- c(
    "Constrained" = "Constrained\nLOEUF < 0.35",
    "Mid-low"     = "Mid-low\n0.35 \u2264 LOEUF < 0.6",
    "Mid-high"    = "Mid-high\n0.6 \u2264 LOEUF < 1.0",
    "Tolerant"    = "Tolerant\nLOEUF \u2265 1.0",
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
cat("  kept subjects =", n_subjects, "; total eligible genes =", n_genes,
    "  (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s)\n")

cat("Loading ", DEF$rds_underlier, " ...\n", sep = "")
z_score_underlier <- read_rds(DEF$rds_underlier)
cat("  ", nrow(z_score_underlier), "underlier (gene, subject) rows after QC filter (",
    DEF$short, ")\n", sep = "")

cat("Loading loeuf_lookup.rds ...\n")
loeuf_lookup <- read_rds(LOEUF_LOOKUP)
cat("  ", nrow(loeuf_lookup), "loeuf rows\n")

# Build per-gene tally over the FULL eligible gene set (so genes with zero
# underliers contribute to per-tier denominators too).
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

# Join LOEUF -- prefer match on full gene_id (gencode-v30 versioned), then
# fall back to gene_id_no_version, then label "No LOEUF".
loeuf_full <- loeuf_lookup |>
    distinct(gene_id, .keep_all = TRUE) |>
    select(gene_id, LOEUF_full = LOEUF, tier_full = LOEUF_tier,
           symbol_full = symbol)
loeuf_nv <- loeuf_lookup |>
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
            ifelse(is.na(LOEUF),
                   "No LOEUF",
                   as.character(coalesce(tier_full, tier_nv))),
            levels = TIER_LEVELS
        )
    ) |>
    select(gene_id, gene_id_no_version, symbol, LOEUF, LOEUF_tier,
           n, pct_underlier)

cat("LOEUF coverage of eligible-gene universe:\n")
print(gene_df |> count(LOEUF_tier))

# Per-tier prevalence table (denominator is full eligible (gene, subject)
# pairs in that tier, not just underlier rows).
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

cat("\nPer-tier prevalence:\n")
print(prev_tbl)

prev_path <- paste0("underlier_loeuf_prevalence_by_tier", DEF$suffix, ".tsv")
write_tsv(prev_tbl, prev_path)
cat("Wrote:", prev_path, "\n")

# Call-outs: top 20 by pct_underlier from Constrained and Tolerant.
callouts <- gene_df |>
    filter(LOEUF_tier %in% c("Constrained", "Tolerant"), n > 0) |>
    group_by(LOEUF_tier) |>
    arrange(desc(pct_underlier), .by_group = TRUE) |>
    slice_head(n = 20) |>
    ungroup() |>
    arrange(LOEUF_tier, desc(pct_underlier))

write_tsv(callouts, DEF$callouts_tsv)
cat("Wrote:", DEF$callouts_tsv, " (", nrow(callouts), "rows )\n")

cat("\nTop 5 Constrained call-outs:\n")
print(callouts |> filter(LOEUF_tier == "Constrained") |> head(5))
cat("\nTop 5 Tolerant call-outs:\n")
print(callouts |> filter(LOEUF_tier == "Tolerant") |> head(5))

# ---- Panel A': per-tier ranked curves --------------------------------------
# Restrict the rank curve to genes that actually have at least one underlier
# (otherwise the y=0 floor dominates and obscures the shape). Per-tier
# denominators in the prevalence table use ALL genes in tier.
# X-axis: rank / n_genes_in_tier so curves are comparable across tiers
# (Constrained vs Tolerant differ greatly in N_genes).
tier_df <- gene_df |>
    filter(n > 0) |>
    group_by(LOEUF_tier) |>
    arrange(desc(pct_underlier), .by_group = TRUE) |>
    mutate(rank_in_tier = row_number()) |>
    ungroup() |>
    left_join(
        prev_tbl |> select(LOEUF_tier, n_genes_in_tier),
        by = "LOEUF_tier"
    ) |>
    mutate(
        frac_tier_rank = rank_in_tier / n_genes_in_tier
    )

# Pretty subtitle / caption strings (wrap so long lines are not clipped at save).
filter_subtitle <- str_wrap(
    sprintf(
        "Underlier per (gene, subject):  %s    %s    \u00b7    LOEUF tiers from gnomAD v4 (lof.oe_ci.upper)\nSubject QC: dropped 71/5,560 with |z|>5 count \u2265 60; analysis on %s kept subjects.",
        DEF$long, DEF$tag, comma(n_subjects)
    ),
    width = 92
)
prev_inline <- prev_tbl |>
    mutate(txt = sprintf("%s: %.1f/10k (%s genes)",
                         LOEUF_tier,
                         prevalence_per_10k,
                         comma(n_genes_in_tier))) |>
    pull(txt) |>
    paste(collapse = "  \u00b7  ")
caption_txt <- sprintf(
    "%s underlier (gene, subject) pairs across %s eligible genes and %s subjects.\nPer-tier prevalence (per 10,000 (gene, subject) pairs):  %s",
    comma(nrow(z_score_underlier)),
    comma(n_genes),
    comma(n_subjects),
    prev_inline
)

# Order tiers so darker constrained colors plot on top.
tier_df$LOEUF_tier <- factor(tier_df$LOEUF_tier, levels = TIER_LEVELS)

# Helper to recolor the legend with the multi-line tier labels.
tier_label_vec <- TIER_LABELS[TIER_LEVELS]

# Pre-compute call-out point coordinates so the geom_repel layer doesn't have
# to do a row-position lookup inside aes() (which broke layer composition).
callouts_xy <- callouts |>
    mutate(
        rank_in_tier = tier_df$rank_in_tier[match(gene_id, tier_df$gene_id)],
        frac_tier_rank = tier_df$frac_tier_rank[match(gene_id, tier_df$gene_id)],
        lbl          = ifelse(!is.na(symbol) & symbol != "",
                              symbol, gene_id_no_version)
    )

# Log y-axis (with a floor at 1/n_subjects since pct_underlier can equal that).
y_floor <- 1 / n_subjects
log_breaks <- c(y_floor, 0.001, 0.002, 0.005, 0.01, 0.02)
log_labels <- percent_format(accuracy = 0.01)(log_breaks)

p_gene <- ggplot(tier_df, aes(x = frac_tier_rank, y = pct_underlier,
                              color = LOEUF_tier)) +
    geom_line(linewidth = 0.55, alpha = 0.95) +
    geom_point(data = callouts_xy,
               aes(x = frac_tier_rank, y = pct_underlier, color = LOEUF_tier),
               size = 1.6, show.legend = FALSE) +
    ggrepel::geom_text_repel(
        data = callouts_xy,
        aes(x = frac_tier_rank, y = pct_underlier, label = lbl,
            color = LOEUF_tier),
        size = 2.7,
        min.segment.length = 0,
        segment.size = 0.2,
        box.padding = 0.4,
        point.padding = 0.25,
        max.overlaps = Inf,
        force = 3,
        seed = 17,
        show.legend = FALSE
    ) +
    scale_color_manual(values = TIER_COLORS, labels = tier_label_vec,
                       drop = FALSE, name = "LOEUF tier") +
    scale_x_continuous(
        labels = label_percent(accuracy = 0.1),
        expand = expansion(mult = c(0.02, 0.06))
    ) +
    scale_y_log10(breaks = log_breaks, labels = log_labels,
                  limits = c(y_floor * 0.9, NA)) +
    labs(
        title    = sprintf("Per-tier gene rank curves  (%s genes \u22651 underlier)",
                           comma(sum(gene_df$n > 0))),
        subtitle = "Within-tier rank / N genes in tier vs % subjects; labelled = top Constrained / Tolerant",
        x        = "Fraction of tier genes (rank among \u22651-underlier genes / N_genes in tier)",
        y        = "% of subjects underlying the gene  (log scale)"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(plot.title.position = "plot",
          legend.position     = "right",
          legend.key.height   = unit(1.0, "lines"),
          legend.text         = element_text(size = 8),
          legend.title        = element_text(size = 9),
          plot.margin         = margin(5.5, 10, 5.5, 8, "pt"))

# ---- inset: boxplot of pct_underlier per tier (underlier-only subset) -------
# Use the same log x-axis as the rank curve so per-tier shapes are comparable.
p_box <- ggplot(tier_df,
                aes(y = LOEUF_tier, x = pct_underlier,
                    fill = LOEUF_tier, color = LOEUF_tier)) +
    geom_boxplot(width = 0.55, outlier.size = 0.4, outlier.alpha = 0.5,
                 alpha = 0.30, linewidth = 0.35) +
    scale_fill_manual(values  = TIER_COLORS, drop = FALSE) +
    scale_color_manual(values = TIER_COLORS, drop = FALSE) +
    scale_x_log10(breaks = log_breaks,
                  labels = percent_format(accuracy = 0.01),
                  limits = c(y_floor * 0.9, NA)) +
    scale_y_discrete(limits = rev(TIER_LEVELS), drop = FALSE) +
    labs(x = "% subjects per gene  (log)", y = NULL,
         title = "Per-tier distribution (genes with \u22651 underlier)") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none",
          plot.title      = element_text(size = 8, face = "bold"),
          axis.text.x     = element_text(size = 7),
          axis.text.y     = element_text(size = 7),
          axis.title.x    = element_text(size = 7),
          plot.margin     = margin(4, 5, 5, 14, "pt"),
          plot.background = element_rect(fill = alpha("white", 0.95),
                                         color = "grey70", linewidth = 0.3))

# Inset shifted right and slightly smaller so tier labels (left) and main
# geom_repel callouts are less likely to overlap or clip.
p_gene_combo <- p_gene +
    inset_element(p_box,
                  left = 0.52, bottom = 0.56,
                  right = 0.99, top = 0.98,
                  align_to = "panel")

# ---- compose + save ---------------------------------------------------------
combined <- p_gene_combo +
    plot_annotation(
        title    = sprintf(
            "TOPMed expression underlier distributions, stratified by LOEUF  %s",
            DEF$tag),
        subtitle = filter_subtitle,
        caption  = caption_txt,
        theme    = theme(
            plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 10, color = "grey25",
                                          lineheight = 1.2),
            plot.caption  = element_text(size = 9,  color = "grey35"),
            plot.margin   = margin(10, 18, 10, 14)
        )
    )

out_png <- paste0("underlier_tallies_loeuf", DEF$suffix, ".png")
out_pdf <- paste0("underlier_tallies_loeuf", DEF$suffix, ".pdf")
ggsave(file.path(OUT_DIR, out_png), combined, width = 11, height = 6.0, dpi = 200)
ggsave(file.path(OUT_DIR, out_pdf), combined, width = 11, height = 6.0, device = cairo_pdf)
cat("\nWrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")
