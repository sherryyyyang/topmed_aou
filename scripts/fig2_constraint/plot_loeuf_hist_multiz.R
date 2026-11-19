#!/usr/bin/env Rscript
# Standalone figure:
# Histograms of per-gene underlier counts (across subjects),
# faceted by LOEUF decile and strict z cutoffs (z<-2, z<-3, z<-4).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG2X


Z_LEVELS <- c("z<-2", "z<-3", "z<-4")
DECILE_LEVELS <- paste0("Decile ", 1:10)
Z_COLORS <- c(
  "z<-2" = "#1B9E77",
  "z<-3" = "#D95F02",
  "z<-4" = "#7570B3"
)

cat("Loading inputs ...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
loeuf_lookup <- read_rds(LOEUF_LOOKUP)
all_eligible_genes <- rownames(expr_cpm_adj)

loeuf_full <- loeuf_lookup |>
  distinct(gene_id, .keep_all = TRUE) |>
  select(gene_id, LOEUF_full = LOEUF)
loeuf_nv <- loeuf_lookup |>
  distinct(gene_id_no_version, .keep_all = TRUE) |>
  select(gene_id_no_version, LOEUF_nv = LOEUF)

gene_df <- tibble(gene_id = all_eligible_genes) |>
  mutate(gene_id_no_version = str_remove(gene_id, "\\.\\d+$")) |>
  left_join(loeuf_full, by = "gene_id") |>
  left_join(loeuf_nv, by = "gene_id_no_version") |>
  mutate(
    LOEUF = coalesce(LOEUF_full, LOEUF_nv)
  ) |>
  select(gene_id, LOEUF)

decile_lookup <- gene_df |>
  filter(!is.na(LOEUF)) |>
  mutate(decile_id = ntile(LOEUF, 10L)) |>
  mutate(
    LOEUF_decile = factor(
      paste0("D", decile_id),
      levels = paste0("D", 1:10),
      labels = DECILE_LEVELS
    )
  ) |>
  select(gene_id, LOEUF_decile)

collect_gene_decile <- function(z_label, rds_path) {
  under <- read_rds(rds_path) |>
    select(gene_id, subjectid) |>
    distinct() |>
    count(gene_id, name = "n_under")

  decile_lookup |>
    left_join(under, by = "gene_id") |>
    mutate(
      n_under = replace_na(n_under, 0L),
      z_label = factor(z_label, levels = Z_LEVELS)
    )
}

hist_df <- bind_rows(
  collect_gene_decile("z<-2", file.path(DATA_UNDERLIER, "z_score_underlier_z_-2_clean.rds")),
  collect_gene_decile("z<-3", file.path(DATA_UNDERLIER, "z_score_underlier_z_-3_clean.rds")),
  collect_gene_decile("z<-4", file.path(DATA_UNDERLIER, "z_score_underlier_z_-4_clean.rds"))
)

facet_stats <- hist_df |>
  group_by(z_label, LOEUF_decile) |>
  summarise(
    median_under = median(n_under),
    max_under = max(n_under),
    .groups = "drop"
  ) |>
  mutate(
    stat_label = paste0("med=", median_under, "  max=", max_under)
  )

binw <- 1L

p <- ggplot(hist_df, aes(x = n_under, fill = z_label, color = z_label)) +
  geom_histogram(
    aes(y = after_stat(count)),
    position = "identity",
    alpha = 0.30,
    linewidth = 0.2,
    binwidth = binw,
    boundary = -0.5 * binw
  ) +
  geom_vline(
    data = facet_stats,
    aes(xintercept = median_under),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "black",
    linewidth = 0.25
  ) +
  geom_vline(
    data = facet_stats,
    aes(xintercept = max_under),
    inherit.aes = FALSE,
    linetype = "dotted",
    color = "grey20",
    linewidth = 0.25
  ) +
  geom_text(
    data = facet_stats,
    aes(x = Inf, y = Inf, label = stat_label),
    inherit.aes = FALSE,
    hjust = 1.02,
    vjust = 1.05,
    size = 2.1,
    color = "grey20"
  ) +
  facet_grid(z_label ~ LOEUF_decile, scales = "free_x") +
  scale_fill_manual(values = Z_COLORS, drop = FALSE) +
  scale_color_manual(values = Z_COLORS, drop = FALSE) +
  scale_x_continuous(
    breaks = function(lims) {
      upper <- max(0, floor(lims[2]))
      if (upper <= 2) {
        c(0, upper)
      } else if (upper <= 8) {
        c(0, floor(upper / 2), upper)
      } else {
        c(0, round(upper / 2), upper)
      }
    },
    labels = label_number(accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 4),
    labels = comma_format(accuracy = 1)
  ) +
  labs(
    title = "Per-gene underlier count histograms by LOEUF decile",
    subtitle = "Each panel uses genes in one LOEUF decile (equal-count bins among genes with LOEUF). Columns are deciles; rows are z-thresholds.",
    x = "Number of underliers",
    y = "Gene count",
    fill = "Strict cutoff",
    color = "Strict cutoff"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )

out_png <- "loeuf_hist_multiz.png"
out_pdf <- "loeuf_hist_multiz.pdf"
ggsave(file.path(OUT_DIR, out_png), p, width = 24, height = 7.5, dpi = 220)
ggsave(file.path(OUT_DIR, out_pdf), p, width = 24, height = 7.5, device = cairo_pdf)
cat("Wrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")

