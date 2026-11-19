#!/usr/bin/env Rscript
# ClinGen HI z<-3 gene-level PTV-vs-underlier scatter (single panel from
# ptv_vs_underlier_gene_level_2x2_clingen_hi) with the 8 candidate genes whose
# per-subject z-score distributions are in candidate_gene_zscore_dist_clingen_z_-3
# highlighted and labelled.
#
# Input : ../ptv_vs_underlier_gene_level_2x2_clingen_hi_genes.tsv
# Output: z_-3_clingen_scatter_highlight_candidates.png / .pdf

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
  library(ggplot2); library(scales); library(ggrepel)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG3C

PW <- ROOT

MODE_Z3 <- "Strict z<-3 + log2CPM"

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_COLORS <- c("Constrained"="#B2182B","Mid-low"="#EF8A62","Mid-high"="#67A9CF",
                 "Tolerant"="#2166AC","No LOEUF"="grey70")

candidates <- tribble(
  ~symbol,   ~corner,
  "NOG",     "High underlier / low PTV",
  "CDKN1C",  "High underlier / low PTV",
  "FGFR1",   "High underlier / low PTV",
  "SLC6A8",  "High underlier / low PTV",
  "DDX3X",   "High PTV / low underlier",
  "RB1",     "High PTV / low underlier",
  "PUF60",   "High PTV / low underlier",
  "CHD7",    "High PTV / low underlier"
)
CORNER_COLORS <- c("High underlier / low PTV"="#1B7837", "High PTV / low underlier"="#762A83")

genes <- read_tsv(file.path(FIG3C, "ptv_vs_underlier_gene_level_2x2_clingen_hi_genes.tsv"),
                  show_col_types = FALSE) |>
  filter(mode == MODE_Z3) |>
  mutate(LOEUF_tier = factor(ifelse(is.na(LOEUF_tier), "No LOEUF", as.character(LOEUF_tier)),
                             levels = TIER_LEVELS))

hi <- genes |>
  inner_join(candidates, by = "symbol") |>
  mutate(corner = factor(corner, levels = c("High underlier / low PTV", "High PTV / low underlier")))
stopifnot(nrow(hi) == nrow(candidates))

# KCNH2: additionally marked (not one of the z-score-distribution candidates),
# shown with a distinct orange ring/label.
KCNH2_COLOR <- "#D95F02"
kcnh2 <- filter(genes, symbol == "KCNH2")
stopifnot(nrow(kcnh2) == 1)

# Some highlighted candidates share identical (pct_ptv, pct_under) coordinates
# (e.g. FGFR1 and CDKN1C: both 0 PTV carriers, 18 underliers) and would plot as a
# single ring. Spread co-located rings vertically by a small amount so both are
# visible; x (PTV frequency) is left exact.
Y_DODGE <- 0.014
hi <- hi |>
  group_by(pct_ptv, pct_under) |>
  mutate(pct_under_plot = pct_under + (row_number() - (n() + 1) / 2) * Y_DODGE) |>
  ungroup()

x_max <- max(genes$pct_ptv)

p <- ggplot(genes, aes(pct_ptv, pct_under)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_point(aes(color = LOEUF_tier), alpha = 0.6, size = 2.2, stroke = 0) +
  geom_point(data = hi, aes(y = pct_under_plot, fill = corner), shape = 21, size = 4.2,
             stroke = 1.2, color = "black") +
  ggrepel::geom_text_repel(
    data = hi, aes(y = pct_under_plot, label = symbol, color = corner), fontface = "bold", size = 3.4,
    min.segment.length = 0, segment.size = 0.3, segment.color = "grey30",
    box.padding = 0.9, point.padding = 0.5, max.overlaps = Inf, seed = 7, show.legend = FALSE) +
  geom_point(data = kcnh2, shape = 21, size = 4.4, stroke = 1.3,
             fill = KCNH2_COLOR, color = "black") +
  ggrepel::geom_text_repel(
    data = kcnh2, aes(label = symbol), color = KCNH2_COLOR, fontface = "bold", size = 3.4,
    min.segment.length = 0, segment.size = 0.3, segment.color = "grey30",
    box.padding = 0.9, point.padding = 0.5, max.overlaps = Inf, seed = 11, show.legend = FALSE) +
  scale_color_manual(values = c(TIER_COLORS, CORNER_COLORS), name = "LOEUF tier",
                     breaks = TIER_LEVELS, drop = FALSE) +
  scale_fill_manual(values = CORNER_COLORS, name = "Highlighted candidate") +
  scale_x_continuous(labels = label_percent(scale = 1, accuracy = 0.01), limits = c(0, x_max),
                     expand = expansion(mult = c(0.03, 0.08))) +
  scale_y_continuous(labels = label_percent(scale = 1, accuracy = 0.01),
                     expand = expansion(mult = c(0.03, 0.08))) +
  labs(
    title = "ClinGen HI genes, strict z<-3 + log2CPM: rare PTV (AF<0.5%) carrier vs RNA underlier frequency",
    subtitle = str_wrap(
      "z<-3 panel from ptv_vs_underlier_gene_level_2x2_clingen_hi, restricted to ClinGen haploinsufficiency 'Sufficient Evidence' genes. Rings = 8 candidate genes whose per-subject z-score distributions are in candidate_gene_zscore_dist_clingen_z_-3. Green = high-underlier / low-PTV; purple = high-PTV / low-underlier; orange = KCNH2 (additionally marked). Background points colored by LOEUF tier.",
      width = 118),
    caption = "Dashed line = equal frequency (y = x). Co-located candidate rings (e.g. FGFR1 & CDKN1C, both 0 PTV / 18 underliers) are offset vertically for visibility; true value is identical.",
    x = "% of subjects carrying a rare PTV (AF<0.5%) in gene",
    y = "% of subjects who are RNA underliers for gene") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, color = "grey25", lineheight = 1.15),
        plot.caption = element_text(size = 7.5, color = "grey35", hjust = 0),
        legend.position = "right", plot.margin = margin(5.5, 14, 5.5, 5.5))

ggsave(file.path(OUT_DIR, "z_-3_clingen_scatter_highlight_candidates.png"), p, width = 9.5, height = 7.2, dpi = 220)
ggsave(file.path(OUT_DIR, "z_-3_clingen_scatter_highlight_candidates.pdf"), p, width = 9.5, height = 7.2, device = cairo_pdf)
cat("Wrote:", normalizePath("z_-3_clingen_scatter_highlight_candidates.png"), "\n")
