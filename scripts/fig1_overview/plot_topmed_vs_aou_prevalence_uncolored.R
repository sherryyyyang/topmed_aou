#!/usr/bin/env Rscript
# Figure 1 prevalence scatter (TOPMed vs AoU, haplo + z < -3).
#
# Canonical R panel for assemble_nature_figures.R (panel e under gene ranks).
# Uses theme_paper() so fonts / axis styling match plot_underlier_tallies_topmed_aou.R.
# The older Python script plot_topmed_vs_aou_prevalence_uncolored.py is kept only
# as a historical renderer; do not use it for the paper composite.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG1C
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_N <- 2L
CUT <- 3
BASE_S <- 10
# Match Python size knots (points^2-ish → ggplot size in mm-ish via scale)
size_area <- function(n) {
  n <- as.numeric(n)
  kn <- c(1, 2, 5, 10, 25, 49)
  ka <- c(BASE_S * (1^0.45), BASE_S * (2^0.45), BASE_S * (5^0.45),
          BASE_S * (10^0.45), 62, 125)
  approx(kn, ka, xout = pmin(pmax(n, 1), 49), rule = 2)$y
}

d <- read_tsv(file.path(OUT_DIR, "topmed_vs_new_z_-3_prevalence_per_gene.tsv"),
              show_col_types = FALSE) |>
  filter(n_topmed >= MIN_N, n_new >= MIN_N) |>
  mutate(x = -log10(prev_topmed), y = -log10(prev_new))

pos <- d |>
  group_by(n_topmed, n_new, x, y) |>
  summarise(n_genes = n(), .groups = "drop") |>
  mutate(
    area = size_area(n_genes),
    in_corner = (x >= CUT & y >= CUT)
  )

lo <- min(c(pos$x, pos$y))
hi <- max(c(pos$x, pos$y))
pad <- 0.06 * (hi - lo)
lim <- c(lo - pad, hi + pad)

# Same theme_paper margins / fonts as gene-rank panel (Fig 1 left column).
p <- ggplot(pos, aes(x, y)) +
  annotate("rect", xmin = CUT, xmax = lim[2], ymin = CUT, ymax = lim[2],
           fill = COL_PTV, alpha = 0.07, color = COL_PTV, linewidth = 0.45) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey55", linewidth = 0.45) +
  geom_vline(xintercept = CUT, linetype = "dotted", color = COL_PTV, linewidth = 0.4) +
  geom_hline(yintercept = CUT, linetype = "dotted", color = COL_PTV, linewidth = 0.4) +
  geom_point(
    data = filter(pos, !in_corner),
    aes(size = n_genes),
    color = "#C8C8C8", alpha = 0.42, stroke = 0.15
  ) +
  geom_point(
    data = filter(pos, in_corner),
    aes(size = n_genes),
    shape = 21, fill = "#C8C8C8", color = COL_PTV, alpha = 0.55, stroke = 0.45
  ) +
  scale_size_area(
    max_size = 8, breaks = {
      mx <- max(pos$n_genes)
      br <- c(1, 2, 5, 10, 25)
      unique(c(br[br <= mx], mx))
    },
    labels = function(x) ifelse(x == 1, "1 gene", paste0(x, " genes")),
    name = "genes stacked\non one point"
  ) +
  # Fill the column width like the gene-rank panel (no coord_fixed shrink),
  # keeping equal data limits so x and y share the same scale.
  coord_cartesian(xlim = lim, ylim = lim, expand = FALSE) +
  labs(
    x = expression(-log[10]~TOPMed~RNA-outlier~prevalence),
    y = expression(-log[10]~AoU~RNA-outlier~prevalence)
  ) +
  theme_paper() +
  theme(
    legend.position = "right",
    legend.justification = "center",
    # Match gene-rank panel left/right margins so y-axes stack-align.
    plot.margin = margin(4, 6, 4, 8)
  )

.paper_plots <- list(e = p)
paper_maybe_ggsave(p, file.path(OUT_DIR, "fig1_topmed_vs_aou_prevalence_uncolored_ggplot"),
                   width = NATURE_DOUBLE_IN / 2, height = 3.6)
if (!is_paper_assemble()) {
  cat("Wrote ggplot prevalence panel; n_genes=", nrow(d),
      " Spearman=", round(cor(d$prev_topmed, d$prev_new, method = "spearman"), 3), "\n", sep = "")
}