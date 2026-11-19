#!/usr/bin/env Rscript
# Figure 2d: TOPMed vs AoU RNA-outlier prevalence, colored by LOEUF tier.
#
# Canonical R panel for assemble_nature_figures.R. Uses theme_paper() so
# axis styling / margins match panels b and c. The older Python renderer
# plot_topmed_vs_aou_prevalence_loeuf.py is kept only as a historical
# reference; do not use it for the paper composite.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG2C
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_N <- 2L
TIER_ORDER <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")

sunflower <- function(k, radius) {
  if (k <= 1L) return(data.frame(dx = 0, dy = 0))
  i <- seq_len(k) - 1
  r <- radius * sqrt((i + 0.5) / k)
  theta <- i * pi * (3 - sqrt(5))
  data.frame(dx = r * cos(theta), dy = r * sin(theta))
}

spread_lattice <- function(d) {
  d <- d |>
    mutate(x_jit = x, y_jit = y, .row = dplyr::row_number())
  pieces <- vector("list", 0)
  for (key in split(d, list(d$n_topmed, d$n_new), drop = TRUE)) {
    k <- nrow(key)
    if (k == 1L) {
      pieces[[length(pieces) + 1L]] <- key
      next
    }
    spacing <- min(
      log10((key$n_topmed[[1]] + 1) / key$n_topmed[[1]]),
      log10((key$n_new[[1]] + 1) / key$n_new[[1]])
    )
    radius <- min(0.45 * spacing, 0.014 * sqrt(k))
    off <- sunflower(k, radius)
    key$x_jit <- key$x + off$dx
    key$y_jit <- key$y + off$dy
    pieces[[length(pieces) + 1L]] <- key
  }
  bind_rows(pieces) |> arrange(.row) |> select(-.row)
}

d <- read_tsv(file.path(FIG1C, "topmed_vs_new_z_-3_prevalence_per_gene.tsv"),
              show_col_types = FALSE) |>
  filter(n_topmed >= MIN_N, n_new >= MIN_N) |>
  mutate(
    LOEUF_tier = ifelse(is.na(LOEUF_tier) | !LOEUF_tier %in% TIER_ORDER,
                        "No LOEUF", as.character(LOEUF_tier)),
    LOEUF_tier = factor(LOEUF_tier, levels = TIER_ORDER),
    x = -log10(prev_topmed),
    y = -log10(prev_new)
  )
d <- spread_lattice(d)

# Shuffle so no tier is systematically painted last.
set.seed(1)
d <- d[sample.int(nrow(d)), , drop = FALSE]

lo <- min(c(d$x_jit, d$y_jit))
hi <- max(c(d$x_jit, d$y_jit))
pad <- 0.06 * (hi - lo)
lim <- c(lo - pad, hi + pad)

sp <- cor(d$prev_topmed, d$prev_new, method = "spearman")
pe <- cor(d$x, d$y, method = "pearson")

# No coord_fixed: fill the cell so x/y axes align with neighbouring ggplots.
p <- ggplot(d, aes(x_jit, y_jit, color = LOEUF_tier)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey55", linewidth = 0.45) +
  geom_point(size = 0.65, alpha = 0.75, stroke = 0.12) +
  scale_color_manual(values = TIER_COLORS, drop = FALSE, name = "LOEUF tier") +
  coord_cartesian(xlim = lim, ylim = lim, expand = FALSE) +
  labs(
    x = expression(-log[10]~TOPMed~RNA-outlier~prevalence),
    y = expression(-log[10]~AoU~RNA-outlier~prevalence)
  ) +
  theme_paper() +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(9, "pt"),
    legend.spacing.y = unit(1, "pt"),
    # Match neighbouring panels so cowplot can align axis lines.
    plot.margin = margin(4, 8, 4, 4)
  )

.paper_plots <- list(d = p)
paper_maybe_ggsave(p, file.path(OUT_DIR, "fig2_topmed_vs_aou_prevalence_loeuf_ggplot"),
                   width = NATURE_DOUBLE_IN / 2, height = NATURE_DOUBLE_IN / 2)
if (!is_paper_assemble()) {
  cat(sprintf(
    "Wrote ggplot Fig2d; n_genes=%s Spearman=%.2f Pearson(-log10)=%.2f\n",
    format(nrow(d), big.mark = ","), sp, pe
  ))
}
