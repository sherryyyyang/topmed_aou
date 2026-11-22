#!/usr/bin/env Rscript
# Figure 2d: TOPMed vs AoU RNA-outlier prevalence, colored by LOEUF tier.
#
# Canonical R panel for assemble_nature_figures.R. Uses theme_paper() so
# axis styling / margins match panels b and c. The older Python renderer
# plot_topmed_vs_aou_prevalence_loeuf.py is kept only as a historical
# reference; do not use it for the paper composite.
#
# Scatter data/geometry match Fig1d exactly (shared helper): same gene set,
# same -log10 prevalences, same axis limits, same Spearman. Only LOEUF
# coloring (plus sunflower spread so stacked genes remain visible) differs.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
source(file.path(ROOT, "scripts/helpers/load_topmed_aou_prevalence_scatter.R"))
OUT_DIR <- FIG2C
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

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

sc <- .load_topmed_aou_prevalence_scatter()
d <- sc$genes |>
  mutate(
    LOEUF_tier = ifelse(is.na(LOEUF_tier) | !LOEUF_tier %in% TIER_ORDER,
                        "No LOEUF", as.character(LOEUF_tier)),
    LOEUF_tier = factor(LOEUF_tier, levels = TIER_ORDER)
  )
d <- spread_lattice(d)

# Shuffle so no tier is systematically painted last.
set.seed(1)
d <- d[sample.int(nrow(d)), , drop = FALSE]

# Axis limits match Fig1d (true -log10 prevalences), not the jitter envelope.
lim <- sc$lim
sp <- sc$sp
pe <- sc$pe

# No coord_fixed: fill the cell so x/y axes align with neighbouring ggplots.
p <- ggplot(d, aes(x_jit, y_jit, color = LOEUF_tier)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey55", linewidth = 0.45) +
  geom_point(size = 0.65, alpha = 0.75, stroke = 0.12) +
  scale_color_manual(values = TIER_COLORS, drop = FALSE, name = "LOEUF tier") +
  coord_cartesian(xlim = lim, ylim = lim, expand = FALSE) +
  # Bottom-left: LOEUF legend already occupies the upper-left (unlike Fig1d).
  annotate(
    "text",
    x = lim[1] + 0.03 * (lim[2] - lim[1]),
    y = lim[1] + 0.03 * (lim[2] - lim[1]),
    label = sprintf("Spearman ρ = %.2f", sp),
    hjust = 0, vjust = 0,
    size = 3.0, family = PAPER_FONT, colour = "grey15"
  ) +
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
    "Wrote ggplot Fig2d; PC=%s n_genes=%s Spearman=%.2f Pearson(-log10)=%.2f\n  data=%s\n",
    PC_SETTING, format(sc$n_genes, big.mark = ","), sp, pe, sc$path
  ))
}
