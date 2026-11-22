#!/usr/bin/env Rscript
# Figure 1 companion panel: how haploinsufficient expression outliers (RNA
# underliers) are selected. Two steps — extreme z-score tails, then retain the
# low-expression (log2CPM) arm.
#
# Canonical: Rscript scripts/fig1_overview/plot_underlier_definition_demo.R
# Assembled next to the cohort overview in assemble_nature_figures.R.
#
# Reference sketch (provenance only; this script redraws the panel):
#   scripts/fig1_overview/assets/fig1_underlier_definition_demo_screenshot.png
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG1A
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COL_TAIL <- "#8B1A1A"
COL_MID  <- "#F5B896"
COL_LINE <- "#4A5568"
COL_STEP <- "#C45C26"

# ---- Step 1: standard normal with |z| >= 3 tails ----
z <- seq(-4.2, 4.2, length.out = 801)
dens <- data.frame(z = z, y = dnorm(z))
mid <- dens |> filter(z > -3, z < 3)
left <- dens |> filter(z < -3)
right <- dens |> filter(z >= 3)

p1 <- ggplot(dens, aes(z, y)) +
  geom_area(data = mid, fill = COL_MID, color = NA, alpha = 0.95) +
  geom_area(data = left, fill = COL_TAIL, color = NA, alpha = 0.95) +
  geom_area(data = right, fill = COL_TAIL, color = NA, alpha = 0.95) +
  geom_line(color = COL_LINE, linewidth = 0.55) +
  geom_vline(xintercept = 0, linetype = "dashed", color = COL_LINE, linewidth = 0.4) +
  annotate("text", x = 0, y = dnorm(0) * 0.55, label = "-3 < z < 3",
           size = 3.1, color = COL_LINE, family = PAPER_FONT) +
  annotate("text", x = -3.55, y = dnorm(-3) * 2.8, label = "z < -3",
           size = 3.0, color = COL_TAIL, fontface = "bold", family = PAPER_FONT) +
  annotate("text", x = 3.55, y = dnorm(3) * 2.8, label = "z >= 3",
           size = 3.0, color = COL_TAIL, fontface = "bold", family = PAPER_FONT) +
  scale_x_continuous(breaks = seq(-4, 4, 2), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Step 1 — extreme z across the gene",
    x = NULL, y = NULL
  ) +
  theme_paper() +
  theme(
    plot.title = element_text(
      family = PAPER_FONT, face = "bold", size = 8.5, colour = COL_LINE,
      hjust = 0, lineheight = 1.05, margin = margin(b = 4)
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(4, 6, 4, 4)
  )

# ---- Step 2: retain low-expression tail (flowchart) ----
# Horizontal flow so panel b reads as one long strip beside cohort overview.
# Inner boxes inset from the grey frame so borders do not overlap.
boxes <- data.frame(
  xmin = c(0.09, 0.385, 0.68),
  xmax = c(0.335, 0.615, 0.91),
  ymin = c(0.34, 0.34, 0.34),
  ymax = c(0.66, 0.66, 0.66),
  fill = c("white", COL_MID, COL_TAIL),
  color = c(COL_TAIL, COL_STEP, COL_TAIL),
  lab = c(
    "All candidates\nwith z < -3",
    "log2CPM <\ngene_mean - 1",
    "RNA\nUnderlier"
  ),
  lab_col = c(COL_TAIL, COL_LINE, "white"),
  stringsAsFactors = FALSE
)
arrows <- data.frame(
  x = c(0.335, 0.615), xend = c(0.385, 0.68),
  y = c(0.50, 0.50), yend = c(0.50, 0.50)
)

p2 <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.02, ymax = 0.98,
           fill = "grey96", color = "grey75", linewidth = 0.4) +
  geom_rect(
    data = boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = boxes$fill, color = boxes$color, linewidth = 0.55
  ) +
  geom_segment(
    data = arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.10, "inches"), type = "closed"),
    color = COL_STEP, linewidth = 0.65
  ) +
  geom_text(
    data = boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = lab),
    color = boxes$lab_col, size = 2.5, fontface = "bold",
    family = PAPER_FONT, lineheight = 0.92
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(title = "Step 2 — retain low-expression tail", x = NULL, y = NULL) +
  theme_void(base_family = PAPER_FONT) +
  theme(
    plot.title = element_text(
      family = PAPER_FONT, face = "bold", size = 8.5, colour = COL_TAIL,
      hjust = 0, margin = margin(b = 4, t = 2)
    ),
    plot.margin = margin(4, 4, 4, 6)
  )

p <- p1 + p2 + plot_layout(widths = c(1.05, 1.35), heights = 1)
.paper_plots <- list(demo = p)

paper_maybe_ggsave(
  p, file.path(OUT_DIR, "fig1_underlier_definition_demo"),
  width = NATURE_DOUBLE_IN, height = 1.70
)
if (!is_paper_assemble()) {
  # Also keep a copy of the reference screenshot beside the redraw.
  ref <- file.path(ROOT, "scripts/fig1_overview/assets/fig1_underlier_definition_demo_screenshot.png")
  if (file.exists(ref)) {
    file.copy(ref, file.path(OUT_DIR, "fig1_underlier_definition_demo_screenshot.png"),
              overwrite = TRUE)
  }
  writeLines(c(
    "Paper Fig1 underlier-definition demo (next to cohort overview).",
    "Script: scripts/fig1_overview/plot_underlier_definition_demo.R",
    "Reference sketch: scripts/fig1_overview/assets/fig1_underlier_definition_demo_screenshot.png"
  ), file.path(OUT_DIR, "SOURCE_underlier_definition_demo.txt"))
  cat("Wrote fig1_underlier_definition_demo\n")
}
