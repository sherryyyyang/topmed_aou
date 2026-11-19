# Shared figure style (Nature final-submission + Word-readable type).
# https://www.nature.com/nature/for-authors/final-submission
#   sans-serif Helvetica or Arial, same family on every panel
#   panel tags: bold upright lowercase a, b, c (Nature 8 pt at print size)
#   other text 5–7 pt at print size; this Word draft uses a larger uniform size
#   double column 183 mm; strokes 0.25–1 pt
#
# Headlines (plot.title / subtitle / caption) are blank: that prose belongs
# in the Word figure legend, not on the graphic. Restore plot.title only for
# short inner identifiers (e.g. a 2 x 3 grid of densities).

THEME_PAPER_LOADED <- TRUE

if (isTRUE(tryCatch(capabilities("cairo"), error = function(e) FALSE))) {
  options(bitmapType = "cairo")
}

PAPER_FONT <- "Arial"
NATURE_DOUBLE_IN <- 183 / 25.4
NATURE_SINGLE_IN <- 89 / 25.4
NATURE_PAGE_IN <- 247 / 25.4
PAPER_DPI <- 300

# Okabe–Ito (colour-blind safe). Sequential colour: viridis.
# https://jfly.uni-koeln.de/color/
OI <- c(
  black      = "#000000",
  orange     = "#E69F00",
  sky        = "#56B4E9",
  green      = "#009E73",
  yellow     = "#F0E442",
  blue       = "#0072B2",
  vermillion = "#D55E00",
  purple     = "#CC79A7"
)
COL_TOPMED <- unname(OI["blue"])
COL_AOU    <- unname(OI["purple"])
COL_PTV    <- unname(OI["vermillion"])
COL_NON    <- "grey55"
COL_ACMG   <- unname(OI["purple"])
TIER_COLORS <- c(
  Constrained = unname(OI["vermillion"]),
  `Mid-low`   = unname(OI["orange"]),
  `Mid-high`  = unname(OI["sky"]),
  Tolerant    = unname(OI["blue"]),
  `No LOEUF`  = unname(OI["black"])
)
CUTOFF_COLORS <- c(
  `z < -2`   = unname(OI["purple"]),
  `z < -3`   = unname(OI["vermillion"]),
  `z < -3.5` = unname(OI["orange"]),
  `z < -4`   = unname(OI["green"]),
  `z < -5`   = unname(OI["blue"])
)

theme_paper <- function(base_size = 11, base_family = PAPER_FONT) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      plot.title          = ggplot2::element_blank(),
      plot.subtitle       = ggplot2::element_blank(),
      plot.caption        = ggplot2::element_blank(),
      axis.title          = ggplot2::element_text(size = 10, color = "black"),
      axis.text           = ggplot2::element_text(size = 8, color = "black"),
      axis.line           = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      axis.ticks          = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      legend.title        = ggplot2::element_text(size = 8, face = "plain"),
      legend.text         = ggplot2::element_text(size = 8),
      legend.background   = ggplot2::element_blank(),
      legend.key          = ggplot2::element_blank(),
      legend.key.size     = ggplot2::unit(10, "pt"),
      strip.background    = ggplot2::element_blank(),
      strip.text          = ggplot2::element_text(size = 9, face = "bold", color = "black"),
      panel.grid.major    = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank(),
      plot.margin         = ggplot2::margin(4, 6, 4, 8)
    )
}

# Short inner identifiers (Fig. 3b grid titles, enrichment class names).
theme_paper_inner_title <- function() {
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      family = PAPER_FONT, face = "bold", size = 9, colour = "black",
      hjust = 0, margin = ggplot2::margin(b = 3)
    )
  )
}

with_phrase <- function(p, title) {
  p + ggplot2::labs(title = title) + theme_paper_inner_title()
}

theme_paper_tag <- function() {
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = PAPER_FONT, face = "bold", size = 14, colour = "black"
    ),
    plot.tag.position = "topleft"
  )
}

is_paper_assemble <- function() isTRUE(getOption("paper.assemble", FALSE))

paper_ggsave <- function(plot, path_stem, width, height, ...) {
  ggplot2::ggsave(paste0(path_stem, ".png"), plot, width = width, height = height,
                  dpi = PAPER_DPI, bg = "white", ...)
  ok_cairo <- isTRUE(tryCatch(capabilities("cairo"), error = function(e) FALSE))
  if (ok_cairo) {
    ggplot2::ggsave(paste0(path_stem, ".pdf"), plot, width = width, height = height,
                    device = grDevices::cairo_pdf, bg = "white", ...)
  }
  invisible(plot)
}

paper_maybe_ggsave <- function(plot, path_stem, width, height, ...) {
  if (is_paper_assemble()) return(invisible(plot))
  paper_ggsave(plot, path_stem, width, height, ...)
}

png_aspect <- function(path) {
  stopifnot(file.exists(path))
  img <- png::readPNG(path)
  dim(img)[1] / dim(img)[2]
}

# Fill the patchwork cell. Size that cell with png_aspect() so the raster
# is neither stretched nor letterboxed in a too-tall row.
wrap_png <- function(path) {
  stopifnot(file.exists(path))
  grob <- grid::rasterGrob(
    png::readPNG(path), interpolate = TRUE,
    width = grid::unit(1, "npc"), height = grid::unit(1, "npc")
  )
  ggplot2::ggplot() +
    ggplot2::annotation_custom(grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 4, 2, 8))
}

# Place a PNG in a vertical band of the cell (top/bottom as npc fractions)
# so it can line up with ggplot axis grobs in a neighbouring column.
wrap_png_vpad <- function(path, top = 0, bottom = 0) {
  stopifnot(file.exists(path), top >= 0, bottom >= 0, top + bottom < 1)
  grob <- grid::rasterGrob(
    png::readPNG(path), interpolate = TRUE,
    width = grid::unit(1, "npc"), height = grid::unit(1, "npc")
  )
  ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob, xmin = 0, xmax = 1, ymin = bottom, ymax = 1 - top
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 4, 0, 8))
}

# Map the PNG's x-axis line (fraction from the top of the PNG) onto y = 0
# of the patchwork cell so it can share a neighbour panel's x-axis line.
# Pixels below that axis (tick labels / axis title) overflow with clip = "off".
wrap_png_pin_xaxis <- function(path, top = 0, xaxis_from_top = 0.919) {
  stopifnot(file.exists(path), top >= 0, xaxis_from_top > 0, xaxis_from_top < 1)
  grob <- grid::rasterGrob(
    png::readPNG(path), interpolate = TRUE,
    width = grid::unit(1, "npc"), height = grid::unit(1, "npc")
  )
  ymax <- 1 - top
  ymin <- ymax * (1 - 1 / xaxis_from_top)
  ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob, xmin = 0, xmax = 1, ymin = ymin, ymax = ymax
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 4, 0, 8))
}

# Keep native aspect; pin to the top of a taller cell (no magick).
wrap_png_top <- function(path) {
  stopifnot(file.exists(path))
  img <- png::readPNG(path)
  ar <- dim(img)[1] / dim(img)[2]
  grob <- grid::rasterGrob(
    img,
    interpolate = TRUE,
    x = grid::unit(0.5, "npc"),
    y = grid::unit(1, "npc"),
    just = c("centre", "top"),
    width = grid::unit(1, "npc"),
    height = grid::unit(ar, "snpc")
  )
  patchwork::wrap_elements(full = grob) +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 4))
}

# Treat a nested ggplot/patchwork (e.g. Fig. 3b 2 x 3) as one tagged panel.
# Optional extra margin (top, right, bottom, left, in pt) so the panel tag
# sits in empty space instead of on an inner title.
as_single_panel <- function(p, margin_pt = c(2, 4, 2, 8)) {
  patchwork::wrap_elements(full = p, clip = FALSE) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        margin_pt[1], margin_pt[2], margin_pt[3], margin_pt[4]
      )
    )
}
