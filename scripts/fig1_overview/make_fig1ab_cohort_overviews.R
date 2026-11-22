#!/usr/bin/env Rscript
# Figure 1a top: TOPMed + AoU cohort overview.
#
# Preferred asset (author-provided):
#   scripts/fig1_overview/assets/fig1a_cohort_overview_provided.png
# Generated fallback:
#   scripts/fig1_overview/assets/publication_figure.png
# Legacy screenshot:
#   scripts/fig1_overview/assets/fig1ab_cohort_overview_screenshot.png
#
# Definition strip nests under this image in assemble_nature_figures.R (still tag "a").
suppressPackageStartupMessages({
  library(png)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
.PKG <- Sys.getenv("TOPMED_SUBMISSION_PACKAGE", unset = "")
if (!nzchar(.PKG)) .PKG <- ROOT
.SCR_ROOT <- if (dir.exists(file.path(.PKG, "scripts"))) .PKG else ROOT
source(file.path(.SCR_ROOT, "scripts/helpers/_paths.R"))
source(file.path(.SCR_ROOT, "scripts/helpers/paper_theme.R"))

FIG1_SCR <- file.path(.SCR_ROOT, "scripts/fig1_overview")
ASSET_DIR <- file.path(FIG1_SCR, "assets")
PROVIDED <- file.path(ASSET_DIR, "fig1a_cohort_overview_provided.png")
PUB_PNG <- file.path(ASSET_DIR, "publication_figure.png")
LEGACY <- file.path(ASSET_DIR, "fig1ab_cohort_overview_screenshot.png")

# Prefer the author-provided overview. Only regenerate via Python when that
# file is absent and FORCE_OVERVIEW_RECREATE=1 (or no provided/publication PNG).
PY <- Sys.getenv("PAPER_PYTHON", unset = "python3")
recreate_py <- file.path(FIG1_SCR, "recreate_publication_figure.py")
force_recreate <- identical(Sys.getenv("FORCE_OVERVIEW_RECREATE", unset = "0"), "1")
use_provided <- file.exists(PROVIDED) && !force_recreate

if (!use_provided) {
  skip <- identical(Sys.getenv("SKIP_PANELS", unset = "0"), "1")
  if ((!skip || !file.exists(PUB_PNG)) && file.exists(recreate_py)) {
    cat("==== python recreate_publication_figure.py ====\n")
    st <- system2(PY, shQuote(recreate_py), stdout = "", stderr = "")
    if (!identical(st, 0L)) stop("recreate_publication_figure.py failed (status ", st, ")")
  }
} else {
  cat("==== using provided Fig1a overview asset ====\n")
  file.copy(PROVIDED, PUB_PNG, overwrite = TRUE)
}

ASSET <- if (file.exists(PROVIDED) && !force_recreate) {
  PROVIDED
} else if (file.exists(PUB_PNG)) {
  PUB_PNG
} else {
  LEGACY
}
if (!file.exists(ASSET)) {
  ASSET <- file.path(ROOT, "scripts/fig1_overview/assets/fig1a_cohort_overview_provided.png")
}
if (!file.exists(ASSET)) {
  ASSET <- file.path(ROOT, "scripts/fig1_overview/assets/publication_figure.png")
}
if (!file.exists(ASSET)) {
  ASSET <- file.path(ROOT, "scripts/fig1_overview/assets/fig1ab_cohort_overview_screenshot.png")
}
stopifnot(file.exists(ASSET))
OUT_DIR <- FIG1A
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

img <- png::readPNG(ASSET)
# Crop outer near-white margins on light canvases.
trim_white <- function(arr, thresh = 0.98, pad = 4L) {
  if (dim(arr)[3] < 3) return(arr)
  ink <- arr[,,1] < thresh | arr[,,2] < thresh | arr[,,3] < thresh
  # If most of the canvas is dark (provided black-bg art), do not trim-white.
  if (mean(ink) > 0.85) return(arr)
  rows <- which(rowSums(ink) > 0)
  cols <- which(colSums(ink) > 0)
  if (!length(rows) || !length(cols)) return(arr)
  r0 <- max(1L, min(rows) - pad); r1 <- min(nrow(arr), max(rows) + pad)
  c0 <- max(1L, min(cols) - pad); c1 <- min(ncol(arr), max(cols) + pad)
  arr[r0:r1, c0:c1, , drop = FALSE]
}
# Crop empty transparent / near-black borders on dark canvases.
# The author-provided overview is RGBA with alpha=0 "black" padding that
# becomes white whitespace in ggsave — trim by opaque-pixel density, then
# flatten onto white so the panel has no hidden transparent margins.
trim_dark_margins <- function(arr, dark_thresh = 0.06, min_frac = 0.02, pad = 3L) {
  if (length(dim(arr)) < 3 || dim(arr)[3] < 3) return(arr)
  if (dim(arr)[3] >= 4L) {
    ink <- arr[,,4] > 0.04
  } else {
    ink <- arr[,,1] > dark_thresh | arr[,,2] > dark_thresh | arr[,,3] > dark_thresh
  }
  # Only for mostly-empty canvases (provided overview is ~30% opaque).
  if (mean(ink) > 0.55) return(arr)
  rows <- which(rowMeans(ink) > min_frac)
  cols <- which(colMeans(ink) > min_frac)
  if (!length(rows) || !length(cols)) return(arr)
  r0 <- max(1L, min(rows) - pad); r1 <- min(nrow(arr), max(rows) + pad)
  c0 <- max(1L, min(cols) - pad); c1 <- min(ncol(arr), max(cols) + pad)
  arr[r0:r1, c0:c1, , drop = FALSE]
}
flatten_rgba_on_white <- function(arr) {
  if (length(dim(arr)) < 3 || dim(arr)[3] < 4L) {
    if (dim(arr)[3] >= 3L) return(arr[,,1:3, drop = FALSE])
    return(arr)
  }
  a <- arr[,,4]
  out <- array(1, dim = c(dim(arr)[1], dim(arr)[2], 3L))
  for (i in 1:3) out[,,i] <- arr[,,i] * a + (1 - a)
  out
}
combined <- trim_white(flatten_rgba_on_white(trim_dark_margins(img)), thresh = 0.985, pad = 2L)
nr <- dim(combined)[1]; nc <- dim(combined)[2]
mid <- as.integer(round(nc * 0.495))
left <- trim_white(combined[, seq_len(mid), , drop = FALSE])
right <- trim_white(combined[, (mid + 1L):nc, , drop = FALSE])

combined_path <- file.path(OUT_DIR, "fig1a_cohort_overview_combined.png")
left_path <- file.path(OUT_DIR, "fig1a_topmed_cohort_sankey.png")
right_path <- file.path(OUT_DIR, "fig1b_aou_cohort_overview.png")

file.copy(ASSET, file.path(OUT_DIR, "publication_figure.png"), overwrite = TRUE)
png::writePNG(combined, combined_path)
png::writePNG(left, left_path)
png::writePNG(right, right_path)

.paper_plots <- list(
  a_path = combined_path,
  left_path = left_path,
  right_path = right_path
)

writeLines(c(
  "Paper 1a top = author-provided cohort overview (or recreate_publication_figure.py fallback).",
  paste0("Asset used: ", ASSET),
  "Script: scripts/fig1_overview/make_fig1ab_cohort_overviews.R",
  "Definition strip nested under this image in assemble (still panel a)."
), file.path(OUT_DIR, "SOURCE.txt"))

if (!is_paper_assemble()) {
  message("wrote ", combined_path)
  message("wrote ", left_path)
  message("wrote ", right_path)
}
