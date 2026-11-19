#!/usr/bin/env Rscript
# Figure 1a top: TOPMed + AoU cohort overview from recreate_publication_figure.py.
#
# Preferred asset (generated):
#   scripts/fig1_overview/assets/publication_figure.png
# Fallback (legacy screenshot):
#   scripts/fig1_overview/assets/fig1ab_cohort_overview_screenshot.png
#
# Definition strip nests under this image in assemble_nature_figures.R (still tag "a").
suppressPackageStartupMessages({
  library(png)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
.PKG <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
.SCR_ROOT <- if (dir.exists(file.path(.PKG, "scripts"))) .PKG else ROOT
source(file.path(.SCR_ROOT, "scripts/helpers/_paths.R"))
source(file.path(.SCR_ROOT, "scripts/helpers/paper_theme.R"))

FIG1_SCR <- file.path(.SCR_ROOT, "scripts/fig1_overview")
ASSET_DIR <- file.path(FIG1_SCR, "assets")
PUB_PNG <- file.path(ASSET_DIR, "publication_figure.png")
LEGACY <- file.path(ASSET_DIR, "fig1ab_cohort_overview_screenshot.png")

# Rebuild publication overview unless SKIP_PANELS=1 and PNG already exists.
PY <- Sys.getenv(
  "PAPER_PYTHON",
  unset = "/net/wonderland/home/sheyang/.pixi/envs/python/bin/python3"
)
recreate_py <- file.path(FIG1_SCR, "recreate_publication_figure.py")
skip <- identical(Sys.getenv("SKIP_PANELS", unset = "0"), "1")
if ((!skip || !file.exists(PUB_PNG)) && file.exists(recreate_py)) {
  cat("==== python recreate_publication_figure.py ====\n")
  st <- system2(PY, shQuote(recreate_py), stdout = "", stderr = "")
  if (!identical(st, 0L)) stop("recreate_publication_figure.py failed (status ", st, ")")
}

ASSET <- if (file.exists(PUB_PNG)) PUB_PNG else LEGACY
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
# Trim outer white margins from the matplotlib export.
trim_white <- function(arr, thresh = 0.98, pad = 4L) {
  ink <- arr[,,1] < thresh | arr[,,2] < thresh | arr[,,3] < thresh
  rows <- which(rowSums(ink) > 0)
  cols <- which(colSums(ink) > 0)
  if (!length(rows) || !length(cols)) return(arr)
  r0 <- max(1L, min(rows) - pad); r1 <- min(nrow(arr), max(rows) + pad)
  c0 <- max(1L, min(cols) - pad); c1 <- min(ncol(arr), max(cols) + pad)
  arr[r0:r1, c0:c1, , drop = FALSE]
}
combined <- trim_white(img)
nr <- dim(combined)[1]; nc <- dim(combined)[2]
# Split near the vertical divider (~mid) for optional single-cohort reuse.
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
  "Paper 1a top = recreate_publication_figure.py (TOPMed + AoU overview).",
  paste0("Asset used: ", ASSET),
  "Script: scripts/fig1_overview/make_fig1ab_cohort_overviews.R",
  "Definition strip nested under this image in assemble (still panel a)."
), file.path(OUT_DIR, "SOURCE.txt"))

if (!is_paper_assemble()) {
  message("wrote ", combined_path)
  message("wrote ", left_path)
  message("wrote ", right_path)
}
