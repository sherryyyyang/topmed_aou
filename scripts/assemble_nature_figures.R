#!/usr/bin/env Rscript
# Orchestrate every paper figure so type, tags, and layout match.
#
# Nature final-submission (https://www.nature.com/nature/for-authors/final-submission):
#   Arial throughout; lowercase bold panel tags a, b, c; double column 183 mm.
# Headlines live in the Word captions, not on the graphics.
#
# Layout: patchwork for tagged ggplot composition; cowplot for wrapping
# Python PNG panels (study design, scatters).
#
# Usage (from anywhere):
#   Rscript scripts/assemble_nature_figures.R
# Optional:
#   SKIP_PANELS=1  reuse existing Python PNGs / skip re-sourcing R panel scripts
#                  (only re-composites)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(grid)
  library(png)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("ANALYSIS_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
# Prefer scripts shipped inside TOPMED_SUBMISSION_PACKAGE when that directory has them.
.PKG <- Sys.getenv("TOPMED_SUBMISSION_PACKAGE", unset = "")
if (!nzchar(.PKG)) .PKG <- ROOT
.SCR_ROOT <- if (dir.exists(file.path(.PKG, "scripts"))) .PKG else ROOT
source(file.path(.SCR_ROOT, "scripts/helpers/_paths.R"))
source(file.path(.SCR_ROOT, "scripts/helpers/paper_theme.R"))
FINAL_DRAFT <- .PKG
# Always stage composites under the package, never under PC55 integrity roots.
PAPER_COMPOSITES <- file.path(FINAL_DRAFT, "paper", "composites")
SCR <- file.path(.SCR_ROOT, "scripts")

options(paper.assemble = TRUE)
Sys.setenv(DEFINITION_LABEL = "z_-3")

dir.create(PAPER_COMPOSITES, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_DRAFT, recursive = TRUE, showWarnings = FALSE)
FIG_MAIN <- file.path(FINAL_DRAFT, "figures", "main")
FIG_SUPP <- file.path(FINAL_DRAFT, "figures", "supplementary")
src_pan <- file.path(FINAL_DRAFT, "source_panels")
as_in <- file.path(src_pan, "as_in_docx")
SUP_FIGS <- FIG_SUPP
dir.create(FIG_MAIN, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_SUPP, recursive = TRUE, showWarnings = FALSE)
dir.create(src_pan, recursive = TRUE, showWarnings = FALSE)
dir.create(as_in, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_FIGS, recursive = TRUE, showWarnings = FALSE)

SKIP <- Sys.getenv("SKIP_PANELS", unset = "0") == "1"
ONLY <- Sys.getenv("ONLY_FIG", unset = "")

source_panel <- function(script) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  e <- new.env(parent = globalenv())
  cat("\n==== sourcing ", basename(script), " ====\n", sep = "")
  source(script, local = e)
  # Nested source(_paths.R) defaults to .GlobalEnv and can clobber these;
  # re-assert package delivery targets after every panel script.
  FINAL_DRAFT <<- .PKG
  PAPER_COMPOSITES <<- file.path(.PKG, "paper", "composites")
  if (is.null(e$.paper_plots)) stop("no .paper_plots from ", script)
  e$.paper_plots
}

PY <- Sys.getenv("PAPER_PYTHON", unset = "python3")

run_py <- function(script) {
  cat("\n==== python ", basename(script), " ====\n", sep = "")
  st <- system2(PY, shQuote(script), stdout = "", stderr = "")
  if (!identical(st, 0L)) stop("python failed: ", script, " (status ", st, ")")
}

tag_it <- function(p) {
  p + plot_annotation(tag_levels = "a") & theme_paper_tag()
}

save_fig <- function(p, name, width = NATURE_DOUBLE_IN, height) {
  stem <- file.path(PAPER_COMPOSITES, name)
  paper_ggsave(p, stem, width = width, height = height)
  png <- paste0(stem, ".png")
  pdf <- paste0(stem, ".pdf")
  dest_dir <- if (startsWith(name, "Supplementary_")) FIG_SUPP else FIG_MAIN
  file.copy(png, file.path(dest_dir, paste0(name, ".png")), overwrite = TRUE)
  file.copy(png, file.path(src_pan, paste0(name, ".png")), overwrite = TRUE)
  file.copy(png, file.path(as_in, paste0(name, ".png")), overwrite = TRUE)
  if (file.exists(pdf)) {
    file.copy(pdf, file.path(dest_dir, paste0(name, ".pdf")), overwrite = TRUE)
  }
  cat("wrote ", name, "\n", sep = "")
}

copy_panel <- function(src, dest_name, dest_dir = src_pan) {
  if (file.exists(src)) {
    file.copy(src, file.path(dest_dir, dest_name), overwrite = TRUE)
  }
}

copy_ggplot <- function(p, dest_name, width, height, dest_dir = src_pan) {
  paper_ggsave(p, file.path(dest_dir, dest_name), width = width, height = height)
}

if (!ONLY %in% c("3", "4")) {
if (!identical(ONLY, "2")) {
# ---------- Figure 1 ----------
# a = cohort overview + underlier-definition strip (one panel)
# b/c = gene / subject ranks
# d = TOPMed vs AoU prevalence (square, right of b|c)
if (!SKIP) {
  f1ab <- source_panel(file.path(SCR, "fig1_overview/make_fig1ab_cohort_overviews.R"))
  f1demo <- source_panel(file.path(SCR, "fig1_overview/plot_underlier_definition_demo.R"))
  f1 <- source_panel(file.path(SCR, "fig1_overview/plot_underlier_tallies_topmed_aou.R"))
  f1e <- source_panel(file.path(SCR, "fig1_overview/plot_topmed_vs_aou_prevalence_uncolored.R"))
} else {
  f1ab <- source_panel(file.path(SCR, "fig1_overview/make_fig1ab_cohort_overviews.R"))
  f1demo <- source_panel(file.path(SCR, "fig1_overview/plot_underlier_definition_demo.R"))
  f1 <- source_panel(file.path(SCR, "fig1_overview/plot_underlier_tallies_topmed_aou.R"))
  f1e <- source_panel(file.path(SCR, "fig1_overview/plot_topmed_vs_aou_prevalence_uncolored.R"))
}

png_1a <- f1ab$a_path
# Prevalence on the right of gene|subject; size legend on top of d.
f1e$e <- f1e$e + theme(
  legend.position = "top",
  legend.direction = "horizontal",
  legend.justification = "center",
  legend.box = "horizontal",
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 2, 0)
)
# Keep provided overview at native aspect (no horizontal stretch).
# Dark empty borders are cropped in make_fig1ab; size the cell to the
# trimmed PNG so Step 1/2 sit tight under the art (no letterbox gap).
h1_overview <- NATURE_DOUBLE_IN * png_aspect(png_1a)
h1_strip <- 1.70
h1_a <- h1_overview + h1_strip
h1_bottom <- 4.80
w_left <- 1.00
w_right <- 1.05
f1$b <- f1$b + theme(
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.box.spacing = unit(0, "pt"),
  plot.margin = margin(2, 4, 0, 4)
)
f1$c <- f1$c + theme(plot.margin = margin(0, 4, 2, 4))
d_plot <- f1e$e + theme(
  legend.position = "none",
  plot.margin = margin(2, 4, 0, 2)
)
d_leg <- cowplot::get_legend(
  f1e$e + theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.key.height = unit(0.28, "lines"),
    legend.key.width = unit(0.9, "lines"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7.5),
    legend.margin = margin(0, 0, 0, 0)
  )
)
# Render overview + strip at the same width, then stack with zero join gap
# (avoids patchwork gutters and mismatched native-PNG vs ggsave widths).
# annotation_raster fills the device exactly (wrap_png/annotation_custom
# can leave a white gutter at the bottom of the overview cell).
ov_rgb <- {
  arr <- png::readPNG(png_1a)
  if (dim(arr)[3] >= 4L) {
    a <- arr[,,4]; out <- arr[,,1:3, drop = FALSE]
    for (i in 1:3) out[,,i] <- out[,,i] * a + (1 - a)
    out
  } else {
    arr[,,1:3, drop = FALSE]
  }
}
ov_plot <- ggplot2::ggplot() +
  ggplot2::annotation_raster(ov_rgb, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
ov_tmp <- tempfile(fileext = ".png")
ggsave(
  ov_tmp, ov_plot,
  width = NATURE_DOUBLE_IN, height = h1_overview, dpi = PAPER_DPI, bg = "white"
)
strip_tmp <- tempfile(fileext = ".png")
ggsave(
  strip_tmp,
  as_single_panel(f1demo$demo, margin_pt = c(0, 2, 2, 2)),
  width = NATURE_DOUBLE_IN, height = h1_strip, dpi = PAPER_DPI, bg = "white"
)
.to_rgb <- function(arr) {
  if (dim(arr)[3] >= 4L) {
    a <- arr[,,4]
    out <- arr[,,1:3, drop = FALSE]
    for (i in 1:3) out[,,i] <- out[,,i] * a + (1 - a)
    out
  } else {
    arr[,,1:3, drop = FALSE]
  }
}
ov <- .to_rgb(png::readPNG(ov_tmp))
st <- .to_rgb(png::readPNG(strip_tmp))
# Drop residual white padding at the overview/strip join.
.trim_white_band <- function(arr, from = c("top", "bottom"), thresh = 0.985, pad = 2L) {
  ink <- arr[,,1] < thresh | arr[,,2] < thresh | arr[,,3] < thresh
  rows <- which(rowSums(ink) > 0)
  if (!length(rows)) return(arr)
  r0 <- 1L; r1 <- nrow(arr)
  if ("top" %in% from) r0 <- max(1L, min(rows) - pad)
  if ("bottom" %in% from) r1 <- min(nrow(arr), max(rows) + pad)
  arr[r0:r1, , , drop = FALSE]
}
ov <- .trim_white_band(ov, from = "bottom", pad = 1L)
st <- .trim_white_band(st, from = c("top", "bottom"), pad = 1L)
wc <- min(dim(ov)[2], dim(st)[2])
ov <- ov[, seq_len(wc), , drop = FALSE]
st <- st[, seq_len(wc), , drop = FALSE]
a_img <- array(1, dim = c(dim(ov)[1] + dim(st)[1], wc, 3L))
a_img[seq_len(dim(ov)[1]), , ] <- ov
a_img[dim(ov)[1] + seq_len(dim(st)[1]), , ] <- st
# Panel-a height from stacked raster (overview + strip, no extra gutter).
h1_a <- NATURE_DOUBLE_IN * (dim(a_img)[1] / dim(a_img)[2])
panel_a <- ggplot() +
  annotation_custom(grid::rasterGrob(a_img, width = unit(1, "npc"), height = unit(1, "npc"))) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
left_bc <- cowplot::plot_grid(
  f1$b, f1$c,
  ncol = 1,
  rel_heights = c(1.1, 1),
  align = "v",
  axis = "lr",
  labels = c("b", "c"),
  label_size = 14,
  label_fontface = "bold",
  label_fontfamily = PAPER_FONT,
  label_x = 0.004,
  label_y = 0.995,
  hjust = 0,
  vjust = 1
)
right_d <- cowplot::plot_grid(
  d_leg,
  cowplot::plot_grid(
    d_plot,
    labels = "d",
    label_size = 14,
    label_fontface = "bold",
    label_fontfamily = PAPER_FONT,
    label_x = 0.004,
    label_y = 0.995,
    hjust = 0,
    vjust = 1
  ),
  ncol = 1,
  rel_heights = c(0.14, 1)
)
bottom_tagged <- cowplot::plot_grid(
  left_bc, right_d,
  nrow = 1,
  rel_widths = c(w_left, w_right)
)
fig1 <- cowplot::plot_grid(
  cowplot::plot_grid(
    panel_a,
    labels = "a",
    label_size = 14,
    label_fontface = "bold",
    label_fontfamily = PAPER_FONT,
    label_x = 0.004,
    label_y = 0.995,
    hjust = 0,
    vjust = 1
  ),
  bottom_tagged,
  ncol = 1,
  rel_heights = c(h1_a, h1_bottom)
)
save_fig(fig1, "Figure1", height = h1_a + h1_bottom)
copy_panel(png_1a, "Fig1a_cohort_overview_combined.png")
copy_panel(f1ab$left_path, "Fig1a_topmed_cohort_sankey.png")
copy_panel(f1ab$right_path, "Fig1a_aou_cohort_overview.png")
copy_ggplot(f1demo$demo, "Fig1a_underlier_definition_demo", NATURE_DOUBLE_IN, h1_strip)
copy_ggplot(f1$b, "Fig1b_gene_ranks", NATURE_DOUBLE_IN * w_left / (w_left + w_right), 1.55)
copy_ggplot(f1$c, "Fig1c_subject_ranks", NATURE_DOUBLE_IN * w_left / (w_left + w_right), 1.40)
copy_ggplot(f1e$e, "Fig1d_topmed_vs_aou_uncolored",
            NATURE_DOUBLE_IN * w_right / (w_left + w_right), 2.85)
copy_panel(file.path(FIG1B, "fig1_underlier_tallies_topmed_aou.png"),
           "Fig1cd_underlier_burden.png")
if (identical(Sys.getenv("ONLY_FIG"), "1")) {
  cat("ONLY_FIG=1: stopped after Figure 1\n")
  quit(save = "no")
}
}

# ---------- Figure 2: a forest, b decile, c constrained/tolerant, d LOEUF scatter ----------
# Equal 2x2: a|b / c|d (same cell size).

f2ab <- source_panel(file.path(SCR, "fig2_constraint/plot_loeuf_calibration.R"))
f2c  <- source_panel(file.path(SCR, "fig2_constraint/plot_underlier_tallies_loeuf_constrained_tolerant.R"))
f2d  <- source_panel(file.path(SCR, "fig2_constraint/plot_topmed_vs_aou_prevalence_loeuf.R"))
w2 <- NATURE_DOUBLE_IN
h2 <- 7.2
# Equal 2x2; all ggplots with theme_paper so axis lines align a↔c / b↔d / c↔d.
axis_margin <- margin(4, 8, 4, 4)
f2ab$a <- f2ab$a + theme(plot.margin = axis_margin)
f2ab$b <- f2ab$b + theme(plot.margin = axis_margin)
f2c$c  <- f2c$c  + theme(plot.margin = axis_margin)
f2d$d  <- f2d$d  + theme(plot.margin = axis_margin)
aligned <- cowplot::align_plots(
  f2ab$a, f2ab$b, f2c$c, f2d$d,
  align = "hv", axis = "tblr"
)
top <- cowplot::plot_grid(
  aligned[[1]], aligned[[2]],
  nrow = 1,
  rel_widths = c(1, 1),
  labels = c("a", "b"),
  label_size = 14,
  label_fontface = "bold",
  label_fontfamily = PAPER_FONT,
  label_x = 0.004,
  label_y = 0.995,
  hjust = 0,
  vjust = 1
)
bot <- cowplot::plot_grid(
  aligned[[3]], aligned[[4]],
  nrow = 1,
  rel_widths = c(1, 1),
  labels = c("c", "d"),
  label_size = 14,
  label_fontface = "bold",
  label_fontfamily = PAPER_FONT,
  label_x = 0.004,
  label_y = 0.995,
  hjust = 0,
  vjust = 1
)
fig2 <- cowplot::plot_grid(
  top, bot,
  ncol = 1,
  rel_heights = c(1, 1)
)
save_fig(fig2, "Figure2", width = w2, height = h2)
copy_ggplot(f2ab$a, "Fig2a_loeuf_forest", NATURE_DOUBLE_IN * 0.55, 3.2)
copy_ggplot(f2ab$b, "Fig2b_loeuf_decile", NATURE_DOUBLE_IN * 0.45, 3.2)
copy_ggplot(as_single_panel(f2c$c), "Fig2c_constrained_vs_tolerant", NATURE_DOUBLE_IN, 3.6)
copy_ggplot(f2d$d, "Fig2d_topmed_vs_aou_loeuf", NATURE_DOUBLE_IN / 2, NATURE_DOUBLE_IN / 2)
copy_panel(file.path(FIG2B, paste0("underlier_tallies_loeuf_z_-3_constrained_tolerant.png")),
           "Fig2c_constrained_vs_tolerant_standalone.png")
if (identical(Sys.getenv("ONLY_FIG"), "2")) {
  cat("ONLY_FIG=2: stopped after Figure 2\n")
  quit(save = "no")
}
}

# ---------- Figure 3: a AF<0.1% z-score row, b gene-level, c ClinGen ----------
# Full-width stacked panels; no plot_spacer() / blank cells.

if (!identical(ONLY, "4")) {
f3a <- source_panel(file.path(SCR, "fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level_2x2.R"))
f3b <- source_panel(file.path(SCR, "fig3_ptv_vs_rna/plot_ptv_vs_nonptv_zscore.R"))
f3c <- source_panel(file.path(SCR, "fig3_ptv_vs_rna/plot_candidate_gene_zscore_dist_clingen_z3.R"))
# Rows a/b swapped vs prior draft: z-score row on top, gene-level second.
fig3 <- cowplot::plot_grid(
  f3b$b,
  f3a$a,
  f3c$c,
  ncol = 1,
  rel_heights = c(0.82, 0.90, 1.18),
  labels = c("a", "b", "c"),
  label_size = 14,
  label_fontface = "bold",
  label_fontfamily = PAPER_FONT,
  label_x = 0.004,
  label_y = 0.995,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "l"
)
save_fig(fig3, "Figure3", height = 8.8)
copy_ggplot(f3b$b, "Fig3a_ptv_vs_nonptv_zscore", NATURE_DOUBLE_IN, 3.4)
copy_ggplot(f3a$a, "Fig3b_ptv_vs_underlier", NATURE_DOUBLE_IN, 2.9)
copy_ggplot(f3c$c, "Fig3c_clingen_hist", NATURE_DOUBLE_IN, 3.4)
if (identical(Sys.getenv("ONLY_FIG"), "3")) {
  cat("ONLY_FIG=3: stopped after Figure 3\n")
  quit(save = "no")
}
}

# ---------- Figure 4: a RBM20 ranks, b 4-dot, c MYH11 ranks, d 4-dot ----------
# Bridge gutters (untagged) link the blue/pink RNA curves at each labelled
# gene to the matching example strip; no plot_spacer() / blank cells.

f4 <- source_panel(file.path(SCR, "fig4_acmg/plot_acmg_4dot_no_ci.R"))
# Letter order is alphabetical in patchwork designs (A..F), so ABC/DEF maps
# ranks | bridge | example on each row. Empty tags keep a–d for the real panels.
fig4 <- (
  f4$a + f4$bridge_rbm20 + f4$b +
    f4$c + f4$bridge_myh11 + f4$d +
    plot_layout(
      design = "ABC\nDEF",
      widths = c(3.45, 0.38, 1.45),
      heights = c(1, 1),
      guides = "collect"
    ) +
    plot_annotation(tag_levels = list(c("a", "", "b", "c", "", "d"))) &
    theme_paper_tag() &
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.text = element_text(size = 8, family = PAPER_FONT),
      legend.box.margin = margin(2, 0, 0, 0),
      plot.margin = margin(6, 8, 4, 4)
    )
)
save_fig(fig4, "Figure4", height = 9.0)
copy_ggplot(f4$a, "Fig4a_RBM20_ranks", NATURE_DOUBLE_IN * 0.71, 4.5)
copy_ggplot(f4$b, "Fig4b_RBM20_4dot", NATURE_DOUBLE_IN * 0.29, 4.5)
copy_ggplot(f4$c, "Fig4c_MYH11_ranks", NATURE_DOUBLE_IN * 0.71, 4.5)
copy_ggplot(f4$d, "Fig4d_MYH11_4dot", NATURE_DOUBLE_IN * 0.29, 4.5)
if (identical(Sys.getenv("ONLY_FIG"), "4")) {
  cat("ONLY_FIG=4: stopped after Figure 4\n")
  quit(save = "no")
}

stale <- c(
  "Fig1A_study_design.png", "Fig1B_rank_curves.png",
  "Fig1C_topmed_vs_aou_uncolored.png", "Fig2A_loeuf_calibration.png",
  "Fig2B_constrained_vs_tolerant.png", "Fig2C_loeuf_scatter_panelB_crop.png",
  "Fig2C_ORIGINAL_two_panel_refined_LOEUF.png",
  "Fig3A_ptv_vs_underlier_z-3_only.png", "Fig3B_ptv_vs_nonptv_zscore.png",
  "Fig3C_clingen_zscore_hist.png", "Fig4A_RBM20.png", "Fig4B_MYH11.png"
)
invisible(file.remove(file.path(src_pan, stale)[file.exists(file.path(src_pan, stale))]))

cat("\nAll figures written to\n  ", PAPER_COMPOSITES, "\n  ", FINAL_DRAFT, "\n  ",
    src_pan, "\n  ", as_in, "\n", sep = "")
