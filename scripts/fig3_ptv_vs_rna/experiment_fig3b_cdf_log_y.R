#!/usr/bin/env Rscript
# EXPERIMENT (not wired into assemble_nature_figures.R): the Fig. 3b centre
# panel with a log10 cumulative axis, so the curves still reach 100% but each
# decade of the left tail gets equal vertical space.
#
# Reads the cached z-score vectors written by plot_ptv_vs_nonptv_zscore.R.
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales); library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

CACHE <- file.path(FIG3B, "ptv_vs_nonptv_zscore_cache_v2.rds")
OUT_DIR <- file.path(ROOT, "fig3_ptv_vs_rna", "experiments")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

r <- readRDS(CACHE)[["AF<0.1%"]]
COL <- c("PTV carrier" = COL_PTV, "non-carrier" = COL_NON)

# Explicit ECDF: stat_ecdf() would emit y = 0 at the left edge, which is
# -Inf once the axis is logged.
ecdf_tbl <- function(z, grp) {
  z <- sort(z)
  tibble(z = z, cum = seq_along(z) / length(z), grp = grp) |>
    group_by(z) |>
    slice_tail(n = 1) |>
    ungroup()
}
cdf <- bind_rows(
  ecdf_tbl(r$ptv_carrier_z, "PTV carrier"),
  ecdf_tbl(r$non_sample_z, "non-carrier")
) |>
  filter(cum > 0)

base <- function(p) {
  p +
    geom_step(linewidth = 0.8, direction = "hv") +
    geom_vline(xintercept = c(-2, -3, -4), linetype = "dotted",
               color = "grey60", linewidth = 0.3) +
    scale_color_manual(values = COL, name = NULL) +
    labs(x = "expression z-score", y = "cumulative \u2264 z") +
    theme_paper() + theme_paper_inner_title() +
    theme(legend.position = "bottom")
}

p_lin <- base(ggplot(cdf, aes(z, cum, color = grp))) +
  coord_cartesian(xlim = c(-6, 4)) +
  scale_y_continuous(labels = label_percent(), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "Linear cumulative axis")

brks <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)
p_log <- base(ggplot(cdf, aes(z, cum, color = grp))) +
  coord_cartesian(xlim = c(-6, 4)) +
  scale_y_log10(breaks = brks,
                labels = c("0.001%", "0.01%", "0.1%", "1%", "10%", "100%"),
                expand = expansion(mult = c(0.03, 0.05))) +
  annotation_logticks(sides = "l", size = 0.25,
                      short = unit(0.04, "cm"), mid = unit(0.07, "cm"),
                      long = unit(0.11, "cm")) +
  labs(title = "Log cumulative axis")

paper_ggsave(p_log, file.path(OUT_DIR, "fig3b_cdf_log_y"),
             width = NATURE_DOUBLE_IN / 2.2, height = 3.4)
paper_ggsave(
  (p_lin | p_log) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom"),
  file.path(OUT_DIR, "fig3b_cdf_linear_vs_log"),
  width = NATURE_DOUBLE_IN, height = 3.6
)

cat("wrote", file.path(OUT_DIR, "fig3b_cdf_log_y.png"), "\n")
cat("wrote", file.path(OUT_DIR, "fig3b_cdf_linear_vs_log.png"), "\n")
cat(sprintf("\nresolution floor: PTV 1/%d = %.4f%%   non-carrier 1/%d = %.4f%%\n",
            length(r$ptv_carrier_z), 100 / length(r$ptv_carrier_z),
            length(r$non_sample_z), 100 / length(r$non_sample_z)))
for (z in c(-4, -3, -2)) {
  cat(sprintf("  z < %-2g  PTV %6.3f%%  non %6.3f%%  fold %.1f\n", z,
              100 * mean(r$ptv_carrier_z < z), 100 * mean(r$non_sample_z < z),
              mean(r$ptv_carrier_z < z) / mean(r$non_sample_z < z)))
}
