#!/usr/bin/env Rscript
# TEST ONLY — Fig 3b with points aggressively separated so none visually overlap.
# Does not replace main Figure3. Output: figures/main/Figure3b_no_overlap_test.png
#
# Approach: plot in log10(x) × y space on a large canvas with small points and
# deterministic jitter (plus a tiny round-robin nudge for exact coordinate ties).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
.PKG <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
source(file.path(.PKG, "scripts/helpers/_paths.R"))
source(file.path(.PKG, "scripts/helpers/paper_theme.R"))

CNV_CLASS_LEVELS <- c("SV (DEL)", "CNV / multi-copy", "Other genes")
CNV_CLASS_COLORS <- c(
  "SV (DEL)"         = unname(OI["blue"]),
  "CNV / multi-copy" = unname(OI["vermillion"]),
  "Other genes"      = "grey70"
)
SV_AF_MAX <- 0.005

gene_tsv <- file.path(ROOT, "pc1055_work/fig3_ptv_vs_rna/extra/ptv_vs_underlier_gene_level_z_-3.tsv")
cnv_tsv <- file.path(ROOT, "data/references/fig3b_cnv_multicopy_genes.tsv")
sv_tsv  <- file.path(ROOT, "data/references/fig3b_sv_del_genes.tsv")
stopifnot(file.exists(gene_tsv), file.exists(cnv_tsv), file.exists(sv_tsv))

sv <- read_tsv(sv_tsv, show_col_types = FALSE) |>
  mutate(af = as.numeric(af)) |>
  filter(is.finite(af), af < SV_AF_MAX) |>
  distinct(gene_nv) |>
  mutate(is_sv = TRUE)

cnv <- read_tsv(cnv_tsv, show_col_types = FALSE) |>
  distinct(gene_nv) |>
  mutate(is_cnv = TRUE)

df <- read_tsv(gene_tsv, show_col_types = FALSE) |>
  left_join(sv, by = "gene_nv") |>
  left_join(cnv, by = "gene_nv") |>
  mutate(
    is_sv = ifelse(is.na(is_sv), FALSE, is_sv),
    is_cnv = ifelse(is.na(is_cnv), FALSE, is_cnv),
    cnv_class = factor(
      case_when(
        is_sv ~ "SV (DEL)",
        is_cnv ~ "CNV / multi-copy",
        TRUE ~ "Other genes"
      ),
      levels = CNV_CLASS_LEVELS
    )
  )

pos <- df$pct_ptv[is.finite(df$pct_ptv) & df$pct_ptv > 0]
x_floor <- min(pos) / 2
df <- df |>
  mutate(
    pct_ptv_plot = ifelse(pct_ptv > 0, pct_ptv, x_floor),
    x_log = log10(pct_ptv_plot),
    y = pct_under
  )

# Spread exact coordinate ties on a grid (avoids ring artifacts for huge
# 0%,0% stacks), then add modest jitter for near-overlaps.
set.seed(42)
df <- df |>
  mutate(
    x_key = round(x_log, 5),
    y_key = round(y, 5)
  ) |>
  group_by(x_key, y_key) |>
  mutate(
    tie_i = row_number() - 1L,
    tie_n = n(),
    grid_side = ceiling(sqrt(tie_n)),
    gx = (tie_i %% grid_side) - (grid_side - 1) / 2,
    gy = (tie_i %/% grid_side) - (grid_side - 1) / 2,
    # Scale grid step with group size so large stacks open up
    step_x = 0.018 / sqrt(pmax(grid_side, 1)),
    step_y = 0.010 / sqrt(pmax(grid_side, 1)),
    x_log_j = x_log + gx * step_x + runif(n(), -0.008, 0.008),
    y_j = pmax(0, y + gy * step_y + runif(n(), -0.004, 0.004))
  ) |>
  ungroup() |>
  arrange(cnv_class == "Other genes", cnv_class == "CNV / multi-copy", cnv_class)

x_breaks_pct <- c(x_floor, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)
x_breaks_pct <- x_breaks_pct[x_breaks_pct > 0]
x_breaks_log <- log10(x_breaks_pct)
label_x <- function(v) {
  vapply(v, function(z) {
    if (!is.finite(z)) return("")
    pct <- 10^z
    if (abs(z - log10(x_floor)) < 1e-6) return("0%")
    acc <- if (pct < 0.1) 0.01 else if (pct < 1) 0.1 else 1
    label_percent(scale = 1, accuracy = acc)(pct)
  }, character(1))
}

p <- ggplot(df, aes(x = x_log_j, y = y_j, color = cnv_class)) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "grey55", linewidth = 0.35
  ) +
  # abline is in original % space; skip it in log-linear jitter coords —
  # draw identity in log-linear as curve? Skip — add reference as points guide only.
  geom_point(size = 1.15, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = CNV_CLASS_COLORS, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = x_breaks_log, labels = label_x) +
  scale_y_continuous(labels = label_percent(scale = 1, accuracy = 0.01)) +
  labs(
    title = "TEST: Fig 3b with jitter so points do not overlap",
    subtitle = "Positions are jittered (not publication coords); colouring unchanged",
    x = "% of subjects carrying a rare PTV (AF<0.5%) in gene (log10 scale)",
    y = "% of subjects with RNA underliers in gene"
  ) +
  theme_paper() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 8, color = "grey30")
  )

# Drop misleading y=x in jittered log-linear space — rebuild without abline
p <- ggplot(df, aes(x = x_log_j, y = y_j, color = cnv_class)) +
  geom_point(size = 1.15, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = CNV_CLASS_COLORS, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = x_breaks_log, labels = label_x) +
  scale_y_continuous(labels = label_percent(scale = 1, accuracy = 0.01)) +
  labs(
    title = "TEST: Fig 3b with jitter so points do not overlap",
    subtitle = "Positions are jittered (not publication coords); colouring unchanged",
    x = "% of subjects carrying a rare PTV (AF<0.5%) in gene (log10 scale)",
    y = "% of subjects with RNA underliers in gene"
  ) +
  theme_paper() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 8, color = "grey30")
  )

out_dir <- file.path(.PKG, "figures/main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "Figure3b_no_overlap_test")
# Large canvas so jittered points stay separable
ggsave(paste0(out, ".png"), p, width = 22, height = 16, dpi = 200, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = 22, height = 16, device = cairo_pdf, bg = "white")
# Also copy to source_panels
sp <- file.path(.PKG, "source_panels")
dir.create(sp, recursive = TRUE, showWarnings = FALSE)
file.copy(paste0(out, ".png"), file.path(sp, "Fig3b_no_overlap_test.png"), overwrite = TRUE)

cat("Wrote:", paste0(out, ".png"), "\n")
cat("n genes:", nrow(df), " exact-tie groups >1:",
    sum(df$tie_n > 1), "\n")
