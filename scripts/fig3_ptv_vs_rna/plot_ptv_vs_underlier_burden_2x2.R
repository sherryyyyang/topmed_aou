#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG3X


ORDER <- c("haplo", "z_-2", "z_-3", "z_-4")
LABELS <- c(
  "haplo" = "Haplo (log2CPM-only)",
  "z_-2" = "Strict z<-2 + log2CPM",
  "z_-3" = "Strict z<-3 + log2CPM",
  "z_-4" = "Strict z<-4 + log2CPM"
)

all_df <- bind_rows(lapply(ORDER, function(m) {
  f <- paste0("ptv_vs_underlier_burden_", m, ".tsv")
  d <- read_tsv(f, show_col_types = FALSE)
  d$mode <- m
  d
})) |>
  mutate(mode = factor(mode, levels = ORDER, labels = LABELS[ORDER]))

x_max <- max(all_df$n_ptv)
y_max <- max(all_df$n_underlier_genes)
overlap_max <- max(all_df$n_ptv_under_overlap_genes)
overlap_levels <- as.character(0:overlap_max)

overlap_colors <- setNames(
  rev(hcl.colors(length(overlap_levels), palette = "viridis")),
  overlap_levels
)

rho_df <- all_df |>
  group_by(mode) |>
  summarise(rho = suppressWarnings(cor(n_ptv, n_underlier_genes, method = "spearman")), .groups = "drop")

mode_labels <- levels(all_df$mode)
rho_lines <- sprintf(
  "%s: rho=%.3f",
  mode_labels,
  rho_df$rho[match(mode_labels, rho_df$mode)]
)
sub_txt <- paste0(
  "Shared x-axis, panel-specific y-axis.\n",
  paste(rho_lines, collapse = "\n")
)
cap_txt <- str_wrap(
  "Color: per-subject count of genes with both a PTV and an RNA underlier.",
  width = 95
)

p <- ggplot(all_df, aes(x = n_ptv, y = n_underlier_genes)) +
  geom_point(
    aes(color = factor(n_ptv_under_overlap_genes, levels = overlap_levels)),
    alpha = 0.45,
    size = 1.3,
    stroke = 0
  ) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.3, color = "grey25") +
  facet_wrap(~mode, ncol = 2, scales = "free_y") +
  scale_color_manual(
    values = overlap_colors,
    breaks = overlap_levels,
    labels = overlap_levels,
    name = "Genes with\nPTV & underlier",
    drop = FALSE
  ) +
  scale_x_continuous(
    labels = comma,
    limits = c(0, x_max),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    title = str_wrap(
      "PTV burden vs RNA underlier burden (haplo, z<-2, z<-3, z<-4)",
      width = 70
    ),
    subtitle = sub_txt,
    caption = cap_txt,
    x = "Number of distinct PTVs carried",
    y = "Number of underlier genes"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey25", lineheight = 1.15),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    plot.caption = element_text(size = 7, color = "grey35", hjust = 0),
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )

out_png <- "ptv_vs_underlier_burden_2x2_haplo_z2_z3_z4.png"
out_pdf <- "ptv_vs_underlier_burden_2x2_haplo_z2_z3_z4.pdf"
ggsave(file.path(OUT_DIR, out_png), p, width = 13.0, height = 9.8, dpi = 220)
ggsave(file.path(OUT_DIR, out_pdf), p, width = 13.0, height = 9.8, device = cairo_pdf)
cat("Wrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")
