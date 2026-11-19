#!/usr/bin/env Rscript
# PTV burden vs RNA underlier burden (per subject).
# Underliers: strict z-score threshold + logCPM vs gene mean (same as pipeline).
# Default z cutoff matches mentor request: z < -3.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG3X


MODE_LABEL <- Sys.getenv("UNDERLIER_LABEL", unset = "z_-3")
if (!MODE_LABEL %in% c("z_-2", "z_-3", "z_-4", "haplo")) {
  stop("UNDERLIER_LABEL must be one of z_-2, z_-3, z_-4, haplo (got ", MODE_LABEL, ")")
}

under_rds <- if (MODE_LABEL == "haplo") {
  file.path(DATA_UNDERLIER, "z_score_underlier_haplo_clean.rds")
} else {
  file.path(DATA_UNDERLIER, paste0("z_score_underlier_", MODE_LABEL, "_clean.rds"))
}
mode_txt <- if (MODE_LABEL == "haplo") {
  "haplo (log2CPM only)"
} else {
  paste0("strict z<", sub("^z_", "", MODE_LABEL), " + log2CPM")
}

PTV_TSV <- Sys.getenv(
  "PTV_TSV",
  unset = "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/ptv_all.tsv"
)

cat("Loading kept_subjects.rds ...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept_subjects)

cat("Loading PTV table:", PTV_TSV, "\n")
ptv <- read_tsv(PTV_TSV, show_col_types = FALSE) |>
  filter(sample %in% kept_subjects)

cat("Loading underliers:", under_rds, "\n")
under <- read_rds(under_rds) |>
  filter(subjectid %in% kept_subjects)

cat("  PTV rows (kept subjects):", nrow(ptv), "\n")
cat("  Underlier rows:", nrow(under), "\n")

ptv_tab <- ptv |>
  mutate(variant_id = paste(chrom, pos, ref, alt, sep = ":")) |>
  group_by(sample) |>
  summarise(n_ptv = n_distinct(variant_id), .groups = "drop") |>
  rename(subjectid = sample)

under_tab <- under |>
  group_by(subjectid) |>
  summarise(n_underlier_genes = n_distinct(gene_id), .groups = "drop")

ptv_genes <- ptv |>
  mutate(gene_nv = str_remove(gene, "\\..*$")) |>
  distinct(subjectid = sample, gene_nv)

under_genes <- under |>
  mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
  distinct(subjectid, gene_nv)

overlap_tab <- ptv_genes |>
  inner_join(under_genes, by = c("subjectid", "gene_nv")) |>
  count(subjectid, name = "n_ptv_under_overlap_genes")

burden <- tibble(subjectid = kept_subjects) |>
  left_join(ptv_tab, by = "subjectid") |>
  left_join(under_tab, by = "subjectid") |>
  left_join(overlap_tab, by = "subjectid") |>
  mutate(
    n_ptv = replace_na(n_ptv, 0L),
    n_underlier_genes = replace_na(n_underlier_genes, 0L),
    n_ptv_under_overlap_genes = replace_na(n_ptv_under_overlap_genes, 0L)
  )

rho <- suppressWarnings(cor(burden$n_ptv, burden$n_underlier_genes, method = "spearman"))

overlap_max <- max(burden$n_ptv_under_overlap_genes)
legend_max <- max(overlap_max, 1L)
overlap_breaks <- seq.int(0L, overlap_max, by = 1L)
overlap_levels <- as.character(overlap_breaks)
overlap_colors <- setNames(
  rev(hcl.colors(length(overlap_levels), palette = "viridis")),
  overlap_levels
)

out_tsv <- paste0("ptv_vs_underlier_burden_", MODE_LABEL, ".tsv")
write_tsv(burden, out_tsv)
cat("Wrote:", normalizePath(out_tsv), "\n")

sub_txt <- str_wrap(
  sprintf(
    "QC cohort as LOEUF figures (n=%s). Spearman rho (PTV vs underlier burden)=%.3f.",
    comma(n_kept),
    rho
  ),
  width = 72
)
cap_txt <- str_wrap(
  paste0(
    "Color: distinct genes with both a PTV and an RNA underlier. ",
    "Discrete colors; darker = higher overlap count. ",
    "Weak correlation / shallow OLS line mean PTV burden and underlier burden are not tightly coupled across people ",
    "(underliers are not redundant with a simple PTV count); ",
    "a slope near 1 would only mean similar linear scaling of counts, not that variants explain all underlier signal."
  ),
  width = 78
)

p <- ggplot(burden, aes(x = n_ptv, y = n_underlier_genes)) +
  geom_point(
    aes(color = factor(n_ptv_under_overlap_genes, levels = overlap_levels)),
    alpha = 0.45,
    size = 1.7,
    stroke = 0
  ) +
  scale_color_manual(
    values = overlap_colors,
    name = "Genes with\nPTV & underlier",
    breaks = overlap_levels,
    labels = overlap_levels,
    drop = FALSE
  ) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.35, color = "grey25") +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.04))) +
  labs(
    title = sprintf("PTV burden vs RNA underlier burden (%s)", mode_txt),
    subtitle = sub_txt,
    caption = cap_txt,
    x = "Number of distinct PTVs carried",
    y = "Number of underlier genes"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    plot.subtitle = element_text(size = 8.5, color = "grey25"),
    plot.caption = element_text(size = 7, color = "grey35", hjust = 0),
    plot.margin = margin(5.5, 14, 5.5, 5.5),
    legend.position = "right",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8)
  )

out_png <- paste0("ptv_vs_underlier_burden_", MODE_LABEL, ".png")
out_pdf <- paste0("ptv_vs_underlier_burden_", MODE_LABEL, ".pdf")
ggsave(file.path(OUT_DIR, out_png), p, width = 8.0, height = 6.8, dpi = 220)
ggsave(file.path(OUT_DIR, out_pdf), p, width = 8.0, height = 6.8, device = cairo_pdf)
cat("Wrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")
