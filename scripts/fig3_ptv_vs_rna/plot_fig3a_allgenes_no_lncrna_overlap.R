#!/usr/bin/env Rscript
# All-gene-HVG Fig 3a scatter (2005 PCs), after dropping GENCODE lncRNAs
# whose exons overlap a protein-coding exon. Does not overwrite the
# published coding-HVG Fig3a_ptv_vs_underlier panel.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

FIG3_TSV <- file.path(
  ROOT, "allgenes_pipeline_out",
  "ptv_vs_underlier_gene_level_z_-3_allgenes_pc2005.tsv"
)
CONS_TSV <- file.path(SHEYANG, "annotation", "gene_constraint_fractions.tsv.gz")
LNC_KEEP <- file.path(
  ROOT, "allgenes_pipeline_out",
  "gencode_v30_lncrna_no_pc_exon_overlap.txt"
)
lnc_keep <- read_lines(LNC_KEEP)
cons <- read_tsv(CONS_TSV, show_col_types = FALSE) |>
  transmute(gene_nv = gene_id_stable, cons_type = gene_type_splitlncrna)

fig3 <- read_tsv(FIG3_TSV, show_col_types = FALSE) |>
  left_join(cons, by = "gene_nv") |>
  mutate(
    map_class = case_when(
      is.na(cons_type) ~ "Unmapped",
      cons_type == "protein_coding" ~ "Protein-coding",
      gene_nv %in% lnc_keep ~ "lncRNA",
      cons_type == "pseudogene" ~ "Pseudogene",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(map_class))

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
fig3 <- fig3 |>
  mutate(
    LOEUF_tier = factor(
      ifelse(is.na(LOEUF_tier), "No LOEUF", as.character(LOEUF_tier)),
      levels = TIER_LEVELS
    )
  )

cat("Fig 3a genes after dropping PC-exon-overlapping lncRNA:\n")
print(fig3 |> count(map_class, name = "n"), n = 20)
cat("n =", nrow(fig3), "\n")

x_max <- max(fig3$pct_ptv)
p <- ggplot(fig3, aes(x = pct_ptv, y = pct_under)) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "grey55", linewidth = 0.35
  ) +
  geom_point(aes(color = LOEUF_tier), alpha = 0.5, size = 1.3, stroke = 0) +
  scale_color_manual(values = TIER_COLORS, name = "LOEUF tier", drop = FALSE) +
  scale_x_continuous(
    labels = label_percent(scale = 1, accuracy = 0.01),
    limits = c(0, x_max),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    labels = label_percent(scale = 1, accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "% of subjects carrying a rare PTV (AF<0.5%) in gene",
    y = "% of subjects with RNA\nunderliers in gene"
  ) +
  theme_paper() +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 4, r = 6, b = 2, l = 16, unit = "mm"),
    axis.title.y = element_text(
      angle = 90, lineheight = 0.92,
      margin = margin(r = 3, unit = "mm")
    )
  )

out <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
src <- file.path(out, "source_panels")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
stem <- "Fig3a_ptv_vs_underlier_allgenes_no_lncrna_overlap"
paper_ggsave(p, file.path(out, stem), width = NATURE_DOUBLE_IN, height = 4.6)
paper_ggsave(p, file.path(src, stem), width = NATURE_DOUBLE_IN, height = 4.6)
write_tsv(
  fig3,
  file.path(out, "supplemental_data", paste0("Source_Data_", stem, ".tsv"))
)
cat("wrote", file.path(out, paste0(stem, ".png")), "\n")
