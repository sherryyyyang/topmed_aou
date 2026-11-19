#!/usr/bin/env Rscript
# Protein-coding vs lncRNA: x = Fig 3a underlier %, y = 0.2-wide bins of
# frac_cons_exon_unique_primate (non-overlapping exon bases only).
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

cons <- read_tsv(CONS_TSV, show_col_types = FALSE) |>
  transmute(
    gene_nv = gene_id_stable,
    cons_type = gene_type_splitlncrna,
    frac_cons_exon_unique_primate
  )

# Protein-coding / lncRNA from constraint type (existing + new lncRNA).
# y uses unique-exon primate fraction so overlapping bases are ignored
# instead of dropping the whole gene.
fig3 <- read_tsv(FIG3_TSV, show_col_types = FALSE) |>
  left_join(cons, by = "gene_nv") |>
  mutate(
    map_class = case_when(
      cons_type == "protein_coding" ~ "Protein-coding",
      cons_type %in% c("existing_lncRNA", "new_lncRNA") ~ "lncRNA",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(map_class), is.finite(frac_cons_exon_unique_primate)) |>
  mutate(
    map_class = factor(map_class, levels = c("Protein-coding", "lncRNA")),
    frac_bin = factor(
      case_when(
        frac_cons_exon_unique_primate < 0.2 ~ "0–0.2",
        frac_cons_exon_unique_primate < 0.4 ~ "0.2–0.4",
        frac_cons_exon_unique_primate < 0.6 ~ "0.4–0.6",
        frac_cons_exon_unique_primate < 0.8 ~ "0.6–0.8",
        TRUE ~ "0.8–1.0"
      ),
      levels = c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
    )
  )

PAL <- c(
  "Protein-coding" = "grey55",
  "lncRNA"         = unname(OI["purple"])
)

n_lab <- fig3 |> count(frac_bin, map_class)
cat("genes per frac_cons_exon_unique_primate bin:\n")
print(n_lab, n = 20)
cat("lncRNA gene_type_splitlncrna:\n")
print(fig3 |> filter(map_class == "lncRNA") |> count(cons_type))
summ <- fig3 |>
  group_by(frac_bin, map_class) |>
  summarise(
    n = n(),
    median_pct_under = median(pct_under),
    .groups = "drop"
  )
cat("median pct_under per bin:\n")
print(summ, n = 20)

x_max <- max(fig3$pct_under)
x_step <- 0.2
x_brks <- seq(0, ceiling(x_max / x_step) * x_step, by = x_step)

set.seed(17)
p <- ggplot(fig3, aes(x = pct_under, y = frac_bin, color = map_class)) +
  geom_jitter(
    data = filter(fig3, map_class == "Protein-coding"),
    height = 0.16, width = 0, size = 0.65, alpha = 0.28, stroke = 0
  ) +
  geom_jitter(
    data = filter(fig3, map_class == "lncRNA"),
    height = 0.14, width = 0, size = 1.05, alpha = 0.75, stroke = 0
  ) +
  geom_point(
    data = summ,
    aes(x = median_pct_under, y = frac_bin, color = map_class),
    shape = 23, size = 2.6, fill = "white", stroke = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_x_continuous(
    breaks = x_brks,
    labels = label_percent(scale = 1, accuracy = 0.1),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    x = "Fig 3a y: % of subjects with RNA underliers in gene",
    y = "frac_cons_exon_unique_primate"
  ) +
  theme_paper() +
  theme(
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "top",
    legend.justification = "left",
    plot.margin = margin(8, 12, 8, 8)
  )

out <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
src <- file.path(out, "source_panels")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
stem <- "Fig3a_constraint_sanity_coding_lncrna"
paper_ggsave(p, file.path(out, stem), width = NATURE_DOUBLE_IN, height = 4.2)
paper_ggsave(p, file.path(src, stem), width = NATURE_DOUBLE_IN, height = 4.2)
paper_ggsave(
  p,
  file.path(ROOT, "allgenes_pipeline_out", stem),
  width = NATURE_DOUBLE_IN, height = 4.2
)
cat("wrote", file.path(out, paste0(stem, ".png")), "\n")
