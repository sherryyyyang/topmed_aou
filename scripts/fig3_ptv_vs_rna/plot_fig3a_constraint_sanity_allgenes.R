#!/usr/bin/env Rscript
# Same layout as plot_fig3a_constraint_sanity.R, but Fig 3a is rebuilt from
# gene_median_tpm_kept.tsv with z < -3 underliers called on noncoding genes
# (see call_underliers_allgenes_pc2005.R). Does not overwrite the coding-HVG
# sanity figure.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
Sys.setenv(PC_SETTING = "2005")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

FIG3_TSV <- file.path(
  PC_WORK, "fig3_ptv_vs_rna", "allgenes",
  "ptv_vs_underlier_gene_level_z_-3_allgenes.tsv"
)
stopifnot(file.exists(FIG3_TSV))

CONS_TSV <- file.path(SHEYANG, "annotation", "gene_constraint_fractions.tsv.gz")
cons <- read_tsv(CONS_TSV, show_col_types = FALSE) |>
  transmute(
    gene_nv = gene_id_stable,
    cons_type = gene_type_splitlncrna
  )
cons_ids <- unique(cons$gene_nv)

fig3 <- read_tsv(FIG3_TSV, show_col_types = FALSE) |>
  left_join(cons, by = "gene_nv") |>
  mutate(
    map_class = case_when(
      is.na(cons_type) ~ "Unmapped",
      cons_type == "protein_coding" ~ "Protein-coding",
      cons_type %in% c("existing_lncRNA", "new_lncRNA") ~ "lncRNA",
      cons_type == "pseudogene" ~ "Pseudogene",
      TRUE ~ cons_type
    ),
    map_class = factor(
      map_class,
      levels = c("Protein-coding", "lncRNA", "Pseudogene", "Unmapped")
    )
  )

PAL <- c(
  "Protein-coding" = "grey70",
  "lncRNA"         = unname(OI["purple"]),
  "Pseudogene"     = unname(OI["green"]),
  "Unmapped"       = unname(OI["vermillion"])
)

cov <- fig3 |>
  count(map_class, name = "n") |>
  mutate(n_tot = sum(n), pct = 100 * n / n_tot)

gencode <- read_tsv(
  GENCODE_SYM,
  col_names = c("gene_nv", "symbol"),
  show_col_types = FALSE
) |>
  mutate(gene_nv = str_remove(gene_nv, "\\..*$")) |>
  distinct(gene_nv)

tpm <- read_tsv(file.path(DATA_CORE, "gene_median_tpm_kept.tsv"), show_col_types = FALSE) |>
  mutate(gene_nv = str_remove(gene_nv, "\\..*$")) |>
  distinct(gene_nv, .keep_all = TRUE)
tpm1 <- tpm |> filter(median_tpm > 1)
loeuf <- readRDS(LOEUF_LOOKUP)
loeuf_ids <- unique(str_remove(loeuf$gene_id_no_version, "\\..*$"))
cons_pc <- unique(cons$gene_nv[cons$cons_type == "protein_coding"])

tally_map <- function(ids, label, onto = cons_ids) {
  ids <- unique(na.omit(ids))
  tibble(
    set = label,
    n_tot = length(ids),
    n_map = sum(ids %in% onto)
  )
}

ref <- bind_rows(
  tally_map(gencode$gene_nv, "GENCODE v30"),
  tally_map(tpm$gene_nv, "TOPMed quantified"),
  tally_map(tpm1$gene_nv, "TPM > 1"),
  tally_map(fig3$gene_nv, "Fig 3a plotted"),
  tally_map(loeuf_ids, "LOEUF"),
  tally_map(cons_pc, "Constraint coding \u2192 LOEUF", onto = loeuf_ids)
) |>
  mutate(
    pct = 100 * n_map / n_tot,
    set = factor(set, levels = c(
      "GENCODE v30", "TOPMed quantified", "TPM > 1", "Fig 3a plotted",
      "LOEUF", "Constraint coding \u2192 LOEUF"
    )),
    lab = sprintf("%s / %s\n%.1f%%", comma(n_map), comma(n_tot), pct),
    direction = ifelse(
      set == "Constraint coding \u2192 LOEUF",
      "onto LOEUF",
      "GENCODE mapped\nonto constraint list"
    )
  )

y_brks <- c(1, 10, 100, 1000, 10000)

pa <- ggplot(cov, aes(x = map_class, y = n, fill = map_class)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = sprintf("%s\n%.1f%%", comma(n), pct)),
    vjust = -0.15, size = 2.4, family = PAPER_FONT, lineheight = 0.9
  ) +
  scale_fill_manual(values = PAL, name = NULL) +
  scale_y_continuous(
    trans = "log10",
    breaks = y_brks,
    labels = comma,
    expand = expansion(mult = c(0.02, 0.28))
  ) +
  labs(x = NULL, y = "Fig 3a genes (log)") +
  theme_paper() +
  theme(axis.text.x = element_text(size = 7.5, angle = 25, hjust = 1))

set.seed(17)
pb <- ggplot(fig3, aes(x = pct_under, y = map_class, color = map_class)) +
  geom_jitter(
    data = filter(fig3, map_class == "Protein-coding"),
    height = 0.22, width = 0, size = 1.2, alpha = 0.35, stroke = 0
  ) +
  geom_jitter(
    data = filter(fig3, map_class != "Protein-coding"),
    height = 0.18, width = 0, size = 2.1, alpha = 0.9, stroke = 0
  ) +
  scale_color_manual(values = PAL, guide = "none") +
  scale_x_continuous(
    labels = label_percent(scale = 1, accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(x = "Fig 3a y: % of subjects with RNA\nunderliers in gene", y = NULL) +
  theme_paper() +
  theme(
    axis.title.x = element_text(size = 7.5, lineheight = 0.95),
    plot.margin = margin(4, 10, 8, 4)
  )

pc <- ggplot(ref, aes(x = set, y = pct, fill = direction)) +
  geom_col(width = 0.66) +
  geom_text(
    aes(label = lab),
    vjust = -0.15, size = 2.35, family = PAPER_FONT, lineheight = 0.95
  ) +
  scale_fill_manual(
    values = c(
      "GENCODE mapped\nonto constraint list" = COL_TOPMED,
      "onto LOEUF" = unname(OI["orange"])
    ),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 118),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(v) paste0(v, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "% of query list that overlaps"
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(size = 7.2, angle = 22, hjust = 1))

fig <- ((pa | pb) / pc) +
  plot_layout(heights = c(1.15, 0.95), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = PAPER_FONT, face = "bold", size = 14, colour = "black"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.text = element_text(size = 7, lineheight = 0.95),
    plot.margin = margin(6, 12, 10, 6)
  )

out <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
src <- file.path(out, "source_panels")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
stem <- "Fig3a_constraint_sanity_allgenes"
paper_ggsave(fig, file.path(out, stem), width = NATURE_DOUBLE_IN, height = 6.4)
paper_ggsave(fig, file.path(src, stem), width = NATURE_DOUBLE_IN, height = 6.4)
file.copy(
  FIG3_TSV,
  file.path(out, "supplemental_data", "Source_Data_Fig3a_constraint_sanity_allgenes.tsv"),
  overwrite = TRUE
)

cat("Fig 3a (all genes / noncoding included) mapped classes:\n")
print(cov, n = 20)
cat("\nReference lists vs constraint:\n")
print(ref, n = 20)
cat("wrote", file.path(out, paste0(stem, ".png")), "\n")
