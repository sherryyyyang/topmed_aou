#!/usr/bin/env Rscript
# Two 1D strips in the Figure 3a point style (PC2005 RNA underliers).
#
# a  Exon mammal phyloP conserved fraction: protein-coding (grey) vs lncRNA.
#    CDS fractions are blank for lncRNA, so exons are the shared score.
# b  TOPMed underlier prevalence for lncRNAs quantified in the cohort.
#
# Universe: genes in gene_median_tpm_kept.tsv (present in TOPMed RNA).
# Underliers: z < -3 AND log2CPM drop, PC_SETTING expression PCs.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

CONS_TSV <- file.path(
  SHEYANG, "annotation", "gene_constraint_fractions.tsv.gz"
)
UNDER_RDS <- file.path(DATA_UNDERLIER, "z_score_underlier_z_-3_clean.rds")
COL_LNCRNA <- unname(OI["purple"])
COL_CODING <- "grey70"

cons <- read_tsv(CONS_TSV, show_col_types = FALSE) |>
  mutate(
    gene_nv = gene_id_stable,
    gene_class = case_when(
      gene_type_splitlncrna == "protein_coding" ~ "Protein-coding",
      gene_type_splitlncrna %in% c("existing_lncRNA", "new_lncRNA") ~ "lncRNA",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(gene_class))

tpm <- read_tsv(file.path(DATA_CORE, "gene_median_tpm_kept.tsv"), show_col_types = FALSE) |>
  distinct(gene_nv, .keep_all = TRUE)
kept <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept)

under_g <- read_rds(UNDER_RDS) |>
  filter(subjectid %in% kept) |>
  mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
  distinct(subjectid, gene_nv) |>
  count(gene_nv, name = "n_underlier_subjects")

g <- cons |>
  inner_join(tpm, by = "gene_nv") |>
  left_join(under_g, by = "gene_nv") |>
  mutate(
    n_underlier_subjects = ifelse(is.na(n_underlier_subjects), 0L, n_underlier_subjects),
    pct_under = 100 * n_underlier_subjects / n_kept,
    cons_exon = frac_cons_exon_mammal
  )

left_df <- g |>
  filter(is.finite(cons_exon)) |>
  mutate(gene_class = factor(gene_class, levels = c("Protein-coding", "lncRNA")))

right_df <- g |>
  filter(gene_class == "lncRNA") |>
  mutate(gene_class = factor("lncRNA", levels = "lncRNA"))

n_coding <- sum(left_df$gene_class == "Protein-coding")
n_lnc_left <- sum(left_df$gene_class == "lncRNA")
n_lnc_right <- nrow(right_df)
n_lnc_under <- sum(right_df$n_underlier_subjects > 0)

cat(sprintf(
  "PC_SETTING=%s  n_kept=%s\n  left: coding=%s  lncRNA=%s (finite exon mammal frac)\n  right: lncRNA in TOPMed=%s  with >=1 underlier=%s\n",
  PC_SETTING, comma(n_kept), comma(n_coding), comma(n_lnc_left),
  comma(n_lnc_right), comma(n_lnc_under)
))

set.seed(17)
pa <- ggplot(left_df, aes(x = cons_exon, y = gene_class, color = gene_class)) +
  geom_jitter(height = 0.22, width = 0, size = 1.3, alpha = 0.45, stroke = 0) +
  scale_color_manual(
    values = c("Protein-coding" = COL_CODING, "lncRNA" = COL_LNCRNA),
    labels = c(
      sprintf("Protein-coding (n=%s)", comma(n_coding)),
      sprintf("lncRNA (n=%s)", comma(n_lnc_left))
    ),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Fraction conserved (exon phyloP, mammal)",
    y = NULL
  ) +
  theme_paper() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.y = element_text(size = 8)
  )

pb <- ggplot(right_df, aes(x = pct_under, y = gene_class)) +
  geom_jitter(
    height = 0.22, width = 0, size = 1.3, alpha = 0.45, stroke = 0,
    color = COL_LNCRNA
  ) +
  scale_x_continuous(
    labels = label_percent(scale = 1, accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = "% of subjects with RNA\nunderliers in gene",
    y = NULL
  ) +
  theme_paper() +
  theme(axis.text.y = element_text(size = 8))

fig <- (pa | pb) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = PAPER_FONT, face = "bold", size = 14, colour = "black"
    )
  )

out_dir <- FINAL_DRAFT
src_pan <- file.path(FINAL_DRAFT, "source_panels")
supp <- file.path(FINAL_DRAFT, "supplemental_data")
dir.create(src_pan, recursive = TRUE, showWarnings = FALSE)
dir.create(supp, recursive = TRUE, showWarnings = FALSE)

paper_ggsave(fig, file.path(out_dir, "Fig_lncRNA_phylop_prevalence"),
             width = NATURE_DOUBLE_IN, height = 3.1)
paper_ggsave(fig, file.path(src_pan, "Fig_lncRNA_phylop_prevalence"),
             width = NATURE_DOUBLE_IN, height = 3.1)
paper_ggsave(pa, file.path(src_pan, "Fig_lncRNA_a_conservation"),
             width = NATURE_DOUBLE_IN * 0.55, height = 3.1)
paper_ggsave(pb, file.path(src_pan, "Fig_lncRNA_b_prevalence"),
             width = NATURE_DOUBLE_IN * 0.45, height = 3.1)

write_tsv(
  g |>
    transmute(
      gene_nv, gene_class, gene_type_splitlncrna, lncrna_category,
      frac_cons_exon_mammal, frac_cons_exon_primate,
      n_underlier_subjects, pct_under, median_tpm
    ) |>
    arrange(gene_class, desc(pct_under)),
  file.path(supp, "Source_Data_Fig_lncRNA_phylop_prevalence.tsv")
)

.paper_plots <- list(a = pa, b = pb, fig = fig)
cat("wrote", file.path(out_dir, "Fig_lncRNA_phylop_prevalence.png"), "\n")
