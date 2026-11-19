#!/usr/bin/env Rscript
# EXPERIMENT (not wired into assemble_nature_figures.R): Figure 1b with a
# log10 y-axis, side by side with the linear version for comparison.
#
# Same data and same shared-top-10% labels as
# fig1_overview/plot_underlier_tallies_topmed_aou.R.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(scales)
  library(patchwork); library(stringr)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

OUT_DIR <- file.path(ROOT, "fig1_overview", "experiments")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

top_g <- read_rds(file.path(DATA_UNDERLIER, "z_score_tally_z_-3_clean.rds"))
aou <- read_tsv(
  file.path(AOU_RNA, "rna_outlier_prevalence_per_gene_haplo_z_-3.tsv"),
  show_col_types = FALSE
) |>
  mutate(
    pct_underlier = n_underlier_subjects / n_subjects_tested,
    gene_nv = str_remove(gene_id, "\\..*$")
  ) |>
  arrange(desc(pct_underlier)) |>
  mutate(rank = row_number(), cohort = "AoU")

sym <- tryCatch(
  read_tsv(GENCODE_SYM, col_names = c("gene_nv", "symbol"),
           show_col_types = FALSE),
  error = function(e) NULL
)
top_gene <- top_g |>
  transmute(rank, pct_underlier, gene_id, cohort = "TOPMed") |>
  mutate(gene_nv = str_remove(gene_id, "\\..*$"))
top_gene <- if (!is.null(sym)) left_join(top_gene, sym, by = "gene_nv") else
  mutate(top_gene, symbol = gene_nv)

gene_long <- bind_rows(
  top_gene |> transmute(rank, pct_underlier, symbol, gene_nv, cohort),
  aou |> transmute(rank, pct_underlier, symbol, gene_nv, cohort)
) |>
  mutate(cohort = factor(cohort, levels = c("TOPMed", "AoU")))

shared <- inner_join(
  top_gene |> transmute(gene_nv, symbol_tm = symbol, rank_tm = rank,
                        pct_tm = pct_underlier),
  aou |> transmute(gene_nv, symbol_aou = symbol, rank_aou = rank,
                   pct_aou = pct_underlier),
  by = "gene_nv"
) |>
  mutate(symbol = coalesce(na_if(symbol_tm, ""), na_if(symbol_aou, ""), gene_nv))
lab_ids <- shared |>
  filter(rank_tm <= floor(0.10 * nrow(top_gene)),
         rank_aou <= floor(0.10 * nrow(aou))) |>
  mutate(geo = sqrt(pct_tm * pct_aou)) |>
  arrange(desc(geo)) |>
  slice_head(n = 9)

lab_pts <- bind_rows(
  lab_ids |> transmute(rank = rank_tm, pct_underlier = pct_tm, symbol,
                       cohort = "TOPMed"),
  lab_ids |> transmute(rank = rank_aou, pct_underlier = pct_aou, symbol,
                       cohort = "AoU")
) |>
  mutate(cohort = factor(cohort, levels = c("TOPMed", "AoU")))
lab_g <- lab_ids |>
  mutate(rank = if_else(pct_tm >= pct_aou, rank_tm, rank_aou),
         pct_underlier = pmax(pct_tm, pct_aou))

med <- gene_long |>
  group_by(cohort) |>
  summarise(med = median(pct_underlier), .groups = "drop")
pal <- c(TOPMed = COL_TOPMED, AoU = COL_AOU)
x_max <- max(gene_long$rank, na.rm = TRUE)
n_lab <- nrow(lab_g)

build <- function(log_y) {
  # Leader targets are laid out in whichever space the axis uses, so the
  # label column stays evenly spaced on the printed page either way.
  lg <- lab_g |>
    mutate(
      x_txt = 0.48 * x_max,
      y_txt = if (log_y) {
        10^seq(log10(0.0265), log10(0.0092), length.out = n_lab)
      } else {
        max(lab_g$pct_underlier) * seq(0.94, 0.36, length.out = n_lab)
      }
    )
  p <- ggplot(gene_long, aes(rank, pct_underlier, color = cohort)) +
    geom_line(linewidth = 0.55) +
    geom_hline(data = med, aes(yintercept = med, color = cohort),
               linetype = "dashed", linewidth = 0.35, show.legend = FALSE) +
    geom_segment(data = lg,
                 aes(x = rank, y = pct_underlier, xend = x_txt, yend = y_txt),
                 color = "grey35", linewidth = 0.25, alpha = 0.7,
                 show.legend = FALSE) +
    geom_point(data = lab_pts, size = 1.35, show.legend = FALSE) +
    geom_label(data = lg, aes(x = x_txt, y = y_txt, label = symbol),
               inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 2.5,
               fontface = "italic", family = PAPER_FONT, linewidth = 0,
               fill = "white", label.padding = unit(0.08, "lines"),
               color = "grey10", show.legend = FALSE) +
    scale_color_manual(values = pal, name = NULL) +
    scale_x_continuous(labels = comma,
                       expand = expansion(mult = c(0.02, 0.06))) +
    labs(x = "Gene rank (descending prevalence)",
         y = if (log_y) "Underlier frequency (log scale)" else
           "Underlier frequency") +
    theme_paper() + theme_paper_inner_title() +
    theme(legend.position = "top")
  if (log_y) {
    brks <- c(0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03)
    p <- p + scale_y_log10(breaks = brks,
                           labels = percent_format(accuracy = 0.01)(brks))
  } else {
    p <- p + scale_y_continuous(labels = percent_format(accuracy = 0.1),
                                expand = expansion(mult = c(0.02, 0.08)))
  }
  p + labs(title = if (log_y) "Genes (log y)" else "Genes (linear y)")
}

p_lin <- build(FALSE)
p_log <- build(TRUE)

paper_ggsave(p_log, file.path(OUT_DIR, "fig1b_log_scale"),
             width = NATURE_DOUBLE_IN / 2, height = 3.4)
paper_ggsave(p_lin | p_log, file.path(OUT_DIR, "fig1b_linear_vs_log"),
             width = NATURE_DOUBLE_IN, height = 3.6)
cat("wrote", file.path(OUT_DIR, "fig1b_log_scale.png"), "\n")
cat("wrote", file.path(OUT_DIR, "fig1b_linear_vs_log.png"), "\n")
cat("min non-zero frequency: TOPMed",
    sprintf("%.4f%%", 100 * min(top_gene$pct_underlier)),
    " AoU", sprintf("%.4f%%", 100 * min(aou$pct_underlier)), "\n")
