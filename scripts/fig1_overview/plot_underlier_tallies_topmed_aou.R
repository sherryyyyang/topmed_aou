#!/usr/bin/env Rscript
# Figure 1 rank curves for the paper:
#   A. Genes ranked by haplo+z<-3 underlier rate in TOPMed and in AoU
#   B. TOPMed subjects ranked by # underlier genes (no sample IDs)
#
# AoU subject-level haplo+z<-3 lists are not in the released gene-prevalence
# table, so panel B is TOPMed-only.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(scales)
  library(patchwork); library(ggrepel); library(stringr)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG1B

PW <- ROOT
OUT_DIR <- OUT_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_TOP <- 5489L
N_AOU <- 8895L

top_g <- read_rds(file.path(DATA_UNDERLIER, "z_score_tally_z_-3_clean.rds"))
top_s <- read_rds(file.path(DATA_UNDERLIER, "z_score_tally_subj_z_-3_clean.rds"))
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

top_gene <- top_g |>
  transmute(rank, pct_underlier, gene_id, cohort = "TOPMed")
# symbol for labels
sym <- tryCatch(
  read_tsv(GENCODE_SYM,
           col_names = c("gene_nv", "symbol"), show_col_types = FALSE),
  error = function(e) NULL
)
if (!is.null(sym)) {
  top_gene <- top_gene |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
    left_join(sym, by = "gene_nv")
} else {
  top_gene$symbol <- top_gene$gene_id
}

gene_long <- bind_rows(
  top_gene |> transmute(rank, pct_underlier, symbol, gene_nv, cohort),
  aou |> transmute(rank, pct_underlier, symbol, gene_nv, cohort)
) |>
  mutate(cohort = factor(cohort, levels = c("TOPMed", "AoU")))

# Labels: genes in the top 10% of BOTH catalogues, then the top 9 by
# geometric mean of the two frequencies (shared hot genes, not cohort-private tops).
shared <- inner_join(
  top_gene |> transmute(gene_nv, symbol_tm = symbol, rank_tm = rank, pct_tm = pct_underlier),
  aou |> transmute(gene_nv, symbol_aou = symbol, rank_aou = rank, pct_aou = pct_underlier),
  by = "gene_nv"
) |>
  mutate(symbol = coalesce(na_if(symbol_tm, ""), na_if(symbol_aou, ""), gene_nv))
cut_tm <- floor(0.10 * nrow(top_gene))
cut_aou <- floor(0.10 * nrow(aou))
lab_ids <- shared |>
  filter(rank_tm <= cut_tm, rank_aou <= cut_aou) |>
  mutate(geo = sqrt(pct_tm * pct_aou)) |>
  arrange(desc(geo)) |>
  slice_head(n = 9)
write_tsv(
  lab_ids |> transmute(gene_nv, symbol, rank_tm, pct_tm, rank_aou, pct_aou, geo),
  file.path(OUT_DIR, "fig1b_shared_top10pct_labels.tsv")
)
lab_pts <- bind_rows(
  lab_ids |> transmute(rank = rank_tm, pct_underlier = pct_tm, symbol, cohort = "TOPMed"),
  lab_ids |> transmute(rank = rank_aou, pct_underlier = pct_aou, symbol, cohort = "AoU")
) |>
  mutate(cohort = factor(cohort, levels = c("TOPMed", "AoU")))
lab_g <- lab_ids |>
  mutate(
    rank = if_else(pct_tm >= pct_aou, rank_tm, rank_aou),
    pct_underlier = pmax(pct_tm, pct_aou)
  )

med <- gene_long |>
  group_by(cohort) |>
  summarise(med = median(pct_underlier), .groups = "drop")

pal <- c(TOPMed = COL_TOPMED, AoU = COL_AOU)
x_max <- max(gene_long$rank, na.rm = TRUE)
y_max <- max(lab_g$pct_underlier, na.rm = TRUE)
n_lab <- nrow(lab_g)
lab_g <- lab_g |>
  mutate(
    x_txt = 0.48 * x_max,
    y_txt = y_max * seq(0.94, 0.36, length.out = n_lab)
  )

p_gene <- ggplot(gene_long, aes(rank, pct_underlier, color = cohort)) +
  geom_line(linewidth = 0.55) +
  geom_hline(data = med, aes(yintercept = med, color = cohort),
             linetype = "dashed", linewidth = 0.35, show.legend = FALSE) +
  geom_segment(
    data = lab_g, aes(x = rank, y = pct_underlier, xend = x_txt, yend = y_txt),
    color = "grey35", linewidth = 0.25, alpha = 0.7, show.legend = FALSE
  ) +
  geom_point(data = lab_pts, size = 1.35, show.legend = FALSE) +
  geom_label(
    data = lab_g, aes(x = x_txt, y = y_txt, label = symbol),
    inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 2.5, fontface = "italic",
    family = PAPER_FONT, linewidth = 0, fill = "white", label.padding = unit(0.08, "lines"),
    color = "grey10", show.legend = FALSE
  ) +
  scale_color_manual(values = pal, name = NULL) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Gene rank (descending prevalence)",
       y = "Underlier frequency") +
  theme_paper() +
  theme_paper_inner_title() +
  labs(title = "Genes") +
  theme(legend.position = "top")

med_s <- median(top_s$n)
p_subj <- ggplot(top_s, aes(rank, n)) +
  geom_line(color = COL_TOPMED, linewidth = 0.5) +
  geom_hline(yintercept = med_s, linetype = "dashed",
             color = COL_TOPMED, linewidth = 0.4) +
  annotate("text", x = nrow(top_s), y = med_s,
           label = sprintf(" median = %s genes", comma(med_s)),
           hjust = 1, vjust = -0.4, color = COL_TOPMED, size = 2.8,
           family = PAPER_FONT) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x = "Subject rank (descending count)",
       y = "Underlier burden per subject") +
  theme_paper() +
  theme_paper_inner_title() +
  labs(title = "Subjects")

.paper_plots <- list(b = p_gene, c = p_subj)
combined <- p_gene | p_subj
out <- file.path(OUT_DIR, "fig1_underlier_tallies_topmed_aou")
paper_ggsave(combined, out, width = NATURE_DOUBLE_IN, height = 3.5)
if (!is_paper_assemble()) cat("Wrote", paste0(out, ".png"), "\n")
