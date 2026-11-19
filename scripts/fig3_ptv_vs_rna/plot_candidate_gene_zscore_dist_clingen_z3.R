#!/usr/bin/env Rscript
# Z-score distribution diagnostics for axis-hugging genes in the ClinGen HI
# z<-3 panel (ptv_vs_underlier_gene_level_2x2_clingen_hi, z<-3 facet).
#
# Mirrors plot_candidate_gene_zscore_dist_clingen.R but for the z<-3 definition.
#
# Inputs (in ../): expr_z_join.rds, kept_subjects.rds,
#   z_score_underlier_z_-3_clean.rds, ../../ptv_all.tsv
# Outputs (this folder):
#   candidate_gene_zscore_dist_clingen_z_-3.png / .pdf
#   candidate_gene_zscore_summary_clingen_z_-3.tsv

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG3C

PW <- ROOT
SH  <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang"

Z_THRESH <- -3

candidates <- tribble(
  ~gene_nv,           ~symbol,   ~corner,
  "ENSG00000183691",  "NOG",     "High underlier / low PTV",
  "ENSG00000129757",  "CDKN1C",  "High underlier / low PTV",
  "ENSG00000077782",  "FGFR1",   "High underlier / low PTV",
  "ENSG00000130821",  "SLC6A8",  "High underlier / low PTV",
  "ENSG00000215301",  "DDX3X",   "High PTV / low underlier",
  "ENSG00000139687",  "RB1",     "High PTV / low underlier",
  "ENSG00000179950",  "PUF60",   "High PTV / low underlier",
  "ENSG00000171316",  "CHD7",    "High PTV / low underlier"
)

cat("Loading kept_subjects.rds ...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept_subjects)

cat("Loading PTV table ...\n")
ptv_carriers <- read_tsv(file.path(SH, "ptv_all.tsv"), show_col_types = FALSE) |>
  filter(sample %in% kept_subjects) |>
  mutate(gene_nv = str_remove(gene, "\\..*$")) |>
  filter(gene_nv %in% candidates$gene_nv) |>
  distinct(subjectid = sample, gene_nv) |>
  mutate(is_ptv = TRUE)

cat("Loading z_-3 underlier definition set ...\n")
underlier_def <- read_rds(file.path(DATA_UNDERLIER, "z_score_underlier_z_-3_clean.rds")) |>
  mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
  filter(gene_nv %in% candidates$gene_nv) |>
  distinct(subjectid, gene_nv) |>
  mutate(is_underlier_def = TRUE)

cat("Loading expr_z_join.rds (large) ...\n")
zj <- read_rds(EXPR_Z_JOIN)
zj <- zj |>
  mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
  filter(gene_nv %in% candidates$gene_nv, subjectid %in% kept_subjects)
gc(verbose = FALSE)

dat <- zj |>
  left_join(ptv_carriers, by = c("subjectid", "gene_nv")) |>
  left_join(underlier_def, by = c("subjectid", "gene_nv")) |>
  mutate(is_ptv = replace_na(is_ptv, FALSE),
         is_underlier_def = replace_na(is_underlier_def, FALSE)) |>
  left_join(candidates, by = "gene_nv") |>
  mutate(is_zc = expression_zscore < Z_THRESH,
         facet = symbol)

ord <- candidates |>
  mutate(facet = symbol) |>
  arrange(factor(corner, levels = c("High underlier / low PTV", "High PTV / low underlier")), symbol)
dat <- dat |> mutate(facet = factor(facet, levels = ord$facet))

summ <- dat |>
  group_by(symbol, corner, gene_nv) |>
  summarise(
    n_subj = n(), n_ptv = sum(is_ptv),
    n_underlier = sum(is_underlier_def), n_zc_raw = sum(is_zc),
    n_ptv_underlier = sum(is_ptv & is_underlier_def),
    z_mean = mean(expression_zscore), z_sd = sd(expression_zscore),
    z_min = min(expression_zscore),
    z_skew = mean(((expression_zscore - mean(expression_zscore)) / sd(expression_zscore))^3),
    z_kurt = mean(((expression_zscore - mean(expression_zscore)) / sd(expression_zscore))^4),
    frac_below_m2 = mean(expression_zscore < -2), logcpm_mean = mean(expression_logcpm),
    .groups = "drop") |>
  arrange(factor(corner, levels = c("High underlier / low PTV", "High PTV / low underlier")), symbol)

write_tsv(summ, file.path(OUT_DIR, "candidate_gene_zscore_summary_clingen_z_-3.tsv"))
cat("\n=== Per-gene z-score summary (ClinGen HI, z<-3) ===\n")
print(summ |> mutate(across(where(is.numeric), ~round(.x, 3))), width = Inf)

# Distinct from Fig 3b (grey / vermillion carrier legend): purple PTV ticks,
# blue underlier ticks + z = -3 line. Gene in the strip; counts on the line below.
COL_3C_PTV <- unname(OI["purple"])
COL_3C_UND <- unname(OI["blue"])
lab_df <- summ |>
  mutate(
    facet_lab = sprintf("%s\n%d PTV · %d underliers", symbol, n_ptv, n_underlier)
  )
facet_levels <- lab_df$facet_lab[match(ord$facet, lab_df$symbol)]
dat <- dat |>
  left_join(lab_df |> select(symbol, facet_lab), by = "symbol") |>
  mutate(facet = factor(facet_lab, levels = facet_levels))
x_lo <- min(dat$expression_zscore)
x_hi <- max(dat$expression_zscore)
rug_ptv <- dplyr::filter(dat, is_ptv) |> mutate(mark = "PTV")
rug_und <- dplyr::filter(dat, is_underlier_def) |> mutate(mark = "underlier")

p <- ggplot(dat, aes(x = expression_zscore)) +
  geom_histogram(bins = 90, fill = "grey85", color = NA) +
  geom_vline(xintercept = Z_THRESH, linetype = "dashed", color = COL_3C_UND, linewidth = 0.4) +
  geom_rug(
    data = rug_ptv, aes(x = expression_zscore, color = mark), sides = "b",
    alpha = 0.75, length = unit(0.06, "npc")
  ) +
  geom_rug(
    data = rug_und, aes(x = expression_zscore, color = mark), sides = "t",
    alpha = 0.80, length = unit(0.06, "npc")
  ) +
  facet_wrap(~facet, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = c(PTV = COL_3C_PTV, underlier = COL_3C_UND),
    breaks = c("PTV", "underlier"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(x_lo, x_hi), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Expression z-score", y = "Number of subjects") +
  theme_paper() +
  theme(
    strip.clip = "off",
    strip.text = element_text(
      face = "bold", size = 6.5, lineheight = 1.10,
      margin = margin(3, 2, 3, 2)
    ),
    legend.position = "bottom",
    legend.justification = "center",
    legend.key.width = unit(10, "pt"),
    legend.key.height = unit(8, "pt"),
    plot.margin = margin(1, 6, 2, 8)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.6)))

.paper_plots <- list(c = p)
paper_ggsave(p, file.path(OUT_DIR, "candidate_gene_zscore_dist_clingen_z_-3"),
             width = 7.2, height = 3.6)
if (!is_paper_assemble()) cat("\nWrote clingen z-score histograms\n")
