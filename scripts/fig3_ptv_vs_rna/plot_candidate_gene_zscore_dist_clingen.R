#!/usr/bin/env Rscript
# Z-score distribution diagnostics for axis-hugging genes in the ClinGen HI
# z<-4 panel (ptv_vs_underlier_gene_level_2x2_clingen_hi, z<-4 facet).
#
# Restricted to ClinGen haploinsufficiency "Sufficient Evidence" disease genes.
# Two corners are examined:
#   * low-x / high-y : more RNA underliers, few/no PTV carriers
#   * high-x / low-y : more PTV carriers, ~no RNA underliers
#
# For each gene we plot the full per-subject expression z-score distribution
# (all kept subjects), mark the z<-4 threshold, the panel's z_-4 underliers
# (red top rug), and rare-PTV carriers (blue bottom rug).
#
# Inputs (in ../ = pre_Watershed_scripts):
#   expr_z_join.rds, kept_subjects.rds, z_score_underlier_z_-4_clean.rds
#   ../../ptv_all.tsv
# Outputs (this folder):
#   candidate_gene_zscore_dist_clingen_z_-4.png / .pdf
#   candidate_gene_zscore_summary_clingen_z_-4.tsv

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
OUT_DIR <- FIG3C

PW <- ROOT
SH  <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang"

Z_THRESH <- -4

# ClinGen HI candidate genes from the z<-4 panel extremes.
candidates <- tribble(
  ~gene_nv,           ~symbol,   ~corner,
  # low-x / high-y : underlier-rich, few/no PTV carriers
  "ENSG00000130821",  "SLC6A8",  "High underlier / low PTV",
  "ENSG00000129757",  "CDKN1C",  "High underlier / low PTV",
  "ENSG00000138326",  "RPS24",   "High underlier / low PTV",
  "ENSG00000055118",  "KCNH2",   "High underlier / low PTV",
  # high-x / low-y : PTV-rich, ~no underliers
  "ENSG00000215301",  "DDX3X",   "High PTV / low underlier",
  "ENSG00000139687",  "RB1",     "High PTV / low underlier",
  "ENSG00000134571",  "MYBPC3",  "High PTV / low underlier",
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

cat("Loading z_-4 underlier definition set ...\n")
underlier_def <- read_rds(file.path(DATA_UNDERLIER, "z_score_underlier_z_-4_clean.rds")) |>
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
  mutate(
    is_ptv           = replace_na(is_ptv, FALSE),
    is_underlier_def = replace_na(is_underlier_def, FALSE)
  ) |>
  left_join(candidates, by = "gene_nv") |>
  mutate(
    is_z4 = expression_zscore < Z_THRESH,
    facet = sprintf("%s\n(%s)", symbol, corner)
  )

ord <- candidates |>
  mutate(facet = sprintf("%s\n(%s)", symbol, corner)) |>
  arrange(factor(corner, levels = c("High underlier / low PTV", "High PTV / low underlier")), symbol)
dat <- dat |> mutate(facet = factor(facet, levels = ord$facet))

summ <- dat |>
  group_by(symbol, corner, gene_nv) |>
  summarise(
    n_subj        = n(),
    n_ptv         = sum(is_ptv),
    n_underlier   = sum(is_underlier_def),
    n_z4_raw      = sum(is_z4),
    n_ptv_underlier = sum(is_ptv & is_underlier_def),
    z_mean        = mean(expression_zscore),
    z_sd          = sd(expression_zscore),
    z_min         = min(expression_zscore),
    z_skew        = mean(((expression_zscore - mean(expression_zscore)) / sd(expression_zscore))^3),
    z_kurt        = mean(((expression_zscore - mean(expression_zscore)) / sd(expression_zscore))^4),
    frac_below_m2 = mean(expression_zscore < -2),
    logcpm_mean   = mean(expression_logcpm),
    .groups = "drop"
  ) |>
  arrange(factor(corner, levels = c("High underlier / low PTV", "High PTV / low underlier")), symbol)

write_tsv(summ, "candidate_gene_zscore_summary_clingen_z_-4.tsv")
cat("\n=== Per-gene z-score summary (ClinGen HI) ===\n")
print(summ |> mutate(across(where(is.numeric), ~round(.x, 3))), width = Inf)

lab_df <- summ |>
  mutate(
    facet = factor(sprintf("%s\n(%s)", symbol, corner), levels = ord$facet),
    lbl = sprintf("PTV carriers: %d\nunderliers (def): %d  (raw z<%g: %d)\nskew: %.2f  kurt: %.1f",
                  n_ptv, n_underlier, Z_THRESH, n_z4_raw, z_skew, z_kurt)
  )

x_lo <- min(dat$expression_zscore); x_hi <- max(dat$expression_zscore)

p <- ggplot(dat, aes(x = expression_zscore)) +
  geom_histogram(bins = 90, fill = "grey70", color = NA) +
  geom_vline(xintercept = Z_THRESH, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  geom_rug(data = dplyr::filter(dat, is_ptv),
           aes(x = expression_zscore), sides = "b",
           color = "#2166AC", alpha = 0.6, length = unit(0.06, "npc")) +
  geom_rug(data = dplyr::filter(dat, is_underlier_def),
           aes(x = expression_zscore), sides = "t",
           color = "#B2182B", alpha = 0.7, length = unit(0.06, "npc")) +
  geom_text(data = lab_df, aes(x = x_hi, y = Inf, label = lbl),
            hjust = 1, vjust = 1.1, size = 2.6, color = "grey20", lineheight = 0.95) +
  facet_wrap(~facet, scales = "free_y", ncol = 4) +
  scale_x_continuous(limits = c(x_lo, x_hi), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Per-subject expression z-score distributions for axis-hugging ClinGen HI genes (z<-4 panel)",
    subtitle = str_wrap(sprintf(
      "Each histogram = all %s kept subjects for one ClinGen HI gene. Red dashed line = z<-4 threshold; red rug (top) = z_-4 underliers (strict z<-4 + log2CPM). Blue rug (bottom) = subjects carrying a rare PTV (AF<0.5%%) in that gene. Top row: high-underlier / low-PTV corner; bottom row: high-PTV / low-underlier corner.",
      comma(n_kept)), width = 130),
    caption = "Bimodality / discrete side-clusters / heavy left shoulders suggest CNV / imprinting / mapping effects driving expression outliers that PTV calls do not capture.",
    x = "Expression z-score", y = "Number of subjects"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "grey25", lineheight = 1.15),
    plot.caption  = element_text(size = 7.5, color = "grey35", hjust = 0),
    strip.text    = element_text(face = "bold", size = 8.3, lineheight = 0.95),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_DIR, "candidate_gene_zscore_dist_clingen_z_-4.png"), p, width = 14, height = 7.5, dpi = 220)
ggsave(file.path(OUT_DIR, "candidate_gene_zscore_dist_clingen_z_-4.pdf"), p, width = 14, height = 7.5, device = cairo_pdf)
cat("\nWrote:", normalizePath("candidate_gene_zscore_dist_clingen_z_-4.png"), "\n")
