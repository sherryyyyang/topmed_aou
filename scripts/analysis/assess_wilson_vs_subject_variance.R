#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG2X


TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_COLORS <- c(
  "Constrained" = "#B2182B",
  "Mid-low"     = "#EF8A62",
  "Mid-high"    = "#67A9CF",
  "Tolerant"    = "#2166AC",
  "No LOEUF"    = "grey60"
)

wilson_ci <- function(k, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- ifelse(n > 0, k / n, NA_real_)
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
  tibble(p = p, lo = pmax(centre - half, 0), hi = pmin(centre + half, 1))
}

set.seed(20260505)

cat("Loading inputs...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_subjects <- length(kept_subjects)
expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
all_eligible_genes <- rownames(expr_cpm_adj)
z_score_underlier <- read_rds(file.path(DATA_UNDERLIER, "z_score_underlier_clean.rds"))  # strict z<-4 clean
loeuf_lookup <- read_rds(LOEUF_LOOKUP)

loeuf_full <- loeuf_lookup %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  select(gene_id, LOEUF_full = LOEUF, tier_full = LOEUF_tier)
loeuf_nv <- loeuf_lookup %>%
  distinct(gene_id_no_version, .keep_all = TRUE) %>%
  select(gene_id_no_version, LOEUF_nv = LOEUF, tier_nv = LOEUF_tier)

gene_df <- tibble(gene_id = all_eligible_genes) %>%
  mutate(gene_id_no_version = str_remove(gene_id, "\\.\\d+$")) %>%
  left_join(loeuf_full, by = "gene_id") %>%
  left_join(loeuf_nv, by = "gene_id_no_version") %>%
  mutate(
    LOEUF = coalesce(LOEUF_full, LOEUF_nv),
    LOEUF_tier = factor(
      ifelse(is.na(LOEUF), "No LOEUF", as.character(coalesce(tier_full, tier_nv))),
      levels = TIER_LEVELS
    )
  ) %>%
  select(gene_id, LOEUF_tier)

# Underlier pairs assigned tier
under_tiered <- z_score_underlier %>%
  ungroup() %>%
  select(gene_id, subjectid) %>%
  inner_join(gene_df, by = "gene_id")

# Tier gene counts for normalization
tier_gene_counts <- gene_df %>%
  count(LOEUF_tier, name = "n_genes_in_tier")

# Subject x tier underlier counts (fill zeros for absent underliers)
subj_tier_counts <- tidyr::crossing(
  subjectid = kept_subjects,
  LOEUF_tier = factor(TIER_LEVELS, levels = TIER_LEVELS)
) %>%
  left_join(
    under_tiered %>% count(subjectid, LOEUF_tier, name = "n_under"),
    by = c("subjectid", "LOEUF_tier")
  ) %>%
  mutate(n_under = replace_na(n_under, 0L)) %>%
  left_join(tier_gene_counts, by = "LOEUF_tier") %>%
  mutate(
    underlier_per_sample_across_genes = n_under / n_genes_in_tier
  )

# Tier-level pooled prevalence (same estimand as calibration plot)
tier_tbl <- under_tiered %>%
  ungroup() %>%
  count(LOEUF_tier, name = "n_underlier_pairs") %>%
  right_join(tier_gene_counts, by = "LOEUF_tier") %>%
  mutate(
    n_underlier_pairs = replace_na(n_underlier_pairs, 0L),
    n_total_pairs = n_genes_in_tier * n_subjects
  ) %>%
  arrange(LOEUF_tier)

ci_w <- wilson_ci(tier_tbl$n_underlier_pairs, tier_tbl$n_total_pairs)
tier_tbl <- bind_cols(tier_tbl, ci_w) %>%
  mutate(
    prevalence_per_10k = 10000 * p,
    wilson_lo_per_10k = 10000 * lo,
    wilson_hi_per_10k = 10000 * hi
  )

# Cluster bootstrap by subject to account for within-subject and within-gene dependence.
B <- 500
subject_ids <- kept_subjects
boot_res <- vector("list", B)
cat("Running subject-bootstrap (B=", B, ")...\n", sep = "")
for (b in seq_len(B)) {
  samp_subj <- sample(subject_ids, size = length(subject_ids), replace = TRUE)
  mult <- as.integer(table(factor(samp_subj, levels = subject_ids)))
  mult_tbl <- tibble(subjectid = subject_ids, m = mult)

  bt <- subj_tier_counts %>%
    left_join(mult_tbl, by = "subjectid") %>%
    group_by(LOEUF_tier) %>%
    summarise(
      n_under_bt = sum(n_under * m),
      n_total_bt = sum(n_genes_in_tier * m),
      .groups = "drop"
    ) %>%
    mutate(prevalence_bt_per_10k = 10000 * (n_under_bt / n_total_bt),
           b = b)
  boot_res[[b]] <- bt
}
boot_df <- bind_rows(boot_res)

boot_ci <- boot_df %>%
  group_by(LOEUF_tier) %>%
  summarise(
    boot_lo_per_10k = quantile(prevalence_bt_per_10k, 0.025, na.rm = TRUE),
    boot_hi_per_10k = quantile(prevalence_bt_per_10k, 0.975, na.rm = TRUE),
    boot_sd_per_10k = sd(prevalence_bt_per_10k, na.rm = TRUE),
    .groups = "drop"
  )

ci_compare <- tier_tbl %>%
  left_join(boot_ci, by = "LOEUF_tier") %>%
  mutate(
    wilson_width = wilson_hi_per_10k - wilson_lo_per_10k,
    boot_width = boot_hi_per_10k - boot_lo_per_10k,
    width_ratio_boot_over_wilson = boot_width / wilson_width
  )

write_tsv(ci_compare, file.path(OUT_DIR, "underlier_loeuf_wilson_vs_subject_bootstrap.tsv"))

# ---------- Plot panel 1: histograms ----------
max_prop <- max(subj_tier_counts$underlier_per_sample_across_genes, na.rm = TRUE)
binw <- max(0.00025, signif(max_prop / 40, 2))

p_hist <- ggplot(
  subj_tier_counts,
  aes(x = underlier_per_sample_across_genes, fill = LOEUF_tier)
) +
  geom_histogram(binwidth = binw, color = "white", linewidth = 0.15) +
  facet_wrap(~LOEUF_tier, scales = "free_y", ncol = 1, as.table = FALSE) +
  scale_fill_manual(values = TIER_COLORS, drop = FALSE) +
  scale_x_continuous(labels = label_percent(accuracy = 0.1)) +
  labs(
    title = "A. Underlier/sample across genes by LOEUF tier",
    subtitle = "Each facet: per-subject underlier burden in that tier (n_underlier_genes_in_tier / n_genes_in_tier).",
    x = "Underlier burden per subject across genes in tier",
    y = "Subjects"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey30")
  )

# ---------- Plot panel 2: CI comparison ----------
ci_for_plot <- ci_compare %>%
  mutate(LOEUF_tier = factor(LOEUF_tier, levels = rev(TIER_LEVELS)))

p_ci <- ggplot(ci_for_plot, aes(y = LOEUF_tier)) +
  geom_errorbar(
    aes(xmin = boot_lo_per_10k, xmax = boot_hi_per_10k),
    width = 0.20, linewidth = 1.1, color = "grey45", orientation = "y"
  ) +
  geom_errorbar(
    aes(xmin = wilson_lo_per_10k, xmax = wilson_hi_per_10k, color = LOEUF_tier),
    width = 0.10, linewidth = 0.8, orientation = "y"
  ) +
  geom_point(aes(x = prevalence_per_10k, color = LOEUF_tier), size = 2.5) +
  scale_color_manual(values = TIER_COLORS, drop = FALSE) +
  labs(
    title = "B. Wilson CI vs subject-bootstrap CI",
    subtitle = "Thin colored = Wilson 95% CI (pair-level binomial). Thick grey = subject-cluster bootstrap 95% CI.",
    x = "Prevalence per 10,000 (gene, subject) pairs",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey30")
  )

combined <- p_hist | p_ci
ggsave(file.path(OUT_DIR, "loeuf_calibration_variance_check.png"), combined, width = 13, height = 9, dpi = 220)
ggsave(file.path(OUT_DIR, "loeuf_calibration_variance_check.pdf"), combined, width = 13, height = 9, device = cairo_pdf)

cat("Wrote:\n")
cat(" - loeuf_calibration_variance_check.png\n")
cat(" - loeuf_calibration_variance_check.pdf\n")
cat(" - underlier_loeuf_wilson_vs_subject_bootstrap.tsv\n")

# Console summary
cat("\nCI width ratio (bootstrap / Wilson):\n")
print(ci_compare %>% select(LOEUF_tier, prevalence_per_10k, wilson_width, boot_width, width_ratio_boot_over_wilson))

