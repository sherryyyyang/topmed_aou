#!/usr/bin/env Rscript
# Regenerate strict (logCPM < gene_mean - 1 AND z < cutoff) LOEUF figures
# for alternative z-score cutoffs.
#
# This workspace currently contains the QC-filtered strict z<-4 artifacts
# (z_score_underlier_clean.rds, loeuf_calibration.png, etc). For this request
# we derive new QC-filtered underlier sets for z<-2 and z<-3 and generate:
#   - underlier_tallies_z_-2/z_-3
#   - underlier_tallies_loeuf_z_-2/z_-3
#   - loeuf_calibration_z_-2/z_-3 (Wilson CI, decile trend)
#   - underlier_loeuf_per_gene_compare_z_-2/z_-3 (enriched callouts; dumbbell)
#
# Output names include the requested "z_-2" label token.

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


wilson_ci <- function(k, n, conf = 0.95) {
  # Wilson 95% CI on a binomial proportion.
  z <- qnorm(1 - (1 - conf) / 2)
  p <- ifelse(n > 0, k / n, NA_real_)
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
  tibble(p = p, lo = pmax(centre - half, 0), hi = pmin(centre + half, 1))
}

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_COLORS <- c(
  "Constrained" = "#B2182B",
  "Mid-low"     = "#EF8A62",
  "Mid-high"    = "#67A9CF",
  "Tolerant"    = "#2166AC",
  "No LOEUF"    = "grey60"
)

resolve_gene_df <- function(loeuf_lookup, all_eligible_genes) {
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
    select(gene_id, LOEUF, LOEUF_tier)

  gene_df
}

derive_strict_underliers <- function(z_cutoff, kept_subjects, expr_z_join) {
  # strict: expression_zscore < z_cutoff AND expression_logcpm < gene_mean(logcpm) - 1
  # gene_mean is computed across all subjects present in expr_z_join.

  gene_means <- expr_z_join %>%
    group_by(gene_id) %>%
    summarise(mean_logcpm = mean(expression_logcpm, na.rm = TRUE), .groups = "drop")

  expr_z_join$mean_logcpm <- gene_means$mean_logcpm[match(expr_z_join$gene_id, gene_means$gene_id)]

  under <- expr_z_join %>%
    filter(expression_zscore < z_cutoff, expression_logcpm < (mean_logcpm - 1)) %>%
    select(gene_id, subjectid, expression_zscore, expression_logcpm)

  # Apply the global "extreme subject" filter after thresholding.
  under_clean <- under %>% filter(subjectid %in% kept_subjects)
  under_clean
}

plot_underlier_tallies <- function(z_tally_gene, z_tally_subj, z_label) {
  # z_tally_gene: gene_id, n, pct_underlier, rank
  # z_tally_subj: subjectid, n, rank
  topk <- z_tally_gene %>% arrange(desc(n)) %>% slice_head(n = 20)

  p_gene <- ggplot(topk, aes(x = reorder(gene_id, n), y = n)) +
    geom_col(fill = "#2166AC") +
    coord_flip() +
    labs(
      title = "Underlier genes (top 20 by count)",
      subtitle = paste0("Strict: logCPM < mean-1 and z < ", sub(".*-","", z_label))
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title.y = element_blank()
    )

  p_subj <- ggplot(z_tally_subj, aes(x = n)) +
    geom_histogram(binwidth = 1, fill = "grey75", color = "white", linewidth = 0.2) +
    labs(
      title = "Underlier genes per subject",
      subtitle = paste0("Strict z < ", sub(".*-","", z_label)),
      x = "Underlier gene count (per subject)",
      y = "Subjects"
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))

  (p_gene | p_subj) + plot_layout(widths = c(1.05, 1))
}

plot_underlier_tallies_loeuf <- function(under_clean_pairs, gene_df, z_label) {
  # under_clean_pairs: gene_id, subjectid
  subj_tier <- under_clean_pairs %>%
    left_join(gene_df, by = "gene_id") %>%
    count(subjectid, LOEUF_tier, name = "n_under")

  subj_tier_full <- tidyr::crossing(
    subjectid = unique(subj_tier$subjectid),
    LOEUF_tier = factor(TIER_LEVELS, levels = TIER_LEVELS)
  ) %>%
    left_join(subj_tier, by = c("subjectid", "LOEUF_tier")) %>%
    mutate(n_under = replace_na(n_under, 0L))

  p <- ggplot(subj_tier_full, aes(x = LOEUF_tier, y = n_under, fill = LOEUF_tier)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, width = 0.6) +
    geom_point(stat = "summary", fun = median, shape = 21, size = 2,
               color = "black", fill = "white", position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = TIER_COLORS, drop = FALSE) +
    labs(
      title = "Underlier genes per subject by LOEUF tier",
      subtitle = paste0("Strict z < ", sub(".*-","", z_label), " (logCPM thresholded)"),
      x = NULL,
      y = "Underlier gene count"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )

  p
}

plot_loeuf_calibration <- function(under_clean, kept_subjects, expr_cpm_adj, loeuf_lookup, z_label) {
  all_eligible_genes <- rownames(expr_cpm_adj)
  n_subjects <- length(kept_subjects)

  gene_df <- resolve_gene_df(loeuf_lookup, all_eligible_genes)

  under_per_gene <- under_clean %>% group_by(gene_id) %>% tally(name = "n")
  gene_df <- gene_df %>%
    left_join(under_per_gene, by = "gene_id") %>%
    mutate(n = replace_na(n, 0L))

  # Tier baseline prevalence + Wilson CI
  tier_tbl <- gene_df %>%
    group_by(LOEUF_tier) %>%
    summarise(
      n_genes_in_tier = n(),
      n_underlier_pairs = sum(n),
      n_total_pairs = n() * n_subjects,
      .groups = "drop"
    ) %>%
    arrange(LOEUF_tier)

  ci_tier <- wilson_ci(tier_tbl$n_underlier_pairs, tier_tbl$n_total_pairs)
  tier_tbl <- bind_cols(tier_tbl, ci_tier) %>%
    mutate(
      prev_per_10k = 10000 * p,
      prev_per_10k_lo = 10000 * lo,
      prev_per_10k_hi = 10000 * hi
    )

  tier_for_plot <- tier_tbl %>%
    mutate(
      LOEUF_tier = factor(LOEUF_tier, levels = rev(TIER_LEVELS)),
      label = sprintf("%.2f (%.2f-%.2f)", prev_per_10k, prev_per_10k_lo, prev_per_10k_hi)
    )

  p_tier <- ggplot(tier_for_plot,
                   aes(y = LOEUF_tier, x = prev_per_10k, color = LOEUF_tier)) +
    geom_errorbar(aes(xmin = prev_per_10k_lo, xmax = prev_per_10k_hi),
                  width = 0.18, linewidth = 0.5, orientation = "y") +
    geom_point(size = 3) +
    geom_text(aes(label = label), nudge_y = 0.32, size = 2.9, color = "grey25", show.legend = FALSE) +
    scale_color_manual(values = TIER_COLORS, drop = FALSE, name = "LOEUF tier") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.10))) +
    labs(
      title = "A. Per-LOEUF-tier underlier prevalence",
      subtitle = paste0("Observed RNA underliers per 10,000 (gene,subject) pairs; Wilson 95% CI (binomial approximation)"),
      x = "Prevalence per 10,000 (gene, subject) pairs",
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      axis.text.y = element_text(face = "bold")
    )

  # Deciles among LOEUF genes only (D1 = most constrained)
  gene_with_loeuf <- gene_df %>% filter(!is.na(LOEUF)) %>% arrange(LOEUF)
  decile_breaks <- quantile(gene_with_loeuf$LOEUF, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  gene_with_loeuf <- gene_with_loeuf %>%
    mutate(LOEUF_decile = cut(LOEUF, breaks = decile_breaks, include.lowest = TRUE,
                              labels = paste0("D", 1:10)))

  decile_tbl <- gene_with_loeuf %>%
    group_by(LOEUF_decile) %>%
    summarise(
      n_genes_in_decile = n(),
      loeuf_min = min(LOEUF),
      loeuf_max = max(LOEUF),
      loeuf_median = median(LOEUF),
      n_underlier_pairs = sum(n),
      n_total_pairs = n() * n_subjects,
      .groups = "drop"
    ) %>%
    arrange(LOEUF_decile)

  ci_decile <- wilson_ci(decile_tbl$n_underlier_pairs, decile_tbl$n_total_pairs)
  decile_tbl <- bind_cols(decile_tbl, ci_decile) %>%
    mutate(
      prev_per_10k = 10000 * p,
      prev_per_10k_lo = 10000 * lo,
      prev_per_10k_hi = 10000 * hi
    )

  p_decile <- ggplot(decile_tbl,
                     aes(x = LOEUF_decile, y = prev_per_10k, group = 1)) +
    geom_ribbon(aes(ymin = prev_per_10k_lo, ymax = prev_per_10k_hi),
                fill = "grey75", alpha = 0.45) +
    geom_line(color = "#2166AC", linewidth = 0.7) +
    geom_point(color = "#2166AC", size = 2.4) +
    geom_text(aes(label = sprintf("LOEUF\nmed=%.2f", loeuf_median)),
              vjust = 1.7, size = 2.4, color = "grey30") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.10))) +
    labs(
      title = "B. Per-LOEUF-decile underlier prevalence",
      subtitle = "D1 = most constrained ... D10 = most tolerant (approx. equal genes/decile). Grey ribbon = Wilson 95% CI.",
      x = "LOEUF decile (gene-level)",
      y = "Prevalence per 10,000 (gene, subject) pairs"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  combined <- (p_tier | p_decile) +
    plot_layout(widths = c(1.18, 1)) +
    plot_annotation(
      title = paste0("Calibration of RNA-underlier prevalence vs LOEUF  [strict ", z_label, "]"),
      subtitle = paste0("Underlier: log2CPM_adj < gene_mean(logCPM) - 1 AND z < ", z_label),
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.subtitle = element_text(size = 10, color = "grey25"))
    )

  list(combined = combined, gene_df = gene_df, tier_tbl = tier_tbl, decile_tbl = decile_tbl)
}

plot_loeuf_per_gene_compare <- function(under_clean, kept_subjects, expr_cpm_adj, loeuf_lookup, z_label) {
  all_eligible_genes <- rownames(expr_cpm_adj)
  n_subjects <- length(kept_subjects)

  gene_df <- resolve_gene_df(loeuf_lookup, all_eligible_genes)

  under_per_gene <- under_clean %>% group_by(gene_id) %>% tally(name = "n")
  gene_df <- gene_df %>%
    left_join(under_per_gene, by = "gene_id") %>%
    mutate(n = replace_na(n, 0L))

  # Tier baseline prevalence
  tier_tbl <- gene_df %>%
    group_by(LOEUF_tier) %>%
    summarise(n_genes_in_tier = n(),
              n_underlier_pairs = sum(n),
              n_total_pairs = n() * n_subjects,
              .groups = "drop") %>%
    mutate(tier_prev_per_10k = 10000 * (n_underlier_pairs / n_total_pairs)) %>%
    select(LOEUF_tier, tier_prev_per_10k)

  tier_prev_map <- tier_tbl %>% {setNames(.$tier_prev_per_10k, .$LOEUF_tier)}

  # Decile baseline prevalence
  gene_with_loeuf <- gene_df %>% filter(!is.na(LOEUF)) %>% arrange(LOEUF)
  decile_breaks <- quantile(gene_with_loeuf$LOEUF, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  gene_with_loeuf <- gene_with_loeuf %>%
    mutate(LOEUF_decile = cut(LOEUF, breaks = decile_breaks, include.lowest = TRUE,
                              labels = paste0("D", 1:10)))

  decile_tbl <- gene_with_loeuf %>%
    group_by(LOEUF_decile) %>%
    summarise(n_genes_in_decile = n(),
              loeuf_median = median(LOEUF),
              n_underlier_pairs = sum(n),
              n_total_pairs = n_genes_in_decile * n_subjects,
              decile_prev_per_10k = 10000 * (n_underlier_pairs / n_total_pairs),
              .groups = "drop")

  # Merge expected into gene_df
  gene_with_decile <- gene_df %>%
    left_join(gene_with_loeuf %>% select(gene_id, LOEUF_decile), by = "gene_id") %>%
    left_join(decile_tbl %>% select(LOEUF_decile, decile_prev_per_10k, loeuf_median), by = "LOEUF_decile") %>%
    left_join(tier_tbl, by = "LOEUF_tier") %>%
    mutate(
      observed_pct = n / n_subjects,
      expected_pct_dec = decile_prev_per_10k / 10000,
      expected_pct_tier = tier_prev_per_10k / 10000,
      fold_decile = ifelse(expected_pct_dec > 0, observed_pct / expected_pct_dec, NA_real_),
      fold_tier = ifelse(expected_pct_tier > 0, observed_pct / expected_pct_tier, NA_real_),
      p_binom = ifelse(is.na(expected_pct_dec), NA_real_, pbinom(n - 1, n_subjects, expected_pct_dec, lower.tail = FALSE))
    )

  # callouts: enriched in decile baseline, significant by p_binom
  callouts <- gene_with_decile %>%
    filter(!is.na(p_binom), fold_decile > 1, p_binom < 0.01) %>%
    arrange(p_binom) %>%
    slice_head(n = 40)

  # Add symbol column from LOEUF lookup (gene_id may have version)
  # Build a best-effort symbol map.
  symbol_map <- loeuf_lookup %>% distinct(gene_id, symbol) %>% filter(!is.na(symbol))
  callouts <- callouts %>%
    left_join(symbol_map, by = "gene_id") %>%
    mutate(symbol = ifelse(is.na(symbol), gene_id, symbol))

  # Plot dumbbell on log x-axis
  callouts <- callouts %>%
    arrange(observed_pct) %>%
    mutate(symbol_f = factor(symbol, levels = unique(symbol)))

  p <- ggplot(callouts, aes(y = symbol_f)) +
    geom_segment(aes(x = expected_pct_dec, xend = observed_pct, y = symbol_f, yend = symbol_f),
                 color = "transparent") +
    geom_segment(aes(x = expected_pct_dec, xend = observed_pct, y = symbol_f, yend = symbol_f),
                 linewidth = 0.35, color = "grey65") +
    geom_point(aes(x = expected_pct_dec), shape = 21, size = 2.3,
               fill = "white", color = "#666666") +
    geom_point(aes(x = observed_pct), shape = 21, size = 2.9,
               fill = "#2166AC", color = "#2166AC") +
    scale_x_continuous(
      trans = log10_trans(),
      labels = label_number(accuracy = 1e-4)
    ) +
    labs(
      title = "Observed vs expected RNA-underlier prevalence per gene",
      subtitle = paste0("Strict ", z_label, " callouts (expected baseline = LOEUF decile prevalence; p_binom = binomial upper-tail)"),
      x = "Prevalence (underlier subjects / total subjects) [log scale]",
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 7)
    )

  # Output TSV with the same core columns as existing strict callouts.
  out_df <- callouts %>%
    transmute(
      gene_id,
      gene_id_no_version = str_remove(gene_id, "\\.\\d+$"),
      symbol,
      LOEUF,
      LOEUF_tier,
      n,
      LOEUF_decile,
      tier_prev_per_10k,
      decile_prev_per_10k = decile_prev_per_10k,
      loeuf_median,
      observed_pct,
      expected_pct_dec,
      expected_pct_tier,
      fold_decile,
      fold_tier,
      p_binom
    )

  list(plot = p, callouts = out_df)
}

main <- function() {
  kept_subjects <- read_rds(KEPT_SUBJECTS)
  N <- length(kept_subjects)

  expr_z_join <- read_rds(EXPR_Z_JOIN)
  # sanity
  stopifnot(all(c("gene_id", "subjectid", "expression_zscore", "expression_logcpm") %in% colnames(expr_z_join)))

  expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
  loeuf_lookup <- read_rds(LOEUF_LOOKUP)

  gene_df_all <- resolve_gene_df(loeuf_lookup, rownames(expr_cpm_adj))

  for (z_cutoff in c(-2, -3)) {
    z_label <- paste0("z_-", abs(z_cutoff))
    message("=== Deriving strict underliers for ", z_label, " ===")

    under_clean <- derive_strict_underliers(z_cutoff, kept_subjects, expr_z_join)

    # Persist underlier pair table.
    saveRDS(under_clean, file.path(DATA_UNDERLIER, paste0("z_score_underlier_", z_label, "_clean.rds")))

    # Gene tally (only genes with >=1 underlier).
    z_tally_gene <- under_clean %>%
      group_by(gene_id) %>%
      tally(name = "n") %>%
      arrange(desc(n)) %>%
      mutate(
        pct_underlier = n / N,
        rank = row_number()
      )
    saveRDS(z_tally_gene, file.path(DATA_UNDERLIER, paste0("z_score_tally_", z_label, "_clean.rds")))

    # Subject tally.
    z_tally_subj <- under_clean %>%
      group_by(subjectid) %>%
      tally(name = "n") %>%
      arrange(desc(n)) %>%
      mutate(rank = row_number())
    saveRDS(z_tally_subj, file.path(DATA_UNDERLIER, paste0("z_score_tally_subj_", z_label, "_clean.rds")))

    # Convenience pairs for LOEUF tier plots.
    under_pairs_simple <- under_clean %>% select(gene_id, subjectid)

    # --------- Figure: underlier tallies ----------
    p_underlier <- plot_underlier_tallies(z_tally_gene, z_tally_subj, z_label)
    ggsave(paste0("underlier_tallies_", z_label, ".png"), p_underlier,
           width = 12, height = 5.0, dpi = 200)
    ggsave(paste0("underlier_tallies_", z_label, ".pdf"), p_underlier,
           width = 12, height = 5.0, device = cairo_pdf())

    p_underlier_loeuf <- plot_underlier_tallies_loeuf(under_pairs_simple, gene_df_all, z_label)
    ggsave(paste0("underlier_tallies_loeuf_", z_label, ".png"), p_underlier_loeuf,
           width = 8.5, height = 5.2, dpi = 200)
    ggsave(paste0("underlier_tallies_loeuf_", z_label, ".pdf"), p_underlier_loeuf,
           width = 8.5, height = 5.2, device = cairo_pdf())

    # --------- Figure: LOEUF calibration ----------
    calib <- plot_loeuf_calibration(
      under_clean = under_clean,
      kept_subjects = kept_subjects,
      expr_cpm_adj = expr_cpm_adj,
      loeuf_lookup = loeuf_lookup,
      z_label = z_label
    )

    ggsave(paste0("loeuf_calibration_", z_label, ".png"), calib$combined,
           width = 12, height = 5.0, dpi = 200)
    ggsave(paste0("loeuf_calibration_", z_label, ".pdf"), calib$combined,
           width = 12, height = 5.0, device = cairo_pdf())

    write_tsv(calib$tier_tbl, paste0("underlier_loeuf_calibration_by_tier_", z_label, ".tsv"))
    write_tsv(calib$decile_tbl, paste0("underlier_loeuf_calibration_by_decile_", z_label, ".tsv"))

    # --------- Figure: per-gene observed vs expected ----------
    pc <- plot_loeuf_per_gene_compare(
      under_clean = under_clean,
      kept_subjects = kept_subjects,
      expr_cpm_adj = expr_cpm_adj,
      loeuf_lookup = loeuf_lookup,
      z_label = z_label
    )

    ggsave(paste0("loeuf_per_gene_compare_", z_label, ".png"), pc$plot,
           width = 10.5, height = 10.5, dpi = 200)
    ggsave(paste0("loeuf_per_gene_compare_", z_label, ".pdf"), pc$plot,
           width = 10.5, height = 10.5, device = cairo_pdf())

    write_tsv(pc$callouts, paste0("underlier_loeuf_callouts_", z_label, ".tsv"))
  }
}

main()

