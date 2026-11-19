#!/usr/bin/env Rscript
# Among Figure 3a genes at 2005 PCs with PTV rate < 0.05%, are genes with
# underlier rate > 0.5% more likely to have a rare nearby SV (AF <= 1%,
# 100 kb) than genes with underlier rate < 0.5%?
suppressPackageStartupMessages({
  library(tidyverse)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
Sys.setenv(PC_SETTING = "2005")
source(file.path(ROOT, "scripts/helpers/_paths.R"))

OUT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
FIG3A_TSV <- file.path(OUT, "supplemental_data", "Source_Data_Fig3a_ptv_vs_underlier_gene.tsv")
SV_RDS <- file.path(DATA_EXPR, "subset_genotypes_100000_0.01_select.rds")

UNDER_CUTS <- c(0.2, 0.3, 0.5)  # percent of subjects
PTV_CUT    <- 0.05              # percent of subjects

fig3a <- read_tsv(FIG3A_TSV, show_col_types = FALSE)
kept <- read_rds(KEPT_SUBJECTS)
sv <- read_rds(SV_RDS) |>
  ungroup() |>
  transmute(gene_nv = str_remove(gene_id, "\\..*$"), sample) |>
  filter(sample %in% kept) |>
  distinct(gene_nv, sample) |>
  count(gene_nv, name = "n_sv_carriers")

d <- fig3a |>
  left_join(sv, by = "gene_nv") |>
  mutate(
    n_sv_carriers = replace_na(n_sv_carriers, 0L),
    has_sv = n_sv_carriers > 0,
    low_ptv = pct_ptv < PTV_CUT
  )

compare_cut <- function(under_cut) {
  n_eq <- sum(d$low_ptv & abs(d$pct_under - under_cut) < 1e-12)
  hi <- d |> filter(low_ptv, pct_under > under_cut)
  lo <- d |> filter(low_ptv, pct_under < under_cut)
  m <- rbind(
    c(sum(hi$has_sv), sum(!hi$has_sv)),
    c(sum(lo$has_sv), sum(!lo$has_sv))
  )
  dimnames(m) <- list(
    group = c(sprintf("under > %s%%, PTV < 0.05%%", under_cut),
              sprintf("under < %s%%, PTV < 0.05%%", under_cut)),
    sv = c("has_rare_SV", "no_rare_SV")
  )
  ft <- fisher.test(m)
  wt <- wilcox.test(hi$n_sv_carriers, lo$n_sv_carriers, alternative = "greater")
  genes <- bind_rows(
    hi |> mutate(group = sprintf("underlier > %s%% and PTV < 0.05%%", under_cut)),
    lo |> mutate(group = sprintf("underlier < %s%% and PTV < 0.05%%", under_cut))
  ) |>
    mutate(under_cut = under_cut)
  summ <- tibble(
    under_cut = under_cut,
    group = c(sprintf("underlier > %s%% and PTV < 0.05%%", under_cut),
              sprintf("underlier < %s%% and PTV < 0.05%%", under_cut)),
    n_genes = c(nrow(hi), nrow(lo)),
    n_with_rare_sv = c(sum(hi$has_sv), sum(lo$has_sv)),
    n_without_rare_sv = c(sum(!hi$has_sv), sum(!lo$has_sv)),
    frac_with_rare_sv = n_with_rare_sv / n_genes,
    mean_sv_carriers = c(mean(hi$n_sv_carriers), mean(lo$n_sv_carriers)),
    median_sv_carriers = c(median(hi$n_sv_carriers), median(lo$n_sv_carriers))
  )
  test <- tibble(
    under_cut = under_cut,
    n_fig3a_genes = nrow(fig3a),
    n_low_ptv = sum(d$low_ptv),
    n_exactly_cut_excluded = n_eq,
    n_high = nrow(hi),
    n_low = nrow(lo),
    high_n_with_sv = sum(hi$has_sv),
    low_n_with_sv = sum(lo$has_sv),
    high_frac = summ$frac_with_rare_sv[1],
    low_frac = summ$frac_with_rare_sv[2],
    fisher_or_has_sv = unname(ft$estimate),
    fisher_ci_lo = ft$conf.int[1],
    fisher_ci_hi = ft$conf.int[2],
    fisher_p = ft$p.value,
    high_mean_sv = mean(hi$n_sv_carriers),
    low_mean_sv = mean(lo$n_sv_carriers),
    high_median_sv = median(hi$n_sv_carriers),
    low_median_sv = median(lo$n_sv_carriers),
    wilcox_greater_p = wt$p.value,
    high_genes = paste(sort(hi$symbol), collapse = ",")
  )
  cat(sprintf("\n======== underlier cutoff %s%% (PTV < 0.05%%) ========\n", under_cut))
  print(as.data.frame(summ), row.names = FALSE)
  cat(sprintf("Exactly at cutoff (excluded): %d\n", n_eq))
  print(m)
  cat(sprintf("Fisher OR (has rare SV, high vs low): %.3f  95%% CI %.3f–%.3f  p = %.3g\n",
              unname(ft$estimate), ft$conf.int[1], ft$conf.int[2], ft$p.value))
  cat(sprintf("Wilcoxon n_sv_carriers high > low: p = %.3g\n", wt$p.value))
  cat(sprintf("High-underlier genes (%d): %s\n", nrow(hi),
              paste(sort(hi$symbol), collapse = ", ")))
  list(genes = genes, summ = summ, test = test)
}

res <- lapply(UNDER_CUTS, compare_cut)
write_tsv(bind_rows(lapply(res, `[[`, "genes")) |>
            arrange(under_cut, desc(pct_under), gene_nv) |>
            dplyr::select(under_cut, group, gene_nv, symbol, LOEUF_tier,
                          n_ptv_carriers, n_underlier_subjects,
                          pct_ptv, pct_under, n_sv_carriers, has_sv),
          file.path(OUT, "fig3a_low_ptv_high_vs_low_underlier_sv_genes.tsv"))
write_tsv(bind_rows(lapply(res, `[[`, "summ")),
          file.path(OUT, "fig3a_low_ptv_high_vs_low_underlier_sv_summary.tsv"))
write_tsv(bind_rows(lapply(res, `[[`, "test")),
          file.path(OUT, "fig3a_low_ptv_high_vs_low_underlier_sv_test.tsv"))
