#!/usr/bin/env Rscript
# Same Figure 3a contrast as the 100 kb test, but rare SVs within 10 kb
# of the gene body (AF <= 1%). Built by restricting the published 100 kb
# gene-SV pairs (already AF-filtered protein-coding overlaps) to those
# whose SV is within 10 kb of the GENCODE v30 gene body.
suppressPackageStartupMessages({
  library(tidyverse)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
Sys.setenv(PC_SETTING = "2005")
source(file.path(ROOT, "scripts/helpers/_paths.R"))

OUT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
FIG3A_TSV <- file.path(OUT, "supplemental_data", "Source_Data_Fig3a_ptv_vs_underlier_gene.tsv")
SV100 <- file.path(DATA_EXPR, "subset_genotypes_100000_0.01_select.rds")
GFF <- file.path(SHEYANG, "annotation", "gencode.v30.annotation.gff3")
CACHE10 <- file.path(OUT, "subset_genotypes_10000_0.01_select.rds")
FLANK <- 10000L
PTV_CUT <- 0.05
UNDER_CUTS <- c(0.2, 0.3)

build_10kb <- function() {
  message("parsing GENCODE v30 gene bodies ...")
  genes <- read_tsv(
    pipe(sprintf("awk '$3==\"gene\"' %s", shQuote(GFF))),
    col_names = c("chr", "src", "type", "start", "end", "score", "strand", "phase", "attr"),
    col_types = "ccciicccc", comment = "#"
  ) |>
    mutate(
      gene_id = str_match(attr, "gene_id=([^;]+)")[, 2],
      gene_nv = str_remove(gene_id, "\\..*$"),
      chr = str_remove(chr, "^chr")
    ) |>
    dplyr::select(gene_id, gene_nv, chr, start, end)
  message("  ", format(nrow(genes), big.mark = ","), " gene rows")

  message("loading 100 kb SV pairs and filtering to 10 kb ...")
  # variant_id can be collapsed with ";", e.g. DEL_8:1-2;DUP_8:3-4
  sv <- read_rds(SV100) |>
    ungroup() |>
    separate_rows(variant_id, sep = ";") |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$"))
  parsed <- str_match(sv$variant_id, "^[^_]+_([^:]+):([0-9]+)-([0-9]+)")
  sv <- sv |>
    mutate(
      sv_chr = str_remove(parsed[, 2], "^chr"),
      sv_start = as.integer(parsed[, 3]),
      sv_end = as.integer(parsed[, 4])
    )
  n_bad <- sum(is.na(sv$sv_start))
  if (n_bad > 0) {
    message("  dropping ", format(n_bad, big.mark = ","),
            " unparsed variant_id rows; examples: ",
            paste(head(unique(sv$variant_id[is.na(sv$sv_start)]), 5), collapse = "; "))
    sv <- filter(sv, !is.na(sv_start), !is.na(sv_end))
  }

  sv10 <- sv |>
    left_join(genes |> dplyr::select(gene_id, g_chr = chr, g_start = start, g_end = end),
              by = "gene_id")
  n_unmatched <- sum(is.na(sv10$g_start))
  if (n_unmatched > 0) {
    message("  ", format(n_unmatched, big.mark = ","), " pairs unmatched to GFF gene_id; trying unversioned")
    genes_nv <- genes |>
      group_by(gene_nv) |>
      slice_min(start, n = 1, with_ties = FALSE) |>
      ungroup() |>
      dplyr::select(gene_nv, g_chr2 = chr, g_start2 = start, g_end2 = end)
    sv10 <- sv10 |>
      left_join(genes_nv, by = "gene_nv") |>
      mutate(
        g_chr = if_else(is.na(g_start), g_chr2, g_chr),
        g_start = if_else(is.na(g_start), g_start2, g_start),
        g_end = if_else(is.na(g_end), g_end2, g_end)
      )
  }
  sv10 <- sv10 |>
    filter(!is.na(g_start), sv_chr == g_chr) |>
    mutate(
      dist = case_when(
        sv_end < g_start ~ g_start - sv_end,
        sv_start > g_end ~ sv_start - g_end,
        TRUE ~ 0L
      )
    ) |>
    filter(dist <= FLANK) |>
    distinct(sample, gene_id, variant_id)

  message(sprintf("  100 kb pairs: %s; 10 kb pairs: %s (%.1f%%)",
                  format(nrow(sv), big.mark = ","),
                  format(nrow(sv10), big.mark = ","),
                  100 * nrow(sv10) / nrow(sv)))
  write_rds(sv10, CACHE10)
  sv10
}

sv10 <- if (file.exists(CACHE10)) read_rds(CACHE10) else build_10kb()

fig3a <- read_tsv(FIG3A_TSV, show_col_types = FALSE)
kept <- read_rds(KEPT_SUBJECTS)
sv_n <- sv10 |>
  filter(sample %in% kept) |>
  mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
  distinct(gene_nv, sample) |>
  count(gene_nv, name = "n_sv_carriers")

d <- fig3a |>
  left_join(sv_n, by = "gene_nv") |>
  mutate(
    n_sv_carriers = replace_na(n_sv_carriers, 0L),
    has_sv = n_sv_carriers > 0,
    low_ptv = pct_ptv < PTV_CUT
  )

compare_cut <- function(under_cut) {
  hi <- d |> filter(low_ptv, pct_under > under_cut)
  lo <- d |> filter(low_ptv, pct_under < under_cut)
  m <- rbind(
    c(sum(hi$has_sv), sum(!hi$has_sv)),
    c(sum(lo$has_sv), sum(!lo$has_sv))
  )
  dimnames(m) <- list(
    group = c(sprintf("under > %s%%, PTV < 0.05%%", under_cut),
              sprintf("under < %s%%, PTV < 0.05%%", under_cut)),
    sv = c("has_rare_SV_10kb", "no_rare_SV_10kb")
  )
  ft <- fisher.test(m)
  wt <- wilcox.test(hi$n_sv_carriers, lo$n_sv_carriers, alternative = "greater")
  summ <- tibble(
    window = "10kb",
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
    window = "10kb",
    under_cut = under_cut,
    n_high = nrow(hi), n_low = nrow(lo),
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
  cat(sprintf("\n======== 10 kb  |  underlier %s%%  |  PTV < 0.05%% ========\n", under_cut))
  print(as.data.frame(summ), row.names = FALSE)
  print(m)
  cat(sprintf("Fisher OR (has rare SV 10kb, high vs low): %.3f  95%% CI %.3f–%.3f  p = %.3g\n",
              unname(ft$estimate), ft$conf.int[1], ft$conf.int[2], ft$p.value))
  cat(sprintf("Wilcoxon n_sv_carriers high > low: p = %.3g\n", wt$p.value))
  list(summ = summ, test = test,
       genes = bind_rows(
         hi |> mutate(group = summ$group[1], under_cut = under_cut),
         lo |> mutate(group = summ$group[2], under_cut = under_cut)
       ))
}

res <- lapply(UNDER_CUTS, compare_cut)
write_tsv(bind_rows(lapply(res, `[[`, "summ")),
          file.path(OUT, "fig3a_low_ptv_sv_10kb_summary.tsv"))
write_tsv(bind_rows(lapply(res, `[[`, "test")),
          file.path(OUT, "fig3a_low_ptv_sv_10kb_test.tsv"))
write_tsv(bind_rows(lapply(res, `[[`, "genes")) |>
            arrange(under_cut, desc(pct_under)) |>
            dplyr::select(under_cut, group, gene_nv, symbol,
                          n_ptv_carriers, n_underlier_subjects,
                          pct_ptv, pct_under, n_sv_carriers, has_sv),
          file.path(OUT, "fig3a_low_ptv_sv_10kb_genes.tsv"))
message("wrote 10 kb gene-SV cache: ", CACHE10)
