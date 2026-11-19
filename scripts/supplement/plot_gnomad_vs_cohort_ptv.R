#!/usr/bin/env Rscript
# gnomAD v4.1 pLoF prevalence vs the cohort PTV-carrier prevalence, on the same
# 624 shared rare-outlier genes as shared_rare_outlier_genes_3dot.png.
#
# The cohort dot counts TOPMed subjects (out of 5,489) carrying a LOFTEE
# high-confidence PTV at AF <= 0.5%. The gnomAD dot is the prevalence the same
# definition predicts from gnomAD v4.1 joint frequencies: cumulative pLoF allele
# frequency q over LOFTEE HC variants, converted to a carrier prevalence under
# Hardy-Weinberg, 1-(1-q)^2. Variants gnomAD has manually curated as not
# loss-of-function are dropped. See fetch_gnomad_plof_af.py.
#
# Two gnomAD flavours are shown because they answer different questions:
#   full       every HC pLoF variant under the same AF <= 0.5% ceiling
#   rare-only  only alleles at AF < 0.01%, the regime a 5,489-subject cohort
#              can actually sample, so it is the fairer like-for-like comparison
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork); library(ggrepel)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- SUP_GN

N_TOP    <- 5489
ZERO_X   <- 6e-5      # where "0 cohort PTV carriers" is drawn on the log axis
ACMG_COL <- "#5E3C99"
GNO_COL  <- "#01665E"
COH_COL  <- "#F0A202"

g3 <- read_tsv("shared_rare_outlier_genes_3dot.tsv", show_col_types = FALSE)
gn <- read_tsv("gnomad_plof_af_per_gene.tsv", show_col_types = FALSE) |>
  filter(freq_source != "absent") |>
  dplyr::select(gene_nv, n_hc_plof, n_used, n_used_rare, max_variant_af,
                n_curated, n_dropped_curated_not_lof,
                gnomad_prev = gnomad_plof_carrier_prev,
                gnomad_prev_rare = gnomad_plof_carrier_prev_rare)

d <- g3 |> left_join(gn, by = "gene_nv")
stopifnot(!any(duplicated(d$gene_nv)))
cat(sprintf("genes: %d | with gnomAD pLoF data: %d | missing: %d\n",
            nrow(d), sum(!is.na(d$gnomad_prev)), sum(is.na(d$gnomad_prev))))
dd <- filter(d, !is.na(gnomad_prev))

# --- how the two genetic estimates compare ---
sp  <- function(a, b) cor(a, b, method = "spearman", use = "complete.obs")
pos <- filter(dd, n_ptv > 0)
cat(sprintf("cohort PTV prevalence: median %s | gnomAD full: %s | gnomAD rare-only: %s\n",
            format(median(dd$prev_ptv), digits = 3), format(median(dd$gnomad_prev), digits = 3),
            format(median(dd$gnomad_prev_rare), digits = 3)))
cat(sprintf("median fold difference gnomAD/cohort (genes with >=1 carrier): full %.1fx, rare-only %.1fx\n",
            median(pos$gnomad_prev / pos$prev_ptv), median(pos$gnomad_prev_rare / pos$prev_ptv)))
cat(sprintf("Spearman vs cohort PTV: full %.3f, rare-only %.3f (n=%d with >=1 carrier)\n",
            sp(pos$gnomad_prev, pos$prev_ptv), sp(pos$gnomad_prev_rare, pos$prev_ptv), nrow(pos)))
cat(sprintf("Spearman vs TOPMed RNA outlier: full %.3f | vs AoU RNA outlier: %.3f\n",
            sp(dd$gnomad_prev, dd$prev_topmed), sp(dd$gnomad_prev, dd$prev_aou)))
cat(sprintf("genes with 0 cohort PTV carriers but a non-zero gnomAD pLoF estimate: %d of %d\n",
            sum(dd$n_ptv == 0 & dd$gnomad_prev > 0), sum(dd$n_ptv == 0)))
cat(sprintf("pLoF variants per gene in gnomAD: median %.0f (AF<0.01%%: %.0f); cohort PTV carriers: median %.0f\n",
            median(dd$n_used), median(dd$n_used_rare), median(dd$n_ptv)))
cat(sprintf("variants dropped by gnomAD curation verdict: %d across %d genes with any curation\n",
            sum(dd$n_dropped_curated_not_lof), sum(dd$n_curated > 0)))

MEAS <- c("TOPMed RNA-outlier prevalence", "AoU RNA-outlier prevalence",
          "TOPMed PTV-carrier prevalence (cohort)", "gnomAD v4.1 pLoF carrier prevalence")
PAL  <- setNames(c("#B2182B", "#2166AC", COH_COL, GNO_COL), MEAS)
SHP  <- setNames(c(16, 16, 17, 18), MEAS)

Y_BREAKS <- c(1e-4, 1e-3, 1e-2, 1e-1)
one_in <- function(v) paste0("1 in ", comma(round(1 / v)))

# ---- Panel A: all genes, four measures, ordered by the gnomAD estimate ----
ord <- dd |> arrange(desc(gnomad_prev)) |> mutate(rank = row_number())
la <- ord |>
  transmute(gene_nv, symbol, rank, is_acmg,
            `TOPMed RNA-outlier prevalence` = prev_topmed,
            `AoU RNA-outlier prevalence`    = prev_aou,
            `TOPMed PTV-carrier prevalence (cohort)` = prev_ptv,
            `gnomAD v4.1 pLoF carrier prevalence`    = gnomad_prev) |>
  pivot_longer(all_of(MEAS), names_to = "measure", values_to = "prev") |>
  mutate(measure = factor(measure, levels = MEAS),
         y = if_else(prev <= 0, ZERO_X, prev))
acmg_lab <- la |> filter(is_acmg) |> group_by(rank, symbol) |>
  summarise(y = max(y), .groups = "drop")

pa <- ggplot(la, aes(rank, y, color = measure, shape = measure)) +
  geom_vline(xintercept = acmg_lab$rank, color = ACMG_COL, alpha = 0.3, linewidth = 0.3) +
  geom_line(aes(group = rank), color = "grey85", linewidth = 0.22) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_text_repel(data = acmg_lab, aes(rank, y, label = symbol), inherit.aes = FALSE,
                  color = ACMG_COL, fontface = "bold", size = 2.7, nudge_y = 0.35,
                  direction = "y", min.segment.length = 0, segment.color = ACMG_COL,
                  segment.size = 0.3, seed = 1) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  scale_y_log10(breaks = c(ZERO_X, Y_BREAKS),
                labels = function(v) ifelse(abs(v - ZERO_X) < 1e-12, "0 carriers", one_in(v))) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(x = sprintf("gene, ranked by gnomAD pLoF carrier prevalence (n = %d)", nrow(ord)),
       y = "prevalence", title = "A. Four prevalence estimates per gene",
       subtitle = str_wrap(paste0(
         "The gnomAD dot is what the PTV definition predicts from 760k reference genomes; the cohort dot is ",
         "what 5,489 TOPMed subjects actually showed. Purple: ACMG SF v3.2 genes."), width = 125)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.8, color = "grey30"),
        legend.position = "top", legend.text = element_text(size = 8),
        plot.margin = margin(5.5, 8, 5.5, 8))

# ---- Panels B and C: gnomAD vs cohort, full and rare-only ----
scatter <- function(ycol, ttl, sub) {
  p <- dd |> mutate(x = if_else(prev_ptv <= 0, ZERO_X, prev_ptv), y = .data[[ycol]])
  rho <- sp(p$y[p$n_ptv > 0], p$prev_ptv[p$n_ptv > 0])
  ggplot(p, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = ZERO_X * 1.6, linetype = "dotted", color = "grey55", linewidth = 0.35) +
    geom_point(aes(color = n_ptv == 0), size = 1.1, alpha = 0.55) +
    scale_color_manual(values = c(`FALSE` = GNO_COL, `TRUE` = "grey60"), guide = "none") +
    geom_point(data = filter(p, is_acmg), size = 2.2, shape = 21, fill = NA,
               color = ACMG_COL, stroke = 0.7) +
    geom_text_repel(data = filter(p, is_acmg), aes(label = symbol), color = ACMG_COL,
                    fontface = "bold", size = 2.6, min.segment.length = 0,
                    segment.color = ACMG_COL, segment.size = 0.3, seed = 2) +
    scale_x_log10(breaks = c(ZERO_X, 2e-4, 5e-4, 1e-3, 5e-3, 2e-2),
                  labels = function(v) ifelse(abs(v - ZERO_X) < 1e-12, "0", one_in(v))) +
    scale_y_log10(labels = one_in) +
    annotate("text", x = ZERO_X, y = Inf, hjust = 0, vjust = 1.4, size = 2.7, color = "grey35",
             label = sprintf("Spearman = %.2f among genes with >=1 cohort carrier\ngrey = 0 cohort carriers (%d genes, plotted at left)",
                             rho, sum(p$n_ptv == 0))) +
    labs(x = "cohort PTV-carrier prevalence (TOPMed, n = 5,489)", y = ycol,
         title = ttl, subtitle = str_wrap(sub, width = 90)) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10.5),
          plot.subtitle = element_text(size = 7.6, color = "grey30"),
          plot.margin = margin(5.5, 8, 5.5, 8))
}

pb <- scatter("gnomad_prev",
              "B. gnomAD pLoF prevalence, all HC variants at AF <= 0.5%",
              paste0("Dashed line is y = x. Points above it are genes where the reference predicts more ",
                     "pLoF carriers than the cohort observed."))
pc <- scatter("gnomad_prev_rare",
              "C. gnomAD pLoF prevalence, rare alleles only (AF < 0.01%)",
              paste0("Restricting to alleles a 5,489-subject cohort could plausibly sample removes the ",
                     "handful of near-cap variants that dominate panel B."))

p <- (pa / (pb | pc)) + plot_layout(heights = c(1, 1.25)) +
  plot_annotation(
    title = "gnomAD v4.1 pLoF prevalence vs observed cohort PTV-carrier prevalence, 624 shared rare-outlier genes",
    subtitle = str_wrap(paste0(
      "Same gene set and same RNA-outlier numbers as shared_rare_outlier_genes_3dot.png. gnomAD prevalence is ",
      "1-(1-q)^2 with q the cumulative LOFTEE high-confidence pLoF allele frequency from the v4.1 joint callset, ",
      "excluding variants gnomAD curated as not loss-of-function. Axes read as 1 carrier in N people."), width = 150),
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 8, color = "grey30")))

W <- 13; H <- 10.4
ggsave(file.path(OUT_DIR, "gnomad_vs_cohort_ptv.pdf"), p, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "gnomad_vs_cohort_ptv.png"), p, width = W, height = H, dpi = 200)
ggsave(file.path(OUT_DIR, "gnomad_vs_cohort_ptv_basefont.pdf"), p, width = W, height = H,
       device = "pdf", useDingbats = FALSE)
write_tsv(arrange(d, desc(gnomad_prev)), "gnomad_vs_cohort_ptv.tsv")
message("wrote gnomad_vs_cohort_ptv.{pdf,png,tsv}")
