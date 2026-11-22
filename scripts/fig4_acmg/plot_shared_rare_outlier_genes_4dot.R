#!/usr/bin/env Rscript
# Same gene set as plot_shared_rare_outlier_genes_3dot.R, with a fourth dot per
# gene: the sum of gnomAD AF_exomes over TOPMed-overlapping PTVs in that gene.
#
# Gene set: -log10 RNA-outlier prevalence >= 3 in BOTH TOPMed and AoU, i.e.
# prevalence <= 1 in 1,000 in both cohorts (624 genes, each resting on 2-5 TOPMed
# and 2-8 AoU underliers).
#
# Four dots:
#   1 TOPMed RNA-outlier prevalence   n_underlier / 5,489
#   2 AoU RNA-outlier prevalence      n_underlier / 8,895
#   3 TOPMed PTV-carrier prevalence   n_PTV_carriers / 5,489   (ptv_all.tsv, kept subjects)
#   4 gnomAD sum AF_exomes            external allele frequency, NO cohort denominator
#
# Definition of dot 4 -----------------------------------------------------------
# Local gnomAD joint v4.1 exomes VCF (INFO/AF_exomes), subset to GENCODE v30 gene
# bodies for these 624 genes (BED from rtracklayer), overlap-joined onto those
# gene GRanges with plyranges::join_overlap_left, then restricted to sites that
# are also a PTV in ptv_all.tsv. Per gene we sum AF_exomes. gnomAD has no
# individual-level genotypes, so this slightly undercounts if one person carries
# more than one PTV in the gene. See sum_af_exome_topmed_overlap.R.
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork); library(ggrepel)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG4D

AF_TSV     <- file.path(FIG4D, "sum_af_exome_topmed_overlap.tsv")
CUT        <- 3          # -log10 prevalence threshold on both axes
TOP_N      <- 40         # genes shown with labels and CIs in panels B and C
# Bands below the axis floor. Smallest summed AF_exomes is 6.84e-7, so every
# band sits below that; checked with stopifnot() further down.
SEP_Y      <- 5.00e-7    # separator: bands below, real frequencies above
ZERO_Y     <- 3.40e-7    # where "0 PTV carriers" is drawn
SEP2_Y     <- 2.40e-7    # separator between the zero band and the no-gnomAD band
NA_Y       <- 1.70e-7    # where "no overlapping gnomAD AF" is drawn

wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  tibble(lo = pmax(ctr - hw, 0), hi = pmin(ctr + hw, 1))
}

pairs_all <- read_tsv(file.path(FIG1C, "topmed_vs_new_z_-3_prevalence_per_gene.tsv"), show_col_types = FALSE)
pairs <- pairs_all |>
  filter(n_topmed >= 2, n_new >= 2,
         -log10(prev_topmed) >= CUT, -log10(prev_new) >= CUT)
# Keep Fig4 example genes even if they fall outside the shared-rare cut
# (MYH11 at PC1055 has only 1 TOPMed underlier).
FOCUS_KEEP <- c("RBM20", "MYH11")
extra <- pairs_all |>
  filter(symbol %in% FOCUS_KEEP, !gene_nv %in% pairs$gene_nv)
if (nrow(extra) > 0) {
  cat(sprintf("force-including Fig4 focus genes outside cut: %s\n",
              paste(unique(extra$symbol), collapse = ", ")))
  pairs <- bind_rows(pairs, extra)
}
N_TOP <- max(pairs$N_topmed); N_AOU <- max(pairs$N_new)
cat(sprintf("genes with -log10 prevalence >= %g in both cohorts (+focus): %d\n", CUT, nrow(pairs)))

kept <- read_rds(KEPT_SUBJECTS)
ptv_carriers <- read_tsv(PTV_TSV, show_col_types = FALSE) |>
  filter(sample %in% kept) |>
  mutate(gene_nv = str_remove(gene, "\\..*$")) |>
  distinct(gene_nv, sample) |> count(gene_nv, name = "n_ptv")

af <- read_tsv(AF_TSV, show_col_types = FALSE) |>
  dplyr::select(gene_nv, n_ptv_overlap, sum_AF_exomes, max_AF_exomes)

g <- pairs |>
  dplyr::select(gene_nv, symbol, LOEUF, LOEUF_tier, n_topmed, n_new) |>
  left_join(ptv_carriers, by = "gene_nv") |>
  left_join(af, by = "gene_nv") |>
  mutate(n_ptv = replace_na(n_ptv, 0L),
         n_ptv_overlap = replace_na(n_ptv_overlap, 0L),
         prev_topmed = n_topmed / N_TOP, prev_aou = n_new / N_AOU,
         prev_ptv = n_ptv / N_TOP)
stopifnot(nrow(g) == nrow(pairs))

N_GN <- sum(!is.na(g$sum_AF_exomes))

# --- small-sample analysis ---
cat(sprintf("cohorts: TOPMed n=%s, AoU n=%s (PTV denominator = TOPMed kept subjects, n=%s)\n",
            comma(N_TOP), comma(N_AOU), comma(length(kept))))
cat(sprintf("underlier counts: TOPMed %d-%d, AoU %d-%d per gene\n",
            min(g$n_topmed), max(g$n_topmed), min(g$n_new), max(g$n_new)))
cat(sprintf("genes with 0 PTV carriers: %d of %d (%.1f%%)\n",
            sum(g$n_ptv == 0), nrow(g), 100 * mean(g$n_ptv == 0)))
cat(sprintf("PTV carriers per gene: median %g, IQR [%g, %g], max %d\n",
            median(g$n_ptv), quantile(g$n_ptv, .25), quantile(g$n_ptv, .75), max(g$n_ptv)))
cat(sprintf("PTV prevalence exceeds TOPMed outlier prevalence in %d of %d genes (%.1f%%)\n",
            sum(g$prev_ptv > g$prev_topmed), nrow(g), 100 * mean(g$prev_ptv > g$prev_topmed)))
cat(sprintf("Spearman, PTV vs TOPMed outlier prevalence: %.3f (p=%.3g)\n",
            cor(g$n_ptv, g$n_topmed, method = "spearman"),
            cor.test(g$n_ptv, g$n_topmed, method = "spearman", exact = FALSE)$p.value))
cat(sprintf("Spearman, TOPMed vs AoU outlier prevalence within this subset: %.3f\n",
            cor(g$n_topmed, g$n_new, method = "spearman")))

cat("\n--- gnomAD sum AF_exomes (dot 4) ---\n")
cat(sprintf("genes with >=1 TOPMed PTV overlapping a gnomAD AF_exomes site: %d of %d (%.1f%%); without: %d (%.1f%%)\n",
            N_GN, nrow(g), 100 * N_GN / nrow(g), nrow(g) - N_GN, 100 * (1 - N_GN / nrow(g))))
gq <- quantile(g$sum_AF_exomes, c(.25, .5, .75), na.rm = TRUE)
cat(sprintf("sum AF_exomes: median %.3g (1 in %s), IQR [%.3g, %.3g], range [%.3g, %.3g]\n",
            gq[2], comma(round(1 / gq[2])), gq[1], gq[3],
            min(g$sum_AF_exomes, na.rm = TRUE), max(g$sum_AF_exomes, na.rm = TRUE)))
for (v in c("prev_topmed", "prev_aou", "prev_ptv")) {
  ct <- cor.test(g$sum_AF_exomes, g[[v]], method = "spearman", exact = FALSE)
  cat(sprintf("Spearman, sum AF_exomes vs %-12s: rho = %+.3f (p = %.3g, n = %d)\n",
              v, unname(ct$estimate), ct$p.value, sum(!is.na(g$sum_AF_exomes))))
}
cat("\nLOEUF tiers:\n"); print(count(g, LOEUF_tier, sort = TRUE))

MEAS <- c("TOPMed RNA-outlier prevalence", "AoU RNA-outlier prevalence",
          "TOPMed PTV-carrier prevalence",
          "gnomAD sum AF_exomes (TOPMed-overlapping PTVs)")
PAL  <- setNames(c("#B2182B", "#2166AC", "#F0A202", "#1B7837"), MEAS)
SHP  <- setNames(c(16, 16, 17, 18), MEAS)
ACMG_COL <- "#5E3C99"

long <- g |>
  transmute(gene_nv, symbol, LOEUF_tier, n_ptv,
            `TOPMed RNA-outlier prevalence` = prev_topmed,
            `AoU RNA-outlier prevalence`    = prev_aou,
            `TOPMed PTV-carrier prevalence` = prev_ptv,
            `gnomAD sum AF_exomes (TOPMed-overlapping PTVs)` = sum_AF_exomes) |>
  pivot_longer(all_of(MEAS), names_to = "measure", values_to = "prevalence") |>
  mutate(measure = factor(measure, levels = MEAS),
         n = case_when(measure == MEAS[2] ~ N_AOU, measure == MEAS[4] ~ NA_real_,
                       TRUE ~ N_TOP),
         x = round(prevalence * n),
         no_data = is.na(prevalence),                       # gnomAD only
         is_zero = !no_data & prevalence == 0,              # 0 PTV carriers
         y_plot  = case_when(no_data ~ NA_Y, is_zero ~ ZERO_Y, TRUE ~ prevalence))
ci <- wilson(long$x, long$n)
long <- bind_cols(long, ci) |>
  mutate(lo = if_else(measure == MEAS[4], NA_real_, lo),     # no CI on summed AF
         hi = if_else(measure == MEAS[4], NA_real_, hi),
         lo_plot = pmax(lo, if_else(is_zero, ZERO_Y, SEP_Y * 1.02)),
         hi_plot = pmax(hi, if_else(is_zero, ZERO_Y, SEP_Y * 1.02)))

stopifnot(min(long$prevalence[long$prevalence > 0], na.rm = TRUE) > SEP_Y)

Y_BREAKS <- c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
y_scale <- function() list(
  scale_y_log10(breaks = Y_BREAKS,
                labels = function(v) paste0("1 in ", comma(round(1 / v))),
                sec.axis = sec_axis(~ ., breaks = Y_BREAKS,
                                    labels = function(v) format(v, scientific = TRUE, digits = 1))),
  geom_hline(yintercept = SEP_Y, linetype = "dotted", color = "grey55", linewidth = 0.35),
  geom_hline(yintercept = SEP2_Y, linetype = "dashed", color = "grey70", linewidth = 0.3))

# Panel A: every gene, ranked by PTV-carrier prevalence
ord_a <- g |> arrange(desc(prev_ptv), desc(prev_topmed)) |> mutate(rank = row_number())
la <- long |> left_join(dplyr::select(ord_a, gene_nv, rank), by = "gene_nv")
pa <- ggplot(filter(la, !no_data), aes(rank, y_plot, color = measure, shape = measure)) +
  geom_line(data = filter(la, !no_data, !is_zero), aes(group = rank),
            color = "grey80", linewidth = 0.25) +
  geom_point(size = 0.9, alpha = 0.85) +
  geom_point(data = filter(la, no_data), shape = 4, size = 0.5, alpha = 0.5,
             color = "grey45", show.legend = FALSE) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.6, alpha = 1)),
         shape = guide_legend(nrow = 2)) +
  y_scale() +
  annotate("text", x = nrow(ord_a), y = ZERO_Y, hjust = 1, vjust = -0.55, size = 2.5,
           color = "grey35", label = sprintf("%d genes with 0 PTV carriers (plotted here)",
                                             sum(g$n_ptv == 0))) +
  annotate("text", x = nrow(ord_a), y = NA_Y, hjust = 1, vjust = -0.55, size = 2.5,
           color = "grey45",
           label = sprintf("%d genes with no TOPMed PTV in gnomAD AF_exomes, drawn as \u00d7",
                           nrow(g) - N_GN)) +
  labs(x = sprintf("gene, ranked by PTV-carrier prevalence (n = %d)", nrow(g)),
       y = "frequency", title = "A. All shared rare-outlier genes",
       subtitle = paste0("One column per gene. Grey line joins the measures for that gene. ",
                         "Dotted line separates the zero band; dashed line the no-gnomAD band.")) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.8, color = "grey30"),
        legend.position = "top", legend.text = element_text(size = 8),
        legend.key.height = unit(10, "pt"),
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank())

labeled_panel <- function(sel, ttl, sub) {
  lb <- long |> inner_join(dplyr::select(sel, gene_nv, lab), by = "gene_nv")
  n_missing <- sum(lb$no_data)
  ggplot(filter(lb, !no_data), aes(lab, y_plot, color = measure, shape = measure)) +
    geom_line(data = filter(lb, !no_data, !is_zero), aes(group = lab),
              color = "grey80", linewidth = 0.3) +
    geom_linerange(aes(ymin = lo_plot, ymax = hi_plot), linewidth = 0.35, alpha = 0.5,
                   position = position_dodge(width = 0.62), na.rm = TRUE) +
    geom_point(size = 1.7, position = position_dodge(width = 0.62)) +
    geom_point(data = filter(lb, no_data), shape = 4, size = 1.1, alpha = 0.6,
               color = "grey45", show.legend = FALSE) +
    scale_color_manual(values = PAL, name = NULL) +
    scale_shape_manual(values = SHP, name = NULL) +
    y_scale() +
    annotate("text", x = 0.6, y = ZERO_Y, hjust = 0, vjust = -0.55, size = 2.4,
             color = "grey35", label = "zero band: 0 PTV carriers") +
    annotate("text", x = 0.6, y = NA_Y, hjust = 0, vjust = -0.55, size = 2.4,
             color = "grey45",
             label = sprintf("\u00d7 = no overlapping gnomAD AF_exomes (%d of %d genes)",
                             n_missing, n_distinct(lb$lab))) +
    labs(x = NULL, y = "frequency", title = ttl, subtitle = sub) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10.5),
          plot.subtitle = element_text(size = 7.8, color = "grey30"),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 6.6),
          legend.position = "none",
          axis.text.y.right = element_text(size = 7, color = "grey45"),
          axis.title.y.right = element_blank())
}

top_b <- ord_a |> slice_head(n = TOP_N) |> mutate(lab = fct_reorder(symbol, rank))
pb <- labeled_panel(
  top_b,
  sprintf("B. Top %d genes by PTV-carrier prevalence, with 95%% binomial CIs", TOP_N),
  sprintf(paste0("CIs are Wilson intervals on the cohort counts (wide: every outlier estimate ",
                 "rests on 2-8 subjects). The green dot has no CI; it is a sum of gnomAD ",
                 "AF_exomes. %d of these %d genes have a gnomAD value."),
          sum(!is.na(top_b$sum_AF_exomes)), nrow(top_b)))

top_c <- g |> filter(!is.na(sum_AF_exomes)) |>
  arrange(desc(sum_AF_exomes)) |> slice_head(n = TOP_N) |>
  mutate(rank_c = row_number(), lab = fct_reorder(symbol, rank_c))
pc <- labeled_panel(
  top_c,
  sprintf("C. Top %d of the %d genes with a gnomAD value, ranked by summed AF_exomes",
          nrow(top_c), N_GN),
  paste0("Green is sum(AF_exomes) over TOPMed-overlapping PTVs in the GENCODE gene body ",
         "(no CI). The other three keep Wilson intervals on cohort counts."))

p <- (pa / pb / pc) + plot_layout(heights = c(1, 1.05, 1.05)) +
  plot_annotation(
    title = sprintf("Shared rare RNA-outlier genes (\u2212log10 prevalence \u2265 %g in both cohorts): outlier, PTV and gnomAD AF_exomes", CUT),
    subtitle = str_wrap(sprintf(paste0(
      "%d genes from the haplo + z<-3 comparison. Left axis reads as 1 in N, right axis as a ",
      "fraction. TOPMed outlier and PTV prevalence share the same denominator (%s subjects); ",
      "AoU outlier prevalence uses %s subjects; there is no AoU PTV call set, so the PTV dot ",
      "is TOPMed-only. The green dot has NO cohort denominator: it is the sum of gnomAD v4.1 ",
      "exome AF_exomes over PTVs that (i) overlap the GENCODE v30 gene body and (ii) are also ",
      "present in the TOPMed PTV table (ptv_all.tsv). Summing allele frequencies slightly ",
      "undercounts if one person carries more than one PTV in the gene. %d of %d genes have ",
      "at least one such overlapping site."),
      nrow(g), comma(N_TOP), comma(N_AOU), N_GN, nrow(g)), width = 128),
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 7.6, color = "grey30")))

W <- 11.5; H <- 14.0
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot.pdf"), p, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot.png"), p, width = W, height = H, dpi = 200)
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot_basefont.pdf"), p, width = W, height = H,
       device = "pdf", useDingbats = FALSE)

out <- g |>
  transmute(gene_nv, symbol, LOEUF, LOEUF_tier,
            n_topmed, N_topmed = N_TOP, prev_topmed,
            n_new, N_new = N_AOU, prev_aou,
            n_ptv, prev_ptv,
            n_ptv_overlap, sum_AF_exomes, max_AF_exomes,
            has_gnomad = !is.na(sum_AF_exomes)) |>
  arrange(desc(prev_ptv), desc(prev_topmed))
write_tsv(out, file.path(FIG4D, "shared_rare_outlier_genes_4dot.tsv"))
message("wrote shared_rare_outlier_genes_4dot.{pdf,png,tsv} and _basefont.pdf")

# --- Panel A only, genes with >=1 TOPMed PTV carrier -------------------------------
g_nz <- g |> filter(n_ptv > 0)
n_dropped <- nrow(g) - nrow(g_nz)
N_GN_nz <- sum(!is.na(g_nz$sum_AF_exomes))
cat(sprintf("panel-A-only plot: %d genes with >=1 PTV carrier; dropped %d with 0 PTV carriers\n",
            nrow(g_nz), n_dropped))

long_nz <- long |> filter(gene_nv %in% g_nz$gene_nv)
ord_nz <- g_nz |> arrange(desc(prev_ptv), desc(prev_topmed)) |> mutate(rank = row_number())
la_nz <- long_nz |> left_join(dplyr::select(ord_nz, gene_nv, rank), by = "gene_nv")

pa_nz <- ggplot(filter(la_nz, !no_data), aes(rank, y_plot, color = measure, shape = measure)) +
  geom_line(data = filter(la_nz, !no_data, !is_zero), aes(group = rank),
            color = "grey80", linewidth = 0.25) +
  geom_point(size = 1.15, alpha = 0.85) +
  geom_point(data = filter(la_nz, no_data), shape = 4, size = 0.6, alpha = 0.5,
             color = "grey45", show.legend = FALSE) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.6, alpha = 1)),
         shape = guide_legend(nrow = 2)) +
  y_scale() +
  annotate("text", x = nrow(ord_nz), y = NA_Y, hjust = 1, vjust = -0.55, size = 2.7,
           color = "grey45",
           label = sprintf("%d genes with no TOPMed PTV in gnomAD AF_exomes, drawn as \u00d7",
                           nrow(g_nz) - N_GN_nz)) +
  labs(x = sprintf("gene, ranked by PTV-carrier prevalence (n = %d)", nrow(g_nz)),
       y = "frequency") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        legend.position = "top", legend.text = element_text(size = 8),
        legend.key.height = unit(10, "pt"),
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank())

p_nz <- pa_nz +
  plot_annotation(
    title = sprintf("Shared rare RNA-outlier genes with \u22651 TOPMed PTV carrier (\u2212log10 prevalence \u2265 %g in both cohorts)", CUT),
    subtitle = str_wrap(sprintf(paste0(
      "Panel A only. The %d genes with 0 PTV carriers among %s TOPMed subjects are omitted ",
      "(original set was %d genes). %d genes remain. Green is sum of gnomAD v4.1 AF_exomes ",
      "over TOPMed-overlapping PTVs in the GENCODE v30 gene body; summing slightly undercounts ",
      "if one person carries more than one PTV in the gene. %d of %d genes have a gnomAD value."),
      n_dropped, comma(N_TOP), nrow(g), nrow(g_nz), N_GN_nz, nrow(g_nz)), width = 128),
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 8, color = "grey30")))

W2 <- 11.5; H2 <- 5.6
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot_panelA_ptv.png"), p_nz, width = W2, height = H2, dpi = 200)
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot_panelA_ptv.pdf"), p_nz, width = W2, height = H2, device = cairo_pdf)
message("wrote shared_rare_outlier_genes_4dot_panelA_ptv.{png,pdf}")

# --- Side-by-side: Panel A + RBM20 (ACMG gene of interest) ---------------------
FOCUS <- "RBM20"
focus_row <- g_nz |> filter(symbol == FOCUS)
if (nrow(focus_row) != 1) stop("expected exactly one ", FOCUS, " row in Panel A")

focus_rank <- ord_nz$rank[ord_nz$symbol == FOCUS]
la_focus <- filter(la_nz, symbol == FOCUS, !no_data)
lab_focus <- la_focus |> group_by(rank, symbol) |> summarise(y = max(y_plot), .groups = "drop")

pa_left <- ggplot(filter(la_nz, !no_data), aes(rank, y_plot, color = measure, shape = measure)) +
  geom_vline(xintercept = focus_rank, color = ACMG_COL, alpha = 0.45, linewidth = 0.45) +
  geom_line(data = filter(la_nz, !no_data, !is_zero), aes(group = rank),
            color = "grey80", linewidth = 0.25) +
  geom_line(data = la_focus, aes(group = rank), color = ACMG_COL, linewidth = 0.85) +
  geom_point(size = 1.05, alpha = 0.8) +
  geom_point(data = la_focus, size = 2.1) +
  geom_point(data = filter(la_nz, no_data), shape = 4, size = 0.55, alpha = 0.5,
             color = "grey45", show.legend = FALSE) +
  geom_text_repel(data = lab_focus, aes(rank, y, label = symbol), inherit.aes = FALSE,
                  color = ACMG_COL, fontface = "bold", size = 3.1, nudge_y = 0.35,
                  direction = "y", min.segment.length = 0, segment.color = ACMG_COL,
                  segment.size = 0.35, seed = 1) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.6, alpha = 1)),
         shape = guide_legend(nrow = 2)) +
  y_scale() +
  labs(x = sprintf("gene, ranked by PTV-carrier prevalence (n = %d)", nrow(g_nz)),
       y = "frequency",
       title = "A. Shared rare-outlier genes with \u22651 TOPMed PTV carrier",
       subtitle = sprintf("Purple marks %s (rank %d of %d), the only ACMG SF v3.2 gene in this panel whose PTV and gnomAD frequencies both sit below both RNA-outlier prevalences.",
                          FOCUS, focus_rank, nrow(g_nz))) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.4, color = "grey30"),
        legend.position = "top", legend.text = element_text(size = 8),
        legend.key.height = unit(10, "pt"),
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank())

# short x labels for the single-gene panel
MEAS_SHORT <- c("TOPMed RNA\noutlier", "AoU RNA\noutlier", "TOPMed PTV\ncarriers", "gnomAD\nsum AF_exomes")
rb <- long_nz |>
  filter(symbol == FOCUS) |>
  mutate(meas_lab = factor(MEAS_SHORT[as.integer(measure)], levels = MEAS_SHORT),
         n_lab = case_when(
           measure == MEAS[1] ~ sprintf("%d / %s", focus_row$n_topmed, comma(N_TOP)),
           measure == MEAS[2] ~ sprintf("%d / %s", focus_row$n_new, comma(N_AOU)),
           measure == MEAS[3] ~ sprintf("%d / %s", focus_row$n_ptv, comma(N_TOP)),
           TRUE ~ "sum of AFs"
         ))

pb_focus <- ggplot(filter(rb, !no_data), aes(meas_lab, y_plot, color = measure, shape = measure)) +
  geom_line(aes(group = 1), color = ACMG_COL, linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 3.2) +
  geom_text(aes(label = n_lab), color = "grey25", size = 2.6, vjust = -1.15,
            show.legend = FALSE) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  y_scale() +
  labs(x = NULL, y = "frequency",
       title = sprintf("B. %s (ACMG SF v3.2)", FOCUS),
       subtitle = paste0("Dilated cardiomyopathy. RNA-outlier prevalence exceeds ",
                         "both the TOPMed PTV-carrier prevalence and the gnomAD ",
                         "summed AF_exomes.")) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5, color = ACMG_COL),
        plot.subtitle = element_text(size = 7.4, color = "grey30"),
        legend.position = "none",
        axis.text.x = element_text(size = 7.4, lineheight = 0.95),
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank())

p_side <- (pa_left | pb_focus) + plot_layout(widths = c(2.2, 1), guides = "collect") +
  plot_annotation(
    title = sprintf("RBM20: RNA-outlier prevalence above PTV and gnomAD frequency (\u2212log10 prevalence \u2265 %g in both cohorts)", CUT),
    subtitle = str_wrap(sprintf(paste0(
      "Left: the %d Panel A genes (\u22651 TOPMed PTV carrier; %d genes with 0 PTV carriers omitted). ",
      "Right: %s, the only ACMG SF v3.2 gene in this set with TOPMed PTV prevalence and gnomAD ",
      "sum AF_exomes both below both RNA-outlier prevalences. MYH11 and TNNI3 are also ACMG genes ",
      "in Panel A but fail that rule; ACTA2 is omitted because it has 0 PTV carriers."),
      nrow(g_nz), n_dropped, FOCUS), width = 145),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8, color = "grey30"),
                  legend.position = "top"))

W3 <- 13.2; H3 <- 5.8
ggsave(file.path(FIG4A, "shared_rare_outlier_genes_4dot_panelA_RBM20.png"), p_side, width = W3, height = H3, dpi = 200)
ggsave(file.path(FIG4A, "shared_rare_outlier_genes_4dot_panelA_RBM20.pdf"), p_side, width = W3, height = H3, device = cairo_pdf)
message("wrote shared_rare_outlier_genes_4dot_panelA_RBM20.{png,pdf}")
