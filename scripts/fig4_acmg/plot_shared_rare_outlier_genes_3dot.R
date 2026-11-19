#!/usr/bin/env Rscript
# Small-sample view of the genes in the top-right corner of the haplo+z<-3 panel:
# -log10 RNA-outlier prevalence >= 3 in BOTH TOPMed and AoU, i.e. prevalence
# <= 1 in 1,000 in both cohorts. These are the rare-outlier genes the two datasets
# agree on, and by construction they rest on very few carriers (2-5 TOPMed
# underliers, 2-8 AoU underliers).
#
# One column per gene, three dots:
#   TOPMed RNA-outlier prevalence   n_underlier / 5,489
#   AoU RNA-outlier prevalence      n_underlier / 8,895
#   TOPMed PTV carrier prevalence   n_PTV_carriers / 5,489   (ptv_all.tsv, kept subjects;
#                                   same convention as pct_ptv in the ClinGen HI tables)
# There is no AoU PTV call set, so the PTV dot is TOPMed-only.
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork); library(ggrepel)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG4D

PTV_TSV  <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/ptv_all.tsv"
ACMG_CSV <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/annotation/variant_annotations/acmg_sf_v3.2.csv"
ACMG_COL <- "#5E3C99"  # fourth colour, for the ACMG SF gene highlight
QUAL_COL <- "#1B7837"  # genes where both RNA-outlier prevalences exceed PTV
CUT      <- 3          # -log10 prevalence threshold on both axes
TOP_N    <- 40         # genes shown with labels and CIs in panel B
ZERO_Y   <- 6e-5       # where "0 PTV carriers" is drawn, below the log axis floor
SEP_Y    <- 1.15e-4    # separator between the zero band and real values

wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  tibble(lo = pmax(ctr - hw, 0), hi = pmin(ctr + hw, 1))
}

pairs <- read_tsv(file.path(FIG1C, "topmed_vs_new_z_-3_prevalence_per_gene.tsv"), show_col_types = FALSE) |>
  filter(n_topmed >= 2, n_new >= 2,
         -log10(prev_topmed) >= CUT, -log10(prev_new) >= CUT)
N_TOP <- max(pairs$N_topmed); N_AOU <- max(pairs$N_new)
cat(sprintf("genes with -log10 prevalence >= %g in both cohorts: %d\n", CUT, nrow(pairs)))

kept <- read_rds(KEPT_SUBJECTS)
ptv_carriers <- read_tsv(PTV_TSV, show_col_types = FALSE) |>
  filter(sample %in% kept) |>
  mutate(gene_nv = str_remove(gene, "\\..*$")) |>
  distinct(gene_nv, sample) |> count(gene_nv, name = "n_ptv")

# ACMG SF v3.2 secondary-finding genes, highlighted in panel A. The first column
# of the csv is the symbol; it repeats genes that have more than one disorder.
acmg_syms <- if (file.exists(ACMG_CSV)) {
  read_csv(ACMG_CSV, show_col_types = FALSE)[[1]] |> str_trim() |> toupper() |>
    (\(x) unique(x[!is.na(x) & str_detect(x, "^[A-Z0-9-]+$")]))()
} else {
  warning("ACMG csv not found, skipping the highlight"); character(0)
}

g <- pairs |>
  dplyr::select(gene_nv, symbol, LOEUF, LOEUF_tier, n_topmed, n_new) |>
  left_join(ptv_carriers, by = "gene_nv") |>
  mutate(n_ptv = replace_na(n_ptv, 0L),
         prev_topmed = n_topmed / N_TOP, prev_aou = n_new / N_AOU,
         prev_ptv = n_ptv / N_TOP,
         is_acmg = toupper(symbol) %in% acmg_syms,
         # genes whose functional signal outruns their coding signal in both
         # cohorts at once; genes with no PTV carrier at all are excluded so the
         # comparison rests on an observed, non-zero PTV rate
         rna_gt_ptv = n_ptv >= 1 & prev_topmed > prev_ptv & prev_aou > prev_ptv)

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
cat("\nLOEUF tiers:\n"); print(count(g, LOEUF_tier, sort = TRUE))
cat(sprintf("\nboth RNA prevalences > PTV prevalence, with >=1 PTV carrier: %d of %d genes\n",
            sum(g$rna_gt_ptv), nrow(g)))
cat(sprintf("  (excluded: %d genes with 0 PTV carriers, %d where PTV >= at least one RNA measure)\n",
            sum(g$n_ptv == 0), sum(g$n_ptv >= 1 & !g$rna_gt_ptv)))
cat(sprintf("\nACMG SF v3.2 genes among these %d: %d of %d on the list\n",
            nrow(g), sum(g$is_acmg), length(acmg_syms)))
print(as.data.frame(g |> filter(is_acmg) |>
  transmute(symbol, n_topmed, n_new, n_ptv,
            ptv_below_both = prev_ptv < prev_topmed & prev_ptv < prev_aou,
            in_panel_b = rna_gt_ptv) |>
  arrange(desc(ptv_below_both), symbol)), row.names = FALSE)

MEAS <- c("TOPMed RNA-outlier prevalence", "AoU RNA-outlier prevalence",
          "TOPMed PTV-carrier prevalence")
PAL  <- setNames(c("#B2182B", "#2166AC", "#F0A202"), MEAS)
SHP  <- setNames(c(16, 16, 17), MEAS)

long <- g |>
  transmute(gene_nv, symbol, LOEUF_tier, n_ptv, is_acmg, rna_gt_ptv,
            `TOPMed RNA-outlier prevalence` = prev_topmed,
            `AoU RNA-outlier prevalence`    = prev_aou,
            `TOPMed PTV-carrier prevalence` = prev_ptv) |>
  pivot_longer(all_of(MEAS), names_to = "measure", values_to = "prevalence") |>
  mutate(measure = factor(measure, levels = MEAS),
         n = case_when(measure == MEAS[2] ~ N_AOU, TRUE ~ N_TOP),
         x = round(prevalence * n),
         is_zero = prevalence == 0,
         y_plot = if_else(is_zero, ZERO_Y, prevalence))
long <- bind_cols(long, wilson(long$x, long$n)) |>
  mutate(lo_plot = pmax(lo, ZERO_Y * 1.05), hi_plot = pmax(hi, ZERO_Y * 1.05))

Y_BREAKS <- c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2)
# show_zero adds a labelled break at the zero band, which otherwise holds dots
# with no tick of their own. Only panel A has zero-PTV genes, so panel B keeps
# the plain breaks.
y_scale <- function(show_zero = FALSE, show_sep = TRUE) {
  brk <- if (show_zero) c(ZERO_Y, Y_BREAKS) else Y_BREAKS
  is_zero_brk <- function(v) !is.na(v) & abs(v - ZERO_Y) < 1e-12
  list(
    scale_y_log10(breaks = brk,
                  labels = function(v) ifelse(is_zero_brk(v), "0 carriers",
                                              paste0("1 in ", comma(round(1 / v)))),
                  sec.axis = sec_axis(~ ., breaks = brk,
                                      labels = function(v) ifelse(is_zero_brk(v), "0",
                                                                  format(v, scientific = TRUE, digits = 1)))),
    # the separator only means something where a zero band is drawn below it
    if (show_sep) geom_hline(yintercept = SEP_Y, linetype = "dotted",
                             color = "grey55", linewidth = 0.35))
}

# Panel A: every gene, ranked by PTV-carrier prevalence
ord_a <- g |> arrange(desc(prev_ptv), desc(prev_topmed)) |> mutate(rank = row_number())
la <- long |> left_join(dplyr::select(ord_a, gene_nv, rank), by = "gene_nv")
acmg_a   <- filter(la, is_acmg)
acmg_lab <- acmg_a |> group_by(rank, symbol) |> summarise(y = max(y_plot), .groups = "drop")
pa <- ggplot(la, aes(rank, y_plot, color = measure, shape = measure)) +
  geom_vline(xintercept = acmg_lab$rank, color = ACMG_COL, alpha = 0.30, linewidth = 0.3) +
  geom_line(aes(group = rank), color = "grey80", linewidth = 0.25) +
  geom_line(data = filter(la, rna_gt_ptv), aes(group = rank), color = QUAL_COL,
            linewidth = 0.4, alpha = 0.9) +
  geom_line(data = acmg_a, aes(group = rank), color = ACMG_COL, linewidth = 0.7) +
  geom_point(size = 0.85, alpha = 0.85) +
  geom_text_repel(data = acmg_lab, aes(rank, y, label = symbol), inherit.aes = FALSE,
                  color = ACMG_COL, fontface = "bold", size = 2.7, nudge_y = 0.3,
                  direction = "y", min.segment.length = 0, segment.color = ACMG_COL,
                  segment.size = 0.3, seed = 1) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  y_scale(show_zero = TRUE) +
  annotate("text", x = nrow(ord_a), y = ZERO_Y, hjust = 1, vjust = -0.8, size = 2.6,
           color = "grey35", label = sprintf("%d genes with 0 PTV carriers (plotted here)",
                                             sum(g$n_ptv == 0))) +
  labs(x = sprintf("gene, ranked by PTV-carrier prevalence (n = %d)", nrow(g)),
       y = "prevalence", title = "A. All shared rare-outlier genes",
       subtitle = str_wrap(sprintf(paste0(
         "One column per gene; the line joins its three measures. Green: the %d genes where both RNA-outlier ",
         "prevalences exceed PTV-carrier prevalence with at least 1 PTV carrier (panel B). Purple: the %d ACMG ",
         "SF v3.2 genes, none of which reaches the panel C top %d."),
         sum(g$rna_gt_ptv), sum(g$is_acmg), TOP_N), width = 125)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.8, color = "grey30"),
        legend.position = "top", legend.text = element_text(size = 8),
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank(),
        plot.margin = margin(5.5, 8, 5.5, 8))

# Panel B: every gene whose two RNA-outlier prevalences both beat its PTV rate,
# ordered by the smaller of the two RNA/PTV ratios so the most discordant lead.
ord_b <- g |> filter(rna_gt_ptv) |>
  mutate(gap = pmin(prev_topmed, prev_aou) / prev_ptv) |>
  arrange(desc(gap), desc(prev_topmed)) |>
  mutate(rank_b = row_number(), lab = fct_reorder(symbol, rank_b))
lb2 <- long |> inner_join(dplyr::select(ord_b, gene_nv, lab, rank_b), by = "gene_nv")
lev_b   <- levels(ord_b$lab)
acmg_b  <- ord_b$symbol[ord_b$is_acmg]
cols_b  <- ifelse(lev_b %in% acmg_b, ACMG_COL, "grey20")
cat(sprintf("panel B: %d genes, ACMG among them: %s\n", nrow(ord_b),
            if (length(acmg_b)) paste(acmg_b, collapse = ", ") else "none"))
pb <- ggplot(lb2, aes(lab, y_plot, color = measure, shape = measure)) +
  geom_vline(xintercept = which(lev_b %in% acmg_b), color = ACMG_COL,
             alpha = 0.35, linewidth = 0.45) +
  geom_line(aes(group = lab), color = QUAL_COL, linewidth = 0.3, alpha = 0.7) +
  geom_point(size = 1.15, alpha = 0.9) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  y_scale(show_sep = FALSE) +
  expand_limits(y = 1e-3) +   # keeps the 1-in-1,000 tick visible
  labs(x = sprintf("gene, most discordant first (n = %d)", nrow(ord_b)), y = "prevalence",
       title = sprintf(paste0("B. The %d genes where BOTH RNA-outlier prevalences exceed ",
                              "PTV-carrier prevalence (>= 1 PTV carrier)"), nrow(ord_b)),
       subtitle = str_wrap(sprintf(paste0(
         "Ordered by the smaller of the two RNA/PTV prevalence ratios, most discordant first. Purple label: ",
         "ACMG SF gene (%s; ACTA2 is excluded here because it has 0 PTV carriers, MYH11 and TNNI3 because ",
         "their PTV rate beats at least one RNA measure). Wilson CIs are omitted at this width, see panel C."),
         if (length(acmg_b)) paste(acmg_b, collapse = ", ") else "none"), width = 145)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.8, color = "grey30"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6,
                                   color = cols_b),
        legend.position = "none",
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank(),
        plot.margin = margin(5.5, 8, 5.5, 8))

# Panel C: the labeled, readable top of the PTV ranking, with binomial CIs
top <- ord_a |> slice_head(n = TOP_N) |> mutate(lab = fct_reorder(symbol, rank))
lb <- long |> inner_join(dplyr::select(top, gene_nv, lab, rank), by = "gene_nv")
pc <- ggplot(lb, aes(lab, y_plot, color = measure, shape = measure)) +
  geom_line(aes(group = lab), color = "grey80", linewidth = 0.3) +
  geom_linerange(aes(ymin = lo_plot, ymax = hi_plot), linewidth = 0.35, alpha = 0.55,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 1.7, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  y_scale() +
  labs(x = NULL, y = "prevalence",
       title = sprintf("C. Top %d of the 624 by PTV-carrier prevalence, with 95%% binomial CIs", TOP_N),
       subtitle = str_wrap(paste0(
         "CIs are Wilson intervals on the underlying counts; they are wide because every outlier estimate ",
         "rests on 2-8 subjects. No gene here meets the panel B criterion: ranking by PTV prevalence puts ",
         "the PTV dot on top by construction."), width = 145)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.8, color = "grey30"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 6.6),
        legend.position = "none",
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank(),
        plot.margin = margin(5.5, 8, 5.5, 8))

p <- (pa / pb / pc) + plot_layout(heights = c(1, 1.2, 1.15)) +
  plot_annotation(
    title = sprintf("Shared rare RNA-outlier genes (\u2212log10 prevalence \u2265 %g in both cohorts): outlier vs PTV prevalence", CUT),
    subtitle = str_wrap(sprintf(paste0(
      "%d genes from the haplo + z<-3 comparison. Left axis reads as 1 carrier in N subjects, ",
      "right axis as a fraction. TOPMed outlier and PTV prevalence share the same denominator ",
      "(%s subjects); AoU outlier prevalence uses %s subjects. There is no AoU PTV call set, so ",
      "the PTV dot is TOPMed-only."), nrow(g), comma(N_TOP), comma(N_AOU)), width = 120),
    theme = theme(plot.title = element_text(face = "bold", size = 11.5),
                  plot.subtitle = element_text(size = 8, color = "grey30")))

# wider than the original two-panel version so 108 gene labels fit in panel B
W <- 13.0; H <- 12.8
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_3dot.pdf"), p, width = W, height = H, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_3dot.png"), p, width = W, height = H, dpi = 200)
# base PDF device as well: some viewers dislike cairo's embedded font subsets
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_3dot_basefont.pdf"), p, width = W, height = H,
       device = "pdf", useDingbats = FALSE)
write_tsv(arrange(g, desc(prev_ptv), desc(prev_topmed)), "shared_rare_outlier_genes_3dot.tsv")
message("wrote shared_rare_outlier_genes_3dot.{pdf,png,tsv}")
