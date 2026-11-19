#!/usr/bin/env Rscript
# Same layout as shared_rare_outlier_genes_4dot_panelA_RBM20.png, marking every
# RBM20-like phenotype gene (RNA-outlier prevalence > PTV in TOPMed and AoU).
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork); library(ggrepel)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG4L

CUT <- 3
FOCUS_TSV <- file.path(FIG4L, "rbm20_like_phenotype_genes.tsv")

wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  tibble(lo = pmax(ctr - hw, 0), hi = pmin(ctr + hw, 1))
}

g_all <- read_tsv(file.path(FIG4D, "shared_rare_outlier_genes_4dot.tsv"), show_col_types = FALSE)
focus_meta <- read_tsv(FOCUS_TSV, show_col_types = FALSE)
N_TOP <- unique(g_all$N_topmed); stopifnot(length(N_TOP) == 1)
N_AOU <- unique(g_all$N_new);    stopifnot(length(N_AOU) == 1)

focus_sym <- focus_meta$symbol
stopifnot(all(focus_sym %in% g_all$symbol))
g_all <- g_all |>
  mutate(is_focus = symbol %in% focus_sym,
         rna_gt_ptv_both = prev_topmed > prev_ptv & prev_aou > prev_ptv)
stopifnot(all(g_all$rna_gt_ptv_both[g_all$is_focus]))

SEP_Y <- 5.00e-7; ZERO_Y <- 3.40e-7; SEP2_Y <- 2.40e-7; NA_Y <- 1.70e-7
MEAS <- c("TOPMed RNA-outlier prevalence", "AoU RNA-outlier prevalence",
          "TOPMed PTV-carrier prevalence",
          "gnomAD sum AF_exomes (TOPMed-overlapping PTVs)")
PAL  <- setNames(c("#B2182B", "#2166AC", "#F0A202", "#1B7837"), MEAS)
SHP  <- setNames(c(16, 16, 17, 18), MEAS)
FOCUS_COL <- "#5E3C99"

long <- g_all |>
  transmute(gene_nv, symbol, n_ptv, is_focus,
            `TOPMed RNA-outlier prevalence` = prev_topmed,
            `AoU RNA-outlier prevalence`    = prev_aou,
            `TOPMed PTV-carrier prevalence` = prev_ptv,
            `gnomAD sum AF_exomes (TOPMed-overlapping PTVs)` = sum_AF_exomes) |>
  pivot_longer(all_of(MEAS), names_to = "measure", values_to = "prevalence") |>
  mutate(measure = factor(measure, levels = MEAS),
         n = case_when(measure == MEAS[2] ~ N_AOU, measure == MEAS[4] ~ NA_real_,
                       TRUE ~ N_TOP),
         x = round(prevalence * n),
         no_data = is.na(prevalence),
         is_zero = !no_data & prevalence == 0,
         y_plot  = case_when(no_data ~ NA_Y, is_zero ~ ZERO_Y, TRUE ~ prevalence))
ci <- wilson(long$x, long$n)
long <- bind_cols(long, ci) |>
  mutate(lo = if_else(measure == MEAS[4], NA_real_, lo),
         hi = if_else(measure == MEAS[4], NA_real_, hi),
         lo_plot = pmax(lo, if_else(is_zero, ZERO_Y, SEP_Y * 1.02)),
         hi_plot = pmax(hi, if_else(is_zero, ZERO_Y, SEP_Y * 1.02)))

Y_BREAKS <- c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
y_scale <- function() list(
  scale_y_log10(breaks = Y_BREAKS,
                labels = function(v) paste0("1 in ", comma(round(1 / v))),
                sec.axis = sec_axis(~ ., breaks = Y_BREAKS,
                                    labels = function(v) format(v, scientific = TRUE, digits = 1))),
  geom_hline(yintercept = SEP_Y, linetype = "dotted", color = "grey55", linewidth = 0.35),
  geom_hline(yintercept = SEP2_Y, linetype = "dashed", color = "grey70", linewidth = 0.3))

# include 0-PTV genes so ACTA2 / KCNJ11 / COL1A2 / PRDM16 / SCN2A appear
ord <- g_all |> arrange(desc(prev_ptv), desc(prev_topmed)) |> mutate(rank = row_number())
la <- long |> left_join(dplyr::select(ord, gene_nv, rank), by = "gene_nv")
la_focus <- filter(la, is_focus, !no_data)
lab_focus <- la_focus |> group_by(rank, symbol) |> summarise(y = max(y_plot), .groups = "drop")
N_GN <- sum(!is.na(g_all$sum_AF_exomes))

pa_left <- ggplot(filter(la, !no_data), aes(rank, y_plot, color = measure, shape = measure)) +
  geom_vline(data = distinct(la_focus, rank), aes(xintercept = rank),
             color = FOCUS_COL, alpha = 0.35, linewidth = 0.4, inherit.aes = FALSE) +
  geom_line(data = filter(la, !no_data, !is_zero), aes(group = rank),
            color = "grey80", linewidth = 0.2) +
  geom_line(data = la_focus, aes(group = rank), color = FOCUS_COL, linewidth = 0.7) +
  geom_point(size = 0.85, alpha = 0.7) +
  geom_point(data = la_focus, size = 1.9) +
  geom_point(data = filter(la, no_data), shape = 4, size = 0.45, alpha = 0.45,
             color = "grey45", show.legend = FALSE) +
  geom_text_repel(data = lab_focus, aes(rank, y, label = symbol), inherit.aes = FALSE,
                  color = FOCUS_COL, fontface = "bold", size = 2.5, nudge_y = 0.28,
                  direction = "both", min.segment.length = 0, segment.color = FOCUS_COL,
                  segment.size = 0.3, max.overlaps = 40, seed = 1) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.6, alpha = 1)),
         shape = guide_legend(nrow = 2)) +
  y_scale() +
  annotate("text", x = nrow(ord), y = ZERO_Y, hjust = 1, vjust = -0.55, size = 2.4,
           color = "grey35",
           label = sprintf("%d genes with 0 PTV carriers (plotted here)", sum(g_all$n_ptv == 0))) +
  labs(x = sprintf("gene, ranked by PTV-carrier prevalence (n = %d)", nrow(g_all)),
       y = "frequency",
       title = "A. Shared rare-outlier genes",
       subtitle = sprintf(paste0("Purple marks %d phenotype genes with RNA-outlier prevalence ",
                                 "above PTV-carrier prevalence in both TOPMed and AoU."),
                          length(focus_sym))) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.4, color = "grey30"),
        legend.position = "top", legend.text = element_text(size = 8),
        legend.key.height = unit(10, "pt"),
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank())

# right panel: one column per focus gene, four measures dodged
ord_f <- focus_meta |>
  arrange(desc(n_topmed_ptv), desc(prev_topmed), symbol) |>
  mutate(lab = factor(symbol, levels = symbol))
lb <- long |>
  inner_join(dplyr::select(ord_f, gene_nv, lab), by = "gene_nv", relationship = "many-to-many") |>
  mutate(lab = factor(as.character(lab), levels = levels(ord_f$lab)))

n_missing <- sum(lb$no_data)
pb <- ggplot(filter(lb, !no_data), aes(lab, y_plot, color = measure, shape = measure)) +
  geom_line(data = filter(lb, !no_data, !is_zero), aes(group = lab),
            color = FOCUS_COL, linewidth = 0.45, alpha = 0.55) +
  geom_linerange(aes(ymin = lo_plot, ymax = hi_plot), linewidth = 0.35, alpha = 0.5,
                 position = position_dodge(width = 0.62), na.rm = TRUE) +
  geom_point(size = 1.9, position = position_dodge(width = 0.62)) +
  geom_point(data = filter(lb, no_data), shape = 4, size = 1.2, alpha = 0.6,
             color = "grey45", show.legend = FALSE) +
  scale_color_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = SHP, name = NULL) +
  y_scale() +
  annotate("text", x = 0.6, y = ZERO_Y, hjust = 0, vjust = -0.55, size = 2.3,
           color = "grey35", label = "zero band: 0 PTV carriers") +
  annotate("text", x = 0.6, y = NA_Y, hjust = 0, vjust = -0.55, size = 2.3,
           color = "grey45",
           label = sprintf("\u00d7 = no overlapping gnomAD AF_exomes (%d of %d genes)",
                           n_missing, n_distinct(lb$lab))) +
  labs(x = NULL, y = "frequency",
       title = "B. RBM20-like phenotype genes",
       subtitle = paste0("ClinGen Definitive/Strong disease genes. Wilson 95% CIs on ",
                         "cohort counts; no CI on gnomAD.")) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5, color = FOCUS_COL),
        plot.subtitle = element_text(size = 7.4, color = "grey30"),
        axis.text.x = element_text(angle = 55, hjust = 1, size = 7.2, face = "bold",
                                   color = FOCUS_COL),
        legend.position = "none",
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank())

n_ptv_pos <- sum(focus_meta$n_topmed_ptv > 0)
n_ptv_zero <- sum(focus_meta$n_topmed_ptv == 0)
p_side <- (pa_left | pb) + plot_layout(widths = c(1.55, 1), guides = "collect") +
  plot_annotation(
    title = sprintf("RBM20-like genes: RNA-outlier prevalence above PTV in TOPMed and AoU (\u2212log10 prevalence \u2265 %g in both cohorts)", CUT),
    subtitle = str_wrap(sprintf(paste0(
      "Left: all %d shared rare-outlier genes (including %d with 0 PTV carriers). ",
      "Right: %d ClinGen phenotype genes with RNA-outlier prevalence > PTV-carrier prevalence ",
      "in both TOPMed (n = %s) and AoU (n = %s). %d have \u22651 TOPMed PTV carrier; %d have 0. ",
      "AoU is RNA only (no AoU PTV call set). Green is gnomAD sum AF_exomes over TOPMed-overlapping PTVs."),
      nrow(g_all), sum(g_all$n_ptv == 0), length(focus_sym),
      comma(N_TOP), comma(N_AOU), n_ptv_pos, n_ptv_zero), width = 150),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8, color = "grey30"),
                  legend.position = "top"))

ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot_panelA_rbm20_like.png"), p_side,
       width = 14.2, height = 6.2, dpi = 200)
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot_panelA_rbm20_like.pdf"), p_side,
       width = 14.2, height = 6.2, device = cairo_pdf)
message("wrote shared_rare_outlier_genes_4dot_panelA_rbm20_like.{png,pdf}")
