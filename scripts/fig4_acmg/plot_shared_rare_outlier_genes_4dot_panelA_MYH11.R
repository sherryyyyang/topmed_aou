#!/usr/bin/env Rscript
# Same layout as shared_rare_outlier_genes_4dot_panelA_RBM20.png, but for the
# opposite rule: TOPMed PTV-carrier prevalence > both RNA-outlier prevalences.
# ACMG SF v3.2 list:
#   /net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/annotation/variant_annotations/acmg_sf_v3.2.csv
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork); library(ggrepel)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG4B

ACMG_CSV <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/annotation/variant_annotations/acmg_sf_v3.2.csv"
CUT <- 3
FOCUS <- "MYH11"

wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  tibble(lo = pmax(ctr - hw, 0), hi = pmin(ctr + hw, 1))
}

g_all <- read_tsv(file.path(FIG4D, "shared_rare_outlier_genes_4dot.tsv"), show_col_types = FALSE)
N_TOP <- unique(g_all$N_topmed); stopifnot(length(N_TOP) == 1)
N_AOU <- unique(g_all$N_new);    stopifnot(length(N_AOU) == 1)

acmg <- read_csv(ACMG_CSV, show_col_types = FALSE)
names(acmg) <- str_trim(names(acmg))
acmg_sym <- unique(toupper(str_trim(acmg$Gene)))
acmg_meta <- acmg |>
  mutate(sym = toupper(str_trim(Gene))) |>
  distinct(sym, .keep_all = TRUE)

g <- g_all |>
  filter(n_ptv > 0) |>
  mutate(is_acmg = toupper(symbol) %in% acmg_sym,
         ptv_gt_both_rna = prev_ptv > prev_topmed & prev_ptv > prev_aou)

n_dropped <- nrow(g_all) - nrow(g)
hits <- g |> filter(ptv_gt_both_rna)
acmg_hits <- hits |> filter(is_acmg)
cat(sprintf("Panel A (>=1 PTV): %d of %d genes\n", nrow(g), nrow(g_all)))
cat(sprintf("PTV prevalence > both RNA-outlier prevalences: %d of %d\n", nrow(hits), nrow(g)))
cat(sprintf("ACMG SF v3.2 among those: %d\n", nrow(acmg_hits)))
print(acmg_hits |> dplyr::select(symbol, n_topmed, n_new, n_ptv, prev_topmed, prev_aou, prev_ptv, sum_AF_exomes))
stopifnot(nrow(acmg_hits) == 1, acmg_hits$symbol == FOCUS)

write_tsv(
  hits |>
    left_join(acmg_meta |> dplyr::select(sym, phenotype = `Disease/Phentyope`,
                                        category = `Phenotype Category`,
                                        inheritance = Inheritance),
              by = c("symbol" = "sym")) |>
    transmute(gene_nv, symbol, is_acmg, phenotype, category, inheritance,
              n_topmed, N_topmed, prev_topmed, n_aou = n_new, N_aou = N_new, prev_aou,
              n_ptv, prev_ptv, sum_AF_exomes) |>
    arrange(desc(prev_ptv), desc(prev_topmed)),
  "panelA_ptv_above_both_rna.tsv")
message("wrote panelA_ptv_above_both_rna.tsv")

SEP_Y <- 5.00e-7; ZERO_Y <- 3.40e-7; SEP2_Y <- 2.40e-7; NA_Y <- 1.70e-7
MEAS <- c("TOPMed RNA-outlier prevalence", "AoU RNA-outlier prevalence",
          "TOPMed PTV-carrier prevalence",
          "gnomAD sum AF_exomes (TOPMed-overlapping PTVs)")
PAL  <- setNames(c("#B2182B", "#2166AC", "#F0A202", "#1B7837"), MEAS)
SHP  <- setNames(c(16, 16, 17, 18), MEAS)
ACMG_COL <- "#5E3C99"

long <- g |>
  transmute(gene_nv, symbol, n_ptv, is_acmg, ptv_gt_both_rna,
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

ord_nz <- g |> arrange(desc(prev_ptv), desc(prev_topmed)) |> mutate(rank = row_number())
la_nz <- long |> left_join(dplyr::select(ord_nz, gene_nv, rank), by = "gene_nv")
N_GN_nz <- sum(!is.na(g$sum_AF_exomes))

focus_row <- g |> filter(symbol == FOCUS)
stopifnot(nrow(focus_row) == 1)
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
  labs(x = sprintf("gene, ranked by PTV-carrier prevalence (n = %d)", nrow(g)),
       y = "frequency",
       title = "A. Shared rare-outlier genes with \u22651 TOPMed PTV carrier",
       subtitle = sprintf(paste0("Purple marks %s (rank %d of %d), the only ACMG SF v3.2 gene in this panel ",
                                 "whose TOPMed PTV-carrier prevalence exceeds both RNA-outlier prevalences. ",
                                 "%d of %d genes meet that rule."),
                          FOCUS, focus_rank, nrow(g), nrow(hits), nrow(g))) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.4, color = "grey30"),
        legend.position = "top", legend.text = element_text(size = 8),
        legend.key.height = unit(10, "pt"),
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank())

MEAS_SHORT <- c("TOPMed RNA\noutlier", "AoU RNA\noutlier", "TOPMed PTV\ncarriers", "gnomAD\nsum AF_exomes")
rb <- long |>
  filter(symbol == FOCUS) |>
  mutate(meas_lab = factor(MEAS_SHORT[as.integer(measure)], levels = MEAS_SHORT),
         n_lab = case_when(
           measure == MEAS[1] ~ sprintf("%d / %s", focus_row$n_topmed, comma(N_TOP)),
           measure == MEAS[2] ~ sprintf("%d / %s", focus_row$n_new, comma(N_AOU)),
           measure == MEAS[3] ~ sprintf("%d / %s", focus_row$n_ptv, comma(N_TOP)),
           TRUE ~ "sum of AFs"
         ))

disease <- acmg_meta$`Disease/Phentyope`[acmg_meta$sym == toupper(FOCUS)]
if (length(disease) != 1) disease <- "Familial thoracic aortic aneurysm"

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
       subtitle = paste0(disease, ". TOPMed PTV-carrier prevalence exceeds both ",
                         "RNA-outlier prevalences.")) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10.5, color = ACMG_COL),
        plot.subtitle = element_text(size = 7.4, color = "grey30"),
        legend.position = "none",
        axis.text.x = element_text(size = 7.4, lineheight = 0.95),
        axis.text.y.right = element_text(size = 7, color = "grey45"),
        axis.title.y.right = element_blank())

p_side <- (pa_left | pb_focus) + plot_layout(widths = c(2.2, 1), guides = "collect") +
  plot_annotation(
    title = sprintf("MYH11: PTV-carrier prevalence above both RNA-outlier prevalences (\u2212log10 prevalence \u2265 %g in both cohorts)", CUT),
    subtitle = str_wrap(sprintf(paste0(
      "Left: the %d Panel A genes (\u22651 TOPMed PTV carrier; %d genes with 0 PTV carriers omitted). ",
      "Right: %s, the only ACMG SF v3.2 gene in this set with TOPMed PTV-carrier prevalence above both ",
      "RNA-outlier prevalences (%d of %d Panel A genes meet that rule). TNNI3 is also ACMG in Panel A ",
      "but PTV exceeds only the AoU RNA prevalence; RBM20 is the opposite pattern; ACTA2 is omitted ",
      "(0 PTV carriers). ACMG list: acmg_sf_v3.2.csv."),
      nrow(g), n_dropped, FOCUS, nrow(hits), nrow(g)), width = 145),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8, color = "grey30"),
                  legend.position = "top"))

W3 <- 13.2; H3 <- 5.8
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot_panelA_MYH11.png"), p_side, width = W3, height = H3, dpi = 200)
ggsave(file.path(OUT_DIR, "shared_rare_outlier_genes_4dot_panelA_MYH11.pdf"), p_side, width = W3, height = H3, device = cairo_pdf)
message("wrote shared_rare_outlier_genes_4dot_panelA_MYH11.{png,pdf}")
