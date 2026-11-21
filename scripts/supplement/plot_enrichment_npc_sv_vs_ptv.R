#!/usr/bin/env Rscript
# Enrichment of underliers for rare variants vs number of expression PCs, with
# the variant side defined two ways:
#   rare SV  -- SV with AF <= 1% within 100 kb of the gene (the original analysis)
#   LoF/PTV  -- LOFTEE high-confidence PTV inside the gene, AF <= 0.5% (ptv_all.tsv)
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork); library(cowplot)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- SUP_PC

COLS <- c("cutoff", "npcs_format", "odds_ratio", "ci_95_lo", "ci_95_hi", "p_value")
load_sweep <- function(f) read_rds(f) |> filter(direction == "underlier") |>
  mutate(npcs_format = as.numeric(npcs_format), cutoff = as.numeric(cutoff)) |>
  dplyr::select(all_of(COLS))

# expr_enrichments_sv_full.rds covers all five cutoffs over the full 5..3955 grid
# and reproduces the three older saved sweeps exactly on their 618 shared points,
# so it supersedes them here (the older files stop at 1055 for z<-2 and z<-3).
sv <- load_sweep(file.path(DATA_ENRICH, "expr_enrichments_sv_full.rds")) |>
  mutate(variant_class = "Rare SV within 100 kb (AF <= 1%)")
ptv <- load_sweep(file.path(DATA_ENRICH, "expr_enrichments_ptv.rds")) |>
  mutate(variant_class = "LoF / PTV in gene (AF <= 0.5%)")

dat <- bind_rows(sv, ptv) |>
  mutate(cutoff_lab = factor(sprintf("z < -%s", cutoff),
                             levels = sprintf("z < -%s", sort(unique(cutoff)))),
         variant_class = factor(variant_class, levels = c(unique(sv$variant_class),
                                                          unique(ptv$variant_class))))
stopifnot(!any(duplicated(dplyr::select(dat, variant_class, npcs_format, cutoff))))

peaks <- dat |> group_by(variant_class, cutoff_lab) |> slice_max(odds_ratio, n = 1) |> ungroup()
summ <- dat |> group_by(variant_class, cutoff_lab) |>
  summarise(npc_min = min(npcs_format), npc_max = max(npcs_format), n_npc = n(),
            or_55 = odds_ratio[match(55, npcs_format)],
            or_2005 = odds_ratio[match(2005, npcs_format)], .groups = "drop") |>
  left_join(dplyr::select(peaks, variant_class, cutoff_lab,
                          peak_npc = npcs_format, peak_or = odds_ratio),
            by = c("variant_class", "cutoff_lab"))
cat("\nsummary by variant class and cutoff:\n"); print(as.data.frame(summ))

PAL <- CUTOFF_COLORS
XBREAKS <- c(0, 1000, 2000, 3000, 4000)

panel <- function(vc) {
  d <- filter(dat, variant_class == vc); pk <- filter(peaks, variant_class == vc)
  short <- if (grepl("Rare SV", vc)) "Rare SVs" else "Rare PTVs"
  ggplot(d, aes(npcs_format, odds_ratio, color = cutoff_lab, fill = cutoff_lab)) +
    geom_hline(yintercept = 1, color = "grey60", linewidth = 0.35) +
    # Vertical guide: 1055 = z<-3 rare-PTV enrichment peak (paper PC).
    geom_vline(xintercept = 1055, linetype = "dashed", color = "grey35", linewidth = 0.45) +
    annotate("text", x = 1055, y = Inf, label = "1055 ", vjust = 1.4, hjust = 1,
             size = 2.6, color = "grey20", fontface = "bold", family = PAPER_FONT) +
    geom_ribbon(aes(ymin = ci_95_lo, ymax = ci_95_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
    geom_point(data = pk, shape = 21, size = 3.4, stroke = 1.1, fill = "white") +
    scale_color_manual(values = PAL, name = "Z-score cutoff") +
    scale_fill_manual(values = PAL, guide = "none") +
    guides(color = guide_legend(
      nrow = 1,
      override.aes = list(size = 1.5, stroke = 0.6, linewidth = 0.5))) +
    scale_x_continuous(labels = comma, breaks = XBREAKS,
                       expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_log10(breaks = c(1, 1.5, 2, 3, 5, 10, 20, 50, 100, 200),
                  labels = c("1", "1.5", "2", "3", "5", "10", "20", "50", "100", "200")) +
    labs(x = "# Expression PCs", y = "Odds ratio",
         title = short) +
    theme_paper() +
    theme_paper_inner_title() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.key.size = unit(7, "pt"),
      legend.key.height = unit(7, "pt"),
      legend.key.width = unit(9, "pt"),
      legend.spacing.x = unit(2, "pt"),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(3, "pt")
    )
}

pa <- panel(levels(dat$variant_class)[1])
pb <- panel(levels(dat$variant_class)[2])
.paper_plots <- list(a = pa, b = pb)

# Side-by-side SV | PTV with shared Z-score legend / caption under the panels.
p <- (pa | pb) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    caption = paste(
      "Enrichment of expression underliers for rare variants as a function of the number of",
      "residualized expression PCs. Left, rare SVs within 100 kb (AF <= 1%).",
      "Right, rare LoF/PTVs in the gene (AF <= 0.5%).",
      "Curves: Fisher odds ratio by Z-score cutoff; ribbons: 95% CI;",
      "open circles: per-cutoff peak; dashed line: 1055 PCs (z<-3 rare-PTV enrichment peak)."
    ),
    theme = theme(
      plot.caption = element_text(
        size = 7.5, hjust = 0, color = "grey20",
        family = PAPER_FONT, margin = margin(t = 6)
      ),
      plot.tag = element_text(face = "bold", size = 14, family = PAPER_FONT)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center"
  )
paper_ggsave(p, file.path(OUT_DIR, "enrichment_npc_sv_vs_ptv"),
             width = NATURE_DOUBLE_IN, height = 4.6)
if (!is_paper_assemble()) {
  # These panels are the no-overlay rendering of S1a/S1b: same sweep, without
  # the cell-type composition line and its right-hand axis. They used to be
  # one-off archived backups, so they drifted from the styling of the overlay
  # versions; regenerating them here keeps the two variants in step.
  S1_NO_OVERLAY_DIR <- file.path(FINAL_DRAFT, "supp_figs")
  dir.create(S1_NO_OVERLAY_DIR, recursive = TRUE, showWarnings = FALSE)
  paper_ggsave(pa, file.path(S1_NO_OVERLAY_DIR, "FigS1a_sv_enrichment_no_overlay"),
               width = NATURE_DOUBLE_IN / 2 + 0.35, height = 4.2)
  paper_ggsave(pb, file.path(S1_NO_OVERLAY_DIR, "FigS1b_ptv_enrichment_no_overlay"),
               width = NATURE_DOUBLE_IN / 2 + 0.35, height = 4.2)

  # Matching composite: a | b over the full-width consequence panel, with one
  # shared legend under the two sweeps.
  s1c_env <- new.env(parent = globalenv())
  options(paper.assemble = TRUE)
  source(file.path(ROOT, "scripts/supplement/plot_ptv_consequence_summary.R"),
         local = s1c_env)
  options(paper.assemble = FALSE)
  s1_panels <- cowplot::plot_grid(
    pa + theme(legend.position = "none"),
    pb + theme(legend.position = "none"),
    ncol = 2, labels = c("a", "b"),
    label_fontfamily = PAPER_FONT, label_size = 14, label_fontface = "bold",
    align = "hv"
  )
  s1_leg <- cowplot::get_legend(pb)
  s1_top <- cowplot::plot_grid(s1_panels, s1_leg, ncol = 1,
                               rel_heights = c(1, 0.07))
  sfig1_no <- cowplot::plot_grid(
    s1_top, s1c_env$.paper_plots$c, ncol = 1, rel_heights = c(1.00, 1.60),
    labels = c("", "c"),
    label_fontfamily = PAPER_FONT, label_size = 14, label_fontface = "bold",
    align = "none"
  )
  paper_ggsave(sfig1_no, file.path(S1_NO_OVERLAY_DIR, "Supplementary_Figure1_no_overlay"),
               width = NATURE_DOUBLE_IN, height = 10.4)

  write_tsv(arrange(dplyr::select(dat, variant_class, cutoff, npcs_format,
                                  odds_ratio, ci_95_lo, ci_95_hi, p_value),
                    variant_class, cutoff, npcs_format),
            file.path(OUT_DIR, "enrichment_npc_sv_vs_ptv_underlier.tsv"))
  message("wrote enrichment_npc_sv_vs_ptv")
}
