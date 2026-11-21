#!/usr/bin/env Rscript
# Paper Figure 4: RBM20 and MYH11 4-dot panels without Wilson CIs on panel B.
# Reads the already-built shared_rare_outlier_genes_4dot.tsv.
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork); library(ggrepel)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG4A

PW <- ROOT
setwd(PW)
OUT_DIR <- OUT_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ACMG_CSV <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/annotation/variant_annotations/acmg_sf_v3.2.csv"
CUT <- 3
SEP_Y <- 5.00e-7; ZERO_Y <- 3.40e-7; SEP2_Y <- 2.40e-7; NA_Y <- 1.70e-7
MEAS <- c("TOPMed RNA-outlier prevalence", "AoU RNA-outlier prevalence",
          "TOPMed PTV-carrier prevalence",
          "gnomAD sum AF_exomes (TOPMed-overlapping PTVs)")
# 4a/4c: trend lines (not 340×4 overlapping dots), four Okabe–Ito colours.
COL_RNA_TOP <- COL_TOPMED             # #0072B2
COL_RNA_AOU <- COL_AOU                # #CC79A7
COL_GNOMAD  <- unname(OI["green"])    # #009E73
PAL  <- setNames(c(COL_RNA_TOP, COL_RNA_AOU, COL_PTV, COL_GNOMAD), MEAS)
SHP  <- setNames(c(16, 16, 17, 18), MEAS)
LTY  <- setNames(c("solid", "solid", "solid", "solid"), MEAS)
ACMG_COL <- COL_ACMG
MEAS_SHORT <- c("TOPMed RNA", "AoU RNA", "TOPMed PTV", "gnomAD AF")
MEAS_LEGEND <- MEAS_SHORT
Y_LIM <- c(5e-7, 0.04)
Y_LIM_RANKS <- c(1e-4, 0.01)
Y_BREAKS <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
Y_BREAKS_RANKS <- c(1e-4, 1e-3, 1e-2)
y_scale <- function(breaks = Y_BREAKS, expand = expansion(mult = c(0.02, 0.08))) list(
  scale_y_log10(
    breaks = breaks,
    labels = function(v) paste0("1 in ", comma(round(1 / v))),
    expand = expand
  )
)
fmt_af <- function(x) {
  if (is.na(x) || !is.finite(x) || x <= 0) return("NA")
  e <- floor(log10(x))
  m <- x / (10^e)
  super <- c("\u2070", "\u00b9", "\u00b2", "\u00b3", "\u2074",
             "\u2075", "\u2076", "\u2077", "\u2078", "\u2079")
  exp_digits <- strsplit(as.character(abs(as.integer(e))), "")[[1]]
  exp_str <- paste0(
    if (e < 0) "\u207b" else "",
    paste0(super[as.integer(exp_digits) + 1L], collapse = "")
  )
  sprintf("%.2f\u00d710%s", m, exp_str)
}

g_all <- read_tsv(file.path(FIG4D, "shared_rare_outlier_genes_4dot.tsv"), show_col_types = FALSE)
N_TOP <- unique(g_all$N_topmed); stopifnot(length(N_TOP) == 1)
N_AOU <- unique(g_all$N_new);    stopifnot(length(N_AOU) == 1)
acmg <- read_csv(ACMG_CSV, show_col_types = FALSE)
names(acmg) <- str_trim(names(acmg))
acmg_sym <- unique(toupper(str_trim(acmg$Gene)))
acmg_meta <- acmg |> mutate(sym = toupper(str_trim(Gene))) |> distinct(sym, .keep_all = TRUE)

make_pair <- function(focus) {
  g <- g_all |> filter(n_ptv > 0) |>
    mutate(is_acmg = toupper(symbol) %in% acmg_sym)
  n_dropped <- nrow(g_all) - nrow(g)
  long <- g |>
    transmute(gene_nv, symbol, n_ptv, is_acmg,
              `TOPMed RNA-outlier prevalence` = prev_topmed,
              `AoU RNA-outlier prevalence`    = prev_aou,
              `TOPMed PTV-carrier prevalence` = prev_ptv,
              `gnomAD sum AF_exomes (TOPMed-overlapping PTVs)` = sum_AF_exomes) |>
    pivot_longer(all_of(MEAS), names_to = "measure", values_to = "prevalence") |>
    mutate(measure = factor(measure, levels = MEAS),
           no_data = is.na(prevalence),
           is_zero = !no_data & prevalence == 0,
           y_plot  = case_when(no_data ~ NA_Y, is_zero ~ ZERO_Y, TRUE ~ prevalence))
  ord_nz <- g |> arrange(desc(prev_ptv), desc(prev_topmed)) |> mutate(rank = row_number())
  la_nz <- long |> left_join(dplyr::select(ord_nz, gene_nv, rank), by = "gene_nv")
  focus_row <- g |> filter(symbol == focus)
  stopifnot(nrow(focus_row) == 1)
  focus_rank <- ord_nz$rank[ord_nz$symbol == focus]
  # Ranks panel: RNA vs PTV only. gnomAD stays on the 4-dot strip.
  MEAS_RANKS <- MEAS[1:3]
  la_focus <- filter(la_nz, symbol == focus, !no_data, measure %in% MEAS_RANKS)
  lab_focus <- la_focus |> group_by(rank, symbol) |> summarise(y = max(y_plot), .groups = "drop")
  la_rna_top <- filter(la_nz, !no_data, measure == MEAS[1])
  la_rna_aou <- filter(la_nz, !no_data, measure == MEAS[2])
  la_ptv     <- filter(la_nz, !no_data, measure == MEAS[3])

  rna_pts <- bind_rows(la_rna_top, la_rna_aou) |> filter(!is_zero)
  ptv_line <- la_ptv |> filter(!is_zero) |> arrange(rank)

  # Same-cohort subtraction: TOPMed RNA − TOPMed PTV.
  y_rna <- focus_row$prev_topmed
  y_ptv <- focus_row$prev_ptv
  delta_lab <- sprintf("%.3f%%", 100 * (y_rna - y_ptv))
  max_rank <- max(la_nz$rank, na.rm = TRUE)
  # Park the bracket in the right margin, past every trend line.
  delta_x <- max_rank + 22
  tick <- 7
  x_lim_ranks <- c(0.5, max_rank + 55)
  y_hi <- max(y_rna, y_ptv)
  y_lo <- min(y_rna, y_ptv)
  y_mid <- exp(mean(log(c(y_hi, y_lo))))
  show_delta <- TRUE  # RNA−PTV bracket on both ranks panels (4a and 4c)
  # Gene-name callout pointing at the focus dots (not the right-margin bridge).
  y_name <- min(
    Y_LIM_RANKS[2] * 0.45,
    10^(log10(max(la_focus$y_plot)) + 0.55)
  )
  x_name <- if (focus_rank > max(la_nz$rank, na.rm = TRUE) * 0.55) {
    focus_rank - 55
  } else {
    focus_rank + 55
  }
  hjust_name <- if (focus_rank > max(la_nz$rank, na.rm = TRUE) * 0.55) 1 else 0

  pa_left <- ggplot(
    filter(la_nz, !no_data, measure %in% MEAS_RANKS),
    aes(rank, y_plot, color = measure, shape = measure)
  ) +
    geom_vline(xintercept = focus_rank, color = ACMG_COL, alpha = 0.45, linewidth = 0.45) +
    geom_point(alpha = 0, size = 0) +
    geom_point(data = rna_pts, size = 0.55, alpha = 0.10, shape = 16, show.legend = FALSE) +
    geom_line(
      data = ptv_line, aes(rank, y_plot, color = measure),
      linewidth = 1.05, linetype = "solid", inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_smooth(
      data = filter(la_rna_top, !is_zero), aes(rank, y_plot, color = measure),
      method = "loess", span = 0.35, se = FALSE, linewidth = 1.05, linetype = "solid",
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_smooth(
      data = filter(la_rna_aou, !is_zero), aes(rank, y_plot, color = measure),
      method = "loess", span = 0.35, se = FALSE, linewidth = 1.05, linetype = "solid",
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    {
      if (show_delta) {
        list(
          geom_segment(
            data = tibble(x = delta_x, xend = delta_x, y = y_hi, yend = y_lo),
            aes(x = x, y = y, xend = xend, yend = yend),
            inherit.aes = FALSE, color = "grey25", linewidth = 0.55,
            arrow = arrow(length = unit(4.2, "pt"), ends = "both",
                          type = "closed", angle = 22)
          ),
          geom_segment(
            data = tibble(
              x = c(delta_x, delta_x),
              xend = c(delta_x - tick, delta_x - tick),
              y = c(y_hi, y_lo),
              yend = c(y_hi, y_lo)
            ),
            aes(x = x, y = y, xend = xend, yend = yend),
            inherit.aes = FALSE, color = "grey25", linewidth = 0.55
          ),
          geom_text(
            data = tibble(x = delta_x + 8, y = y_mid),
            aes(x = x, y = y, label = delta_lab),
            inherit.aes = FALSE, color = "grey20", fontface = "bold",
            size = 3.1, family = PAPER_FONT, hjust = 0, vjust = 0.5
          )
        )
      } else {
        list()
      }
    } +
    geom_point(data = la_focus, size = 2.6, alpha = 1) +
    geom_segment(
      data = tibble(
        x = x_name, y = y_name,
        xend = focus_rank, yend = max(la_focus$y_plot)
      ),
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, color = "grey15", linewidth = 0.45,
      arrow = arrow(length = unit(6, "pt"), type = "closed", angle = 22)
    ) +
    geom_text(
      data = tibble(x = x_name, y = y_name, hjust = hjust_name),
      aes(x = x, y = y, label = focus, hjust = hjust),
      inherit.aes = FALSE, color = ACMG_COL, fontface = "bold",
      size = 3.2, family = PAPER_FONT, vjust = 0.35
    ) +
    # Zoom cue: vertical focus line + gene callout at the dots; bridge arrow at right
    # has no gene-name text (names only as strip titles on b/d).
    scale_color_manual(
      values = PAL[MEAS_RANKS], name = NULL,
      labels = MEAS_SHORT[seq_along(MEAS_RANKS)],
      breaks = MEAS_RANKS, drop = TRUE
    ) +
    scale_shape_manual(
      values = SHP[MEAS_RANKS], name = NULL,
      labels = MEAS_SHORT[seq_along(MEAS_RANKS)],
      breaks = MEAS_RANKS, drop = TRUE
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        override.aes = list(
          size = 2.4, alpha = 1,
          shape = unname(SHP[MEAS_RANKS]),
          linetype = unname(LTY[MEAS_RANKS]),
          linewidth = 0.9
        )
      ),
      shape = "none"
    ) +
    y_scale(Y_BREAKS_RANKS, expand = expansion(mult = c(0.04, 0.05))) +
    labs(x = "Gene, ranked by PTV-carrier prevalence",
         y = "Frequency") +
    coord_cartesian(xlim = x_lim_ranks, ylim = Y_LIM_RANKS, clip = "off") +
    theme_paper() +
    theme(legend.position = "top",
          legend.key.height = unit(9, "pt"),
          legend.text = element_text(size = 8, family = PAPER_FONT),
          plot.margin = margin(8, 4, 4, 6))

  # Compact strip. Drop gnomAD so example strips match the ranks panels (no green AF point).
  drop_gnomad <- TRUE
  x_labs <- c("TOPMed", "AoU", "PTV", if (!drop_gnomad) "gnomAD")
  y_vals <- c(
    focus_row$prev_topmed,
    focus_row$prev_aou,
    focus_row$prev_ptv,
    if (!drop_gnomad) focus_row$sum_AF_exomes
  )
  n_labs <- c(
    sprintf("%s/%s", comma(focus_row$n_topmed), comma(N_TOP)),
    sprintf("%s/%s", comma(focus_row$n_new), comma(N_AOU)),
    sprintf("%s/%s", comma(focus_row$n_ptv), comma(N_TOP)),
    if (!drop_gnomad) fmt_af(focus_row$sum_AF_exomes)
  )
  if (drop_gnomad) {
    y_lim_dot <- Y_LIM_RANKS
    y_breaks_dot <- Y_BREAKS_RANKS
    y_floor <- Y_LIM_RANKS[1]
    # One shared row above the highest point so AoU is not pulled down with its dot.
    y_lab_vals <- rep(10^(log10(max(y_vals)) + 0.72), length(y_vals))
  } else {
    y_dex <- c(-0.52, -0.52, 0.52, -0.52)
    y_lim_dot <- Y_LIM
    y_breaks_dot <- Y_BREAKS
    y_floor <- Y_LIM[1]
    y_lab_vals <- 10^(log10(pmax(y_vals, y_floor)) + y_dex)
  }
  rb4 <- tibble(
    x = seq_along(x_labs),
    y_plot = y_vals,
    n_lab = n_labs,
    col = c(COL_RNA_TOP, COL_RNA_AOU, COL_PTV, if (!drop_gnomad) COL_GNOMAD),
    shp = c(16, 16, 17, if (!drop_gnomad) 18),
    y_lab = y_lab_vals
  )

  pb_focus <- ggplot(rb4, aes(x, y_plot)) +
    geom_line(aes(group = 1), color = ACMG_COL, linewidth = 0.7, alpha = 0.7) +
    geom_point(aes(colour = col, shape = shp), size = 2.8, stroke = 0.9) +
    geom_text(
      aes(y = y_lab, label = n_lab),
      color = "grey20", size = 2.15, family = PAPER_FONT, show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_shape_identity() +
    scale_x_continuous(
      breaks = seq_along(x_labs),
      labels = x_labs,
      expand = expansion(mult = c(0.22, 0.22))
    ) +
    y_scale(y_breaks_dot) +
    coord_cartesian(ylim = y_lim_dot, clip = "off") +
    labs(x = NULL, y = "Frequency", title = focus) +
    theme_paper() +
    theme_paper_inner_title() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 5.5, angle = 40, hjust = 1,
                                 margin = margin(t = 1)),
      axis.text.y = element_text(size = 6.5),
      axis.title.y = element_text(size = 9),
      plot.margin = margin(8, 6, 2, 2)
    )

  # Gutter bridge: right-pointing zoom arrow only (gene name is on strip b/d).
  bridge <- ggplot() +
    annotate(
      "segment", x = 0.12, xend = 0.92, y = 0.50, yend = 0.50,
      colour = "grey15", linewidth = 0.85,
      arrow = arrow(length = unit(8, "pt"), type = "closed", angle = 22)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_family = PAPER_FONT) +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.tag = element_blank()
    )

  list(
    left = pa_left, right = pb_focus, bridge = bridge,
    n_genes = nrow(g), n_dropped = n_dropped
  )
}

pair_r <- make_pair("RBM20")
pair_m <- make_pair("MYH11")
.paper_plots <- list(
  a = pair_r$left, b = pair_r$right,
  c = pair_m$left, d = pair_m$right,
  bridge_rbm20 = pair_r$bridge, bridge_myh11 = pair_m$bridge
)
cat(sprintf("Fig 4 caption: RBM20/MYH11 ranks n_genes=%d (dropped zero-PTV=%d)\n",
            pair_r$n_genes, pair_r$n_dropped))

pair_layout <- function(left, bridge, right) {
  (left | bridge | right) +
    plot_layout(widths = c(3.45, 0.38, 1.45))
}
p_r <- pair_layout(pair_r$left, pair_r$bridge, pair_r$right)
p_m <- pair_layout(pair_m$left, pair_m$bridge, pair_m$right)
paper_ggsave(p_r, file.path(FIG4A, "fig4_RBM20_4dot_no_ci"),
             width = NATURE_DOUBLE_IN, height = 3.8)
paper_ggsave(p_m, file.path(FIG4B, "fig4_MYH11_4dot_no_ci"),
             width = NATURE_DOUBLE_IN, height = 3.8)
