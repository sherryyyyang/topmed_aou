#!/usr/bin/env Rscript
# Redraw enrichment_npc_all_cutoffs.pdf with (a) the lenient cutoffs z<-2 and
# z<-3 added and (b) the sweep extended below 55 PCs, so it is possible to see
# whether 55 is actually the optimum at z<-2.
#
# Three saved sweeps are combined:
#   expr_enrichments_low_npc.rds          cutoffs 2/3/3.5/4/5, npc 5..55  (step 5)
#   expr_enrichments_all_cutoffs_2_5.rds  cutoffs 2/3/3.5/4/5, npc 55..1055 (step 100)
#   expr_enrichments_all_cutoffs.rds      cutoffs 3.5/4/5,     npc 55..3955 (step 50)
# Where they overlap they agree exactly (checked below), so on the shared points
# the older saved sweeps are kept.
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- SUP_PC

load_sweep <- function(f) read_rds(f) |> filter(direction == "underlier") |>
  mutate(npcs_format = as.numeric(npcs_format), cutoff = as.numeric(cutoff))
wide   <- load_sweep("expr_enrichments_all_cutoffs.rds")
narrow <- load_sweep("expr_enrichments_all_cutoffs_2_5.rds")
low    <- load_sweep("expr_enrichments_low_npc.rds")

# --- sanity check: do the sweeps agree where they overlap? ---
check_overlap <- function(a, b, lab) {
  ov <- inner_join(dplyr::select(a, npcs_format, cutoff, or_a = odds_ratio),
                   dplyr::select(b, npcs_format, cutoff, or_b = odds_ratio),
                   by = c("npcs_format", "cutoff"))
  cat(sprintf("overlap %-22s %3d shared points, max |diff| = %.3g\n",
              lab, nrow(ov), if (nrow(ov)) max(abs(ov$or_a - ov$or_b)) else NA))
}
check_overlap(wide, narrow, "wide vs narrow:")
check_overlap(narrow, low, "narrow vs low-npc:")
check_overlap(wide, low, "wide vs low-npc:")

dat <- bind_rows(wide,
                 filter(narrow, cutoff %in% c(2, 3)),
                 anti_join(low, bind_rows(wide, narrow), by = c("npcs_format", "cutoff"))) |>
  mutate(cutoff_lab = factor(sprintf("z < -%s", cutoff),
                             levels = sprintf("z < -%s", sort(unique(cutoff)))))
stopifnot(!any(duplicated(dplyr::select(dat, npcs_format, cutoff))))

peaks <- dat |> group_by(cutoff_lab) |> slice_max(odds_ratio, n = 1) |> ungroup()
cat("\nnpc range per cutoff, and each curve's optimum:\n")
print(as.data.frame(dat |> group_by(cutoff_lab) |>
  summarise(npc_min = min(npcs_format), npc_max = max(npcs_format), n = n()) |>
  left_join(dplyr::select(peaks, cutoff_lab, peak_npc = npcs_format,
                          peak_or = odds_ratio), by = "cutoff_lab")))

z2 <- filter(dat, cutoff == 2) |> arrange(npcs_format)
cat("\nz<-2 curve (the cutoff that motivated 55 PCs):\n")
print(as.data.frame(dplyr::select(z2, npcs_format, odds_ratio, ci_95_lo, ci_95_hi)))
at55  <- filter(z2, npcs_format == 55)
best2 <- slice_max(z2, odds_ratio, n = 1)
cat(sprintf("\nz<-2: 55 PCs OR = %.4f [%.4f, %.4f]; best is %g PCs OR = %.4f -> %s\n",
            at55$odds_ratio, at55$ci_95_lo, at55$ci_95_hi,
            best2$npcs_format, best2$odds_ratio,
            if (best2$odds_ratio <= at55$ci_95_hi)
              "best value lies inside the 95% CI at 55 (indistinguishable)" else
              "best value lies outside the 95% CI at 55"))
cat(sprintf("z<-2: 55 PCs is rank %d of %d PC values tested\n",
            sum(z2$odds_ratio > at55$odds_ratio) + 1, nrow(z2)))

PAL <- c("z < -2" = "#7570B3", "z < -3" = "#E7298A", "z < -3.5" = "#D95F02",
         "z < -4" = "#1B9E77", "z < -5" = "#2171B5")

base_layers <- function(p) p +
  geom_hline(yintercept = 1.5, color = "red", linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(aes(ymin = ci_95_lo, ymax = ci_95_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.9) +
  scale_color_manual(values = PAL, name = "Z-score cutoff") +
  scale_fill_manual(values = PAL, guide = "none") +
  theme_classic(base_size = 11)

# Panel A: everything, log x so the 5-55 region is visible next to the 4000 tail
pa <- base_layers(ggplot(dat, aes(npcs_format, odds_ratio,
                                  color = cutoff_lab, fill = cutoff_lab))) +
  geom_vline(xintercept = c(55, 2005), linetype = "dashed", color = "grey35", linewidth = 0.45) +
  annotate("text", x = 55, y = Inf, label = "55 (used)", vjust = 1.4, hjust = 1.05,
           size = 2.8, color = "grey20", fontface = "bold") +
  annotate("text", x = 2005, y = Inf, label = "2005 (Gavish-Donoho)", vjust = 1.4,
           hjust = 1.05, size = 2.8, color = "grey20", fontface = "bold") +
  geom_point(data = peaks, shape = 21, size = 3.1, stroke = 1.1, fill = "white") +
  scale_x_log10(labels = comma, breaks = c(5, 10, 25, 55, 100, 250, 500, 1000, 2000, 4000)) +
  labs(x = "# Expression PCs (log scale)", y = "RV/underlier odds ratio",
       title = "A. Full sweep, 5 to 3,955 expression PCs",
       subtitle = "Open circles = each cutoff's optimum") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.3, color = "grey30"),
        legend.position = "right")

# Panel B: z<-2 only, linear x over the region that was searched originally
pb <- base_layers(ggplot(z2, aes(npcs_format, odds_ratio,
                                 color = cutoff_lab, fill = cutoff_lab))) +
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  geom_hline(yintercept = at55$ci_95_hi, linetype = "dotted", color = "grey45", linewidth = 0.4) +
  geom_point(data = best2, shape = 21, size = 3.4, stroke = 1.1, fill = "white") +
  annotate("text", x = 55, y = min(z2$ci_95_lo), label = " 55 PCs", hjust = 0,
           vjust = 0, size = 2.9, color = "grey20", fontface = "bold") +
  annotate("text", x = best2$npcs_format, y = best2$odds_ratio, label = "optimum ",
           hjust = 1, vjust = -1.1, size = 2.9, color = "grey20") +
  annotate("text", x = max(z2$npcs_format), y = at55$ci_95_hi,
           label = "upper 95% CI at 55 PCs ", hjust = 1, vjust = -0.6,
           size = 2.7, color = "grey45") +
  scale_x_continuous(labels = comma, breaks = c(5, 55, seq(155, 1055, 200))) +
  coord_cartesian(ylim = range(c(z2$ci_95_lo, z2$ci_95_hi))) +
  labs(x = "# Expression PCs", y = "RV/underlier odds ratio",
       title = "B. z < -2 only: is 55 PCs the maximum?",
       subtitle = str_wrap(paste0(
         "Enrichment rises steeply from 5 to ~55 PCs and then plateaus. The optimum sits at ",
         best2$npcs_format, " PCs, but every value from ~55 upward falls inside the 95% CI at ",
         "55 PCs, so 55 is the start of the plateau rather than a distinct peak."), width = 78)) +
  guides(color = "none") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8, color = "grey30"))

p <- (pa / pb) + plot_layout(heights = c(1, 0.95)) +
  plot_annotation(
    title = "Rare-variant / underlier enrichment vs number of expression PCs",
    subtitle = str_wrap(paste0(
      "The sweep now extends below 55 PCs. At the lenient z<-2 cutoff, 55 PCs is the knee of a ",
      "flat plateau, which is how it came to be used. At the strict cutoffs actually applied ",
      "downstream, enrichment keeps improving well past 55 and peaks in the low thousands ",
      "(z<-4: OR 1.95 at 5 PCs, 2.88 at 55, 5.29 at 2005, 5.45 at its 2205 optimum)."), width = 112),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.5, color = "grey30")))

ggsave(file.path(OUT_DIR, "enrichment_npc_all_cutoffs_with_z2_z3.pdf"), p, width = 8.4, height = 9.2, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "enrichment_npc_all_cutoffs_with_z2_z3.png"), p, width = 8.4, height = 9.2, dpi = 190)
write_tsv(arrange(dplyr::select(dat, cutoff, npcs_format, odds_ratio, ci_95_lo, ci_95_hi, p_value),
                  cutoff, npcs_format), "enrichment_npc_all_cutoffs_underlier.tsv")
message("wrote enrichment_npc_all_cutoffs_with_z2_z3.{pdf,png} and enrichment_npc_all_cutoffs_underlier.tsv")
