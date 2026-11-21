#!/usr/bin/env Rscript
# TEST only: Fig3b as boxplots of gene-level prevalences (PTV vs RNA underlier),
# overall and split by NMD status. Does not replace the main scatter panel.
#
# Inputs: pc1055 gene-level TSV (matches underliers_pc1055 / current Fig3b) +
#         nmd_escape_vep/fig3b_ptv_nmd_escape_by_variant.tsv
# Outputs:
#   figures/main/Figure3b_boxplot_test.png|pdf
#   figures/main/Figure3_boxplot_test.png|pdf  (optional composite with a/c)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
PKG  <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
# Prefer package helpers if present; fall back to ROOT.
.helpers <- if (file.exists(file.path(PKG, "scripts/helpers/paper_theme.R"))) {
  file.path(PKG, "scripts/helpers")
} else {
  file.path(ROOT, "scripts/helpers")
}
source(file.path(.helpers, "paper_theme.R"))

PC_SETTING <- Sys.getenv("PC_SETTING", unset = "1055")
GENE_TSV <- file.path(
  ROOT, paste0("pc", PC_SETTING, "_work"),
  "fig3_ptv_vs_rna/extra/ptv_vs_underlier_gene_level_z_-3.tsv"
)
if (!file.exists(GENE_TSV)) {
  stop("Missing gene-level TSV for PC_SETTING=", PC_SETTING, ":\n  ", GENE_TSV)
}
NMD_BY_VARIANT <- file.path(ROOT, "nmd_escape_vep", "fig3b_ptv_nmd_escape_by_variant.tsv")
NMD_AF_MAX <- 0.001
NMD_CLASS_LEVELS <- c("NMD-escaping", "NMD-triggering", "Unknown")
NMD_CLASS_COLORS <- c(
  "NMD-escaping"   = unname(OI["vermillion"]),
  "NMD-triggering" = unname(OI["blue"]),
  "Unknown"        = "grey55"
)
METRIC_LEVELS <- c("PTV carrier %", "RNA underlier %")
METRIC_COLORS <- c(
  "PTV carrier %"    = unname(OI["vermillion"]),
  "RNA underlier %"  = unname(OI["blue"])
)

OUT_DIR <- file.path(PKG, "figures", "main")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("PC_SETTING =", PC_SETTING, "\n")
cat("Loading", GENE_TSV, "\n")
genes <- read_tsv(GENE_TSV, show_col_types = FALSE)

cat("Loading NMD annotations:", NMD_BY_VARIANT, "\n")
nmd_map <- read_tsv(NMD_BY_VARIANT, show_col_types = FALSE) |>
  mutate(
    gene_nv = str_remove(as.character(gene), "\\..*$"),
    af = as.numeric(af),
    nmd_status = as.character(nmd_status)
  ) |>
  filter(!is.na(af), af < NMD_AF_MAX) |>
  group_by(gene_nv) |>
  summarise(
    n_escape = sum(nmd_status == "NMD_escaping", na.rm = TRUE),
    n_trigger = sum(nmd_status == "NMD_triggering", na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    nmd_class = case_when(
      n_escape > 0 ~ "NMD-escaping",
      n_trigger > 0 ~ "NMD-triggering",
      TRUE ~ "Unknown"
    ),
    nmd_class = factor(nmd_class, levels = NMD_CLASS_LEVELS)
  ) |>
  select(gene_nv, nmd_class)

genes <- genes |>
  left_join(nmd_map, by = "gene_nv") |>
  mutate(
    nmd_class = factor(
      ifelse(is.na(nmd_class), "Unknown", as.character(nmd_class)),
      levels = NMD_CLASS_LEVELS
    )
  )

long <- genes |>
  transmute(
    gene_nv, nmd_class,
    `PTV carrier %` = pct_ptv,
    `RNA underlier %` = pct_under
  ) |>
  pivot_longer(
    cols = all_of(METRIC_LEVELS),
    names_to = "metric",
    values_to = "pct"
  ) |>
  mutate(metric = factor(metric, levels = METRIC_LEVELS))

cat("Genes:", nrow(genes), "  long rows:", nrow(long), "\n")
print(as.data.frame(table(genes$nmd_class)))

# --- Panel 1: overall PTV vs underlier prevalence ---
p_overall <- ggplot(long, aes(x = metric, y = pct, fill = metric)) +
  geom_boxplot(
    outlier.size = 0.55, outlier.alpha = 0.35, outlier.stroke = 0,
    width = 0.55, linewidth = 0.35, show.legend = FALSE
  ) +
  scale_fill_manual(values = METRIC_COLORS) +
  scale_y_continuous(
    labels = label_percent(scale = 1, accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = NULL,
    y = "Gene-level prevalence\n(% of subjects)",
    title = "All genes (Fig. 3b set)"
  ) +
  theme_paper() +
  theme(
    plot.title = element_text(size = 10, face = "plain", hjust = 0),
    axis.text.x = element_text(size = 9)
  )

# --- Panel 2: same metrics split by NMD (escaping / triggering / unknown) ---
p_nmd <- ggplot(long, aes(x = nmd_class, y = pct, fill = metric)) +
  geom_boxplot(
    outlier.size = 0.45, outlier.alpha = 0.30, outlier.stroke = 0,
    width = 0.65, linewidth = 0.35,
    position = position_dodge(width = 0.72)
  ) +
  scale_fill_manual(values = METRIC_COLORS, name = NULL) +
  scale_y_continuous(
    labels = label_percent(scale = 1, accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = "Gene-level NMD status",
    y = "Gene-level prevalence\n(% of subjects)",
    title = "By NMD status"
  ) +
  theme_paper() +
  theme(
    plot.title = element_text(size = 10, face = "plain", hjust = 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 8.5, angle = 15, hjust = 1)
  )

# --- Panel 3 (optional clarity): underlier % only among escaping vs triggering ---
known <- long |> filter(nmd_class %in% c("NMD-escaping", "NMD-triggering"))
p_known <- ggplot(
  known |> filter(metric == "RNA underlier %"),
  aes(x = nmd_class, y = pct, fill = nmd_class)
) +
  geom_boxplot(
    outlier.size = 0.55, outlier.alpha = 0.35, outlier.stroke = 0,
    width = 0.55, linewidth = 0.35, show.legend = FALSE
  ) +
  scale_fill_manual(values = NMD_CLASS_COLORS) +
  scale_y_continuous(
    labels = label_percent(scale = 1, accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = NULL,
    y = "RNA underlier prevalence\n(% of subjects)",
    title = "Underlier % (known NMD only)"
  ) +
  theme_paper() +
  theme(
    plot.title = element_text(size = 10, face = "plain", hjust = 0),
    axis.text.x = element_text(size = 8.5, angle = 15, hjust = 1)
  )

fig <- (p_overall | p_nmd | p_known) +
  plot_annotation(
    title = paste0(
      "Fig. 3b test — gene-level prevalence boxplots (PC_SETTING=",
      PC_SETTING, ", z < -3)"
    ),
    theme = theme(
      plot.title = element_text(family = PAPER_FONT, size = 11, face = "bold")
    )
  )

stem <- file.path(OUT_DIR, "Figure3b_boxplot_test")
ggsave(paste0(stem, ".png"), fig, width = NATURE_DOUBLE_IN, height = 3.6, dpi = 300)
ggsave(paste0(stem, ".pdf"), fig, width = NATURE_DOUBLE_IN, height = 3.6, device = cairo_pdf)
cat("Wrote:", paste0(stem, ".png"), "\n")

# --- Optional: Figure3 composite with boxplot replacing panel b ---
# Built separately with Python/PIL when magick is unavailable (see script run notes).
# Expected path: figures/main/Figure3_boxplot_test.png
cat("Boxplot panel ready; composite assembled externally if needed.\n")

# Median summary for the log
cat("\nMedian prevalences (%):\n")
print(
  long |>
    group_by(nmd_class, metric) |>
    summarise(n = n(), median = median(pct), q75 = quantile(pct, 0.75), .groups = "drop")
)
