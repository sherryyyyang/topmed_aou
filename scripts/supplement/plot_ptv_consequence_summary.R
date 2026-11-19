#!/usr/bin/env Rscript
# Bar chart of VEP consequence strings (column "consequence" in ptv_all.tsv).
# Each distinct string is its own category (including compound "&" annotations).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- SUP_PTV


PTV_TSV <- Sys.getenv(
  "PTV_TSV",
  unset = "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/ptv_all.tsv"
)
USE_KEPT_ONLY <- Sys.getenv("PTV_CONSEQUENCE_USE_KEPT", unset = "1") == "1"

cat("Loading:", PTV_TSV, "\n")
ptv <- read_tsv(PTV_TSV, show_col_types = FALSE)

if (USE_KEPT_ONLY) {
  kept <- read_rds(KEPT_SUBJECTS)
  ptv <- ptv |> filter(sample %in% kept)
  cohort_lbl <- sprintf("QC cohort (n subjects in PTV table = %s)", comma(n_distinct(ptv$sample)))
} else {
  cohort_lbl <- sprintf("All rows in PTV file (n samples = %s)", comma(n_distinct(ptv$sample)))
}

if (!"consequence" %in% names(ptv)) {
  stop("Expected column 'consequence' in PTV table")
}

# Truncating to a fixed width collapsed distinct annotations onto a shared
# label (e.g. four splice_acceptor + coding_sequence... compounds), which
# silently stacked their bars and stacked their count labels on top of each
# other. Drop the uninformative "_variant" suffix that every SO term carries
# and wrap instead, which keeps all 31 annotations on their own bar.
summ <- ptv |>
  count(consequence, name = "n_variants", sort = TRUE) |>
  mutate(
    lab_flat = str_replace_all(consequence, "&", " + "),
    lab_flat = str_replace_all(lab_flat, "_variant\\b", ""),
    lab = factor(str_wrap(lab_flat, width = 40),
                 levels = rev(unique(str_wrap(lab_flat, width = 40))))
  )
stopifnot(!any(duplicated(levels(summ$lab))), nlevels(summ$lab) == nrow(summ))

out_tsv <- if (USE_KEPT_ONLY) {
  "ptv_consequence_counts_qc_cohort.tsv"
} else {
  "ptv_consequence_counts_all.tsv"
}
write_tsv(
  summ |> transmute(consequence, display_label = lab_flat, n_variants) |>
    arrange(desc(n_variants)),
  file.path(OUT_DIR, out_tsv)
)
if (!is_paper_assemble()) cat("Wrote:", file.path(OUT_DIR, out_tsv), "\n")

# Log10 count axis so rare compound annotations remain visible next to
# frameshift / stop_gained. Map counts to log10 space first: geom_col baselines
# at 0, which is invalid under scale_x_log10(count) but equals count=1 here.
stopifnot(all(summ$n_variants >= 1))
x_breaks <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000)
p <- ggplot(summ, aes(x = log10(n_variants), y = lab)) +
  geom_col(fill = COL_TOPMED, width = 0.85) +
  geom_text(aes(label = comma(n_variants)), hjust = -0.12, size = 2.3,
            family = PAPER_FONT, color = "grey20") +
  scale_x_continuous(
    breaks = log10(x_breaks),
    labels = comma(x_breaks),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(x = "Variant rows (count, log scale)", y = NULL) +
  theme_paper() +
  theme(axis.text.y = element_text(size = 7, lineheight = 0.85))

.paper_plots <- list(c = p)
out_png <- if (USE_KEPT_ONLY) "ptv_consequence_summary_qc_cohort" else "ptv_consequence_summary_all"
paper_ggsave(p, file.path(OUT_DIR, out_png), width = 7.2, height = 8.4)
if (!is_paper_assemble()) cat("Wrote:", out_png, "\n")
