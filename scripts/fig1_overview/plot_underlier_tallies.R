#!/usr/bin/env Rscript
# Ranked-distribution plots of underlier tallies.
#
# Reproduces lines 329, 331, 334 of outlier_processing.R:
#   z_score_underlier  <- group_by(expr_z_join, gene_id) |> group_modify(get_underlier)
#   z_score_tally      <- group_by(z_score_underlier, gene_id)   |> tally() |> arrange(desc(n))
#       z_score_tally$pct_underlier <- z_score_tally$n / ncol(expr_cpm)
#   z_score_tally_subj <- group_by(z_score_underlier, subjectid) |> tally() |> arrange(desc(n))
#
# Underlier filter (lines 318-322):
#   per gene g: subject s is an underlier if
#       z-score(g, s) < -4   AND   logCPM_adj(g, s) < mean_g(logCPM_adj) - 1
#
# Inputs (in working dir):
#   z_score_underlier_clean.rds   -- underliers after subject QC (same definition as
#                                    z_score_underlier.rds; rows from bad subjects removed)
#   kept_subjects.rds             -- QC-pass subject ids (from apply_extreme_subject_filter.R)
#   expr_cpm_adj.rds              -- matrix; nrow() = gene universe; ncol() unused for n_subjects
# Output:
#   underlier_tallies.png / .pdf  -- two-panel ranked distribution figure
#   z_score_tally.rds / z_score_tally_subj.rds  -- tallies on the QC-filtered cohort

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(scales)
    library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG1B

source(file.path(ROOT, "scripts/helpers/_underlier_def.R"))
DEF <- resolve_underlier_def()
z_cutoff <- ifelse(is.null(DEF$z_cutoff), -4, DEF$z_cutoff)
logcpm_off <- 1L
cat("Underlier definition: ", DEF$short, " ", DEF$tag, "\n", sep = "")

cat("Loading kept_subjects.rds ...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_subjects <- length(kept_subjects)
cat("  kept subjects (QC) =", n_subjects, "\n")

cat("Loading ", DEF$rds_underlier, " ...\n", sep = "")
t0 <- Sys.time()
z_score_underlier <- read_rds(DEF$rds_underlier)
cat("  ", nrow(z_score_underlier), "underlier rows (QC-filtered) in",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")

cat("Loading expr_cpm_adj.rds for gene count ...\n")
t0 <- Sys.time()
expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
n_genes <- nrow(expr_cpm_adj)
rm(expr_cpm_adj); gc(verbose = FALSE)
cat("  eligible genes =", n_genes,
    "  (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s)\n")

# ---- gene-level tally (lines 331-332) ----
z_score_tally <- z_score_underlier |>
    group_by(gene_id) |>
    tally() |>
    arrange(desc(n)) |>
    mutate(pct_underlier = n / n_subjects,
           rank          = row_number())

# ---- subject-level tally (line 334) ----
z_score_tally_subj <- z_score_underlier |>
    group_by(subjectid) |>
    tally() |>
    arrange(desc(n)) |>
    mutate(pct_genes = n / n_genes,
           rank      = row_number())

cat("Genes with >=1 underlier:    ", nrow(z_score_tally),     "/", n_genes,    "\n")
cat("Subjects with >=1 underlier: ", nrow(z_score_tally_subj), "/", n_subjects, "\n")
cat("Top genes by pct_underlier:\n");           print(head(z_score_tally, 5))
cat("Top subjects by # underlier genes:\n");    print(head(z_score_tally_subj, 5))

out_tally_g <- ifelse(DEF$label == "strict", "z_score_tally.rds",
                      paste0("z_score_tally", DEF$suffix, ".rds"))
out_tally_s <- ifelse(DEF$label == "strict", "z_score_tally_subj.rds",
                      paste0("z_score_tally_subj", DEF$suffix, ".rds"))
write_rds(z_score_tally, file.path(DATA_UNDERLIER, out_tally_g))
write_rds(z_score_tally_subj, file.path(DATA_UNDERLIER, out_tally_s))

# ---- annotation strings ----
filter_subtitle <- sprintf(
    paste0(
        "Underlier per (gene, subject):  z-score < %g  AND  log\u2082CPM_adj < gene_mean \u2212 %g\n",
        "Subject QC: dropped 71/5,560 with |z|>5 count \u2265 60; analysis on %s kept subjects."
    ),
    z_cutoff, logcpm_off, comma(n_subjects)
)
caption_txt <- sprintf(
    "%s underlier (gene, subject) pairs across %s genes and %s kept subjects (QC-filtered cohort).",
    comma(nrow(z_score_underlier)), comma(n_genes), comma(n_subjects)
)

# ---- helper: log10 ticks that include 1 ----
log10_breaks <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)

# ---- Panel A: genes ranked by pct_underlier ----
top_genes <- head(z_score_tally, 5)
med_pct_g <- median(z_score_tally$pct_underlier)

p_gene <- ggplot(z_score_tally, aes(x = rank, y = pct_underlier)) +
    geom_line(color = "#2C7FB8", linewidth = 0.4) +
    geom_hline(yintercept = med_pct_g, linetype = "dashed",
               color = "firebrick", linewidth = 0.4) +
    annotate("text",
             x = nrow(z_score_tally), y = med_pct_g,
             label = sprintf(" median = %s", percent(med_pct_g, accuracy = 0.01)),
             hjust = 1, vjust = -0.4, color = "firebrick", size = 3.1) +
    geom_point(data = top_genes, color = "#08306B", size = 1.6) +
    ggrepel::geom_text_repel(
        data = top_genes, aes(label = gene_id),
        size = 3, color = "grey20",
        nudge_x = nrow(z_score_tally) * 0.02,
        nudge_y = 0.005,
        segment.size = 0.2, min.segment.length = 0,
        max.overlaps = Inf
    ) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(title    = sprintf("A. Genes ranked by underlier rate  (n = %s)",
                            comma(nrow(z_score_tally))),
         subtitle = "Sorted by pct_underlier = #underlier_subjects / #subjects",
         x        = "Gene rank (descending pct_underlier)",
         y        = "% of subjects underlying the gene") +
    theme_classic(base_size = 11) +
    theme(plot.title.position = "plot")

# ---- Panel B: subjects ranked by # underlier genes ----
top_subj <- head(z_score_tally_subj, 5)
med_n_s  <- median(z_score_tally_subj$n)

p_subj <- ggplot(z_score_tally_subj, aes(x = rank, y = n)) +
    geom_line(color = "#D95F0E", linewidth = 0.4) +
    geom_hline(yintercept = med_n_s, linetype = "dashed",
               color = "firebrick", linewidth = 0.4) +
    annotate("text",
             x = nrow(z_score_tally_subj), y = med_n_s,
             label = sprintf(" median = %s genes", comma(med_n_s)),
             hjust = 1, vjust = -0.4, color = "firebrick", size = 3.1) +
    geom_point(data = top_subj, color = "#7F2704", size = 1.6) +
    ggrepel::geom_text_repel(
        data = top_subj, aes(label = subjectid),
        size = 3, color = "grey20",
        nudge_x = nrow(z_score_tally_subj) * 0.02,
        nudge_y = max(top_subj$n) * 0.02,
        segment.size = 0.2, min.segment.length = 0,
        max.overlaps = Inf
    ) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(
        labels = comma,
        sec.axis = sec_axis(~ . / n_genes,
                            name   = "% of all genes",
                            labels = percent_format(accuracy = 0.1))
    ) +
    labs(title    = sprintf("B. Subjects ranked by # underlier genes  (n = %s)",
                            comma(nrow(z_score_tally_subj))),
         subtitle = "Sorted by # genes for which the subject is an underlier",
         x        = "Subject rank (descending count)",
         y        = "# underlier genes per subject") +
    theme_classic(base_size = 11) +
    theme(plot.title.position = "plot")

combined <- (p_gene | p_subj) +
    plot_annotation(
        title    = "TOPMed expression underlier distributions (QC-filtered subjects)",
        subtitle = filter_subtitle,
        caption  = caption_txt,
        theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                         plot.subtitle = element_text(size = 10, color = "grey25"),
                         plot.caption  = element_text(size = 9,  color = "grey35"))
    )

out_png <- "underlier_tallies.png"
out_pdf <- "underlier_tallies.pdf"
if (DEF$label != "strict") {
    out_png <- paste0("underlier_tallies", DEF$suffix, ".png")
    out_pdf <- paste0("underlier_tallies", DEF$suffix, ".pdf")
}
ggsave(file.path(OUT_DIR, out_png), combined, width = 11, height = 4.5, dpi = 200)
ggsave(file.path(OUT_DIR, out_pdf), combined, width = 11, height = 4.5, device = cairo_pdf)
cat("Wrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")
