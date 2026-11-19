#!/usr/bin/env Rscript
# de Goede 2021 known-function lncRNAs vs the pink lncRNAs on
# Fig3a_constraint_sanity_coding_lncrna.png.
#
# de Goede supplies only the gene list (Table S1 of_known_function).
# Underlier rate is TOPMed, not de Goede:
#   n kept_subjects with z_2005 < -3 AND log2CPM < gene_mean-1,
#   denominator 5,489 (same x as that figure).
#
# de Goede genes in the all-gene HVG z matrix are scored, including rate 0.
# Genes not in that matrix were not tested and are left out of the plot.
# The 539 plot lncRNAs are the control; protein-coding is not used.
#
# Does not overwrite Fig3a_constraint_sanity_coding_lncrna.png or the
# published Source_Data_Fig3a_ptv_vs_underlier_gene.tsv.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

FIG3_TSV <- file.path(
  ROOT, "allgenes_pipeline_out",
  "ptv_vs_underlier_gene_level_z_-3_allgenes_pc2005.tsv"
)
CONS_TSV <- file.path(SHEYANG, "annotation", "gene_constraint_fractions.tsv.gz")
KNOWN_TSV <- file.path(SHEYANG, "annotation", "goede_2021_known_function_lncrnas.tsv")
HVG_TSV <- file.path(ROOT, "allgenes_pipeline_out", "hvg_allgenes_gene_annotation.tsv")
Z_CACHE <- file.path(ROOT, "allgenes_pipeline_out", "expr_zscore_pc2005_allgenes.rds")
CPM_RDS <- file.path(ROOT, "allgenes_pipeline_out", "expr_cpm_adj.rds")
UNDER_COUNT_TSV <- file.path(
  ROOT, "allgenes_pipeline_out",
  "underlier_gene_counts_z_-3_allgenes_pc2005_kept.tsv"
)
TPM_TSV <- file.path(DATA_CORE, "gene_median_tpm_kept.tsv")
Z_CUT <- -3
LOGCPM_DROP <- 1

DRAFT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
SRC <- file.path(DRAFT, "source_panels")
SUPP <- file.path(DRAFT, "supplemental_data")
PUB_FIG3A_TSV <- file.path(SUPP, "Source_Data_Fig3a_ptv_vs_underlier_gene.tsv")
STEM <- "Fig3a_degoede_vs_plot_lncrna"

COL_DEGOEDE <- unname(OI["orange"])
COL_PLOT <- unname(OI["purple"])
PAL <- c(
  "de Goede known function" = COL_DEGOEDE,
  "Fig 3a plot lncRNA"      = COL_PLOT
)

nv <- function(x) str_remove(x, "\\..*$")

kept <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept)
stopifnot(n_kept == 5489L)

known <- read_tsv(KNOWN_TSV, show_col_types = FALSE) |>
  mutate(
    gene_nv = gene_stable,
    function_source = coalesce(function_source, ""),
    literature = str_detect(
      function_source, regex("hugo|lncrnadb|manual", ignore_case = TRUE)
    )
  )

cons <- read_tsv(CONS_TSV, show_col_types = FALSE) |>
  transmute(
    gene_nv = gene_id_stable,
    cons_type = gene_type_splitlncrna,
    frac_cons_exon_unique_primate
  )

hvg <- read_tsv(HVG_TSV, show_col_types = FALSE) |>
  mutate(gene_nv = nv(gene_nv)) |>
  distinct(gene_nv, .keep_all = TRUE)

tpm <- read_tsv(TPM_TSV, show_col_types = FALSE) |>
  mutate(gene_nv = nv(gene_nv)) |>
  distinct(gene_nv, .keep_all = TRUE)

# Pair files under allgenes_pipeline_out/underliers_z_-3.* are a different
# call set (not kept-subject Fig 3a 2005-PC). Recount from the same z/CPM
# matrices as build_ptv_vs_underlier_gene_level_allgenes_pc2005.R.
if (file.exists(UNDER_COUNT_TSV)) {
  under_n <- read_tsv(UNDER_COUNT_TSV, show_col_types = FALSE)
} else {
  cat("Recounting kept-subject underliers from all-gene 2005-PC z ...\n")
  z_mat <- read_rds(Z_CACHE)
  cpm <- as.matrix(read_rds(CPM_RDS))
  if (!setequal(rownames(z_mat), rownames(cpm))) {
    stop("z and expr_cpm_adj gene IDs do not match")
  }
  if (!setequal(colnames(z_mat), colnames(cpm))) {
    stop("z and expr_cpm_adj sample IDs do not match")
  }
  cpm <- cpm[rownames(z_mat), colnames(z_mat), drop = FALSE]
  gene_mean <- rowMeans(cpm)
  kept_in_z <- intersect(kept, colnames(z_mat))
  if (length(kept_in_z) != n_kept) {
    stop("kept_subjects missing from z grid: ", n_kept - length(kept_in_z))
  }
  n_hit <- rowSums(
    (z_mat[, kept_in_z, drop = FALSE] < Z_CUT) &
      (cpm[, kept_in_z, drop = FALSE] < (gene_mean - LOGCPM_DROP))
  )
  under_n <- tibble::tibble(
    gene_id = names(n_hit),
    gene_nv = nv(names(n_hit)),
    n_underlier_subjects = as.integer(n_hit)
  )
  rm(z_mat, cpm)
  gc(verbose = FALSE)
  write_tsv(under_n, UNDER_COUNT_TSV)
  cat("wrote ", UNDER_COUNT_TSV, "  genes=", nrow(under_n), "\n", sep = "")
}
under_n <- under_n |>
  mutate(gene_nv = nv(gene_nv)) |>
  group_by(gene_nv) |>
  summarise(n_underlier_subjects = sum(n_underlier_subjects), .groups = "drop")

fig3 <- read_tsv(FIG3_TSV, show_col_types = FALSE) |>
  left_join(cons, by = "gene_nv") |>
  mutate(
    map_class = case_when(
      cons_type == "protein_coding" ~ "Protein-coding",
      cons_type %in% c("existing_lncRNA", "new_lncRNA") ~ "lncRNA",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(map_class), is.finite(frac_cons_exon_unique_primate))

plot_lnc <- fig3 |>
  filter(map_class == "lncRNA") |>
  transmute(
    gene_nv, symbol, n_underlier_subjects, pct_under,
    cons_type, frac_cons_exon_unique_primate
  )
stopifnot(nrow(plot_lnc) == 539L)

# Same-assay rates for de Goede genes in the all-gene HVG matrix (zeros kept).
degoede <- known |>
  mutate(
    in_hvg = gene_nv %in% hvg$gene_nv,
    in_plot_lncrna = gene_nv %in% plot_lnc$gene_nv
  ) |>
  left_join(cons, by = "gene_nv") |>
  left_join(tpm |> select(gene_nv, median_tpm), by = "gene_nv") |>
  left_join(under_n, by = "gene_nv") |>
  mutate(
    n_underlier_subjects = if_else(
      in_hvg, as.integer(replace(n_underlier_subjects, is.na(n_underlier_subjects), 0L)),
      NA_integer_
    ),
    pct_under = if_else(in_hvg, 100 * n_underlier_subjects / n_kept, NA_real_),
    tpm_gt1 = is.finite(median_tpm) & median_tpm > 1
  )

# Pink-gene rates from the Fig 3a table must match kept-subject underlier counts.
chk <- plot_lnc |>
  left_join(under_n, by = "gene_nv", suffix = c("_fig3", "_call")) |>
  mutate(n_call = replace(n_underlier_subjects_call, is.na(n_underlier_subjects_call), 0L))
n_mismatch <- sum(chk$n_underlier_subjects_fig3 != chk$n_call)
if (n_mismatch > 0L) {
  stop("Fig 3a n_underlier_subjects disagrees with underliers_z_-3.tsv.gz for ",
       n_mismatch, " plot lncRNAs")
}

scored <- degoede |> filter(in_hvg)

lab_degoede <- sprintf("de Goede known function (n=%s)", comma(nrow(scored)))
lab_plot <- sprintf("Fig 3a plot lncRNA (n=%s)", comma(nrow(plot_lnc)))

summ_stats <- function(df, grp) {
  tibble(
    group = grp,
    n = nrow(df),
    n_ge1 = sum(df$pct_under > 0),
    frac_ge1 = mean(df$pct_under > 0),
    mean_pct_under = mean(df$pct_under),
    median_pct_under = median(df$pct_under)
  )
}

stats_a <- bind_rows(
  summ_stats(scored, "de Goede known function (HVG-scored)"),
  summ_stats(plot_lnc, "Fig 3a plot lncRNA")
)
scored_tpm <- scored |> filter(tpm_gt1)
stats_tpm <- bind_rows(
  summ_stats(scored_tpm, "de Goede HVG TPM>1"),
  summ_stats(plot_lnc, "Fig 3a plot lncRNA")
)

wilcox_p <- function(a, b) {
  wilcox.test(a$pct_under, b$pct_under, alternative = "two.sided")$p.value
}

p_a <- wilcox_p(scored, plot_lnc)
p_tpm <- wilcox_p(scored_tpm, plot_lnc)

cat("n_kept=", n_kept, "\n", sep = "")
cat("de Goede list=", nrow(known),
    "  in all-gene HVG=", nrow(scored),
    "  overlap with plot lncRNA=", sum(scored$in_plot_lncrna),
    "  not tested (not HVG)=", sum(!degoede$in_hvg), "\n", sep = "")
cat("plot lncRNA (pink)=", nrow(plot_lnc), "\n", sep = "")
cat("\nPanel a (HVG-scored de Goede vs plot lncRNA)\n")
print(stats_a, n = 10)
cat("Wilcoxon p=", signif(p_a, 3), "\n", sep = "")
cat("\nSensitivity: de Goede HVG TPM>1 vs plot lncRNA\n")
print(stats_tpm, n = 10)
cat("Wilcoxon p=", signif(p_tpm, 3),
    "  overlap=", sum(scored_tpm$in_plot_lncrna), "\n", sep = "")

strip_plot <- function(df, y_levels) {
  med <- df |>
    group_by(grp) |>
    summarise(median_pct_under = median(pct_under), .groups = "drop")
  pal <- c(COL_DEGOEDE, COL_PLOT)
  names(pal) <- y_levels
  ggplot(df, aes(x = pct_under, y = grp, color = grp)) +
    geom_jitter(height = 0.18, width = 0, size = 1.05, alpha = 0.55, stroke = 0) +
    geom_point(
      data = med,
      aes(x = median_pct_under, y = grp),
      shape = 23, size = 2.6, fill = "white", stroke = 0.8,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_continuous(
      labels = label_percent(scale = 1, accuracy = 0.1),
      expand = expansion(mult = c(0.01, 0.04))
    ) +
    labs(
      x = "Fig 3a y: % of subjects with RNA underliers in gene",
      y = NULL
    ) +
    theme_paper() +
    theme(
      axis.title.x = element_text(size = 8),
      axis.text.y = element_text(size = 8)
    )
}

df_a <- bind_rows(
  scored |>
    transmute(
      gene_nv, pct_under,
      grp = factor(lab_degoede, levels = c(lab_degoede, lab_plot))
    ),
  plot_lnc |>
    transmute(
      gene_nv, pct_under,
      grp = factor(lab_plot, levels = c(lab_degoede, lab_plot))
    )
)
fig <- strip_plot(df_a, c(lab_degoede, lab_plot))

dir.create(SRC, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)

if (normalizePath(file.path(SUPP, paste0("Source_Data_", STEM, ".tsv")),
                  mustWork = FALSE) ==
    normalizePath(PUB_FIG3A_TSV, mustWork = FALSE)) {
  stop("Refusing to overwrite published Source_Data_Fig3a")
}

paper_ggsave(fig, file.path(DRAFT, STEM), width = NATURE_DOUBLE_IN, height = 2.6)
paper_ggsave(fig, file.path(SRC, STEM), width = NATURE_DOUBLE_IN, height = 2.6)
paper_ggsave(
  fig,
  file.path(ROOT, "allgenes_pipeline_out", STEM),
  width = NATURE_DOUBLE_IN, height = 2.6
)

gene_tbl <- bind_rows(
  degoede |>
    transmute(
      gene_nv,
      symbol,
      in_degoede = TRUE,
      in_plot_lncrna,
      in_hvg,
      tpm_gt1,
      cons_type,
      function_source,
      literature,
      median_tpm,
      frac_cons_exon_unique_primate,
      n_underlier_subjects,
      pct_under
    ),
  plot_lnc |>
    filter(!gene_nv %in% known$gene_nv) |>
    left_join(tpm |> select(gene_nv, median_tpm), by = "gene_nv") |>
    transmute(
      gene_nv,
      symbol,
      in_degoede = FALSE,
      in_plot_lncrna = TRUE,
      in_hvg = TRUE,
      tpm_gt1 = is.finite(median_tpm) & median_tpm > 1,
      cons_type,
      function_source = NA_character_,
      literature = NA,
      median_tpm,
      frac_cons_exon_unique_primate,
      n_underlier_subjects,
      pct_under
    )
) |>
  arrange(desc(in_degoede), desc(in_plot_lncrna), desc(pct_under))

write_tsv(gene_tbl, file.path(SUPP, paste0("Source_Data_", STEM, ".tsv")))
write_tsv(
  bind_rows(
    stats_a |> mutate(comparison = "a_hvg_scored", wilcox_p = p_a),
    stats_tpm |> mutate(comparison = "sensitivity_hvg_tpmgt1", wilcox_p = p_tpm)
  ),
  file.path(SUPP, paste0("Source_Data_", STEM, "_summary.tsv"))
)

cat("wrote ", file.path(DRAFT, paste0(STEM, ".png")), "\n", sep = "")
cat("published Source_Data_Fig3a left alone: ", PUB_FIG3A_TSV, "\n", sep = "")
