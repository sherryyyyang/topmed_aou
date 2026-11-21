# Shared TOPMed vs AoU RNA-outlier prevalence scatter inputs for Fig1d / Fig2d.
#
# Both panels must use identical gene rows, -log10 coordinates, n, and
# Spearman/Pearson. Only the renderer (bubble stack vs LOEUF coloring) differs.
#
# Requires _paths.R already sourced (FIG1C, PC_SETTING).
# Returns a list:
#   genes  — one row per gene (n_topmed/n_new >= MIN_N), with x/y = -log10(prev)
#   lim    — shared axis limits from true (non-jittered) x/y
#   sp, pe — Spearman on prevalence; Pearson on -log10 axes
#   path   — TSV path used

.load_topmed_aou_prevalence_scatter <- function(min_n = 2L, pad_frac = 0.06) {
  path <- file.path(FIG1C, "topmed_vs_new_z_-3_prevalence_per_gene.tsv")
  if (!file.exists(path)) {
    stop("missing prevalence table for Fig1d/Fig2d: ", path,
         " (PC_SETTING=", PC_SETTING, ")")
  }
  genes <- readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::filter(n_topmed >= min_n, n_new >= min_n) |>
    dplyr::mutate(
      x = -log10(prev_topmed),
      y = -log10(prev_new)
    )
  if (!nrow(genes)) stop("no genes left after n>=", min_n, " filter in ", path)

  lo <- min(c(genes$x, genes$y))
  hi <- max(c(genes$x, genes$y))
  pad <- pad_frac * (hi - lo)
  lim <- c(lo - pad, hi + pad)

  list(
    genes = genes,
    lim = lim,
    sp = cor(genes$prev_topmed, genes$prev_new, method = "spearman"),
    pe = cor(genes$x, genes$y, method = "pearson"),
    path = path,
    min_n = as.integer(min_n),
    n_genes = nrow(genes)
  )
}
