#!/usr/bin/env Rscript
# Per-gene Pearson (and Spearman) correlation of residualized expression z
# (55 expression PCs + 5 genotype PCs) with NNLS cell-type proportions.
# Also correlates uncorrected log2CPM (expr_cpm.rds) for the same genes.
#
# Does NOT load expr_z_join.rds. Does not rebuild S2–S7 compute.
# Writes Supplementary Figure 8 and reassembles Supplementary Figure 1 as
# a 2 x 2: a | b (OR + composition R² overlay) over c (gene-r) | d (VEP).
#
# Samples: kept_subjects (5,489) ∩ NNLS ∩ z. Report n vs the 5,560 z grid.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(cowplot)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

OUT <- file.path(FINAL_DRAFT, "supp_figs")
TSV <- file.path(OUT, "tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TSV, recursive = TRUE, showWarnings = FALSE)

# Follows PC_SETTING (see scripts/helpers/_paths.R); 55-PC table by default.
Z_LONG <- Z_LONGFORM
PROP_PATH <- file.path(TSV, "S7_celltype_nnls_proportions.tsv")
CUM_PATH <- file.path(TSV, "S7_cumulative_composition_r2.tsv")

CT_ORDER <- c(
  "Erythroid", "Platelet", "Neutrophil", "Monocyte_macrophage",
  "CD4_T", "CD8_T", "NK", "B", "Plasma"
)
CT_LABEL <- c(
  Erythroid = "Erythroid",
  Platelet = "Platelet",
  Neutrophil = "Neutrophil",
  Monocyte_macrophage = "Monocyte / macrophage",
  CD4_T = "CD4 T",
  CD8_T = "CD8 T",
  NK = "NK",
  B = "B",
  Plasma = "Plasma"
)
CT_COLS <- c(
  Erythroid            = unname(OI["vermillion"]),
  Platelet             = unname(OI["orange"]),
  Neutrophil           = unname(OI["yellow"]),
  Monocyte_macrophage  = unname(OI["blue"]),
  CD4_T                = unname(OI["sky"]),
  CD8_T                = unname(OI["green"]),
  NK                   = unname(OI["purple"]),
  B                    = unname(OI["black"]),
  Plasma               = "grey55"
)

SKIP_COMPUTE <- Sys.getenv("SKIP_COMPUTE", unset = "0") == "1"
ASSEMBLE_S1 <- Sys.getenv("ASSEMBLE_S1", unset = "1") == "1"
PY <- Sys.getenv("PAPER_PYTHON", unset = "/net/wonderland/home/sheyang/.pixi/envs/python/bin/python3")

write_tsv_comment <- function(df, path, comments) {
  tmp <- tempfile(fileext = ".tsv")
  data.table::fwrite(df, tmp, sep = "\t", quote = FALSE, na = "")
  writeLines(c(comments, readLines(tmp, warn = FALSE)), path)
  unlink(tmp)
  cat("  wrote ", basename(path), " (", nrow(df), " rows)\n", sep = "")
  invisible(path)
}

pearson_xy <- function(X, Y) {
  # X: samples x genes; Y: samples x types. Both numeric, no NA.
  n <- nrow(X)
  Xs <- scale(X)
  Ys <- scale(Y)
  crossprod(Xs, Ys) / (n - 1)
}

p_from_r <- function(r, n) {
  df <- n - 2
  den <- pmax(1 - r * r, .Machine$double.eps)
  tstat <- r * sqrt(df / den)
  2 * pt(-abs(tstat), df = df)
}

rank_cols <- function(M) {
  apply(M, 2, rank, ties.method = "average", na.last = "keep")
}

long_cor <- function(r_mat, p_mat, n, genes, symbols, method_note = "pearson") {
  ct <- colnames(r_mat)
  out <- vector("list", length(ct))
  for (j in seq_along(ct)) {
    out[[j]] <- data.table(
      gene_id = genes,
      gene_nv = sub("\\.\\d+$", "", genes),
      symbol = unname(symbols[genes]),
      cell_type = ct[j],
      r = as.numeric(r_mat[, j]),
      p = as.numeric(p_mat[, j]),
      n = n,
      method = method_note
    )
  }
  rbindlist(out)
}

summarise_r <- function(long_dt, cell_types = CT_ORDER) {
  long_dt <- as.data.table(long_dt)
  long_dt[, abs_r := abs(r)]
  out <- long_dt[, .(
    n_genes = .N,
    n_samples = as.integer(median(n)),
    median_abs_r = median(abs_r, na.rm = TRUE),
    p95_abs_r = as.numeric(quantile(abs_r, 0.95, na.rm = TRUE, names = FALSE)),
    n_abs_r_gt_0.2 = sum(abs_r > 0.2, na.rm = TRUE),
    n_abs_r_gt_0.3 = sum(abs_r > 0.3, na.rm = TRUE),
    n_abs_r_gt_0.5 = sum(abs_r > 0.5, na.rm = TRUE)
  ), by = cell_type]
  out[, cell_type := factor(cell_type, levels = cell_types)]
  setorder(out, cell_type)
  out[, cell_type := as.character(cell_type)]
  out
}

ct_lab_factor <- function(x, rev = FALSE) {
  labs <- unname(CT_LABEL[as.character(x)])
  levs <- unname(CT_LABEL[CT_ORDER])
  if (rev) levs <- rev(levs)
  factor(labs, levels = levs)
}

INTERP <- paste(
  "If residualization worked, gene-z vs composition |r| should be small",
  "(near 0) compared with uncorrected log2CPM."
)

# ---------------------------------------------------------------------------
# A. Correlations
# ---------------------------------------------------------------------------
if (SKIP_COMPUTE) {
  cat("SKIP_COMPUTE=1: reusing S8 TSVs\n")
  pear <- as.data.table(read_tsv(
    file.path(TSV, "S8_gene_z_vs_celltype_pearson.tsv"),
    comment = "#", show_col_types = FALSE
  ))
  summ <- as.data.table(read_tsv(
    file.path(TSV, "S8_gene_z_vs_celltype_summary.tsv"),
    comment = "#", show_col_types = FALSE
  ))
  cpm_pear <- as.data.table(read_tsv(
    file.path(TSV, "S8_gene_log2cpm_vs_celltype_pearson.tsv"),
    comment = "#", show_col_types = FALSE
  ))
  n_used <- as.integer(median(pear$n, na.rm = TRUE))
  n_genes <- uniqueN(pear$gene_id)
  n_z_grid <- 5560L
  n_kept <- 5489L
} else {
  stopifnot(file.exists(Z_LONG), file.exists(PROP_PATH))

  cat("Loading kept subjects, NNLS proportions, symbols...\n")
  kept <- readRDS(KEPT_SUBJECTS)
  prop_df <- fread(PROP_PATH)
  ct_keep <- intersect(CT_ORDER, names(prop_df))
  stopifnot(length(ct_keep) == 9L)

  genc <- fread(
    GENCODE_SYM, header = FALSE, col.names = c("gene_nv", "symbol")
  )
  genc <- unique(genc, by = "gene_nv")

  cat("Loading uncorrected log2CPM (expr_cpm.rds; not expr_cpm_adj)...\n")
  cpm <- readRDS(file.path(DATA_EXPR, "expr_cpm.rds"))
  cpm_genes <- rownames(cpm)
  cpm_samples <- colnames(cpm)

  nnls_samples <- prop_df$sample
  z_samples_expected <- unique(c(nnls_samples, cpm_samples))
  cat(
    "NNLS samples: ", uniqueN(nnls_samples),
    "  CPM samples: ", length(cpm_samples),
    "  kept: ", length(kept), "\n", sep = ""
  )

  cat("fread residualized z longer-form (55 expr PCs + 5 geno PCs)...\n")
  zlong <- fread(
    Z_LONG,
    colClasses = c(
      gene_id = "character",
      subjectid = "character",
      expression_zscore = "numeric"
    ),
    showProgress = TRUE
  )
  n_z_grid <- uniqueN(zlong$subjectid)
  n_z_genes <- uniqueN(zlong$gene_id)
  cat("z longer-form: ", n_z_genes, " genes x ", n_z_grid, " subjects\n", sep = "")

  genes <- intersect(unique(zlong$gene_id), cpm_genes)
  cat("gene universe z ∩ expr_cpm: ", length(genes), "\n", sep = "")
  stopifnot(length(genes) >= 8000)

  samples_5560 <- intersect(intersect(unique(zlong$subjectid), nnls_samples), cpm_samples)
  samples_used <- intersect(samples_5560, kept)
  n_used <- length(samples_used)
  n_kept <- length(kept)
  cat(
    "samples z ∩ NNLS ∩ CPM: ", length(samples_5560),
    "  ∩ kept: ", n_used, "\n", sep = ""
  )
  stopifnot(n_used >= 5000)

  zlong <- zlong[gene_id %in% genes & subjectid %in% samples_used]
  n_per <- zlong[, .N, by = gene_id]
  stopifnot(all(n_per$N == n_used))

  setorder(zlong, gene_id, subjectid)
  gene_levels <- unique(zlong$gene_id)
  subj_ref <- zlong[gene_id == gene_levels[1L], subjectid]
  last_g <- gene_levels[length(gene_levels)]
  subj_last <- zlong[gene_id == last_g, subjectid]
  balanced <- identical(subj_ref, subj_last) &&
    nrow(zlong) == length(gene_levels) * n_used

  if (balanced) {
    cat("Reshaping z with matrix() (balanced gene-major grid)...\n")
    zmat <- matrix(
      zlong$expression_zscore,
      nrow = length(gene_levels),
      ncol = n_used,
      byrow = TRUE,
      dimnames = list(gene_levels, subj_ref)
    )
  } else {
    cat("Reshaping z with dcast...\n")
    zw <- dcast(zlong, gene_id ~ subjectid, value.var = "expression_zscore")
    zmat <- as.matrix(zw[, -1, with = FALSE])
    rownames(zmat) <- zw$gene_id
    rm(zw)
  }
  rm(zlong)
  gc(verbose = FALSE)

  samples_used <- intersect(colnames(zmat), samples_used)
  samples_used <- samples_used[order(match(samples_used, subj_ref))]
  zmat <- zmat[, samples_used, drop = FALSE]
  n_used <- ncol(zmat)
  n_genes <- nrow(zmat)

  prop <- as.matrix(prop_df[match(samples_used, sample), ..ct_keep])
  rownames(prop) <- samples_used
  stopifnot(!anyNA(prop), nrow(prop) == n_used)

  cpm <- cpm[rownames(zmat), samples_used, drop = FALSE]
  stopifnot(identical(rownames(cpm), rownames(zmat)))
  stopifnot(identical(colnames(cpm), samples_used))

  n_na_z <- sum(is.na(zmat))
  n_na_cpm <- sum(is.na(cpm))
  cat("NA in z matrix: ", n_na_z, "  NA in CPM: ", n_na_cpm, "\n", sep = "")
  if (n_na_z > 0 || n_na_cpm > 0) {
    stop("NA present; pairwise-complete path not implemented. Inspect inputs.")
  }

  cat("Pearson + Spearman (residualized z and uncorrected log2CPM)...\n")
  Xz <- t(zmat)
  Xc <- t(cpm)
  r_z <- pearson_xy(Xz, prop)
  r_cpm <- pearson_xy(Xc, prop)
  p_z <- p_from_r(r_z, n_used)
  p_cpm <- p_from_r(r_cpm, n_used)

  r_z_sp <- pearson_xy(rank_cols(Xz), rank_cols(prop))
  r_cpm_sp <- pearson_xy(rank_cols(Xc), rank_cols(prop))
  p_z_sp <- p_from_r(r_z_sp, n_used)
  p_cpm_sp <- p_from_r(r_cpm_sp, n_used)

  symbols <- genc$symbol[match(sub("\\.\\d+$", "", rownames(zmat)), genc$gene_nv)]
  names(symbols) <- rownames(zmat)

  pear <- long_cor(r_z, p_z, n_used, rownames(zmat), symbols, "pearson")
  pear[, q_BH := p.adjust(p, method = "BH")]
  pear[, spearman_r := as.numeric(r_z_sp[cbind(gene_id, cell_type)])]
  pear[, spearman_p := as.numeric(p_z_sp[cbind(gene_id, cell_type)])]
  pear[, spearman_q_BH := p.adjust(spearman_p, method = "BH")]
  pear[, cell_type := factor(cell_type, levels = ct_keep)]
  setorder(pear, cell_type, gene_id)
  pear[, cell_type := as.character(cell_type)]

  cpm_pear <- long_cor(r_cpm, p_cpm, n_used, rownames(zmat), symbols, "pearson")
  cpm_pear[, q_BH := p.adjust(p, method = "BH")]
  cpm_pear[, spearman_r := as.numeric(r_cpm_sp[cbind(gene_id, cell_type)])]
  cpm_pear[, spearman_p := as.numeric(p_cpm_sp[cbind(gene_id, cell_type)])]
  cpm_pear[, spearman_q_BH := p.adjust(spearman_p, method = "BH")]
  cpm_pear[, cell_type := factor(cell_type, levels = ct_keep)]
  setorder(cpm_pear, cell_type, gene_id)
  cpm_pear[, cell_type := as.character(cell_type)]

  summ <- summarise_r(pear, ct_keep)
  summ_cpm <- summarise_r(cpm_pear, ct_keep)
  setnames(
    summ_cpm,
    setdiff(names(summ_cpm), "cell_type"),
    paste0(setdiff(names(summ_cpm), "cell_type"), "_uncorrected")
  )
  summ <- merge(summ, summ_cpm, by = "cell_type", all.x = TRUE, sort = FALSE)
  summ[, cell_type := factor(cell_type, levels = ct_keep)]
  setorder(summ, cell_type)
  summ[, cell_type := as.character(cell_type)]

  cmp <- merge(
    pear[, .(gene_id, gene_nv, symbol, cell_type, r_z = r, p_z = p, q_BH_z = q_BH,
             spearman_r_z = spearman_r)],
    cpm_pear[, .(gene_id, cell_type, r_uncorrected = r, p_uncorrected = p,
                 q_BH_uncorrected = q_BH, spearman_r_uncorrected = spearman_r)],
    by = c("gene_id", "cell_type")
  )
  cmp[, cell_type := factor(cell_type, levels = ct_keep)]
  setorder(cmp, cell_type, gene_id)
  cmp[, cell_type := as.character(cell_type)]

  hdr <- c(
    paste0("# ", INTERP),
    paste0(
      "# residualized z = 55 expression PCs + 5 genotype PCs, per-gene scaled; ",
      "source=", basename(Z_LONG)
    ),
    paste0(
      "# n_genes=", n_genes,
      "; n_samples_used=", n_used,
      " (kept ∩ NNLS ∩ z ∩ CPM); n_z_grid=", n_z_grid,
      "; n_kept=", n_kept,
      "; n_tests=", nrow(pear),
      "; BH-FDR across all gene x cell-type Pearson tests"
    ),
    "# proportions: S7_celltype_nnls_proportions.tsv (uncorrected CPM NNLS)"
  )
  hdr_cpm <- c(
    "# Uncorrected log2CPM (expr_cpm.rds, NOT expr_cpm_adj) vs NNLS proportions.",
    paste0(
      "# Same genes and samples as S8_gene_z_vs_celltype_pearson.tsv; ",
      "n_genes=", n_genes, "; n_samples_used=", n_used
    )
  )
  hdr_sum <- c(
    paste0("# ", INTERP),
    paste0(
      "# residualized-z summary; n_genes=", n_genes,
      "; n_samples_used=", n_used,
      "; uncorrected columns are |r| from expr_cpm.rds"
    )
  )

  pear_out <- pear[, .(
    gene_id, gene_nv, symbol, cell_type, r, p, n, q_BH,
    spearman_r, spearman_p, spearman_q_BH
  )]
  cpm_out <- cpm_pear[, .(
    gene_id, gene_nv, symbol, cell_type, r, p, n, q_BH,
    spearman_r, spearman_p, spearman_q_BH
  )]

  write_tsv_comment(pear_out, file.path(TSV, "S8_gene_z_vs_celltype_pearson.tsv"), hdr)
  write_tsv_comment(summ, file.path(TSV, "S8_gene_z_vs_celltype_summary.tsv"), hdr_sum)
  write_tsv_comment(
    cpm_out, file.path(TSV, "S8_gene_log2cpm_vs_celltype_pearson.tsv"), hdr_cpm
  )
  write_tsv_comment(
    cmp, file.path(TSV, "S8_gene_z_vs_uncorrected_r.tsv"),
    c(
      paste0("# ", INTERP),
      paste0("# paired r_z (residualized) vs r_uncorrected (log2CPM); n=", n_used)
    )
  )

  cat("\n========== KEY STATS (not on graphics) ==========\n")
  cat("n genes x types: ", n_genes, " x ", length(ct_keep),
      " = ", nrow(pear), "\n", sep = "")
  cat("n_samples_used (kept): ", n_used,
      "  n_z_grid: ", n_z_grid,
      "  n_kept_file: ", n_kept, "\n", sep = "")
  cat("median |r| residualized z per cell type:\n")
  print(as.data.frame(summ[, .(cell_type, median_abs_r, p95_abs_r,
                               n_abs_r_gt_0.2, n_abs_r_gt_0.3, n_abs_r_gt_0.5)]))
  cat("median |r| uncorrected log2CPM per cell type:\n")
  print(as.data.frame(summ[, .(
    cell_type,
    median_abs_r_uncorrected,
    n_abs_r_gt_0.3_uncorrected
  )]))
  n_gt03 <- sum(abs(pear$r) > 0.3, na.rm = TRUE)
  n_gt03_cpm <- sum(abs(cpm_pear$r) > 0.3, na.rm = TRUE)
  cat(
    "gene x type tests with |r|>0.3: residualized z ", n_gt03,
    " / ", nrow(pear),
    " ; uncorrected ", n_gt03_cpm, " / ", nrow(cpm_pear), "\n", sep = ""
  )

  globins <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  gchk <- cmp[symbol %in% globins & cell_type == "Erythroid",
              .(symbol, r_z, r_uncorrected)]
  if (nrow(gchk)) {
    cat("Sanity (globin vs Erythroid):\n")
    print(as.data.frame(gchk))
  }
  rm(Xz, Xc, zmat, cpm, prop)
  gc(verbose = FALSE)
}

# ---------------------------------------------------------------------------
# B. Plots — S8 and the S1 gene-r panel
# ---------------------------------------------------------------------------
pear[, ct_lab := ct_lab_factor(cell_type)]
cpm_pear[, ct_lab := ct_lab_factor(cell_type)]

# Residualized |r| is ~0.01; ±1 hides the violins. ±0.5 still shows they are small
# and keeps the single Neutrophil tail (S100A8 r = 0.36) on scale.
p_violin <- ggplot(pear, aes(ct_lab, r, fill = cell_type)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_violin(
    color = NA, alpha = 0.9, width = 0.92, trim = TRUE, scale = "width"
  ) +
  geom_boxplot(
    width = 0.18, outlier.shape = 16, outlier.size = 0.55,
    fill = "white", linewidth = 0.35, color = "grey20"
  ) +
  scale_fill_manual(values = CT_COLS, guide = "none") +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, 0.25)) +
  labs(x = NULL, y = "Pearson r") +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 8),
    plot.margin = margin(4, 6, 8, 8)
  )

# Top 20 genes by max |r| across types (union of top 10/type is ~90 near-zero
# genes and unreadable). Neutrophil residual signal dominates this list.
gene_ord <- pear[, .(max_abs = max(abs(r))), by = .(gene_id, symbol)]
setorder(gene_ord, -max_abs)
gene_ord <- gene_ord[seq_len(min(20L, .N))]
gene_ord[, ylab := fifelse(
  !is.na(symbol) & symbol != "",
  symbol,
  sub("\\.\\d+$", "", gene_id)
)]
gene_ord[, ylab := ifelse(
  duplicated(ylab) | duplicated(ylab, fromLast = TRUE),
  paste0(ylab, " (", sub("\\.\\d+$", "", gene_id), ")"),
  ylab
)]
heat <- pear[gene_id %in% gene_ord$gene_id]
heat[, ylab := factor(
  gene_ord$ylab[match(gene_id, gene_ord$gene_id)],
  levels = rev(gene_ord$ylab)
)]
heat[, ct_lab := ct_lab_factor(cell_type)]

p_heat <- ggplot(heat, aes(ct_lab, ylab, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(
    low = unname(OI["blue"]), mid = "white", high = unname(OI["vermillion"]),
    midpoint = 0, limits = c(-0.5, 0.5), oob = squish, name = "Pearson r",
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

cmp_long <- rbindlist(list(
  pear[, .(cell_type, ct_lab, abs_r = abs(r), matrix = "Residualized z")],
  cpm_pear[, .(cell_type, ct_lab, abs_r = abs(r), matrix = "Uncorrected log2CPM")]
))
cmp_long[, matrix := factor(
  matrix, levels = c("Uncorrected log2CPM", "Residualized z")
)]

p_cmp <- ggplot(cmp_long, aes(ct_lab, abs_r, fill = matrix)) +
  geom_violin(
    color = NA, alpha = 0.88, width = 0.92, trim = TRUE,
    position = position_dodge(width = 0.78)
  ) +
  geom_boxplot(
    width = 0.14, outlier.shape = NA, color = "grey20", linewidth = 0.3,
    position = position_dodge(width = 0.78)
  ) +
  scale_fill_manual(
    values = c(
      "Uncorrected log2CPM" = unname(OI["orange"]),
      "Residualized z" = unname(OI["blue"])
    ),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "|Pearson r|") +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 8),
    legend.position = "bottom",
    legend.justification = "center"
  )

tag_s <- function(p) p + plot_annotation(tag_levels = "a") & theme_paper_tag()
sfig8 <- tag_s(
  p_violin / p_heat / p_cmp + plot_layout(heights = c(1.20, 1.70, 1.25))
)
paper_ggsave(sfig8, file.path(OUT, "Supplementary_Figure8"),
             width = NATURE_DOUBLE_IN, height = 8.8)
paper_ggsave(p_violin, file.path(OUT, "FigS8a_gene_z_celltype_r"),
             width = NATURE_DOUBLE_IN, height = 3.0)
paper_ggsave(p_heat, file.path(OUT, "FigS8b_top_gene_r_heatmap"),
             width = NATURE_DOUBLE_IN, height = 4.2)
paper_ggsave(p_cmp, file.path(OUT, "FigS8c_uncorrected_vs_residualized_r"),
             width = NATURE_DOUBLE_IN, height = 3.2)
# Extra single-panel copy used as S1c.
paper_ggsave(p_violin, file.path(OUT, "FigS1c_gene_z_celltype_r"),
             width = NATURE_DOUBLE_IN / 2 + 0.15, height = 3.6)

cat("Wrote Supplementary_Figure8 and FigS8a–c\n")

# ---------------------------------------------------------------------------
# C. Reassemble Supplementary Figure 1 as a | b only
# Gene-z vs cell-type r stays on Supplementary Figure 8; PTV consequences is
# Supplementary Figure 10 (written below). Do not stack c|d under S1.
# ---------------------------------------------------------------------------
if (ASSEMBLE_S1) {
  cat("Assembling Supplementary Figure 1 (a | b enrichment only)...\n")
  COLS <- c("cutoff", "npcs_format", "odds_ratio", "ci_95_lo", "ci_95_hi", "p_value")
  load_sweep <- function(f) {
    readRDS(f) |>
      filter(direction == "underlier") |>
      mutate(npcs_format = as.numeric(npcs_format), cutoff = as.numeric(cutoff)) |>
      dplyr::select(all_of(COLS))
  }
  sv <- load_sweep(file.path(DATA_ENRICH, "expr_enrichments_sv_full.rds")) |>
    mutate(variant_class = "Rare SV within 100 kb (AF <= 1%)")
  ptv <- load_sweep(file.path(DATA_ENRICH, "expr_enrichments_ptv.rds")) |>
    mutate(variant_class = "LoF / PTV in gene (AF <= 0.5%)")
  dat <- bind_rows(sv, ptv) |>
    mutate(
      cutoff_lab = factor(
        sprintf("z < -%s", cutoff),
        levels = sprintf("z < -%s", sort(unique(cutoff)))
      ),
      variant_class = factor(
        variant_class,
        levels = c(unique(sv$variant_class), unique(ptv$variant_class))
      )
    )
  peaks <- dat |>
    group_by(variant_class, cutoff_lab) |>
    slice_max(odds_ratio, n = 1) |>
    ungroup()
  cum <- read_tsv(CUM_PATH, show_col_types = FALSE)
  overlay <- cum |> transmute(npcs_format = n_pcs, comp_r2 = mean_adj_r2)
  XBREAKS <- c(0, 1000, 2000, 3000, 4000)
  r2_to_or <- function(r2, or_lo, or_hi) {
    10^(log10(or_lo) + pmin(pmax(r2, 0), 1) * (log10(or_hi) - log10(or_lo)))
  }
  panel_overlay <- function(vc) {
    d <- filter(dat, variant_class == vc)
    pk <- filter(peaks, variant_class == vc)
    short <- if (grepl("Rare SV", vc)) "Rare SVs" else "Rare PTVs"
    or_lo <- min(c(d$ci_95_lo, d$odds_ratio, 1), na.rm = TRUE)
    or_hi <- max(c(d$ci_95_hi, d$odds_ratio), na.rm = TRUE)
    ov <- overlay |>
      filter(npcs_format >= min(d$npcs_format), npcs_format <= max(d$npcs_format)) |>
      mutate(or_y = r2_to_or(comp_r2, or_lo, or_hi))
    r2_breaks <- c(0, 0.25, 0.5, 0.75, 1)
    ggplot(d, aes(npcs_format, odds_ratio, color = cutoff_lab, fill = cutoff_lab)) +
      geom_hline(yintercept = 1, color = "grey60", linewidth = 0.35) +
      geom_vline(
        xintercept = c(55, 2005), linetype = "dashed",
        color = "grey35", linewidth = 0.45
      ) +
      annotate(
        "text", x = 55, y = Inf, label = " 55", vjust = 1.4, hjust = 0,
        size = 2.6, color = "grey20", fontface = "bold", family = PAPER_FONT
      ) +
      annotate(
        "text", x = 2005, y = Inf, label = "2005 ", vjust = 1.4, hjust = 1,
        size = 2.6, color = "grey20", fontface = "bold", family = PAPER_FONT
      ) +
      geom_ribbon(aes(ymin = ci_95_lo, ymax = ci_95_hi), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.2) +
      geom_point(data = pk, shape = 21, size = 3.4, stroke = 1.1, fill = "white") +
      geom_line(
        data = ov, aes(npcs_format, or_y),
        inherit.aes = FALSE, color = unname(OI["black"]),
        linetype = "longdash", linewidth = 0.7
      ) +
      scale_color_manual(values = CUTOFF_COLORS, name = "Z-score cutoff") +
      scale_fill_manual(values = CUTOFF_COLORS, guide = "none") +
      guides(color = guide_legend(
        nrow = 1,
        override.aes = list(size = 1.5, stroke = 0.6, linewidth = 0.5)
      )) +
      scale_x_continuous(
        labels = comma, breaks = XBREAKS,
        expand = expansion(mult = c(0.02, 0.03))
      ) +
      scale_y_log10(
        breaks = c(1, 1.5, 2, 3, 5, 10, 20, 50, 100, 200),
        labels = c("1", "1.5", "2", "3", "5", "10", "20", "50", "100", "200"),
        sec.axis = sec_axis(
          transform = ~ (log10(.) - log10(or_lo)) / (log10(or_hi) - log10(or_lo)),
          name = expression("Composition " * R^2),
          breaks = r2_breaks,
          labels = c("0", "0.25", "0.50", "0.75", "1")
        )
      ) +
      labs(
        x = "# Expression PCs",
        y = "Odds ratio",
        title = short
      ) +
      theme_paper() +
      theme_paper_inner_title() +
      theme(
        axis.title.y.right = element_text(size = 10, color = "black"),
        axis.text.y.right = element_text(size = 8, color = "black"),
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
  pa <- panel_overlay(levels(dat$variant_class)[1])
  pb <- panel_overlay(levels(dat$variant_class)[2])
  paper_ggsave(pa, file.path(OUT, "FigS1a_sv_enrichment"),
               width = NATURE_DOUBLE_IN / 2 + 0.35, height = 4.2)
  paper_ggsave(pb, file.path(OUT, "FigS1b_ptv_enrichment"),
               width = NATURE_DOUBLE_IN / 2 + 0.35, height = 4.2)

  s1c_env <- new.env(parent = globalenv())
  options(paper.assemble = TRUE)
  source(
    file.path(ROOT, "scripts/supplement/plot_ptv_consequence_summary.R"),
    local = s1c_env
  )
  options(paper.assemble = FALSE)
  paper_ggsave(s1c_env$.paper_plots$c, file.path(OUT, "FigS1c_ptv_consequences"),
               width = NATURE_DOUBLE_IN, height = 5.6)
  paper_ggsave(s1c_env$.paper_plots$c, file.path(OUT, "FigS10a_ptv_consequences"),
               width = NATURE_DOUBLE_IN, height = 5.6)
  paper_ggsave(s1c_env$.paper_plots$c, file.path(OUT, "Supplementary_Figure10"),
               width = NATURE_DOUBLE_IN, height = 5.6)

  s1_panels <- cowplot::plot_grid(
    pa + theme(legend.position = "none"),
    pb + theme(legend.position = "none"),
    ncol = 2, labels = c("a", "b"),
    label_fontfamily = PAPER_FONT, label_size = 14, label_fontface = "bold",
    align = "hv"
  )
  s1_leg <- cowplot::get_legend(
    pb + theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.justification = "center"
    )
  )
  sfig1 <- cowplot::plot_grid(s1_panels, s1_leg, ncol = 1, rel_heights = c(1, 0.10))
  paper_ggsave(sfig1, file.path(OUT, "Supplementary_Figure1"),
               width = NATURE_DOUBLE_IN, height = 4.8)
  cat("Wrote Supplementary_Figure1 (a|b) and Supplementary_Figure10 (PTV consequences)\n")
}

# ---------------------------------------------------------------------------
# D. Excel tabs + MANIFEST (do not rebuild S2–S7 figures)
# ---------------------------------------------------------------------------
pack <- file.path(ROOT, "scripts/supplement/pack_supplementary_tables.py")
if (file.exists(pack)) {
  st <- system2(PY, shQuote(pack), stdout = "", stderr = "")
  if (!identical(st, 0L)) warning("Excel packager failed: ", st)
}

manifest <- tibble(
  file = list.files(OUT, pattern = "Supplementary_Figure|FigS", full.names = FALSE),
  kind = "figure"
) |>
  bind_rows(tibble(
    file = file.path("supplemental_tables", list.files(TSV)),
    kind = "table"
  ))
write_tsv(manifest, file.path(TSV, "MANIFEST.tsv"))
cat("Done. Figures in ", OUT, "\n", sep = "")
