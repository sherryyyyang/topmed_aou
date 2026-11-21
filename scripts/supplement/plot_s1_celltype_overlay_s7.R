#!/usr/bin/env Rscript
# Cell-type composition vs expression PCs, S1a/b overlay, and Supplementary
# Figure 7. Does not load expr_z_join.rds. Does not overwrite the pre-overlay
# S1 backups (*_no_overlay.*).
#
# Samples: n = 5,560 in expr_pca.rds ∩ expr_cpm.rds (SV-present PCA object).
# The paper QC cohort is 5,489 after dropping 71 extreme subjects; enrichment
# S1 also uses 5,560, so the overlay matches the OR sweep.

suppressPackageStartupMessages({
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


PY <- Sys.getenv("PAPER_PYTHON", unset = "/net/wonderland/home/sheyang/.pixi/envs/python/bin/python3")
SC_PATH <- file.path(
  ROOT, "pc_optimization",
  "single_cell_reference_600_per_donor_all_genes_including_globins.tsv.gz"
)
K_HEAT <- 40L
K_PC <- 80L
MAX_CUM_PC <- 3955L
NPC_MARKS <- c(55L, 2005L)
XBREAKS <- c(0, 1000, 2000, 3000, 4000)
GLOBINS <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

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

write_out <- function(df, name) {
  path <- file.path(TSV, name)
  write_tsv(df, path)
  cat("  wrote ", name, " (", nrow(df), " rows)\n", sep = "")
  invisible(path)
}

# Lawson–Hanson NNLS on the Gram matrix (k is small).
nnls_gram <- function(AtA, Atb, tol = 1e-10) {
  k <- length(Atb)
  P <- integer(0)
  R <- seq_len(k)
  x <- numeric(k)
  w <- as.numeric(Atb)
  max_iter <- 3L * k
  iter <- 0L
  while (length(R) && any(w[R] > tol) && iter < max_iter) {
    iter <- iter + 1L
    j <- R[which.max(w[R])]
    P <- sort.int(c(P, j))
    R <- setdiff(R, j)
    s <- numeric(k)
    s[P] <- as.numeric(solve(AtA[P, P, drop = FALSE], Atb[P], tol = 1e-12))
    inner <- 0L
    while (length(P) && any(s[P] <= 0) && inner < k) {
      inner <- inner + 1L
      Q <- P[s[P] <= 0]
      alpha <- min(x[Q] / pmax(x[Q] - s[Q], .Machine$double.eps))
      x <- x + alpha * (s - x)
      hit0 <- P[x[P] <= tol]
      R <- sort.int(unique(c(R, hit0)))
      P <- setdiff(P, R)
      s <- numeric(k)
      if (length(P)) {
        s[P] <- as.numeric(solve(AtA[P, P, drop = FALSE], Atb[P], tol = 1e-12))
      }
    }
    x <- s
    w <- as.numeric(Atb - AtA %*% x)
  }
  pmax(x, 0)
}

nnls_matrix <- function(A, B) {
  # A: genes x K; B: genes x samples. Returns samples x K, rows sum to 1.
  AtA <- crossprod(A)
  AtB <- crossprod(A, B)
  k <- ncol(A)
  n <- ncol(B)
  out <- matrix(0, n, k, dimnames = list(colnames(B), colnames(A)))
  for (i in seq_len(n)) {
    x <- nnls_gram(AtA, AtB[, i])
    s <- sum(x)
    out[i, ] <- if (s > 0) x / s else rep(1 / k, k)
  }
  out
}

r2_lm <- function(y, X) {
  y <- as.numeric(y)
  yc <- y - mean(y)
  sst <- sum(yc^2)
  if (sst <= 0) return(c(r2 = NA_real_, adj_r2 = NA_real_))
  fit <- .lm.fit(cbind(1, X), y)
  sse <- sum(fit$residuals^2)
  r2 <- 1 - sse / sst
  n <- length(y)
  p <- ncol(X)
  adj <- if (n > p + 1) 1 - (1 - r2) * (n - 1) / (n - p - 1) else NA_real_
  c(r2 = r2, adj_r2 = adj)
}

SKIP_COMPUTE <- Sys.getenv("SKIP_COMPUTE", unset = "0") == "1"

# ---------------------------------------------------------------------------
# A. Gene universe = PCA HVGs (rownames of expr_cpm / expr_cpm_adj)
# ---------------------------------------------------------------------------
if (SKIP_COMPUTE) {
  cat("SKIP_COMPUTE=1: reusing TSVs for plots\n")
  cum <- read_tsv(file.path(TSV, "S7_cumulative_composition_r2.tsv"), show_col_types = FALSE)
  per_pc <- read_tsv(file.path(TSV, "S7_pc_multivar_r2.tsv"), show_col_types = FALSE)
  r_long <- read_tsv(file.path(TSV, "S7_pc_celltype_pearson_r.tsv"), show_col_types = FALSE)
  r_mat <- as.matrix(r_long[, setdiff(names(r_long), "pc"), drop = FALSE])
  rownames(r_mat) <- paste0("PC", r_long$pc)
  prop_df <- read_tsv(file.path(TSV, "S7_celltype_nnls_proportions.tsv"), show_col_types = FALSE)
  ct_names <- intersect(CT_ORDER, names(prop_df))
  mean_prop <- colMeans(as.matrix(prop_df[, ct_names, drop = FALSE]))
  names(mean_prop) <- ct_names
  k_use <- nrow(per_pc)
} else {
cat("Loading HVG IDs and gencode map...\n")
cpm_ids <- rownames(readRDS(file.path(DATA_EXPR, "expr_cpm_adj.rds")))
genc <- read.delim(
  file.path(DATA_CORE, "gencode_v30_gene_symbols.tsv"),
  header = FALSE, stringsAsFactors = FALSE
)
names(genc) <- c("ensg", "symbol")
map <- tibble(
  ensg_version = cpm_ids,
  ensg = str_remove(cpm_ids, "\\.\\d+$"),
  symbol = genc$symbol[match(str_remove(cpm_ids, "\\.\\d+$"), genc$ensg)]
)
stopifnot(!anyNA(map$symbol))

dup_sym <- map |>
  count(symbol) |>
  filter(n > 1)
ens <- readRDS(file.path(DATA_EXPR, "ensembl_96_df.rds"))
gene_rows <- ens |>
  filter(type == "gene") |>
  transmute(ensg = str_remove(gene_id, "\\.\\d+$"), gene_type) |>
  distinct(ensg, .keep_all = TRUE)
map <- map |>
  left_join(gene_rows, by = "ensg")
stopifnot(nrow(map) == length(cpm_ids))
n_pc <- sum(map$gene_type == "protein_coding", na.rm = TRUE)

sym_path <- file.path(TSV, "S7_hvg_symbols.tsv")
write_out(distinct(map, symbol), "S7_hvg_symbols.tsv")

# Broader protein-coding check: every HVG is already protein-coding; the
# available bulk matrices (expr_cpm, expr_cpm_adj, expr_qc_normalized) are
# the same 8,124-gene PCA universe. A larger SC gene set exists but would
# not describe THESE PCs.
cat(
  "HVGs: ", nrow(map), " Ensembl IDs / ", n_distinct(map$symbol),
  " symbols; protein-coding: ", n_pc, "\n", sep = ""
)

# ---------------------------------------------------------------------------
# B. Cell-type signatures (stream SC; mean raw counts)
# ---------------------------------------------------------------------------
sig_path <- file.path(TSV, "S7_celltype_signature_mean_counts.tsv")
overlap_sc_path <- file.path(TSV, "S7_sc_overlap_raw.tsv")
cells_path <- file.path(TSV, "S7_sc_cell_counts.tsv")
need_sig <- !all(file.exists(c(sig_path, overlap_sc_path, cells_path)))
if (need_sig) {
  cat("Building cell-type signatures from single-cell reference...\n")
  py_script <- file.path(ROOT, "scripts/supplement/build_celltype_signatures.py")
  st <- system2(
    PY,
    c(
      py_script,
      "--sc", SC_PATH,
      "--symbols", sym_path,
      "--out-signature", sig_path,
      "--out-overlap", overlap_sc_path,
      "--out-cells", cells_path
    ),
    stdout = "", stderr = ""
  )
  if (!identical(st, 0L)) stop("build_celltype_signatures.py failed: ", st)
} else {
  cat("Reusing existing cell-type signatures\n")
}

sig_df <- read_tsv(sig_path, show_col_types = FALSE)
ct_keep <- intersect(CT_ORDER, setdiff(names(sig_df), "symbol"))
sig_mat <- as.matrix(sig_df[, ct_keep, drop = FALSE])
rownames(sig_mat) <- sig_df$symbol
sc_overlap <- read_tsv(overlap_sc_path, show_col_types = FALSE)
sc_cells <- read_tsv(cells_path, show_col_types = FALSE)

# ---------------------------------------------------------------------------
# C. Uncorrected bulk (log2 CPM) + PCs; match samples
# ---------------------------------------------------------------------------
cat("Loading expr_cpm.rds and expr_pca.rds...\n")
bulk <- readRDS(file.path(DATA_EXPR, "expr_cpm.rds"))
pca <- readRDS(file.path(DATA_EXPR, "expr_pca.rds"))
pcs_all <- as.matrix(pca$rotated)
rm(pca)
gc(verbose = FALSE)

common_samples <- intersect(colnames(bulk), rownames(pcs_all))
common_samples <- common_samples[order(match(common_samples, colnames(bulk)))]
n_samples <- length(common_samples)
cat("Samples in PCA ∩ uncorrected CPM: ", n_samples, "\n", sep = "")
stopifnot(n_samples >= 5000)
bulk <- bulk[, common_samples, drop = FALSE]
pcs_all <- pcs_all[common_samples, , drop = FALSE]

bulk_map <- map |>
  mutate(in_sc = symbol %in% rownames(sig_mat))
# Average duplicate symbols in bulk (IGF2, MATR3, SCO2).
bulk_sym <- rowsum(bulk, group = bulk_map$symbol, reorder = FALSE)
n_iso <- as.numeric(table(bulk_map$symbol)[rownames(bulk_sym)])
bulk_sym <- bulk_sym / n_iso
overlap_genes <- intersect(rownames(bulk_sym), rownames(sig_mat))
overlap_genes <- sort(overlap_genes)
bulk_o <- bulk_sym[overlap_genes, , drop = FALSE]
sig_o <- sig_mat[overlap_genes, , drop = FALSE]
n_overlap <- length(overlap_genes)
globin_in <- intersect(overlap_genes, GLOBINS)
cat("Gene overlap for deconvolution: ", n_overlap, " (globins: ",
    paste(globin_in, collapse = ","), ")\n", sep = "")

# Linear-scale NNLS: signature = mean counts, column-normalized to 1e6;
# bulk = 2^log2CPM (positive CPM-like). Same HVG ∩ SC gene universe.
sig_lin <- sweep(sig_o, 2, colSums(sig_o), "/") * 1e6
sig_lin[is.na(sig_lin)] <- 0
bulk_lin <- 2^bulk_o

cat("NNLS deconvolution (", n_samples, " samples x ", ncol(sig_lin),
    " cell types)...\n", sep = "")
prop <- nnls_matrix(sig_lin, bulk_lin)
prop <- prop[, CT_ORDER[CT_ORDER %in% colnames(prop)], drop = FALSE]

# ---------------------------------------------------------------------------
# D. PC–composition correlations and cumulative R²
# ---------------------------------------------------------------------------
n_pc_avail <- ncol(pcs_all)
k_use <- min(K_PC, n_pc_avail)
k_cum_max <- min(MAX_CUM_PC, n_pc_avail)
pcs_k <- pcs_all[, seq_len(k_cum_max), drop = FALSE]
colnames(pcs_k) <- paste0("PC", seq_len(k_cum_max))

ct_names <- colnames(prop)
r_mat <- matrix(
  NA_real_, k_use, length(ct_names),
  dimnames = list(paste0("PC", seq_len(k_use)), ct_names)
)
for (j in seq_along(ct_names)) {
  r_mat[, j] <- cor(pcs_k[, seq_len(k_use), drop = FALSE], prop[, j])[, 1]
}

per_pc <- lapply(seq_len(k_use), function(k) {
  rr <- r2_lm(pcs_k[, k], prop)
  tibble(
    pc = k,
    max_abs_r = max(abs(r_mat[k, ])),
    top_cell_type = ct_names[which.max(abs(r_mat[k, ]))],
    top_r = r_mat[k, which.max(abs(r_mat[k, ]))],
    multivar_r2 = unname(rr["r2"]),
    multivar_adj_r2 = unname(rr["adj_r2"])
  )
}) |> bind_rows()

enr <- read_tsv(
  file.path(SUP_PC, "enrichment_npc_sv_vs_ptv_underlier.tsv"),
  show_col_types = FALSE
)
npc_grid <- sort(unique(c(
  as.integer(enr$npcs_format),
  NPC_MARKS,
  1:k_use
)))
npc_grid <- npc_grid[npc_grid >= 1 & npc_grid <= k_cum_max]

cat("QR of centered PCs for cumulative composition R² (", k_cum_max,
    " columns)...\n", sep = "")
P_c <- scale(prop, center = TRUE, scale = FALSE)
X_c <- scale(pcs_k, center = TRUE, scale = FALSE)
sst_ct <- colSums(P_c^2)
sst_tot <- sum(P_c^2)
nobs <- nrow(P_c)
# One QR; Q[, 1:n] spans PC1..n plus the implicit intercept (columns centered).
Q <- qr.Q(qr(X_c, LAPACK = TRUE))
cat("Cumulative composition R² on ", length(npc_grid), " n-PC values...\n", sep = "")
cum <- lapply(npc_grid, function(n) {
  Qn <- Q[, seq_len(n), drop = FALSE]
  fitted <- Qn %*% (crossprod(Qn, P_c))
  sse_ct <- colSums((P_c - fitted)^2)
  r2_ct <- 1 - sse_ct / sst_ct
  r2_tot <- 1 - sum((P_c - fitted)^2) / sst_tot
  adj <- function(r2, p) {
    if (nobs > p + 1) 1 - (1 - r2) * (nobs - 1) / (nobs - p - 1) else NA_real_
  }
  tibble(
    n_pcs = n,
    mean_r2 = mean(r2_ct),
    mean_adj_r2 = mean(vapply(r2_ct, adj, numeric(1), p = n)),
    total_r2 = r2_tot,
    total_adj_r2 = adj(r2_tot, n)
  )
}) |> bind_rows()

# ---------------------------------------------------------------------------
# E. Stats TSV (not drawn on graphics)
# ---------------------------------------------------------------------------
r_at <- function(n) {
  row <- cum[match(n, cum$n_pcs), ]
  if (nrow(row) == 0 || is.na(row$n_pcs[1])) {
    return(tibble(n_pcs = n, mean_r2 = NA, mean_adj_r2 = NA,
                  total_r2 = NA, total_adj_r2 = NA))
  }
  row
}
s55 <- r_at(55L)
s2005 <- r_at(2005L)
# Has composition R² already saturated before 55?
pre55 <- filter(cum, n_pcs <= 55)
sat_n <- pre55$n_pcs[which.max(pre55$mean_adj_r2)]
delta_after <- s2005$mean_adj_r2 - s55$mean_adj_r2
peak_early <- max(per_pc$multivar_r2[per_pc$pc <= 10], na.rm = TRUE)
past_peak <- (s55$mean_adj_r2 >= 0.9 * max(pre55$mean_adj_r2, na.rm = TRUE)) &&
  (delta_after < 0.10)

pc1_5 <- per_pc |>
  filter(pc <= 5) |>
  mutate(
    note = sprintf("PC%d ~ %s r=%+.3f; max|r|=%.3f; R2(PC~types)=%.3f",
                   pc, CT_LABEL[top_cell_type], top_r, max_abs_r, multivar_r2)
  )

overlap_tab <- tibble(
  item = c(
    "n_hvg_ensembl",
    "n_hvg_unique_symbols",
    "n_hvg_protein_coding",
    "n_hvg_duplicate_symbols",
    "duplicate_symbols",
    "n_sc_rows",
    "n_sc_unique_symbols",
    "n_overlap_hvg_sc",
    "n_overlap_globins",
    "overlap_globins",
    "n_samples_pca_and_bulk",
    "sample_set",
    "n_kept_paper_qc",
    "n_cell_types",
    "cell_types",
    "broader_proteincoding_note",
    "bulk_matrix",
    "signature",
    "deconvolution",
    "cum_r2_at_55_mean_adj",
    "cum_r2_at_55_total_adj",
    "cum_r2_at_2005_mean_adj",
    "cum_r2_at_2005_total_adj",
    "cum_r2_saturates_by_55",
    "n_pcs_max_mean_adj_r2_among_le_55",
    "delta_mean_adj_r2_2005_minus_55",
    "pc1_to_5"
  ),
  value = c(
    nrow(map),
    n_distinct(map$symbol),
    n_pc,
    nrow(dup_sym),
    paste(dup_sym$symbol, collapse = ","),
    sc_overlap$value[sc_overlap$item == "n_sc_rows"],
    sc_overlap$value[sc_overlap$item == "n_sc_unique_symbols"],
    n_overlap,
    length(globin_in),
    paste(globin_in, collapse = ","),
    n_samples,
    "5560 PCA / expr_cpm intersection (SV-present); not the 5489 kept set",
    5489,
    length(ct_names),
    paste(ct_names, collapse = ","),
    paste0(
      "All 8124 HVGs are protein-coding. expr_cpm, expr_cpm_adj and ",
      "expr_qc_normalized are the same 8124-gene PCA universe; a broader ",
      "bulk protein-coding matrix is not on disk. Restricted to these HVGs ",
      "so the correlations describe the PCs that were actually computed."
    ),
    "data/expression/expr_cpm.rds (log2 CPM, not 55-PC residualized)",
    "mean raw counts per cell type, columns scaled to 1e6; genes = HVG ∩ SC",
    "NNLS (Lawson-Hanson) of linear 2^log2CPM on the signature; proportions sum to 1",
    sprintf("%.4f", s55$mean_adj_r2),
    sprintf("%.4f", s55$total_adj_r2),
    sprintf("%.4f", s2005$mean_adj_r2),
    sprintf("%.4f", s2005$total_adj_r2),
    if (isTRUE(past_peak)) "yes" else "no",
    sat_n,
    sprintf("%.4f", delta_after),
    paste(pc1_5$note, collapse = " | ")
  )
)

write_out(overlap_tab, "S7_gene_overlap_and_stats.tsv")
write_out(
  as_tibble(prop, rownames = "sample") |>
    mutate(in_kept_5489 = sample %in% readRDS(KEPT_SUBJECTS)),
  "S7_celltype_nnls_proportions.tsv"
)
write_out(
  as_tibble(r_mat, rownames = "pc") |>
    mutate(pc = as.integer(str_remove(pc, "^PC"))),
  "S7_pc_celltype_pearson_r.tsv"
)
write_out(per_pc, "S7_pc_multivar_r2.tsv")
write_out(cum, "S7_cumulative_composition_r2.tsv")
write_out(
  sc_cells |>
    mutate(display = unname(CT_LABEL[cell_type])),
  "S7_sc_cell_counts.tsv"
)
write_out(
  tibble(symbol = overlap_genes, in_globin = symbol %in% GLOBINS),
  "S7_deconvolution_genes.tsv"
)

mean_prop <- colMeans(prop)
cat("\n========== KEY STATS (TSV only; not on graphics) ==========\n")
cat("overlap genes:", n_overlap, " samples:", n_samples, "\n")
cat("mean NNLS proportions:\n")
print(round(mean_prop, 3))
cat("PC1-5:\n")
print(as.data.frame(pc1_5[, c("pc", "top_cell_type", "top_r", "max_abs_r", "multivar_r2")]))
cat(sprintf("cumulative mean adj R² at 55:   %.4f (total adj %.4f)\n",
            s55$mean_adj_r2, s55$total_adj_r2))
cat(sprintf("cumulative mean adj R² at 2005: %.4f (total adj %.4f)\n",
            s2005$mean_adj_r2, s2005$total_adj_r2))
cat("saturates by 55:", past_peak, " peak n<=55:", sat_n, "\n")
}

# ---------------------------------------------------------------------------
# F. S1a / S1b with composition overlay
# ---------------------------------------------------------------------------
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

# Overlay uses adjusted mean R² (each cell-type proportion ~ PC1..n).
overlay <- cum |>
  transmute(npcs_format = n_pcs, comp_r2 = mean_adj_r2)

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
    # Vertical guide: 1055 = z<-3 rare-PTV enrichment peak (paper PC).
    geom_vline(
      xintercept = 1055, linetype = "dashed",
      color = "grey35", linewidth = 0.45
    ) +
    annotate(
      "text", x = 1055, y = Inf, label = "1055 ", vjust = 1.4, hjust = 1,
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

# Panel c from the existing consequence script (unchanged).
s1c_env <- new.env(parent = globalenv())
options(paper.assemble = TRUE)
source(
  file.path(ROOT, "scripts/supplement/plot_ptv_consequence_summary.R"),
  local = s1c_env
)
options(paper.assemble = FALSE)
pc <- s1c_env$.paper_plots$c
paper_ggsave(pc, file.path(OUT, "FigS1c_ptv_consequences"),
             width = NATURE_DOUBLE_IN, height = 4.2)

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
s1_top <- cowplot::plot_grid(s1_panels, s1_leg, ncol = 1, rel_heights = c(1, 0.07))

# S1 is enrichment a|b only. Gene-z vs cell-type r stays on S8; PTV
# consequences is Supplementary Figure 10.
cat("Assembling Supplementary Figure 1 as a|b only; consequences -> S10\n")
paper_ggsave(pc, file.path(OUT, "FigS1c_ptv_consequences"),
             width = NATURE_DOUBLE_IN, height = 5.6)
paper_ggsave(pc, file.path(OUT, "FigS10a_ptv_consequences"),
             width = NATURE_DOUBLE_IN, height = 5.6)
paper_ggsave(pc, file.path(OUT, "Supplementary_Figure10"),
             width = NATURE_DOUBLE_IN, height = 5.6)
sfig1 <- s1_top
paper_ggsave(sfig1, file.path(OUT, "Supplementary_Figure1"),
             width = NATURE_DOUBLE_IN, height = 4.8)

# ---------------------------------------------------------------------------
# G. Supplementary Figure 7 — evaluation
# ---------------------------------------------------------------------------
n_heat <- min(K_HEAT, if (exists("k_use")) k_use else nrow(r_mat))
heat_df <- as_tibble(r_mat[seq_len(n_heat), , drop = FALSE], rownames = "pc") |>
  mutate(pc = as.integer(str_remove(pc, "^PC"))) |>
  pivot_longer(-pc, names_to = "cell_type", values_to = "r") |>
  mutate(
    cell_type = factor(cell_type, levels = ct_names),
    ct_lab = factor(
      unname(CT_LABEL[as.character(cell_type)]),
      levels = rev(unname(CT_LABEL[ct_names]))
    )
  )

p7a <- ggplot(heat_df, aes(pc, ct_lab, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(
    low = unname(OI["blue"]), mid = "white", high = unname(OI["vermillion"]),
    midpoint = 0, limits = c(-1, 1), name = "Pearson r",
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20, 25, 30, 35, 40),
    limits = c(0.5, n_heat + 0.5),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Expression PC", y = NULL) +
  theme_paper() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

p7b <- ggplot(per_pc, aes(pc, multivar_r2)) +
  geom_col(fill = unname(OI["blue"]), width = 0.85) +
  geom_vline(xintercept = 55, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  scale_x_continuous(breaks = c(1, 10, 20, 40, 55, 80), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Expression PC", y = expression(R^2 * " (PC ~ cell types)")) +
  theme_paper()

stack_df <- tibble(
  cell_type = factor(names(mean_prop), levels = rev(ct_names)),
  ct_lab = factor(unname(CT_LABEL[names(mean_prop)]),
                  levels = rev(unname(CT_LABEL[ct_names]))),
  proportion = as.numeric(mean_prop),
  dummy = "mean"
)
p7c <- ggplot(stack_df, aes(dummy, proportion, fill = ct_lab)) +
  geom_col(width = 0.45, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(
    values = setNames(unname(CT_COLS[ct_names]), unname(CT_LABEL[ct_names])),
    name = NULL, guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean NNLS proportion") +
  theme_paper() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

tag_s <- function(p) p + plot_annotation(tag_levels = "a") & theme_paper_tag()
sfig7 <- tag_s(
  p7a / p7b / p7c + plot_layout(heights = c(2.15, 1.45, 1.15))
)
paper_ggsave(sfig7, file.path(OUT, "Supplementary_Figure7"),
             width = NATURE_DOUBLE_IN, height = 7.6)
paper_ggsave(sfig7, file.path(OUT, "FigS_celltype_pc_eval"),
             width = NATURE_DOUBLE_IN, height = 7.6)
paper_ggsave(p7a, file.path(OUT, "FigS7a_pc_celltype_r_heatmap"),
             width = NATURE_DOUBLE_IN, height = 3.2)
paper_ggsave(p7b, file.path(OUT, "FigS7b_pc_multivar_r2"),
             width = NATURE_DOUBLE_IN, height = 2.4)
paper_ggsave(p7c, file.path(OUT, "FigS7c_mean_nnls_proportions"),
             width = NATURE_DOUBLE_IN, height = 2.0)

# Excel tabs
pack <- file.path(ROOT, "scripts/supplement/pack_supplementary_tables.py")
st <- system2(PY, shQuote(pack), stdout = "", stderr = "")
if (!identical(st, 0L)) warning("Excel packager failed: ", st)

manifest <- tibble(
  file = list.files(OUT, pattern = "Supplementary_Figure|FigS", full.names = FALSE),
  kind = "figure"
) |>
  bind_rows(tibble(
    file = file.path("supplemental_tables", list.files(TSV)),
    kind = "table"
  ))
write_out(manifest, "MANIFEST.tsv")

cat("\nDone. Figures in ", OUT, "\n", sep = "")
