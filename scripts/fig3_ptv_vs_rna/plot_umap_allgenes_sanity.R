#!/usr/bin/env Rscript
# Sample-scattering UMAP sanity check for the all-gene (coding+noncoding) HVG
# pipeline, matching Supplementary Figure 3 layout:
#   a  sequencing centre, before vs after expression-PC residualization
#   b  study cohort, same facets
#
# BEFORE PCs: first 30 of expr_pc_df.rds (PCAtools::pca on qc_normalized,
# scale=TRUE — same object the pipeline already cached).
# AFTER PCs: removeBatchEffect(qc_normalized, first 55 expr PCs) then
# PCAtools::pca(scale=TRUE, center=TRUE), as in outlier_processing.R /
# make_supplementary_s2_s5.R. Cached under allgenes_pipeline_out/.
#
# Does not overwrite published S3 or data/expression/expr_cpm_adj.rds.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(uwot)
  library(limma)
  library(PCAtools)
  library(tibble)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

OUT <- file.path(ROOT, "allgenes_pipeline_out")
DRAFT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
STEM <- "Fig_hvg_allgenes_umap"
N_PCS_UMAP <- 30L
N_EXPR_PCS <- 55L
set.seed(1)

COHORT_COLS <- c(
  GALAII     = unname(OI["sky"]),
  SPIROMICS  = unname(OI["blue"]),
  WHI        = unname(OI["vermillion"]),
  FHS        = unname(OI["green"]),
  SAGE       = unname(OI["orange"]),
  COPDGene   = unname(OI["purple"])
)
CENTER_COLS <- c(Broad = unname(OI["blue"]), NWGC = unname(OI["orange"]))
tag_s <- function(p) p + plot_annotation(tag_levels = "a") & theme_paper_tag()

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

stopifnot(dir.exists(OUT))
meta <- read_tsv(METADATA, show_col_types = FALSE)

# ---------------------------------------------------------------------------
# BEFORE: sample scores already stored by the all-gene pipeline
# ---------------------------------------------------------------------------
msg("Loading expr_pc_df (before residualization) ...")
expr_pc_df <- readRDS(file.path(OUT, "expr_pc_df.rds"))
pc_names <- paste0("PC", seq_len(max(N_PCS_UMAP, N_EXPR_PCS)))
missing_pc <- setdiff(pc_names, names(expr_pc_df))
if (length(missing_pc)) stop("expr_pc_df missing columns: ", paste(missing_pc, collapse = ", "))
id_col <- if ("ID" %in% names(expr_pc_df)) "ID" else names(expr_pc_df)[1]
rot_b <- as.matrix(expr_pc_df[, pc_names, drop = FALSE])
rownames(rot_b) <- as.character(expr_pc_df[[id_col]])
msg("  before PCs: ", nrow(rot_b), " samples x ", ncol(rot_b), " PCs")

# ---------------------------------------------------------------------------
# AFTER: residualize Freeman-Tukey HVG on first 55 expr PCs, then PCA
# ---------------------------------------------------------------------------
after_cache <- file.path(OUT, "umap_sanity_pca_after_rotated30.rds")
if (file.exists(after_cache)) {
  msg("Loading cached after-residualization PCs ...")
  rot_a <- readRDS(after_cache)
} else {
  msg("Loading qc_normalized ...")
  qc_normalized <- readRDS(file.path(OUT, "qc_normalized.rds"))
  shared <- intersect(colnames(qc_normalized), rownames(rot_b))
  if (length(shared) < 2L) stop("No shared samples between qc_normalized and expr_pc_df")
  msg(
    "  qc_normalized: ", nrow(qc_normalized), " genes x ", ncol(qc_normalized),
    " samples; shared=", length(shared)
  )
  msg("removeBatchEffect on first ", N_EXPR_PCS, " expr PCs ...")
  covar <- rot_b[shared, paste0("PC", seq_len(N_EXPR_PCS)), drop = FALSE]
  storage.mode(covar) <- "double"
  adj <- removeBatchEffect(as.matrix(qc_normalized[, shared, drop = FALSE]), covariates = covar)
  rm(qc_normalized)
  gc(verbose = FALSE)
  msg("PCAtools::pca on residualized matrix (scale=TRUE, center=TRUE) ...")
  pca_after <- pca(adj, scale = TRUE, center = TRUE)
  rot_a <- as.matrix(pca_after$rotated[, seq_len(N_PCS_UMAP), drop = FALSE])
  if (is.null(rownames(rot_a)) || !length(rownames(rot_a))) {
    rownames(rot_a) <- colnames(adj)
  }
  rm(adj, pca_after, covar)
  gc(verbose = FALSE)
  saveRDS(rot_a, after_cache)
  msg("  cached ", after_cache)
}

align_meta <- function(rot) {
  ids <- rownames(rot)
  md <- meta |>
    filter(WGS %in% ids) |>
    distinct(WGS, .keep_all = TRUE)
  keep <- ids[ids %in% md$WGS]
  list(rot = rot[keep, , drop = FALSE], meta = md[match(keep, md$WGS), ])
}
ab <- align_meta(rot_b)
aa <- align_meta(rot_a)
msg("  aligned before=", nrow(ab$rot), " after=", nrow(aa$rot))

# ---------------------------------------------------------------------------
# UMAP (same uwot settings / seed as S3)
# ---------------------------------------------------------------------------
msg("UMAP on first ", N_PCS_UMAP, " PCs ...")
umap_one <- function(rot, md, stage) {
  coords <- uwot::umap(
    as.matrix(rot[, seq_len(N_PCS_UMAP), drop = FALSE]),
    n_neighbors = 15, min_dist = 0.3, metric = "euclidean",
    verbose = FALSE, seed = 1
  )
  tibble(
    WGS = rownames(rot),
    UMAP_1 = coords[, 1],
    UMAP_2 = coords[, 2],
    stage = stage
  ) |>
    left_join(md |> select(WGS, COHORT, SEQUENCING_CENTER, SEX), by = "WGS")
}
umap_df <- bind_rows(
  umap_one(ab$rot, ab$meta, "Before residualization"),
  umap_one(aa$rot, aa$meta, "After residualization")
)
umap_df$stage <- factor(umap_df$stage, levels = c("Before residualization", "After residualization"))
umap_df$COHORT <- factor(umap_df$COHORT, levels = names(COHORT_COLS))
write_tsv(umap_df, file.path(OUT, "Fig_hvg_allgenes_umap_coords.tsv"))

# ---------------------------------------------------------------------------
# Plot (S3 layout)
# ---------------------------------------------------------------------------
umap_pt <- function(df, color_aes, cols, legend_name) {
  ggplot(df, aes(UMAP_1, UMAP_2, color = .data[[color_aes]])) +
    geom_point(size = 0.25, alpha = 0.45, stroke = 0) +
    scale_color_manual(values = cols, name = legend_name, drop = FALSE) +
    facet_wrap(~ stage, nrow = 1) +
    coord_equal() +
    theme_paper() +
    theme(
      legend.position = "bottom",
      legend.key.size = grid::unit(8, "pt"),
      strip.text = element_text(size = 8, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1)) +
    labs(x = "UMAP 1", y = "UMAP 2")
}

p_center <- umap_pt(umap_df, "SEQUENCING_CENTER", CENTER_COLS, "Centre")
p_coh <- umap_pt(umap_df, "COHORT", COHORT_COLS, "Study")
fig <- tag_s(p_center / p_coh)
paper_ggsave(fig, file.path(OUT, STEM), width = NATURE_DOUBLE_IN, height = 7.6)

if (dir.exists(DRAFT)) {
  file.copy(paste0(file.path(OUT, STEM), ".png"), file.path(DRAFT, paste0(STEM, ".png")), overwrite = TRUE)
  file.copy(paste0(file.path(OUT, STEM), ".pdf"), file.path(DRAFT, paste0(STEM, ".pdf")), overwrite = TRUE)
}

# ---------------------------------------------------------------------------
# One-line structure diagnostic (UMAP-coordinate R2)
# ---------------------------------------------------------------------------
diag_df <- bind_rows(lapply(levels(umap_df$stage), function(st) {
  d <- umap_df[umap_df$stage == st, , drop = FALSE]
  tibble(
    stage = st,
    n = nrow(d),
    r2_umap1_center = summary(lm(UMAP_1 ~ SEQUENCING_CENTER, d))$r.squared,
    r2_umap2_center = summary(lm(UMAP_2 ~ SEQUENCING_CENTER, d))$r.squared,
    r2_umap1_cohort = summary(lm(UMAP_1 ~ COHORT, d))$r.squared,
    r2_umap2_cohort = summary(lm(UMAP_2 ~ COHORT, d))$r.squared
  )
}))
write_tsv(diag_df, file.path(OUT, "Fig_hvg_allgenes_umap_structure_r2.tsv"))
print(diag_df, n = 10)
msg("wrote ", file.path(OUT, paste0(STEM, ".png")))
