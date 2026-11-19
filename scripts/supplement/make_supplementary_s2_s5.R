#!/usr/bin/env Rscript
# Supplementary Figures 2–6 + intermediate TSVs in
# final_paper_draft/supp_figs/
#
#   S2  RNA-seq QC (connectivity, extreme-subject filter, scree, cohort tests)
#   S3  UMAP before vs after expression-PC residualization (batch)
#   S4  TOPMed vs AoU underlier-count correlation at each cutoff
#   S5  z-cutoff landscape (n genes/pairs + LOEUF calibration)
#
# Also writes gene lists and every intermediate as TSV under
# supp_figs/tables/.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(uwot)
  library(vegan)
  library(Matrix)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

OUT <- file.path(FINAL_DRAFT, "supp_figs")
TSV <- file.path(OUT, "tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TSV, recursive = TRUE, showWarnings = FALSE)

N_TOPMED <- 5489L
N_AOU <- 8895L
N_PCS_UMAP <- 30L
N_PCS_TEST <- 10L
set.seed(1)

wilson_ci <- function(k, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- ifelse(n > 0, k / n, NA_real_)
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
  tibble(p = p, lo = pmax(centre - half, 0), hi = pmin(centre + half, 1))
}

save_panel <- function(plot, stem, width, height) {
  paper_ggsave(plot, file.path(OUT, stem), width = width, height = height)
}

save_comp <- function(plot, name, width, height) {
  paper_ggsave(plot, file.path(OUT, name), width = width, height = height)
}

write_out <- function(df, name) {
  path <- file.path(TSV, name)
  write_tsv(df, path)
  cat("  wrote ", name, " (", nrow(df), " rows)\n", sep = "")
  invisible(path)
}

nv <- function(x) str_remove(x, "\\.\\d+$")
tag_s <- function(p) p + plot_annotation(tag_levels = "a") & theme_paper_tag()

COHORT_COLS <- c(
  GALAII     = unname(OI["sky"]),
  SPIROMICS  = unname(OI["blue"]),
  WHI        = unname(OI["vermillion"]),
  FHS        = unname(OI["green"]),
  SAGE       = unname(OI["orange"]),
  COPDGene   = unname(OI["purple"])
)
CENTER_COLS <- c(Broad = unname(OI["blue"]), NWGC = unname(OI["orange"]))
CUT_LEVELS <- c("haplo", "z < -2", "z < -3", "z < -4")
CUT_COLS <- c(
  haplo      = unname(OI["sky"]),
  `z < -2`   = unname(OI["purple"]),
  `z < -3`   = unname(OI["vermillion"]),
  `z < -4`   = unname(OI["green"])
)

# ---------------------------------------------------------------------------
# Move Supplementary Figure 1 into supp_figs (canonical home)
# ---------------------------------------------------------------------------
cat("Moving Supplementary Figure 1 into supp_figs ...\n")
for (ext in c("png", "pdf")) {
  for (src in c(
    file.path(FINAL_DRAFT, paste0("Supplementary_Figure1.", ext)),
    file.path(FINAL_DRAFT, "source_panels", paste0("Supplementary_Figure1.", ext))
  )) {
    dst <- file.path(OUT, paste0("Supplementary_Figure1.", ext))
    if (file.exists(src)) {
      if (!file.exists(dst) || !identical(unname(tools::md5sum(src)), unname(tools::md5sum(dst)))) {
        file.copy(src, dst, overwrite = TRUE)
      }
      if (!identical(normalizePath(src, mustWork = FALSE), normalizePath(dst, mustWork = TRUE))) {
        unlink(src)
      }
    }
  }
}
# S1 source tables
s1_map <- c(
  "Source_Data_S1ab_enrichment_or.tsv" = file.path(SUP_PC, "enrichment_npc_sv_vs_ptv_underlier.tsv"),
  "Source_Data_S1ab_enrichment_bootstrap.tsv" = file.path(SUP_PC, "enrichment_peak_vs_gd_bootstrap.tsv"),
  "Source_Data_S1c_ptv_consequences.tsv" = file.path(SUP_PTV, "ptv_consequence_counts_qc_cohort.tsv")
)
for (nm in names(s1_map)) {
  if (file.exists(s1_map[[nm]])) file.copy(s1_map[[nm]], file.path(TSV, nm), overwrite = TRUE)
}

# ---------------------------------------------------------------------------
# Shared annotations
# ---------------------------------------------------------------------------
cat("Loading annotations ...\n")
kept_subjects <- readRDS(KEPT_SUBJECTS)
meta <- read_tsv(METADATA, show_col_types = FALSE)
loeuf <- readRDS(LOEUF_LOOKUP) |>
  distinct(gene_id, .keep_all = TRUE) |>
  transmute(gene_id, gene_nv = gene_id_no_version, symbol, LOEUF, LOEUF_tier)
sym <- read_tsv(GENCODE_SYM, col_names = c("gene_nv", "symbol_gencode"), show_col_types = FALSE)

extreme <- read_tsv(file.path(DATA_CORE, "subject_extreme_outlier_counts.tsv"), show_col_types = FALSE)
write_out(extreme, "S2_subject_extreme_outlier_counts.tsv")

n_drop_ext <- sum(!extreme$keep)
n_pre_ext <- nrow(extreme)
prop_ci <- wilson_ci(n_drop_ext, n_pre_ext)

SKIP_S23 <- file.exists(file.path(OUT, "Supplementary_Figure2.png")) &&
  file.exists(file.path(OUT, "Supplementary_Figure3.png"))
if (SKIP_S23) {
  cat("S2 and S3 already in supp_figs; skipping rebuild\n")
} else {
# ---------------------------------------------------------------------------
# S2a connectivity Z (recompute from expr_normalized, same rule as pipeline)
# ---------------------------------------------------------------------------
cat("S2: connectivity Z-scores ...\n")
conn_path <- file.path(TSV, "S2_connectivity_z.tsv")
if (file.exists(conn_path)) {
  conn_df <- read_tsv(conn_path, show_col_types = FALSE)
} else {
  expr_normalized <- readRDS(file.path(DATA_EXPR, "expr_normalized.rds"))
  mat <- as.matrix(expr_normalized)
  rm(expr_normalized); gc(verbose = FALSE)
  cat("  cor() on", ncol(mat), "samples x", nrow(mat), "genes ...\n")
  cormat <- cor(mat)
  rm(mat); gc(verbose = FALSE)
  adj <- 0.5 + 0.5 * cormat
  rm(cormat); gc(verbose = FALSE)
  conn <- WGCNA::fundamentalNetworkConcepts(adj)$Connectivity
  rm(adj); gc(verbose = FALSE)
  conn_z <- as.numeric(scale(conn))
  conn_df <- tibble(
    subjectid = names(conn),
    connectivity = as.numeric(conn),
    connectivity_z = conn_z,
    drop_connectivity = conn_z < -3
  )
  write_out(conn_df, "S2_connectivity_z.tsv")
}
n_drop_conn <- sum(conn_df$drop_connectivity, na.rm = TRUE)
conn_keep <- conn_df$connectivity_z[conn_df$connectivity_z >= -3]
# Two-sided Shapiro is too sensitive at n~5k; report KS vs N(0,1) on retained samples
ks_conn <- ks.test(as.numeric(scale(conn_keep)), "pnorm", exact = FALSE)

p_s2a <- ggplot(conn_df, aes(connectivity_z, fill = drop_connectivity)) +
  geom_histogram(binwidth = 0.2, color = "white", linewidth = 0.1) +
  geom_vline(xintercept = -3, linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(
    values = c(`TRUE` = unname(OI["vermillion"]), `FALSE` = "grey70"),
    labels = c(`TRUE` = "Dropped (Z < -3)", `FALSE` = "Retained"),
    name = NULL
  ) +
  labs(x = "Sample connectivity Z", y = "Subjects") +
  theme_paper() +
  theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# S2b extreme |Z|>5 subject filter
# ---------------------------------------------------------------------------
p_s2b <- ggplot(extreme, aes(n_extreme_outliers, fill = keep)) +
  geom_histogram(binwidth = 5, color = "white", linewidth = 0.1) +
  geom_vline(xintercept = 60, linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(
    values = c(`TRUE` = "grey70", `FALSE` = unname(OI["vermillion"])),
    labels = c(`TRUE` = "Kept", `FALSE` = "Dropped (>= 60)"),
    name = NULL
  ) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 5, 20, 60, 200)) +
  labs(x = "Events at |Z| > 5 per subject", y = "Subjects") +
  theme_paper() +
  theme(legend.position = "bottom")

# ---------------------------------------------------------------------------
# S2c scree + S3 UMAP + MANOVA/PERMANOVA from PCA objects
# ---------------------------------------------------------------------------
cat("Loading expr_pca (before residualization) ...\n")
pca_before <- readRDS(file.path(DATA_EXPR, "expr_pca.rds"))
rot_b <- pca_before$rotated[, seq_len(max(N_PCS_UMAP, N_PCS_TEST)), drop = FALSE]
var_b <- as.numeric(pca_before$variance)
rm(pca_before); gc(verbose = FALSE)

cat("Loading expr_norm_adj_pca (after expr-PC residualization) ...\n")
pca_after <- readRDS(file.path(DATA_EXPR, "expr_norm_adj_pca.rds"))
rot_a <- pca_after$rotated[, seq_len(max(N_PCS_UMAP, N_PCS_TEST)), drop = FALSE]
var_a <- as.numeric(pca_after$variance)
rm(pca_after); gc(verbose = FALSE)

n_comp_show <- 80L
scree_df <- bind_rows(
  tibble(pc = seq_len(n_comp_show), variance_pct = var_b[seq_len(n_comp_show)], stage = "Before residualization"),
  tibble(pc = seq_len(n_comp_show), variance_pct = var_a[seq_len(n_comp_show)], stage = "After residualization")
)
scree_df$stage <- factor(scree_df$stage, levels = c("Before residualization", "After residualization"))
write_out(scree_df, "S2_pca_scree_variance.tsv")
cum_b_55 <- sum(var_b[1:55], na.rm = TRUE)
cum_a_55 <- sum(var_a[1:55], na.rm = TRUE)

p_s2c <- ggplot(scree_df, aes(pc, variance_pct, color = stage)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = 55, linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(values = c(
    "Before residualization" = unname(OI["blue"]),
    "After residualization" = unname(OI["vermillion"])
  ), name = NULL) +
  labs(x = "Expression PC", y = "Variance explained (%)") +
  theme_paper() +
  theme(legend.position = "bottom")

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

run_tests <- function(rot, md, stage) {
  pcs <- as.matrix(rot[, seq_len(N_PCS_TEST), drop = FALSE])
  storage.mode(pcs) <- "double"
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  md$COHORT <- factor(md$COHORT)
  md$SEQUENCING_CENTER <- factor(md$SEQUENCING_CENTER)
  man_coh <- summary(manova(pcs ~ md$COHORT), test = "Pillai")
  man_cen <- summary(manova(pcs ~ md$SEQUENCING_CENTER), test = "Pillai")
  pillai <- function(s) {
    st <- s$stats
    tibble(pillai = unname(st[1, "Pillai"]), p_value = unname(st[1, "Pr(>F)"]))
  }
  set.seed(1)
  ad_coh <- adonis2(pcs ~ COHORT, data = md, permutations = 199, method = "euclidean")
  set.seed(1)
  ad_cen <- adonis2(pcs ~ SEQUENCING_CENTER, data = md, permutations = 199, method = "euclidean")
  r2_pc1_coh <- summary(lm(pcs[, 1] ~ md$COHORT))$r.squared
  r2_pc1_cen <- summary(lm(pcs[, 1] ~ md$SEQUENCING_CENTER))$r.squared
  bind_rows(
    tibble(stage = stage, grouping = "COHORT", test = "MANOVA Pillai (PC1-10)",
           statistic = pillai(man_coh)$pillai, r2 = NA_real_, p_value = pillai(man_coh)$p_value),
    tibble(stage = stage, grouping = "SEQUENCING_CENTER", test = "MANOVA Pillai (PC1-10)",
           statistic = pillai(man_cen)$pillai, r2 = NA_real_, p_value = pillai(man_cen)$p_value),
    tibble(stage = stage, grouping = "COHORT", test = "PERMANOVA (PC1-10, 199 perm)",
           statistic = NA_real_, r2 = ad_coh$R2[1], p_value = ad_coh$`Pr(>F)`[1]),
    tibble(stage = stage, grouping = "SEQUENCING_CENTER", test = "PERMANOVA (PC1-10, 199 perm)",
           statistic = NA_real_, r2 = ad_cen$R2[1], p_value = ad_cen$`Pr(>F)`[1]),
    tibble(stage = stage, grouping = "COHORT", test = "lm R2 of PC1",
           statistic = NA_real_, r2 = r2_pc1_coh, p_value = anova(lm(pcs[, 1] ~ md$COHORT))$`Pr(>F)`[1]),
    tibble(stage = stage, grouping = "SEQUENCING_CENTER", test = "lm R2 of PC1",
           statistic = NA_real_, r2 = r2_pc1_cen, p_value = anova(lm(pcs[, 1] ~ md$SEQUENCING_CENTER))$`Pr(>F)`[1])
  )
}

cat("S2: MANOVA / PERMANOVA ...\n")
qc_tests <- bind_rows(
  run_tests(ab$rot, ab$meta, "Before residualization"),
  run_tests(aa$rot, aa$meta, "After residualization")
) |>
  mutate(
    extreme_subjects_dropped = n_drop_ext,
    extreme_subjects_tested = n_pre_ext,
    extreme_drop_prop = prop_ci$p[1],
    extreme_drop_lo = prop_ci$lo[1],
    extreme_drop_hi = prop_ci$hi[1],
    connectivity_dropped = n_drop_conn,
    ks_connectivity_D = unname(ks_conn$statistic),
    ks_connectivity_p = ks_conn$p.value,
    gavish_donoho_k = 55L,
    variance_pct_pc1_55_before = cum_b_55,
    variance_pct_pc1_55_after = cum_a_55
  )
write_out(qc_tests, "S2_qc_statistical_tests.tsv")

perm_plot_df <- qc_tests |>
  filter(test == "PERMANOVA (PC1-10, 199 perm)") |>
  mutate(
    stage = factor(stage, levels = c("Before residualization", "After residualization")),
    grouping = recode(grouping, SEQUENCING_CENTER = "Sequencing centre", COHORT = "Study")
  )
p_s2d <- ggplot(perm_plot_df, aes(r2, grouping, color = stage)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
  scale_color_manual(values = c(
    "Before residualization" = unname(OI["blue"]),
    "After residualization" = unname(OI["vermillion"])
  ), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(x = "PERMANOVA R\u00b2 (first 10 PCs)", y = NULL) +
  theme_paper() +
  theme(legend.position = "bottom")

fig_s2 <- tag_s(
  (p_s2a | p_s2b) / (p_s2c | p_s2d) + plot_layout(heights = c(1, 1))
)
save_comp(fig_s2, "Supplementary_Figure2", NATURE_DOUBLE_IN, 7.2)
save_panel(p_s2a, "FigS2a_connectivity_z", NATURE_DOUBLE_IN / 2, 3.4)
save_panel(p_s2b, "FigS2b_extreme_subject_filter", NATURE_DOUBLE_IN / 2, 3.4)
save_panel(p_s2c, "FigS2c_pca_scree", NATURE_DOUBLE_IN / 2, 3.4)
save_panel(p_s2d, "FigS2d_permanova_r2", NATURE_DOUBLE_IN / 2, 3.4)
cat("wrote Supplementary_Figure2\n")

# ---------------------------------------------------------------------------
# S3 UMAP
# ---------------------------------------------------------------------------
cat("S3: UMAP ...\n")
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
write_out(umap_df, "S3_umap_coords_before_after.tsv")

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

p_s3_center <- umap_pt(umap_df, "SEQUENCING_CENTER", CENTER_COLS, "Centre")
p_s3_coh <- umap_pt(umap_df, "COHORT", COHORT_COLS, "Study")
fig_s3 <- tag_s(p_s3_center / p_s3_coh)
save_comp(fig_s3, "Supplementary_Figure3", NATURE_DOUBLE_IN, 7.6)
save_panel(p_s3_center, "FigS3a_umap_sequencing_center", NATURE_DOUBLE_IN, 3.6)
save_panel(p_s3_coh, "FigS3b_umap_cohort", NATURE_DOUBLE_IN, 3.6)
cat("wrote Supplementary_Figure3\n")
rm(rot_b, rot_a, ab, aa); gc(verbose = FALSE)
}  # end SKIP_S23 else

# ---------------------------------------------------------------------------
# Gene-level counts at each cutoff (TOPMed + AoU) + pair tables
# ---------------------------------------------------------------------------
cat("Building gene-level underlier counts ...\n")
tally_of <- function(label, rds) {
  readRDS(rds) |>
    transmute(gene_id, gene_nv = nv(gene_id), n_topmed = as.integer(n), cutoff = label)
}
tallies <- bind_rows(
  tally_of("haplo", file.path(DATA_UNDERLIER, "z_score_tally_haplo_clean.rds")),
  tally_of("z < -2", file.path(DATA_UNDERLIER, "z_score_tally_z_-2_clean.rds")),
  tally_of("z < -3", file.path(DATA_UNDERLIER, "z_score_tally_z_-3_clean.rds")),
  tally_of("z < -4", file.path(DATA_UNDERLIER, "z_score_tally_z_-4_clean.rds"))
)

AOU_PREV <- file.path(ROOT, "call-RunRnaUnderlier", "glob-833705bb57df602e29e5b680691e557c")
aou_of <- function(label, f) {
  x <- read_tsv(file.path(AOU_PREV, f), show_col_types = FALSE)
  if (!"gene_nv" %in% names(x)) x$gene_nv <- nv(x$gene_id)
  x |>
    transmute(
      gene_nv,
      n_aou = as.integer(n_underlier_subjects),
      cutoff = label
    )
}
aou_tallies <- bind_rows(
  aou_of("haplo", "rna_outlier_prevalence_per_gene_haplo.tsv"),
  aou_of("z < -2", "rna_outlier_prevalence_per_gene_z_-2.tsv"),
  aou_of("z < -3", "rna_outlier_prevalence_per_gene_z_-3.tsv"),
  aou_of("z < -4", "rna_outlier_prevalence_per_gene_z_-4.tsv")
)

# Eligible TOPMed gene universe
cat("  eligible gene universe from expr_cpm_adj ...\n")
expr_cpm_adj <- readRDS(EXPR_CPM_ADJ)
eligible <- tibble(gene_id = rownames(expr_cpm_adj), gene_nv = nv(rownames(expr_cpm_adj)))
rm(expr_cpm_adj); gc(verbose = FALSE)

ann <- eligible |>
  left_join(loeuf |> select(gene_id, symbol, LOEUF, LOEUF_tier), by = "gene_id") |>
  left_join(sym, by = "gene_nv") |>
  mutate(symbol = coalesce(symbol, symbol_gencode)) |>
  select(gene_id, gene_nv, symbol, LOEUF, LOEUF_tier)

gene_counts <- tidyr::expand_grid(gene_nv = unique(c(ann$gene_nv, aou_tallies$gene_nv)), cutoff = CUT_LEVELS) |>
  left_join(ann |> distinct(gene_nv, .keep_all = TRUE), by = "gene_nv") |>
  left_join(tallies |> select(gene_nv, cutoff, n_topmed), by = c("gene_nv", "cutoff")) |>
  left_join(aou_tallies, by = c("gene_nv", "cutoff")) |>
  mutate(
    n_topmed = replace_na(n_topmed, 0L),
    n_aou = replace_na(n_aou, 0L),
    cutoff = factor(cutoff, levels = CUT_LEVELS)
  )
write_out(gene_counts, "S_gene_underlier_counts_by_cutoff.tsv")

# Wide gene list (one row per gene)
gene_wide <- gene_counts |>
  select(gene_id, gene_nv, symbol, LOEUF, LOEUF_tier, cutoff, n_topmed, n_aou) |>
  pivot_wider(
    names_from = cutoff,
    values_from = c(n_topmed, n_aou),
    names_glue = "{.value}_{cutoff}"
  )
# clean names
names(gene_wide) <- names(gene_wide) |>
  str_replace_all("n_topmed_", "n_topmed_") |>
  str_replace_all("z < -", "z_") |>
  str_replace_all(" ", "")
write_out(gene_wide, "S_gene_list_all_cutoffs.tsv")

# Pair tables (TOPMed)
export_pairs <- function(label, rds, out_name) {
  cat("  pairs ", label, " ...\n", sep = "")
  x <- readRDS(rds) |>
    transmute(
      gene_id, gene_nv = nv(gene_id), subjectid,
      expression_zscore, expression_logcpm, cutoff = label
    )
  write_tsv(x, file.path(TSV, out_name))
  cat("  wrote ", out_name, " (", nrow(x), " rows)\n", sep = "")
  invisible(nrow(x))
}
n_pairs <- tibble(
  cutoff = CUT_LEVELS,
  n_pairs_topmed = c(
    export_pairs("haplo", file.path(DATA_UNDERLIER, "z_score_underlier_haplo_clean.rds"), "S_underlier_pairs_topmed_haplo.tsv.gz"),
    export_pairs("z < -2", file.path(DATA_UNDERLIER, "z_score_underlier_z_-2_clean.rds"), "S_underlier_pairs_topmed_z-2.tsv.gz"),
    export_pairs("z < -3", file.path(DATA_UNDERLIER, "z_score_underlier_z_-3_clean.rds"), "S_underlier_pairs_topmed_z-3.tsv.gz"),
    export_pairs("z < -4", file.path(DATA_UNDERLIER, "z_score_underlier_z_-4_clean.rds"), "S_underlier_pairs_topmed_z-4.tsv.gz")
  )
)

# AoU pair lists (already gzipped) — copy
aou_glob_pairs <- file.path(ROOT, "call-RunRnaUnderlier", "glob-6e53194078465641cb773324937722ad")
for (lab in c("haplo", "z_-2", "z_-3", "z_-4")) {
  src <- file.path(aou_glob_pairs, paste0("underliers_", lab, ".tsv.gz"))
  if (!file.exists(src) && lab != "haplo") src <- file.path(aou_glob_pairs, paste0("underliers_", lab, ".tsv.gz"))
  if (file.exists(src)) {
    file.copy(src, file.path(TSV, paste0("S_underlier_pairs_aou_", lab, ".tsv.gz")), overwrite = TRUE)
  }
}

# ---------------------------------------------------------------------------
# S4 TOPMed vs AoU count correlation
# ---------------------------------------------------------------------------
cat("S4: count correlation ...\n")
cor_one <- function(df) {
  both1 <- df |> filter(n_topmed >= 1, n_aou >= 1)
  both2 <- df |> filter(n_topmed >= 2, n_aou >= 2)
  sp1 <- if (nrow(both1) >= 3) cor.test(both1$n_topmed, both1$n_aou, method = "spearman", exact = FALSE) else NULL
  pe1 <- if (nrow(both1) >= 3) cor.test(log10(both1$n_topmed), log10(both1$n_aou), method = "pearson") else NULL
  sp2 <- if (nrow(both2) >= 3) cor.test(both2$n_topmed, both2$n_aou, method = "spearman", exact = FALSE) else NULL
  tibble(
    n_genes_union = nrow(df),
    n_genes_ge1_both = nrow(both1),
    n_genes_ge2_both = nrow(both2),
    spearman_ge1 = if (is.null(sp1)) NA_real_ else unname(sp1$estimate),
    spearman_ge1_p = if (is.null(sp1)) NA_real_ else sp1$p.value,
    pearson_log10_ge1 = if (is.null(pe1)) NA_real_ else unname(pe1$estimate),
    pearson_log10_ge1_p = if (is.null(pe1)) NA_real_ else pe1$p.value,
    spearman_ge2 = if (is.null(sp2)) NA_real_ else unname(sp2$estimate),
    spearman_ge2_p = if (is.null(sp2)) NA_real_ else sp2$p.value
  )
}
cor_stats <- gene_counts |>
  group_by(cutoff) |>
  group_modify(~ cor_one(.x)) |>
  ungroup()
write_out(cor_stats, "S4_topmed_vs_aou_count_correlation.tsv")
write_out(gene_counts, "S4_topmed_vs_aou_counts_long.tsv")

p_s4 <- ggplot(gene_counts, aes(n_topmed + 1, n_aou + 1)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_point(size = 0.45, alpha = 0.35, color = unname(OI["blue"])) +
  scale_x_log10(breaks = c(1, 2, 11, 101), labels = c("0", "1", "10", "100")) +
  scale_y_log10(breaks = c(1, 2, 11, 101), labels = c("0", "1", "10", "100")) +
  facet_wrap(~ cutoff, nrow = 2) +
  labs(x = "TOPMed underlier count + 1", y = "AoU underlier count + 1") +
  theme_paper() +
  theme(strip.text = element_text(size = 8, face = "bold"))
fig_s4 <- tag_s(p_s4)
save_comp(fig_s4, "Supplementary_Figure4", NATURE_DOUBLE_IN, 6.8)
save_panel(p_s4, "FigS4_topmed_vs_aou_counts", NATURE_DOUBLE_IN, 6.6)
cat("wrote Supplementary_Figure4\n")

# ---------------------------------------------------------------------------
# S5 z-cutoff landscape
# ---------------------------------------------------------------------------
cat("S5: cutoff landscape ...\n")
subj_tally <- function(label, rds) {
  x <- readRDS(rds)
  tibble(
    cutoff = label,
    n_genes_with_ge1 = nrow(x),
    n_subjects_with_ge1 = NA_integer_,
    median_per_subject = NA_real_
  )
}
subj_files <- c(
  haplo = "z_score_tally_subj_haplo_clean.rds",
  `z < -2` = "z_score_tally_subj_z_-2_clean.rds",
  `z < -3` = "z_score_tally_subj_z_-3_clean.rds",
  `z < -4` = "z_score_tally_subj_z_-4_clean.rds"
)
subj_sum <- lapply(names(subj_files), function(lab) {
  x <- readRDS(file.path(DATA_UNDERLIER, subj_files[[lab]]))
  tibble(
    cutoff = lab,
    n_subjects_with_ge1 = nrow(x),
    median_underliers_per_subject = median(x$n),
    mean_underliers_per_subject = mean(x$n)
  )
}) |> bind_rows()

cutoff_summary <- gene_counts |>
  group_by(cutoff) |>
  summarise(
    n_genes_topmed_ge1 = sum(n_topmed >= 1),
    n_genes_aou_ge1 = sum(n_aou >= 1),
    n_pairs_topmed = sum(n_topmed),
    n_pairs_aou = sum(n_aou),
    .groups = "drop"
  ) |>
  left_join(subj_sum, by = "cutoff") |>
  left_join(n_pairs |> rename(n_pairs_topmed_from_calls = n_pairs_topmed), by = "cutoff") |>
  left_join(cor_stats |> select(cutoff, spearman_ge1, spearman_ge1_p, spearman_ge2), by = "cutoff")
write_out(cutoff_summary, "S5_cutoff_landscape_summary.tsv")

count_long <- cutoff_summary |>
  select(cutoff, n_genes_topmed_ge1, n_pairs_topmed) |>
  pivot_longer(-cutoff, names_to = "metric", values_to = "value") |>
  mutate(
    metric = recode(metric,
      n_genes_topmed_ge1 = "Genes with >= 1 underlier",
      n_pairs_topmed = "Subject-gene pairs"
    ),
    cutoff = factor(cutoff, levels = CUT_LEVELS)
  )
p_s5a <- ggplot(count_long, aes(cutoff, value, fill = cutoff)) +
  geom_col(width = 0.7) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = CUT_COLS, guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = NULL) +
  theme_paper() +
  theme(strip.text = element_text(size = 8, face = "bold"))

# LOEUF forests already computed
loeuf_files <- c(
  haplo = file.path(FIG2A, "underlier_loeuf_calibration_by_tier_haplo.tsv"),
  `z < -2` = file.path(FIG2A, "underlier_loeuf_calibration_by_tier_z_-2.tsv"),
  `z < -3` = file.path(FIG2A, "underlier_loeuf_calibration_by_tier_z_-3.tsv"),
  `z < -4` = file.path(FIG2A, "underlier_loeuf_calibration_by_tier_z_-4.tsv")
)
loeuf_all <- lapply(names(loeuf_files), function(lab) {
  read_tsv(loeuf_files[[lab]], show_col_types = FALSE) |>
    mutate(cutoff = lab)
}) |> bind_rows()
loeuf_all$cutoff <- factor(loeuf_all$cutoff, levels = CUT_LEVELS)
loeuf_all$LOEUF_tier <- factor(
  loeuf_all$LOEUF_tier,
  levels = c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
)
write_out(loeuf_all, "S5_loeuf_tier_by_cutoff.tsv")

p_s5b <- ggplot(loeuf_all, aes(prev_per_10k, LOEUF_tier, color = cutoff)) +
  geom_pointrange(
    aes(xmin = prev_per_10k_lo, xmax = prev_per_10k_hi),
    position = position_dodge(width = 0.55),
    linewidth = 0.4, fatten = 1.6
  ) +
  scale_color_manual(values = CUT_COLS, name = NULL) +
  labs(x = "Prevalence per 10,000 pairs", y = NULL) +
  theme_paper() +
  theme(legend.position = "bottom")

fig_s5 <- tag_s(p_s5a / p_s5b + plot_layout(heights = c(0.9, 1.1)))
save_comp(fig_s5, "Supplementary_Figure5", NATURE_DOUBLE_IN, 7.4)
save_panel(p_s5a, "FigS5a_cutoff_counts", NATURE_DOUBLE_IN, 3.2)
save_panel(p_s5b, "FigS5b_loeuf_by_cutoff", NATURE_DOUBLE_IN, 3.8)
cat("wrote Supplementary_Figure5\n")

# Copy remaining existing robustness TSVs
extra <- c(
  file.path(FIG2X, "underlier_loeuf_wilson_vs_subject_bootstrap.tsv"),
  file.path(SUP_GN, "cohort_vs_gnomad_incidence_per_gene.tsv")
)
file.copy(extra[file.exists(extra)], TSV, overwrite = TRUE)

# ---------------------------------------------------------------------------
# S6 enrichment OR at the paper setting (55 expression PCs)
# ---------------------------------------------------------------------------
cat("S6: enrichment odds ratios at ", PC_SETTING, " PCs ...\n", sep = "")
enr <- read_tsv(file.path(SUP_PC, "enrichment_npc_sv_vs_ptv_underlier.tsv"), show_col_types = FALSE)
s6_df <- enr |>
  filter(as.numeric(npcs_format) == as.numeric(PC_SETTING),
         as.numeric(cutoff) %in% c(2, 3, 4)) |>
  mutate(
    cutoff_lab = factor(
      sprintf("z < -%s", as.integer(cutoff)),
      levels = c("z < -2", "z < -3", "z < -4")
    ),
    variant_short = if_else(
      grepl("PTV", variant_class), "Rare PTV (AF \u2264 0.5%)", "Rare SV (AF \u2264 1%, 100 kb)"
    )
  )
write_out(s6_df, sprintf("S6_enrichment_or_%spcs.tsv", PC_SETTING))

p_s6a <- ggplot(s6_df, aes(odds_ratio, cutoff_lab, color = variant_short)) +
  geom_vline(xintercept = 1, linewidth = 0.35, colour = "grey55") +
  geom_pointrange(
    aes(xmin = ci_95_lo, xmax = ci_95_hi),
    position = position_dodge(width = 0.5),
    linewidth = 0.45, fatten = 1.8
  ) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 15)) +
  scale_color_manual(
    values = c(
      "Rare SV (AF \u2264 1%, 100 kb)" = unname(OI["blue"]),
      "Rare PTV (AF \u2264 0.5%)" = unname(OI["vermillion"])
    ),
    name = NULL
  ) +
  labs(x = "Fisher odds ratio (95% CI)", y = NULL) +
  theme_paper() +
  theme(legend.position = "bottom")

p_s6b <- ggplot(s6_df, aes(cutoff_lab, odds_ratio, color = variant_short, group = variant_short)) +
  geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey55") +
  geom_pointrange(aes(ymin = ci_95_lo, ymax = ci_95_hi), linewidth = 0.45, fatten = 1.8) +
  geom_line(linewidth = 0.4) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 15)) +
  scale_color_manual(
    values = c(
      "Rare SV (AF \u2264 1%, 100 kb)" = unname(OI["blue"]),
      "Rare PTV (AF \u2264 0.5%)" = unname(OI["vermillion"])
    ),
    name = NULL
  ) +
  labs(x = NULL, y = "Fisher odds ratio (95% CI)") +
  theme_paper() +
  theme(legend.position = "bottom")

fig_s6 <- tag_s(p_s6a | p_s6b)
save_comp(fig_s6, "Supplementary_Figure6", NATURE_DOUBLE_IN, 3.8)
save_panel(p_s6a, "FigS6a_enrichment_or_forest", NATURE_DOUBLE_IN / 2, 3.6)
save_panel(p_s6b, "FigS6b_enrichment_or_by_cutoff", NATURE_DOUBLE_IN / 2, 3.6)
cat("wrote Supplementary_Figure6\n")

# ---------------------------------------------------------------------------
# Participant characteristics (TOPMed QC cohort). No ethnicity on disk.
# ---------------------------------------------------------------------------
cat("Participant characteristics ...\n")
kept_meta <- meta |>
  filter(WGS %in% kept_subjects) |>
  distinct(WGS, .keep_all = TRUE)
stopifnot(nrow(kept_meta) == N_TOPMED)

pct <- function(n, N = nrow(kept_meta)) round(100 * n / N, 2)
n_f <- sum(kept_meta$SEX == "Female", na.rm = TRUE)
n_m <- sum(kept_meta$SEX == "Male", na.rm = TRUE)
age_ok <- kept_meta$AGE[is.finite(kept_meta$AGE)]
n_age_miss <- sum(!is.finite(kept_meta$AGE))

cohort_tab <- kept_meta |>
  count(COHORT, name = "n") |>
  mutate(
    section = "Study",
    level = COHORT,
    percent = pct(n),
    extra = NA_character_
  ) |>
  select(section, level, n, percent, extra)

center_tab <- kept_meta |>
  count(SEQUENCING_CENTER, name = "n") |>
  mutate(
    section = "Sequencing centre",
    level = SEQUENCING_CENTER,
    percent = pct(n),
    extra = NA_character_
  ) |>
  select(section, level, n, percent, extra)

participants <- bind_rows(
  tibble(section = "Cohort", level = "TOPMed QC kept subjects", n = N_TOPMED, percent = 100, extra = "kept_subjects.rds"),
  tibble(section = "Sex", level = "Female", n = n_f, percent = pct(n_f), extra = NA_character_),
  tibble(section = "Sex", level = "Male", n = n_m, percent = pct(n_m), extra = NA_character_),
  tibble(section = "Age (years)", level = "Median (IQR)", n = length(age_ok), percent = NA_real_,
         extra = sprintf("%.1f (%.1f-%.1f); n_with_age=%d; n_missing=%d (all WHI)",
                         median(age_ok), quantile(age_ok, 0.25), quantile(age_ok, 0.75),
                         length(age_ok), n_age_miss)),
  tibble(section = "Age (years)", level = "Mean (SD)", n = length(age_ok), percent = NA_real_,
         extra = sprintf("%.2f (%.2f)", mean(age_ok), sd(age_ok))),
  tibble(section = "Ethnicity / race", level = "Not available", n = NA_integer_, percent = NA_real_,
         extra = "Not in metadata.filtered_with_sv.txt or 2025-06-11/data/metadata.filtered.txt"),
  tibble(section = "AoU", level = "Participant characteristics", n = N_AOU, percent = NA_real_,
         extra = "RNA n=8895 gene-level only; no AoU sex/age/ethnicity table on disk"),
  cohort_tab,
  center_tab
)
write_out(participants, "S_participant_characteristics.tsv")

# ---------------------------------------------------------------------------
# Excel workbook (Python / openpyxl)
# ---------------------------------------------------------------------------
cat("Packing Supplementary_Tables.xlsx ...\n")
py <- Sys.getenv("PAPER_PYTHON", unset = "/net/wonderland/home/sheyang/.pixi/envs/python/bin/python3")
xlsx_script <- file.path(ROOT, "scripts/supplement/pack_supplementary_tables.py")
st <- system2(py, shQuote(xlsx_script), stdout = "", stderr = "")
if (!identical(st, 0L)) warning("Excel packager failed with status ", st)

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
manifest <- tibble(
  file = list.files(OUT, pattern = "Supplementary_Figure|FigS", full.names = FALSE),
  kind = "figure"
) |>
  bind_rows(tibble(
    file = file.path("supplemental_tables", list.files(TSV)),
    kind = "table"
  ))
write_out(manifest, "MANIFEST.tsv")

cat("\nDone. Figures in\n  ", OUT, "\nTables in\n  ", TSV, "\n", sep = "")
