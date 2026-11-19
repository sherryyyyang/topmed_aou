#!/usr/bin/env Rscript
# Compare residualized z-scores for the SAME GENCODE v30 protein-coding genes:
#   published 55-PC coding-only HVG pipeline
#     vs
#   all-gene HVG pipeline (HVG included noncoding; n_expr_pcs=55, n_geno_pcs=5)
#
# Fair PC match: published 55-PC z, NOT the 2005-PC replicate.
# Subjects are intersected (published 5,560 vs all-gene 6,024).
# Does not load expr_scaled_flat_list.rds. Loads only shared coding genes
# from the long-form z tables (data.table subset immediately after read).
#
# Writes:
#   allgenes_pipeline_out/Fig_coding_zscore_hvg_coding_vs_allgenes.{png,pdf}
#   allgenes_pipeline_out/coding_z_hvg_vs_allgenes_summary.tsv
#   allgenes_pipeline_out/coding_z_hvg_vs_allgenes_per_gene.tsv
# and copies the figure (new name) into submission_package_nov_2026/.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

OUT <- file.path(ROOT, "allgenes_pipeline_out")
DRAFT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
STEM <- "Fig_coding_zscore_hvg_coding_vs_allgenes"
PUB_Z <- file.path(ROOT, "data/expression/expr_z_join.rds")
PUB_CPM <- file.path(ROOT, "data/expression/expr_cpm_adj.rds")
PUB_UNDER <- file.path(ROOT, "data/underliers/z_score_underlier_z_-3_clean.rds")
NEW_Z <- file.path(OUT, "expr_zscore_long.rds")
NEW_CPM <- file.path(OUT, "expr_cpm_adj.rds")
NEW_UNDER <- file.path(OUT, "underliers_z_-3.rds")
GENE_META <- file.path(OUT, "gene_meta.rds")
N_SCATTER <- 2e6L
DELTA_LOTS <- c(5L, 10L, 20L)
CACHE <- file.path(OUT, "coding_z_hvg_vs_allgenes_cache.rds")
set.seed(1)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

nv <- function(x) str_remove(x, "\\..*$")

# --------------------------------------------------------------------------
if (file.exists(CACHE)) {
  msg("Loading cached pairwise comparison ...")
  .cache <- read_rds(CACHE)
  merged <- .cache$merged
  per_gene <- .cache$per_gene
  gm <- .cache$gm
  env_stats <- .cache$env_stats
  list2env(env_stats, envir = environment())
  rm(.cache)
  gc(verbose = FALSE)
  msg("  cache: ", comma(nrow(merged)), " pairs, ", comma(nrow(per_gene)), " genes")
} else {
  msg("Gene / subject axes (matrices + GENCODE types; not the z tables) ...")
  pub_cpm <- read_rds(PUB_CPM)
  new_cpm <- read_rds(NEW_CPM)
  gm <- as.data.table(read_rds(GENE_META))
  pub_ids <- rownames(pub_cpm)
  new_ids <- rownames(new_cpm)
  pub_subj <- colnames(pub_cpm)
  new_subj <- colnames(new_cpm)
  rm(pub_cpm, new_cpm)
  gc(verbose = FALSE)

  coding_nv <- unique(gm$gene_nv[gm$gene_type == "protein_coding"])
  pub_nv <- nv(pub_ids)
  new_nv <- nv(new_ids)
  stopifnot(!anyDuplicated(pub_nv), !anyDuplicated(new_nv))

  new_coding_nv <- intersect(new_nv, coding_nv)
  shared_nv <- intersect(pub_nv, new_coding_nv)
  pub_only_nv <- setdiff(pub_nv, new_nv)
  new_only_coding_nv <- setdiff(new_coding_nv, pub_nv)
  shared_ids <- pub_ids[pub_nv %in% shared_nv]
  shared_subj <- intersect(pub_subj, new_subj)

  # versioned IDs must match for shared genes (same GENCODE v30)
  new_shared_ids <- new_ids[new_nv %in% shared_nv]
  if (!setequal(shared_ids, new_shared_ids)) {
    stop("Shared coding genes have versioned ID mismatches between matrices")
  }

  n_shared_genes <- length(shared_nv)
  n_shared_subj <- length(shared_subj)
  n_pairs <- as.numeric(n_shared_genes) * n_shared_subj
  msg(sprintf(
    "  coding with z: shared=%s  published-only=%s  allgene-only=%s",
    comma(n_shared_genes), comma(length(pub_only_nv)), comma(length(new_only_coding_nv))
  ))
  msg(sprintf(
    "  subjects: shared=%s  published-only=%s  allgene-only=%s  pairs=%s",
    comma(n_shared_subj),
    comma(length(setdiff(pub_subj, new_subj))),
    comma(length(setdiff(new_subj, pub_subj))),
    comma(n_pairs)
  ))
  if (length(pub_only_nv)) {
    pub_only_lab <- gm[gene_nv %in% pub_only_nv, paste(symbol, collapse = ", ")]
    msg("  published-only coding HVG: ", pub_only_lab)
  }

  # Published z: 8,124 coding genes x 5,560 samples. Keep shared genes only.
  msg("Loading published 55-PC expr_z_join (then subset) ...")
  old <- as.data.table(read_rds(PUB_Z))
  old <- old[gene_id %chin% shared_ids,
             .(gene_id, subjectid, z_old = expression_zscore)]
  gc(verbose = FALSE)
  msg("  published rows kept: ", comma(nrow(old)))

  # All-gene z long form (no logcpm) — subset immediately to shared gene x subject.
  msg("Loading all-gene expr_zscore_long (then subset; not the 6G join) ...")
  new <- as.data.table(read_rds(NEW_Z))
  new <- new[gene_id %chin% shared_ids & subjectid %chin% shared_subj,
             .(gene_id, subjectid, z_new = expression_zscore)]
  gc(verbose = FALSE)
  msg("  all-gene rows kept: ", comma(nrow(new)))

  msg("Inner-joining shared gene x subject pairs ...")
  setkey(old, gene_id, subjectid)
  setkey(new, gene_id, subjectid)
  merged <- old[new, nomatch = 0]
  rm(old, new)
  gc(verbose = FALSE)
  merged <- merged[is.finite(z_old) & is.finite(z_new)]
  msg("  finite shared pairs: ", comma(nrow(merged)))
  if (nrow(merged) == 0L) stop("No shared finite gene x subject pairs")

  msg("Pairwise and per-gene correlations ...")
  pearson_all <- cor(merged$z_old, merged$z_new, method = "pearson")
  spearman_all <- cor(merged$z_old, merged$z_new, method = "spearman")
  abs_dz <- abs(merged$z_old - merged$z_new)
  mae <- mean(abs_dz)
  med_abs <- median(abs_dz)
  rmse <- sqrt(mean((merged$z_old - merged$z_new)^2))
  frac_gt_0.5 <- mean(abs_dz > 0.5)
  frac_gt_1 <- mean(abs_dz > 1)
  frac_gt_2 <- mean(abs_dz > 2)
  rm(abs_dz)
  gc(verbose = FALSE)

  per_gene <- merged[, {
    n <- .N
    if (n < 10L) {
      .(n = n, pearson = NA_real_, spearman = NA_real_)
    } else {
      .(n = n,
        pearson = cor(z_old, z_new, method = "pearson"),
        spearman = cor(z_old, z_new, method = "spearman"))
    }
  }, by = gene_id]
  per_gene[, gene_nv := nv(gene_id)]
  r_ok <- per_gene$pearson[is.finite(per_gene$pearson)]
  r_q <- quantile(r_ok, probs = c(0.25, 0.5, 0.75), names = FALSE)
  msg(sprintf(
    "  overall Pearson=%.4f  Spearman=%.4f  MAE=%.4f  median|dz|=%.4f",
    pearson_all, spearman_all, mae, med_abs
  ))
  msg(sprintf(
    "  per-gene Pearson median=%.4f  IQR=[%.4f, %.4f]  n_genes=%s",
    r_q[2], r_q[1], r_q[3], comma(length(r_ok))
  ))

  msg("Underlier concordance on shared pairs (z<-3 AND log2CPM drop) ...")
  u_old <- as.data.table(read_rds(PUB_UNDER))
  u_new <- as.data.table(read_rds(NEW_UNDER))
  u_old <- unique(u_old[gene_id %chin% shared_ids & subjectid %chin% shared_subj,
                        .(gene_id, subjectid)])
  u_new <- unique(u_new[gene_id %chin% shared_ids & subjectid %chin% shared_subj,
                        .(gene_id, subjectid)])
  u_old[, old_under := TRUE]
  u_new[, new_under := TRUE]
  setkey(u_old, gene_id, subjectid)
  setkey(u_new, gene_id, subjectid)
  u_join <- merge(u_old, u_new, all = TRUE)
  u_join[is.na(old_under), old_under := FALSE]
  u_join[is.na(new_under), new_under := FALSE]
  n_both <- as.integer(u_join[old_under == TRUE & new_under == TRUE, .N])
  n_old_only <- as.integer(u_join[old_under == TRUE & new_under == FALSE, .N])
  n_new_only <- as.integer(u_join[old_under == FALSE & new_under == TRUE, .N])
  n_old_u <- n_both + n_old_only
  n_new_u <- n_both + n_new_only
  n_neither <- nrow(merged) - n_both - n_old_only - n_new_only
  jaccard <- if ((n_both + n_old_only + n_new_only) > 0L) {
    n_both / (n_both + n_old_only + n_new_only)
  } else {
    NA_real_
  }

  n_old_g <- u_old[, .(n_under_old = .N), by = gene_id]
  n_new_g <- u_new[, .(n_under_new = .N), by = gene_id]
  per_gene <- merge(per_gene, n_old_g, by = "gene_id", all.x = TRUE)
  per_gene <- merge(per_gene, n_new_g, by = "gene_id", all.x = TRUE)
  per_gene[is.na(n_under_old), n_under_old := 0L]
  per_gene[is.na(n_under_new), n_under_new := 0L]
  per_gene[, delta_n := n_under_new - n_under_old]
  n_cor <- cor(per_gene$n_under_old, per_gene$n_under_new, method = "pearson")
  n_cor_s <- cor(per_gene$n_under_old, per_gene$n_under_new, method = "spearman")
  n_chg <- sapply(DELTA_LOTS, function(k) as.integer(sum(abs(per_gene$delta_n) >= k)))
  names(n_chg) <- paste0("n_genes_abs_delta_ge_", DELTA_LOTS)
  n_merged <- nrow(merged)

  msg(sprintf(
    "  2x2  both=%s  old-only=%s  new-only=%s  neither=%s  Jaccard=%.4f",
    comma(n_both), comma(n_old_only), comma(n_new_only), comma(n_neither), jaccard
  ))
  msg(sprintf(
    "  n_underlier per gene: Pearson=%.4f  Spearman=%.4f  |d|>=5: %s  |d|>=10: %s  |d|>=20: %s",
    n_cor, n_cor_s, comma(n_chg[1]), comma(n_chg[2]), comma(n_chg[3])
  ))

  env_stats <- list(
    pub_nv = pub_nv, new_coding_nv = new_coding_nv,
    n_shared_genes = n_shared_genes, pub_only_nv = pub_only_nv,
    new_only_coding_nv = new_only_coding_nv, pub_subj = pub_subj,
    new_subj = new_subj, n_shared_subj = n_shared_subj,
    pearson_all = pearson_all, spearman_all = spearman_all,
    mae = mae, med_abs = med_abs, rmse = rmse,
    frac_gt_0.5 = frac_gt_0.5, frac_gt_1 = frac_gt_1, frac_gt_2 = frac_gt_2,
    r_ok = r_ok, r_q = r_q, n_both = n_both, n_old_only = n_old_only,
    n_new_only = n_new_only, n_neither = n_neither, jaccard = jaccard,
    n_old_u = n_old_u, n_new_u = n_new_u, n_cor = n_cor, n_cor_s = n_cor_s,
    n_chg = n_chg, n_merged = nrow(merged)
  )
  # Keep a 2e6-row scatter sample in the cache, not all 45M pairs.
  sc_keep <- merged[sample.int(nrow(merged), min(N_SCATTER, nrow(merged)))]
  write_rds(
    list(merged = sc_keep, per_gene = per_gene, gm = gm, env_stats = env_stats,
         n_merged = nrow(merged)),
    CACHE, compress = "gz"
  )
  msg("  wrote cache (scatter subsample + per-gene + stats)")
}

# --------------------------------------------------------------------------
msg("Writing summary / per-gene tables ...")
pub_only_sym <- gm[gene_nv %in% pub_only_nv, paste(sort(symbol), collapse = ";")]
summary_dt <- data.table(
  metric = c(
    "n_published_coding_hvg",
    "n_allgene_coding_hvg",
    "n_shared_coding_genes",
    "n_published_only_coding",
    "n_allgene_only_coding",
    "published_only_symbols",
    "n_published_subjects",
    "n_allgene_subjects",
    "n_shared_subjects",
    "n_allgene_only_subjects",
    "n_shared_finite_pairs",
    "pearson_all_pairs",
    "spearman_all_pairs",
    "mae_abs_dz",
    "median_abs_dz",
    "rmse_dz",
    "frac_pairs_abs_dz_gt_0.5",
    "frac_pairs_abs_dz_gt_1",
    "frac_pairs_abs_dz_gt_2",
    "per_gene_pearson_median",
    "per_gene_pearson_q25",
    "per_gene_pearson_q75",
    "per_gene_spearman_median",
    "per_gene_spearman_q25",
    "per_gene_spearman_q75",
    "n_genes_per_gene_cor",
    "underlier_both",
    "underlier_old_only",
    "underlier_new_only",
    "underlier_neither",
    "underlier_jaccard",
    "n_underlier_old_shared",
    "n_underlier_new_shared",
    "per_gene_n_underlier_pearson",
    "per_gene_n_underlier_spearman",
    "median_delta_n_underlier",
    "iqr_delta_n_underlier_q25",
    "iqr_delta_n_underlier_q75",
    names(n_chg),
    "pc_match",
    "underlier_definition"
  ),
  value = c(
    as.character(length(pub_nv)),
    as.character(length(new_coding_nv)),
    as.character(n_shared_genes),
    as.character(length(pub_only_nv)),
    as.character(length(new_only_coding_nv)),
    pub_only_sym,
    as.character(length(pub_subj)),
    as.character(length(new_subj)),
    as.character(n_shared_subj),
    as.character(length(setdiff(new_subj, pub_subj))),
    as.character(n_merged),
    sprintf("%.6f", pearson_all),
    sprintf("%.6f", spearman_all),
    sprintf("%.6f", mae),
    sprintf("%.6f", med_abs),
    sprintf("%.6f", rmse),
    sprintf("%.6f", frac_gt_0.5),
    sprintf("%.6f", frac_gt_1),
    sprintf("%.6f", frac_gt_2),
    sprintf("%.6f", r_q[2]),
    sprintf("%.6f", r_q[1]),
    sprintf("%.6f", r_q[3]),
    sprintf("%.6f", median(per_gene$spearman, na.rm = TRUE)),
    sprintf("%.6f", quantile(per_gene$spearman, 0.25, na.rm = TRUE)),
    sprintf("%.6f", quantile(per_gene$spearman, 0.75, na.rm = TRUE)),
    as.character(length(r_ok)),
    as.character(n_both),
    as.character(n_old_only),
    as.character(n_new_only),
    as.character(n_neither),
    sprintf("%.6f", jaccard),
    as.character(n_old_u),
    as.character(n_new_u),
    sprintf("%.6f", n_cor),
    sprintf("%.6f", n_cor_s),
    sprintf("%.6f", median(per_gene$delta_n)),
    sprintf("%.6f", quantile(per_gene$delta_n, 0.25)),
    sprintf("%.6f", quantile(per_gene$delta_n, 0.75)),
    as.character(n_chg),
    "published_55PC_vs_allgene_55PC",
    "z < -3 AND log2CPM < gene_mean - 1"
  )
)
fwrite(summary_dt, file.path(OUT, "coding_z_hvg_vs_allgenes_summary.tsv"), sep = "\t")
if (!"symbol" %in% names(per_gene)) {
  per_gene <- merge(per_gene, gm[, .(gene_nv, symbol, gene_type)], by = "gene_nv", all.x = TRUE)
}
per_gene[, abs_delta := abs(delta_n)]
setorder(per_gene, -abs_delta, gene_nv)
per_gene[, abs_delta := NULL]
fwrite(
  per_gene[, .(gene_id, gene_nv, symbol, n, pearson, spearman,
               n_under_old, n_under_new, delta_n)],
  file.path(OUT, "coding_z_hvg_vs_allgenes_per_gene.tsv"),
  sep = "\t"
)

# --------------------------------------------------------------------------
msg("Plotting sanity figure ...")
n_plot <- min(N_SCATTER, nrow(merged))
sc <- merged[sample.int(nrow(merged), n_plot)]
lim <- c(-8, 8)
ann <- sprintf(
  "Pearson r = %.3f\nSpearman ρ = %.3f\n%s pairs (plot: %s)",
  pearson_all, spearman_all, comma(n_merged), comma(n_plot)
)

p_sc <- ggplot(sc, aes(z_old, z_new)) +
  geom_bin2d(bins = 80, boundary = 0) +
  geom_abline(slope = 1, intercept = 0, colour = "white", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, colour = "grey20", linewidth = 0.35,
              linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10", option = "C", name = "pairs") +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  labs(
    x = "Published z (coding HVG, 55 PCs)",
    y = "All-gene z (55 PCs)"
  ) +
  annotate(
    "text", x = -7.6, y = 7.4, label = ann, hjust = 0, vjust = 1,
    size = 2.6, family = PAPER_FONT, lineheight = 0.95
  ) +
  theme_paper() +
  theme(legend.position = "right")

p_hist <- ggplot(per_gene[is.finite(pearson)], aes(x = pearson)) +
  geom_histogram(bins = 40, fill = unname(OI["blue"]), colour = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = r_q[2], colour = unname(OI["vermillion"]),
             linewidth = 0.6) +
  scale_x_continuous(limits = c(min(0.5, min(r_ok) - 0.02), 1.001),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Per-gene Pearson r (old z vs new z)",
    y = "Coding genes"
  ) +
  annotate(
    "text",
    x = r_q[2], y = Inf,
    label = sprintf("median = %.3f\nIQR [%s, %s]", r_q[2],
                    format(round(r_q[1], 3), nsmall = 3),
                    format(round(r_q[3], 3), nsmall = 3)),
    hjust = 1.05, vjust = 1.2, size = 2.6, family = PAPER_FONT,
    colour = unname(OI["vermillion"]), lineheight = 0.95
  ) +
  theme_paper()

p_n <- ggplot(per_gene, aes(n_under_old, n_under_new)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_point(alpha = 0.25, size = 0.7, colour = unname(OI["blue"]),
             stroke = 0) +
  labs(
    x = "n underliers, published (shared subjects)",
    y = "n underliers, all-gene"
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf(
      "r = %.3f   Jaccard (pairs) = %.3f\n|Δn|≥10: %s / %s genes",
      n_cor, jaccard, comma(n_chg[2]), comma(nrow(per_gene))
    ),
    hjust = -0.05, vjust = 1.2, size = 2.6, family = PAPER_FONT, lineheight = 0.95
  ) +
  theme_paper()

fig <- (p_sc | p_hist) / p_n +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = PAPER_FONT, face = "bold", size = 14, colour = "black"
    )
  )

paper_ggsave(fig, file.path(OUT, STEM), width = NATURE_DOUBLE_IN, height = 6.4)
if (dir.exists(DRAFT)) {
  file.copy(paste0(file.path(OUT, STEM), ".png"),
            file.path(DRAFT, paste0(STEM, ".png")), overwrite = TRUE)
  file.copy(paste0(file.path(OUT, STEM), ".pdf"),
            file.path(DRAFT, paste0(STEM, ".pdf")), overwrite = TRUE)
  src <- file.path(DRAFT, "source_panels")
  if (dir.exists(src)) {
    file.copy(paste0(file.path(OUT, STEM), ".png"),
              file.path(src, paste0(STEM, ".png")), overwrite = TRUE)
    file.copy(paste0(file.path(OUT, STEM), ".pdf"),
              file.path(src, paste0(STEM, ".pdf")), overwrite = TRUE)
  }
}

cat("\n==== SUMMARY ====\n")
print(summary_dt, nrows = 80)
cat("\nTop |Δn_underlier| genes:\n")
print(head(per_gene[, .(symbol, gene_nv, n_under_old, n_under_new, delta_n, pearson)], 15))
msg("DONE  figure ", file.path(OUT, paste0(STEM, ".png")))
