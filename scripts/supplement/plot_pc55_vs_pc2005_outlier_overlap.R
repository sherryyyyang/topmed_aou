#!/usr/bin/env Rscript
# Unique vs shared underlier calls at 55 vs 2005 expression PCs.
#
# 55 PCs is the setting used throughout the paper; 2005 is the Gavish-Donoho
# choice that the PC sweep (scripts/supplement/sweep_enrichment_*.R) shows to
# peak for rare-variant enrichment. This script asks how much the two settings
# actually disagree at the level of individual gene x subject underlier calls,
# and whether the calls unique to one setting look more like true regulatory
# events (rare PTV / rare SV carriers) than the calls unique to the other.
#
# Data source: the two long-form z-score tables
#   data/expression/expr_scaled_genotype_pc_5_expression_pc_55_longer_form.tsv
#   data/expression/expr_scaled_genotype_pc_5_expression_pc_2005_longer_form.tsv
# are byte-identical in their gene_id/subjectid columns (verified below), i.e.
# both are the native row order of the 8,124 gene x 5,560 subject expression
# matrix. Nothing has to be re-residualized: only the z column of each file is
# read, so the 28 GB expr_scaled_flat_list.rds is never touched.
#
# Underlier definition is the paper definition, held fixed across PC settings:
#   expression_zscore < -Z  AND  log2CPM_adj < gene_mean(log2CPM_adj) - 1
# followed by the global extreme-subject filter (kept_subjects.rds). The logCPM
# side comes from expr_cpm_adj.rds and is identical for both settings, so the
# only thing that varies is the z-score matrix.
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(scales); library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

OUT_DIR    <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
PANEL_DIR  <- file.path(OUT_DIR, "source_panels")
CACHE      <- file.path(OUT_DIR, "pc55_vs_pc2005_overlap_cache.rds")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

Z55_TSV   <- file.path(DATA_EXPR, "expr_scaled_genotype_pc_5_expression_pc_55_longer_form.tsv")
Z2005_TSV <- file.path(DATA_EXPR, "expr_scaled_genotype_pc_5_expression_pc_2005_longer_form.tsv")
SV_RDS    <- file.path(DATA_EXPR, "subset_genotypes_100000_0.01_select.rds")

CUTOFFS      <- c(2, 3, 3.5, 4, 5)
PRIMARY_CUT  <- 3
LAB55        <- "55 PCs only"
LAB20        <- "2005 PCs only"
LABSH        <- "Shared"
CAT_LEVELS   <- c(LAB55, LABSH, LAB20)
CAT_COLORS   <- setNames(c(unname(OI["blue"]), "grey60", unname(OI["orange"])), CAT_LEVELS)

# ---------------------------------------------------------------- 2x2 helpers
# Haldane-corrected odds ratio with a Woolf 95% CI. fisher.test() is not usable
# here: the background cell counts run to tens of millions of gene x subject
# pairs.
or_woolf <- function(n_cat_var, n_cat_novar, n_bg_var, n_bg_novar) {
  a <- n_cat_var + 0.5; b <- n_cat_novar + 0.5
  c <- n_bg_var  + 0.5; d <- n_bg_novar  + 0.5
  lor <- log((a * d) / (b * c))
  se  <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  tibble(odds_ratio = exp(lor), ci_95_lo = exp(lor - 1.96 * se),
         ci_95_hi = exp(lor + 1.96 * se), p_value = 2 * pnorm(-abs(lor / se)))
}

# --------------------------------------------------------------------- stage 1
# Everything expensive happens once and is cached; re-runs only replot.
build_cache <- function() {
  message("[", format(Sys.time(), "%H:%M:%S"), "] loading expr_cpm_adj.rds ...")
  cpm <- read_rds(EXPR_CPM_ADJ)
  nsub <- ncol(cpm); ngene <- nrow(cpm)
  message(sprintf("  expression grid: %s genes x %s subjects = %s pairs",
                  comma(ngene), comma(nsub), comma(as.numeric(ngene) * nsub)))

  # The long-form tables are one block of all subjects per gene. The subject
  # order matches expr_cpm_adj columns, but the gene blocks are a permutation of
  # its rows, so recover the file's own axes and align the matrix to them.
  axes_rds <- file.path(OUT_DIR, "longform_axes.rds")
  axes <- if (file.exists(axes_rds)) read_rds(axes_rds) else {
    message("[", format(Sys.time(), "%H:%M:%S"), "] recovering long-form gene/subject axes ...")
    a <- list(
      subjects = fread(cmd = sprintf("head -n %d %s | cut -f2 | tail -n +2", nsub + 1, Z55_TSV),
                       header = FALSE)[[1]],
      genes = fread(cmd = sprintf("awk -F'\\t' 'NR>1 && (NR-2)%%%d==0 {print $1}' %s", nsub, Z55_TSV),
                    header = FALSE)[[1]])
    write_rds(a, axes_rds)
    a
  }
  genes <- axes$genes; subjects <- axes$subjects
  stopifnot(length(genes) == ngene, length(subjects) == nsub,
            !any(duplicated(genes)), !any(duplicated(subjects)),
            setequal(genes, rownames(cpm)), identical(subjects, colnames(cpm)))
  cpm <- cpm[genes, , drop = FALSE]
  message("  aligned expr_cpm_adj to the long-form gene/subject axes")

  # logCPM side of the underlier definition: identical for both PC settings.
  # Recycling a length-ngene vector down the rows of a genes x subjects matrix
  # is the per-gene comparison we want.
  pass_logcpm <- t(cpm < (rowMeans(cpm) - 1))          # subjects x genes
  rm(cpm); gc()
  kept <- subjects %in% read_rds(KEPT_SUBJECTS)        # recycles down rows
  message(sprintf("  logCPM-eligible pairs: %s; kept subjects: %s / %s",
                  comma(sum(pass_logcpm)), comma(sum(kept)), comma(nsub)))

  read_z <- function(f, tag) {
    message("[", format(Sys.time(), "%H:%M:%S"), "] reading z column: ", tag)
    t0 <- Sys.time()
    z <- fread(f, select = 3, showProgress = FALSE)[[1]]
    stopifnot(length(z) == as.numeric(ngene) * nsub, !anyNA(z))
    dim(z) <- c(nsub, ngene)                           # subject varies fastest
    message(sprintf("  %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    z
  }
  z55 <- read_z(Z55_TSV, "55 PCs")
  z20 <- read_z(Z2005_TSV, "2005 PCs")

  # --- rare-variant carrier matrices, same shape as the call matrices --------
  message("[", format(Sys.time(), "%H:%M:%S"), "] building carrier matrices ...")
  carrier_matrix <- function(gene_key, sample_key, gene_axis) {
    ri <- match(sample_key, subjects)
    ci <- match(gene_key, gene_axis)
    ok <- !is.na(ri) & !is.na(ci)
    m <- matrix(FALSE, nsub, ngene)
    m[cbind(ri[ok], ci[ok])] <- TRUE
    m
  }
  ptv <- fread(PTV_TSV, select = c("gene", "sample"))
  ptv_m <- carrier_matrix(ptv$gene, ptv$sample, str_remove(genes, "\\..*$"))
  sv <- read_rds(SV_RDS) |> ungroup() |> distinct(gene_id, sample)
  sv_m <- carrier_matrix(sv$gene_id, sv$sample, genes)
  message(sprintf("  PTV carrier pairs on grid: %s; SV carrier pairs on grid: %s",
                  comma(sum(ptv_m)), comma(sum(sv_m))))

  # --- per-cutoff overlap and variant enrichment -----------------------------
  overlap <- list(); enrich <- list(); per_gene <- NULL; per_subject <- NULL
  for (zc in CUTOFFS) {
    u55 <- (z55 < -zc) & pass_logcpm & kept
    u20 <- (z20 < -zc) & pass_logcpm & kept
    shared <- u55 & u20
    n_shared <- sum(shared); n_55 <- sum(u55); n_20 <- sum(u20)
    n_only55 <- n_55 - n_shared; n_only20 <- n_20 - n_shared

    # Background: eligible (kept-subject) pairs called by neither setting.
    n_eligible <- sum(kept) * as.numeric(ngene)
    n_bg <- n_eligible - (n_shared + n_only55 + n_only20)

    overlap[[as.character(zc)]] <- tibble(
      cutoff = zc, n_underlier_pc55 = n_55, n_underlier_pc2005 = n_20,
      n_shared = n_shared, n_only_pc55 = n_only55, n_only_pc2005 = n_only20,
      n_eligible_pairs = n_eligible,
      jaccard = n_shared / (n_shared + n_only55 + n_only20),
      frac_pc55_shared = n_shared / n_55, frac_pc2005_shared = n_shared / n_20)

    for (vc in c("PTV", "SV")) {
      vm <- if (vc == "PTV") ptv_m else sv_m
      n_var_shared <- sum(shared & vm)
      n_var_55 <- sum(u55 & vm) - n_var_shared
      n_var_20 <- sum(u20 & vm) - n_var_shared
      n_var_elig <- sum(vm & kept)
      n_var_bg <- n_var_elig - (n_var_shared + n_var_55 + n_var_20)
      n_pairs_cat   <- c(n_only55, n_shared, n_only20)
      n_carrier_cat <- c(n_var_55, n_var_shared, n_var_20)
      enrich[[paste(zc, vc)]] <- tibble(
        cutoff = zc, variant_class = vc,
        category = factor(CAT_LEVELS, levels = CAT_LEVELS),
        n_pairs = n_pairs_cat, n_carrier = n_carrier_cat) |>
        mutate(carrier_rate_per_10k = 1e4 * n_carrier / n_pairs,
               bg_rate_per_10k = 1e4 * n_var_bg / n_bg) |>
        bind_cols(or_woolf(n_carrier_cat, n_pairs_cat - n_carrier_cat,
                           n_var_bg, n_bg - n_var_bg))
    }

    if (zc == PRIMARY_CUT) {
      # The 55-PC calls at the paper cutoff already exist as a saved table built
      # by a completely different code path (run_strict_zcutoff_figures.R via
      # expr_z_join.rds). Require an exact set match before trusting the 2005-PC
      # calls produced here.
      idx <- which(u55)
      got  <- paste(genes[(idx - 1) %/% nsub + 1], subjects[(idx - 1) %% nsub + 1])
      saved <- read_rds(file.path(DATA_UNDERLIER,
                                  sprintf("z_score_underlier_z_-%s_clean.rds", zc)))
      want <- paste(saved$gene_id, saved$subjectid)
      message(sprintf("  validation vs saved 55-PC calls: %s here, %s saved, %s shared",
                      comma(length(got)), comma(length(want)),
                      comma(length(intersect(got, want)))))
      stopifnot(setequal(got, want))
      message("  55-PC underlier calls reproduce the saved table exactly")

      per_gene <- tibble(gene_id = genes,
                         n_pc55 = colSums(u55), n_pc2005 = colSums(u20),
                         n_shared = colSums(shared)) |>
        mutate(n_only_pc55 = n_pc55 - n_shared, n_only_pc2005 = n_pc2005 - n_shared)
      per_subject <- tibble(subjectid = subjects,
                            n_pc55 = rowSums(u55), n_pc2005 = rowSums(u20),
                            n_shared = rowSums(shared)) |>
        filter(subjectid %in% subjects[kept]) |>
        mutate(n_only_pc55 = n_pc55 - n_shared, n_only_pc2005 = n_pc2005 - n_shared)
      # How extreme are the unique calls under the *other* PC setting? This
      # separates "genuinely discordant" from "just across the threshold".
      z_unique <- tibble(
        category = factor(c(rep(LAB55, n_only55), rep(LAB20, n_only20)), levels = CAT_LEVELS),
        z_own   = c(z55[u55 & !u20], z20[u20 & !u55]),
        z_other = c(z20[u55 & !u20], z55[u20 & !u55]))
    }
    message(sprintf("[%s] z<-%-3s  pc55=%s  pc2005=%s  shared=%s  Jaccard=%.3f",
                    format(Sys.time(), "%H:%M:%S"), zc, comma(n_55), comma(n_20),
                    comma(n_shared), n_shared / (n_55 + n_20 - n_shared)))
    rm(u55, u20, shared); gc()
  }

  out <- list(overlap = bind_rows(overlap), enrich = bind_rows(enrich),
              per_gene = per_gene, per_subject = per_subject, z_unique = z_unique,
              ngene = ngene, nsub = nsub, n_kept = sum(kept))
  write_rds(out, CACHE, compress = "gz")
  message("wrote cache: ", CACHE)
  out
}

D <- if (file.exists(CACHE)) {
  message("using cache: ", CACHE); read_rds(CACHE)
} else build_cache()

ov  <- D$overlap
ovp <- filter(ov, cutoff == PRIMARY_CUT)
cat("\noverlap of underlier calls by z cutoff:\n"); print(as.data.frame(ov))
cat("\nrare-variant enrichment of unique vs shared calls:\n")
print(as.data.frame(filter(D$enrich, cutoff == PRIMARY_CUT)))

cut_lab <- function(z) factor(sprintf("z < -%s", z), levels = sprintf("z < -%s", CUTOFFS))
# Short forms of the category names for crowded discrete axes.
CAT_SHORT <- setNames(c("55 PCs\nonly", "Shared", "2005 PCs\nonly"), CAT_LEVELS)

# "x of y lose / z gain" caption for the paired-count scatters.
gain_lose_label <- function(d, unit) {
  n_lose <- sum(d$n_pc2005 < d$n_pc55); n_gain <- sum(d$n_pc2005 > d$n_pc55)
  sprintf("%s of %s %s lose calls at 2005 PCs\n%s gain, %s unchanged",
          comma(n_lose), comma(nrow(d)), unit, comma(n_gain),
          comma(nrow(d) - n_lose - n_gain))
}

# ------------------------------------------------------------------- panel a
# Absolute breakdown at the paper cutoff.
pa_dat <- tibble(category = factor(CAT_LEVELS, levels = CAT_LEVELS),
                 n = c(ovp$n_only_pc55, ovp$n_shared, ovp$n_only_pc2005))
pa <- ggplot(pa_dat, aes(category, n, fill = category)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = comma(n)), vjust = -0.45, size = 2.8,
            family = PAPER_FONT, fontface = "bold") +
  scale_fill_manual(values = CAT_COLORS, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Underlier gene x subject pairs") +
  theme_paper() +
  theme(axis.text.x = element_text(size = 8))

# ------------------------------------------------------------------- panel b
# Same breakdown as a proportion, across z cutoffs.
pb_dat <- ov |>
  dplyr::select(cutoff, !!LAB55 := n_only_pc55, !!LABSH := n_shared, !!LAB20 := n_only_pc2005) |>
  pivot_longer(-cutoff, names_to = "category", values_to = "n") |>
  mutate(category = factor(category, levels = CAT_LEVELS), cutoff_lab = cut_lab(cutoff)) |>
  group_by(cutoff) |> mutate(frac = n / sum(n)) |> ungroup()
pb_lab <- ov |> mutate(cutoff_lab = cut_lab(cutoff),
                       lab = sprintf("%.0f%%", 100 * jaccard))
pb <- ggplot(pb_dat, aes(cutoff_lab, frac, fill = category)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.25) +
  geom_text(data = pb_lab, aes(x = cutoff_lab, y = 1.04, label = lab), inherit.aes = FALSE,
            size = 2.6, family = PAPER_FONT, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = CAT_COLORS, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.09),
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Share of union of calls") +
  theme_paper() +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, -4, 0))

# ------------------------------------------------------------------- panel c
ps <- D$per_subject
sub_lim <- c(0, max(ps$n_pc55, ps$n_pc2005))
pc <- ggplot(ps, aes(n_pc55, n_pc2005)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 0.5, alpha = 0.35, colour = unname(OI["blue"])) +
  annotate("text", x = sub_lim[2] * 0.02, y = sub_lim[2], hjust = 0, vjust = 1,
           label = gain_lose_label(ps, "subjects"), size = 2.4,
           family = PAPER_FONT, colour = "grey20", lineheight = 1.05) +
  coord_cartesian(xlim = sub_lim, ylim = sub_lim) +
  scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) +
  labs(x = "Underliers per subject, 55 PCs", y = "Underliers per subject, 2005 PCs") +
  theme_paper()

# ------------------------------------------------------------------- panel d
pg <- D$per_gene
gene_max <- max(pg$n_pc55, pg$n_pc2005) + 1
pd <- ggplot(pg, aes(n_pc55 + 1, n_pc2005 + 1)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 0.5, alpha = 0.3, colour = unname(OI["vermillion"])) +
  annotate("text", x = 1, y = gene_max, hjust = 0, vjust = 1,
           label = gain_lose_label(pg, "genes"), size = 2.4,
           family = PAPER_FONT, colour = "grey20", lineheight = 1.05) +
  scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
  coord_cartesian(xlim = c(1, gene_max), ylim = c(1, gene_max)) +
  labs(x = "Underliers per gene + 1, 55 PCs", y = "Underliers per gene + 1, 2005 PCs") +
  theme_paper()

# ------------------------------------------------------------------- panel e
# Do the calls unique to one setting carry rare variants at a higher rate?
pe_dat <- filter(D$enrich, cutoff == PRIMARY_CUT) |>
  mutate(variant_lab = factor(ifelse(variant_class == "PTV",
                                     "Rare PTV in gene", "Rare SV within 100 kb"),
                              levels = c("Rare PTV in gene", "Rare SV within 100 kb")))
pe <- ggplot(pe_dat, aes(category, odds_ratio, colour = category)) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.35) +
  geom_errorbar(aes(ymin = ci_95_lo, ymax = ci_95_hi), width = 0.14, linewidth = 0.55) +
  geom_point(size = 2.2) +
  facet_wrap(~variant_lab, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = CAT_COLORS, guide = "none") +
  scale_x_discrete(labels = CAT_SHORT) +
  scale_y_log10() +
  labs(x = NULL, y = "Odds ratio vs non-underlier pairs") +
  theme_paper() +
  theme(axis.text.x = element_text(size = 7.5))

# ------------------------------------------------------------------- panel f
# Unique calls scored under the other PC setting.
pf <- ggplot(D$z_unique, aes(z_other, fill = category, colour = category)) +
  geom_vline(xintercept = -PRIMARY_CUT, linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  scale_fill_manual(values = CAT_COLORS[c(LAB55, LAB20)], name = NULL) +
  scale_colour_manual(values = CAT_COLORS[c(LAB55, LAB20)], name = NULL) +
  coord_cartesian(xlim = c(-6, 2)) +
  labs(x = "Z-score under the other PC setting", y = "Density") +
  theme_paper() +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

# ---------------------------------------------------------------- assemble
.paper_plots <- list(a = pa, b = pb, c = pc, d = pd, e = pe, f = pf)
for (nm in names(.paper_plots)) {
  paper_ggsave(.paper_plots[[nm]],
               file.path(PANEL_DIR, paste0("pc55_vs_pc2005_panel_", nm)),
               width = NATURE_DOUBLE_IN / 2, height = 2.7)
}

# Explicit 12-column design rather than (pa | pb) / ...: panel e carries two
# facets plus three discrete labels and needs more than half the row.
LAYOUT <- "
AAAAAABBBBBB
CCCCCCDDDDDD
EEEEEEEFFFFF
"
fig <- (pa + pb + pc + pd + pe + pf) +
  plot_layout(design = LAYOUT) +
  plot_annotation(tag_levels = "a") &
  theme_paper_tag()
paper_ggsave(fig, file.path(OUT_DIR, "Supplementary_Figure_PC55_vs_PC2005_outlier_overlap"),
             width = NATURE_DOUBLE_IN, height = 8.6)

# ---------------------------------------------------------------- source data
write_tsv(ov, file.path(OUT_DIR, "pc55_vs_pc2005_overlap_by_cutoff.tsv"))
write_tsv(D$enrich, file.path(OUT_DIR, "pc55_vs_pc2005_unique_call_variant_enrichment.tsv"))
write_tsv(D$per_gene, file.path(OUT_DIR, "pc55_vs_pc2005_underliers_per_gene.tsv"))
write_tsv(D$per_subject, file.path(OUT_DIR, "pc55_vs_pc2005_underliers_per_subject.tsv"))
write_tsv(
  D$z_unique |> group_by(category) |>
    summarise(n = n(), median_z_own = median(z_own), median_z_other = median(z_other),
              frac_z_other_lt_2 = mean(z_other < -2), .groups = "drop"),
  file.path(OUT_DIR, "pc55_vs_pc2005_unique_call_z_summary.tsv"))

cat(sprintf("\nat z < -%s: %s shared, %s unique to 55 PCs, %s unique to 2005 PCs (Jaccard %.3f)\n",
            PRIMARY_CUT, comma(ovp$n_shared), comma(ovp$n_only_pc55),
            comma(ovp$n_only_pc2005), ovp$jaccard))
message("wrote ", OUT_DIR)
