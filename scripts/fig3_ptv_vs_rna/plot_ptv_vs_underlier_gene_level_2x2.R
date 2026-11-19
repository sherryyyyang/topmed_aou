#!/usr/bin/env Rscript
# Gene-level PTV carrier frequency vs RNA underlier frequency.
# Primary underlier definition only: strict z < -3 + log2CPM.
#
# Companion to plot_ptv_vs_underlier_gene_level.R (single-definition sweep)
# and the per-subject burden plots.
#
# Each point = a gene. X = % of subjects carrying a PTV in that gene.
# Y = % of subjects with RNA underliers in that gene. Genes high on the
# y-axis with near-zero x carry expression-outlier signal PTVs do not reveal.
# Points are coloured by gene-level NMD escape status from Ensembl VEP
# NMD plugin annotations (see nmd_escape_vep/).
#
# Outputs:
#   ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4.png / .pdf
#   (filename kept so Figure 3 assembly still finds this panel)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})
ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG3A


ORDER <- c("z_-3")
LABELS <- c(
  "z_-3"  = "Strict z<-3 + log2CPM"
)
under_rds_for <- function(m) {
  file.path(DATA_UNDERLIER, if (m == "haplo") "z_score_underlier_haplo_clean.rds"
            else paste0("z_score_underlier_", m, "_clean.rds"))
}

PTV_TSV <- Sys.getenv(
  "PTV_TSV",
  unset = "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/ptv_all.tsv"
)

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")

# Gene-level NMD escape status for Fig. 3b colouring.
# Variant-level VEP NMD plugin calls (selected transcript: MANE > CANONICAL)
# are aggregated per gene:
#   NMD-escaping   = >=1 AF<0.1% HC PTV annotated NMD_escaping
#   NMD-triggering = >=1 AF<0.1% HC PTV annotated NMD_triggering, none escaping
#   Unknown        = no AF<0.1% HC PTV with a conclusive NMD call
#                   (includes underlier-only genes and inconclusive-only)
NMD_CLASS_LEVELS <- c("NMD-escaping", "NMD-triggering", "Unknown")
# NMD-escaping: Okabe–Ito vermillion (#D55E00). Avoid semi-transparent points
# with a solid legend key (alpha 0.55 on white washed vermillion to yellow-peach
# while the swatch stayed true orange). One opaque vermillion for geom + legend.
NMD_CLASS_COLORS <- c(
  "NMD-escaping"   = unname(OI["vermillion"]),  # #D55E00
  "NMD-triggering" = unname(OI["blue"]),        # #0072B2
  "Unknown"        = "grey55"
)
NMD_POINT_ALPHA <- 1
NMD_BY_VARIANT <- file.path(
  ROOT, "nmd_escape_vep", "fig3b_ptv_nmd_escape_by_variant.tsv"
)
# Match the Fig3b AF panel cutoff used in NMD_ESCAPE_SUMMARY.txt.
NMD_AF_MAX <- 0.001

load_nmd_gene_map <- function() {
  cat("Loading NMD escape annotations:", NMD_BY_VARIANT, "\n")
  nmd <- read_tsv(NMD_BY_VARIANT, show_col_types = FALSE) |>
    mutate(
      gene_nv = str_remove(as.character(gene), "\\..*$"),
      af = as.numeric(af),
      nmd_status = as.character(nmd_status)
    ) |>
    filter(!is.na(af), af < NMD_AF_MAX)
  cat("  AF<", NMD_AF_MAX, " variants:", nrow(nmd), "\n")
  print(as.data.frame(sort(table(nmd$nmd_status), decreasing = TRUE)))
  out <- nmd |>
    group_by(gene_nv) |>
    summarise(
      n_escape = sum(nmd_status == "NMD_escaping", na.rm = TRUE),
      n_trigger = sum(nmd_status == "NMD_triggering", na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      nmd_class = case_when(
        n_escape > 0 ~ "NMD-escaping",
        n_trigger > 0 ~ "NMD-triggering",
        TRUE ~ "Unknown"
      ),
      nmd_class = factor(nmd_class, levels = NMD_CLASS_LEVELS)
    ) |>
    select(gene_nv, nmd_class)
  cat("  gene-level NMD class counts:\n")
  print(as.data.frame(sort(table(out$nmd_class), decreasing = TRUE)))
  out
}

cat("Loading kept_subjects.rds ...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept_subjects)

cat("Loading PTV table:", PTV_TSV, "\n")
ptv <- read_tsv(PTV_TSV, show_col_types = FALSE) |>
  filter(sample %in% kept_subjects)

cat("Loading loeuf_lookup.rds ...\n")
loeuf_lookup <- read_rds(LOEUF_LOOKUP)
loeuf_nv <- loeuf_lookup |>
  distinct(gene_id_no_version, .keep_all = TRUE) |>
  select(gene_nv = gene_id_no_version, symbol, LOEUF, LOEUF_tier)

# Eligible gene universe = expression-tested genes (rownames of expr_cpm_adj).
# PTVs are called genome-wide, but a gene that was never expression-tested
# cannot be an RNA underlier; restricting to this set avoids PTV-only genes
# (~5,645 not expression-tested) being plotted as false y=0 zeros.
cat("Loading expr_cpm_adj.rds (eligible gene universe) ...\n")
.expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
eligible_nv <- unique(str_remove(rownames(.expr_cpm_adj), "\\..*$"))
rm(.expr_cpm_adj); gc(verbose = FALSE)
cat("  eligible (expression-tested) genes:", length(eligible_nv), "\n")

# Expression filter: keep genes with median TPM > 1 across kept subjects.
# Per-gene median TPM precomputed from gene-tpm.parquet -> gene_median_tpm_kept.tsv.
TPM_MIN <- 1
tpm_med <- read_tsv(file.path(DATA_CORE, "gene_median_tpm_kept.tsv"), show_col_types = FALSE) |>
  distinct(gene_nv, .keep_all = TRUE)
tpm_ok_nv <- tpm_med |> filter(median_tpm > TPM_MIN) |> pull(gene_nv)
cat("  genes with median TPM >", TPM_MIN, ":", length(tpm_ok_nv), "\n")

# PTV carrier count per gene (independent of underlier definition).
ptv_gene <- ptv |>
  mutate(gene_nv = str_remove(gene, "\\..*$")) |>
  distinct(subjectid = sample, gene_nv) |>
  count(gene_nv, name = "n_ptv_carriers")

gene_table_for_mode <- function(m) {
  under <- read_rds(under_rds_for(m)) |>
    filter(subjectid %in% kept_subjects)
  under_gene <- under |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
    distinct(subjectid, gene_nv) |>
    count(gene_nv, name = "n_underlier_subjects")

  full_join(ptv_gene, under_gene, by = "gene_nv") |>
    mutate(
      n_ptv_carriers       = replace_na(n_ptv_carriers, 0L),
      n_underlier_subjects = replace_na(n_underlier_subjects, 0L)
    ) |>
    filter(gene_nv %in% eligible_nv) |>
    filter(gene_nv %in% tpm_ok_nv) |>
    filter(n_ptv_carriers > 0 | n_underlier_subjects > 0) |>
    left_join(loeuf_nv, by = "gene_nv") |>
    left_join(nmd_map, by = "gene_nv") |>
    mutate(
      pct_ptv    = 100 * n_ptv_carriers / n_kept,
      pct_under  = 100 * n_underlier_subjects / n_kept,
      LOEUF_tier = factor(
        ifelse(is.na(LOEUF_tier), "No LOEUF", as.character(LOEUF_tier)),
        levels = TIER_LEVELS
      ),
      nmd_class = factor(
        ifelse(is.na(nmd_class), "Unknown", as.character(nmd_class)),
        levels = NMD_CLASS_LEVELS
      ),
      mode = m
    )
}

# Gene-level NMD map from VEP NMD plugin (AF<0.1% HC PTVs).
nmd_map <- load_nmd_gene_map()

all_df <- bind_rows(lapply(ORDER, function(m) {
  cat("Computing gene-level table for:", m, "\n")
  gene_table_for_mode(m)
})) |>
  mutate(mode = factor(mode, levels = ORDER, labels = LABELS[ORDER])) |>
  # Draw Unknown underneath; NMD-escaping on top where points overlap.
  arrange(nmd_class == "Unknown", nmd_class == "NMD-triggering", nmd_class)

# Per-facet summary: Spearman rho + p-value + underlier-only gene count.
# Large-sample (exact = FALSE) Spearman p-value via the t approximation:
#   t = rho * sqrt((n-2)/(1-rho^2)),  p = 2 * P(T_{n-2} < -|t|).
# With ~12k genes the p-value underflows double precision (cor.test would
# return 0), so compute it on the natural-log scale and report the true value
# in scientific notation. This matches cor.test(method="spearman", exact=FALSE).
spearman_fit <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  r <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  tstat <- r * sqrt((n - 2) / (1 - r^2))
  logp  <- log(2) + stats::pt(-abs(tstat), df = n - 2, log.p = TRUE)  # natural log
  tibble(rho = r, n_pairs = n, logp = logp, p_value = exp(logp))
}
# Format a natural-log p-value as "m.me-E" without underflowing.
fmt_logp <- function(logp) {
  log10p <- logp / log(10)
  e <- floor(log10p)
  m <- 10^(log10p - e)
  sprintf("%.1fe%d", m, as.integer(e))
}

# rho is computed on genes with >=1 PTV carrier only (see note below); Pearson
# over all genes is reported as a no-trend reference.
summ <- all_df |>
  group_by(mode) |>
  summarise(
    spearman_fit(pct_ptv[n_ptv_carriers > 0], pct_under[n_ptv_carriers > 0]),
    pearson_all  = suppressWarnings(stats::cor(pct_ptv, pct_under, method = "pearson")),
    n_genes      = n(),
    n_ptv_genes  = sum(n_ptv_carriers > 0),
    n_under_only = sum(n_ptv_carriers == 0 & n_underlier_subjects > 0),
    .groups = "drop"
  ) |>
  mutate(p_label = fmt_logp(logp))

cat("\nPer-definition correlation (gene-level):\n")
print(summ |> select(mode, rho, p_label, pearson_all, n_ptv_genes, n_genes, n_under_only))

x_max <- max(all_df$pct_ptv)

p <- ggplot(all_df, aes(x = pct_ptv, y = pct_under)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey55", linewidth = 0.35) +
  geom_point(aes(color = nmd_class), alpha = NMD_POINT_ALPHA, size = 1.3, stroke = 0) +
  scale_color_manual(values = NMD_CLASS_COLORS, name = "NMD status", drop = FALSE) +
  guides(color = guide_legend(override.aes = list(alpha = NMD_POINT_ALPHA, size = 2.4))) +
  scale_x_continuous(labels = label_percent(scale = 1, accuracy = 0.01),
                     limits = c(0, x_max),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = label_percent(scale = 1, accuracy = 0.01),
                     expand = expansion(mult = c(0.02, 0.10))) +
  coord_cartesian(clip = "off") +
  labs(
    x = "% of subjects carrying a rare PTV (AF<0.5%) in gene",
    y = "% of subjects with RNA\nunderliers in gene"
  ) +
  theme_paper() +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 4, r = 6, b = 2, l = 16, unit = "mm"),
    axis.title.y = element_text(
      angle = 90, lineheight = 0.92,
      margin = margin(r = 3, unit = "mm")
    )
  )

.paper_plots <- list(a = p)
out_stem <- file.path(OUT_DIR, "ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4")
paper_ggsave(p, out_stem, width = NATURE_DOUBLE_IN, height = 4.6)
if (!is_paper_assemble()) cat("Wrote:", out_stem, "\n")
