#!/usr/bin/env Rscript
# Gene-level PTV carrier frequency vs RNA underlier frequency.
#
# Message: underliers reveal signal that PTV genotypes do not. Each point is a
# gene; genes high on the y-axis with near-zero x-axis (many underlier subjects,
# few/no PTV carriers) are underlier signal that variant calls cannot reveal.
#
# X-axis: % of subjects carrying a PTV in that gene.
# Y-axis: % of subjects who are RNA underliers for that gene.
#
# Underlier definition picked via env UNDERLIER_LABEL: z_-2 | z_-3 | z_-4 | haplo
#   (default haplo, matching ptv_vs_underlier_burden_haplo.png).
#
# Outputs:
#   ptv_vs_underlier_gene_level_{LABEL}.tsv
#   ptv_vs_underlier_gene_level_{LABEL}.png / .pdf

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG3X


MODE_LABEL <- Sys.getenv("UNDERLIER_LABEL", unset = "haplo")
if (!MODE_LABEL %in% c("z_-2", "z_-3", "z_-4", "haplo")) {
  stop("UNDERLIER_LABEL must be one of z_-2, z_-3, z_-4, haplo (got ", MODE_LABEL, ")")
}

under_rds <- if (MODE_LABEL == "haplo") {
  file.path(DATA_UNDERLIER, "z_score_underlier_haplo_clean.rds")
} else {
  file.path(DATA_UNDERLIER, paste0("z_score_underlier_", MODE_LABEL, "_clean.rds"))
}
mode_txt <- if (MODE_LABEL == "haplo") {
  "haplo (log2CPM only)"
} else {
  paste0("strict z<", sub("^z_", "", MODE_LABEL), " + log2CPM")
}

PTV_TSV <- Sys.getenv(
  "PTV_TSV",
  unset = "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/ptv_all.tsv"
)

TIER_LEVELS <- c("Constrained", "Mid-low", "Mid-high", "Tolerant", "No LOEUF")
TIER_COLORS <- c(
  "Constrained" = "#B2182B",
  "Mid-low"     = "#EF8A62",
  "Mid-high"    = "#67A9CF",
  "Tolerant"    = "#2166AC",
  "No LOEUF"    = "grey70"
)

cat("Loading kept_subjects.rds ...\n")
kept_subjects <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept_subjects)

cat("Loading PTV table:", PTV_TSV, "\n")
ptv <- read_tsv(PTV_TSV, show_col_types = FALSE) |>
  filter(sample %in% kept_subjects)

cat("Loading underliers:", under_rds, "\n")
under <- read_rds(under_rds) |>
  filter(subjectid %in% kept_subjects)

cat("Loading loeuf_lookup.rds ...\n")
loeuf_lookup <- read_rds(LOEUF_LOOKUP)

# Eligible gene universe = expression-tested genes (rownames of expr_cpm_adj),
# the SAME set used by the LOEUF/underlier figures. PTVs occur genome-wide, but
# a gene that was never expression-tested cannot be an RNA underlier, so points
# are restricted to this universe (otherwise PTV-only genes inflate the count).
cat("Loading expr_cpm_adj.rds (eligible gene universe) ...\n")
.expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
eligible_nv <- unique(str_remove(rownames(.expr_cpm_adj), "\\..*$"))
rm(.expr_cpm_adj); gc(verbose = FALSE)
cat("  eligible (expression-tested) genes:", length(eligible_nv), "\n")

# Expression filter: keep genes with median TPM > 1 across kept subjects.
# Per-gene median TPM precomputed from gene-tpm.parquet (TOR ids mapped to NWD,
# restricted to kept_subjects) -> gene_median_tpm_kept.tsv.
TPM_MIN <- 1
tpm_med <- read_tsv(file.path(DATA_CORE, "gene_median_tpm_kept.tsv"), show_col_types = FALSE) |>
  distinct(gene_nv, .keep_all = TRUE)
tpm_ok_nv <- tpm_med |> filter(median_tpm > TPM_MIN) |> pull(gene_nv)
cat("  genes with median TPM >", TPM_MIN, ":", length(tpm_ok_nv), "\n")

# --- Per-gene PTV carrier count (distinct subjects with >=1 PTV in gene) ---
ptv_gene <- ptv |>
  mutate(gene_nv = str_remove(gene, "\\..*$")) |>
  distinct(subjectid = sample, gene_nv) |>
  count(gene_nv, name = "n_ptv_carriers")

# --- Per-gene underlier count (distinct subjects who are underliers in gene) ---
under_gene <- under |>
  mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
  distinct(subjectid, gene_nv) |>
  count(gene_nv, name = "n_underlier_subjects")

loeuf_nv <- loeuf_lookup |>
  distinct(gene_id_no_version, .keep_all = TRUE) |>
  select(gene_nv = gene_id_no_version, symbol, LOEUF, LOEUF_tier)

# HGNC symbol map (GENCODE v30, same annotation the gene IDs come from).
# Built once via: zcat gencode.v30.annotation.gff3.gz | awk ... > gencode_v30_gene_symbols.tsv
GENCODE_SYMS <- GENCODE_SYM
gencode_nv <- read_tsv(
  GENCODE_SYMS, col_names = c("gene_nv", "gencode_symbol"), show_col_types = FALSE
) |>
  distinct(gene_nv, .keep_all = TRUE)

genes <- full_join(ptv_gene, under_gene, by = "gene_nv") |>
  mutate(
    n_ptv_carriers       = replace_na(n_ptv_carriers, 0L),
    n_underlier_subjects = replace_na(n_underlier_subjects, 0L)
  ) |>
  filter(gene_nv %in% eligible_nv) |>
  filter(gene_nv %in% tpm_ok_nv) |>
  filter(n_ptv_carriers > 0 | n_underlier_subjects > 0) |>
  left_join(loeuf_nv, by = "gene_nv") |>
  left_join(gencode_nv, by = "gene_nv") |>
  mutate(
    # Prefer GENCODE HGNC symbol; fall back to LOEUF-table symbol, then Ensembl ID.
    symbol = coalesce(
      na_if(gencode_symbol, ""),
      na_if(symbol, ""),
      gene_nv
    ),
    pct_ptv   = 100 * n_ptv_carriers / n_kept,
    pct_under = 100 * n_underlier_subjects / n_kept,
    LOEUF_tier = factor(
      ifelse(is.na(LOEUF_tier), "No LOEUF", as.character(LOEUF_tier)),
      levels = TIER_LEVELS
    ),
    gene_label = symbol
  )

out_tsv <- file.path(OUT_DIR, paste0("ptv_vs_underlier_gene_level_", MODE_LABEL, ".tsv"))
write_tsv(
  genes |>
    select(
      gene_nv, gene_label, symbol, LOEUF, LOEUF_tier,
      n_ptv_carriers, n_underlier_subjects, pct_ptv, pct_under
    ) |>
    arrange(desc(pct_under)),
  out_tsv
)
cat("Wrote:", normalizePath(out_tsv), "(", nrow(genes), "genes )\n")

# Spearman rho restricted to genes with >=1 PTV carrier. Over the full table
# both axes are integer carrier counts / n_kept, so a huge block of 0-underlier
# genes (many of which carry PTVs) ties at the y-rank floor and drags a
# whole-table rank correlation spuriously negative at strict z cutoffs. The
# question of interest -- among genes that DO carry PTVs, do more carriers track
# more underliers? -- is answered on the PTV>0 subset (weakly positive). Pearson
# over all genes (~0) is reported alongside as a no-trend reference.
# Large-sample t approximation (exact = FALSE); p computed on the natural-log
# scale to avoid double-precision underflow, then shown in scientific notation.
.ptvpos  <- genes$n_ptv_carriers > 0
.n_pairs <- sum(.ptvpos)
rho      <- suppressWarnings(
  stats::cor(genes$pct_ptv[.ptvpos], genes$pct_under[.ptvpos], method = "spearman")
)
.tstat   <- rho * sqrt((.n_pairs - 2) / (1 - rho^2))
.logp    <- log(2) + stats::pt(-abs(.tstat), df = .n_pairs - 2, log.p = TRUE)
p_value  <- exp(.logp)
p_txt <- {
  log10p <- .logp / log(10)
  e <- floor(log10p)
  m <- 10^(log10p - e)
  sprintf("%.1fe%d", m, as.integer(e))
}
pearson_all <- suppressWarnings(
  stats::cor(genes$pct_ptv, genes$pct_under, method = "pearson")
)
cat(sprintf("Spearman rho|PTV>0 = %.4f (n=%d), p = %s; Pearson(all) = %.4f\n",
            rho, .n_pairs, p_txt, pearson_all))

# Genes where the underlier signal is "invisible" to PTVs: high underlier rate,
# zero PTV carriers. Label the most extreme of these.
underlier_only <- genes |>
  filter(n_ptv_carriers == 0, n_underlier_subjects > 0) |>
  arrange(desc(pct_under)) |>
  slice_head(n = 15)

n_under_only  <- sum(genes$n_ptv_carriers == 0 & genes$n_underlier_subjects > 0)
n_ptv_only    <- sum(genes$n_ptv_carriers > 0 & genes$n_underlier_subjects == 0)
n_both        <- sum(genes$n_ptv_carriers > 0 & genes$n_underlier_subjects > 0)

sub_txt <- str_wrap(
  sprintf(
    "Each point = an expression-tested gene, median TPM>1 (n=%s with >=1 PTV carrier or >=1 underlier; cohort n=%s). Spearman rho among genes with >=1 PTV carrier (n=%s) = %.3f (p %s); Pearson over all genes = %.3f (~0, no monotonic trend). Underlier-only genes (0 PTV carriers): %s; PTV-only: %s; both: %s.",
    comma(nrow(genes)), comma(n_kept), comma(.n_pairs), rho, p_txt, pearson_all,
    comma(n_under_only), comma(n_ptv_only), comma(n_both)
  ),
  width = 95
)
cap_txt <- str_wrap(
  paste0(
    "Genes high on the y-axis with near-zero x (many RNA underlier subjects, few/no PTV carriers) ",
    "carry expression-outlier signal that PTV genotypes do not reveal. ",
    "Dashed line = equal frequency (y = x)."
  ),
  width = 110
)

p <- ggplot(genes, aes(x = pct_ptv, y = pct_under)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_point(aes(color = LOEUF_tier), alpha = 0.55, size = 1.6, stroke = 0) +
  ggrepel::geom_text_repel(
    data = underlier_only,
    aes(label = gene_label),
    size = 2.6, color = "grey20",
    min.segment.length = 0, segment.size = 0.2,
    box.padding = 0.4, max.overlaps = Inf, seed = 17
  ) +
  scale_color_manual(values = TIER_COLORS, name = "LOEUF tier", drop = FALSE) +
  scale_x_continuous(labels = label_percent(scale = 1, accuracy = 0.01),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(labels = label_percent(scale = 1, accuracy = 0.01),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = sprintf("Gene-level rare PTV (AF<0.5%%) carrier vs RNA underlier frequency (%s)", mode_txt),
    subtitle = sub_txt,
    caption = cap_txt,
    x = "% of subjects carrying a rare PTV (AF<0.5%) in gene",
    y = "% of subjects who are RNA underliers for gene"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "grey25"),
    plot.caption  = element_text(size = 7.5, color = "grey35", hjust = 0),
    legend.position = "right",
    plot.margin = margin(5.5, 14, 5.5, 5.5)
  )

out_png <- paste0("ptv_vs_underlier_gene_level_", MODE_LABEL, ".png")
out_pdf <- paste0("ptv_vs_underlier_gene_level_", MODE_LABEL, ".pdf")
ggsave(file.path(OUT_DIR, out_png), p, width = 8.5, height = 7.0, dpi = 220)
ggsave(file.path(OUT_DIR, out_pdf), p, width = 8.5, height = 7.0, device = cairo_pdf)
cat("Wrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")
