#!/usr/bin/env Rscript
# 2x2 gene-level PTV carrier frequency vs RNA underlier frequency,
# RESTRICTED to ClinGen haploinsufficiency (HI) "Sufficient Evidence" genes.
#
# Identical logic to plot_ptv_vs_underlier_gene_level_2x2.R, but the gene set is
# limited to the ClinGen HI dosage-sensitivity gene list so we can see whether the
# PTV vs underlier correlation differs among established dominant/HI disease genes.
#
# Outputs:
#   ptv_vs_underlier_gene_level_2x2_clingen_hi_haplo_z2_z3_z4.png / .pdf
#   ptv_vs_underlier_gene_level_2x2_clingen_hi_correlations.tsv

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
OUT_DIR <- FIG3C


CLINGEN_CSV <- Sys.getenv(
  "CLINGEN_CSV",
  unset = CLINGEN_HI
)
# gencode v30 symbol map (matches the expression-data gene IDs); 2 cols, no header:
#   <ensembl_id_no_version>\t<symbol>
GENCODE_MAP <- Sys.getenv("GENCODE_MAP", unset = "gencode_v30_gene_symbols.tsv")

ORDER <- c("haplo", "z_-2", "z_-3", "z_-4")
LABELS <- c(
  "haplo" = "Haplo (log2CPM-only)",
  "z_-2"  = "Strict z<-2 + log2CPM",
  "z_-3"  = "Strict z<-3 + log2CPM",
  "z_-4"  = "Strict z<-4 + log2CPM"
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
TIER_COLORS <- c(
  "Constrained" = "#B2182B",
  "Mid-low"     = "#EF8A62",
  "Mid-high"    = "#67A9CF",
  "Tolerant"    = "#2166AC",
  "No LOEUF"    = "grey70"
)

# ---- ClinGen HI gene list --------------------------------------------------
# The "Gene Symbol /Region Name" column encodes symbol + HGNC id, e.g.
# "AAGABHGNC:25662"; take everything before "HGNC". Keep only Type == Gene rows
# (drop "Region" CNV intervals, which have no single symbol).
cat("Loading ClinGen HI list:", CLINGEN_CSV, "\n")
clingen_raw <- read_csv(CLINGEN_CSV, show_col_types = FALSE)
sym_col <- names(clingen_raw)[grepl("Gene Symbol", names(clingen_raw))][1]
clingen_genes <- clingen_raw |>
  filter(.data[["Type"]] == "Gene") |>
  mutate(symbol = str_trim(str_remove(.data[[sym_col]], "HGNC.*$"))) |>
  filter(symbol != "") |>
  distinct(symbol) |>
  pull(symbol)
clingen_genes <- toupper(clingen_genes)
cat("  ClinGen HI 'Sufficient Evidence' genes:", length(clingen_genes), "\n")

# Map ClinGen symbols -> Ensembl gene_nv via the gencode v30 map used by the
# expression data (recovers more genes than the loeuf_lookup symbol column).
cat("Loading gencode v30 symbol map:", GENCODE_MAP, "\n")
gencode_map <- read_tsv(GENCODE_MAP, show_col_types = FALSE,
                        col_names = c("gene_nv", "symbol")) |>
  mutate(gene_nv = str_remove(gene_nv, "\\..*$"), symbol_uc = toupper(symbol))
clingen_nv <- gencode_map |>
  filter(symbol_uc %in% clingen_genes) |>
  distinct(gene_nv) |>
  pull(gene_nv)
cat("  ClinGen HI symbols mapped to Ensembl IDs:", length(clingen_nv), "\n")

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

cat("Loading expr_cpm_adj.rds (eligible gene universe) ...\n")
.expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
eligible_nv <- unique(str_remove(rownames(.expr_cpm_adj), "\\..*$"))
rm(.expr_cpm_adj); gc(verbose = FALSE)
cat("  eligible (expression-tested) genes:", length(eligible_nv), "\n")

TPM_MIN <- 1
tpm_med <- read_tsv("gene_median_tpm_kept.tsv", show_col_types = FALSE) |>
  distinct(gene_nv, .keep_all = TRUE)
tpm_ok_nv <- tpm_med |> filter(median_tpm > TPM_MIN) |> pull(gene_nv)
cat("  genes with median TPM >", TPM_MIN, ":", length(tpm_ok_nv), "\n")

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
    mutate(
      pct_ptv    = 100 * n_ptv_carriers / n_kept,
      pct_under  = 100 * n_underlier_subjects / n_kept,
      LOEUF_tier = factor(
        ifelse(is.na(LOEUF_tier), "No LOEUF", as.character(LOEUF_tier)),
        levels = TIER_LEVELS
      ),
      mode = m
    )
}

all_df <- bind_rows(lapply(ORDER, function(m) {
  cat("Computing gene-level table for:", m, "\n")
  gene_table_for_mode(m)
})) |>
  mutate(mode = factor(mode, levels = ORDER, labels = LABELS[ORDER]))

# ---- Restrict to ClinGen HI genes (by Ensembl gene_nv) ---------------------
all_df_full <- all_df
all_df <- all_df |>
  filter(gene_nv %in% clingen_nv)

matched_nv <- sort(unique(all_df$gene_nv))
cat("  ClinGen HI genes present in plotted table (any mode):", length(matched_nv),
    "of", length(clingen_nv), " mapped (", length(clingen_genes), " symbols)\n")

# Per-gene table export (all modes), with a resolved HGNC symbol (gencode v30,
# falling back to loeuf symbol then Ensembl id) so downstream figures can pick
# and label ClinGen HI axis-hugging genes.
gencode_sym <- gencode_map |>
  distinct(gene_nv, .keep_all = TRUE) |>
  select(gene_nv, gencode_symbol = symbol)
clingen_gene_table <- all_df |>
  left_join(gencode_sym, by = "gene_nv") |>
  mutate(symbol = coalesce(na_if(gencode_symbol, ""), na_if(symbol, ""), gene_nv)) |>
  select(mode, gene_nv, symbol, LOEUF, LOEUF_tier,
         n_ptv_carriers, n_underlier_subjects, pct_ptv, pct_under) |>
  arrange(mode, desc(pct_under))
write_tsv(clingen_gene_table, "ptv_vs_underlier_gene_level_2x2_clingen_hi_genes.tsv")
cat("  Wrote per-gene table: ptv_vs_underlier_gene_level_2x2_clingen_hi_genes.tsv\n")

spearman_fit <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(rho = NA_real_, n_pairs = n, logp = NA_real_, p_value = NA_real_))
  }
  r <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  tstat <- r * sqrt((n - 2) / (1 - r^2))
  logp  <- log(2) + stats::pt(-abs(tstat), df = n - 2, log.p = TRUE)
  tibble(rho = r, n_pairs = n, logp = logp, p_value = exp(logp))
}
fmt_logp <- function(logp) {
  if (is.na(logp)) return("NA")
  log10p <- logp / log(10)
  e <- floor(log10p)
  m <- 10^(log10p - e)
  sprintf("%.1fe%d", m, as.integer(e))
}

summ <- all_df |>
  group_by(mode) |>
  summarise(
    spearman_fit(pct_ptv[n_ptv_carriers > 0], pct_under[n_ptv_carriers > 0]),
    pearson_all  = suppressWarnings(stats::cor(pct_ptv, pct_under, method = "pearson")),
    n_genes      = n(),
    n_ptv_genes  = sum(n_ptv_carriers > 0),
    n_under_only = sum(n_ptv_carriers == 0 & n_underlier_subjects > 0),
    y_top        = max(pct_under),
    .groups = "drop"
  ) |>
  mutate(
    p_label = vapply(logp, fmt_logp, character(1)),
    lbl = sprintf("rho|PTV>0 = %s\np = %s\nPearson(all) = %s\n%s PTV genes / %s total\n%s underlier-only (0 PTV)",
                  ifelse(is.na(rho), "NA", sprintf("%.3f", rho)),
                  p_label,
                  ifelse(is.na(pearson_all), "NA", sprintf("%.3f", pearson_all)),
                  comma(n_ptv_genes), comma(n_genes), comma(n_under_only)),
    x   = max(all_df$pct_ptv) * 0.98,
    y   = y_top
  )

cat("\nPer-definition correlation (gene-level, ClinGen HI only):\n")
print(summ |> select(mode, rho, p_label, pearson_all, n_ptv_genes, n_genes, n_under_only))

write_tsv(
  summ |> select(mode, rho, logp, p_value, pearson_all, n_pairs, n_genes, n_ptv_genes, n_under_only),
  "ptv_vs_underlier_gene_level_2x2_clingen_hi_correlations.tsv"
)

x_max <- max(all_df$pct_ptv)

p <- ggplot(all_df, aes(x = pct_ptv, y = pct_under)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey55", linewidth = 0.35) +
  geom_point(aes(color = LOEUF_tier), alpha = 0.6, size = 1.5, stroke = 0) +
  geom_text(
    data = summ, aes(x = x, y = y, label = lbl),
    hjust = 1, vjust = 1, size = 2.7, color = "grey20", lineheight = 0.95
  ) +
  facet_wrap(~mode, ncol = 2, scales = "free_y") +
  scale_color_manual(values = TIER_COLORS, name = "LOEUF tier", drop = FALSE) +
  scale_x_continuous(labels = label_percent(scale = 1, accuracy = 0.01),
                     limits = c(0, x_max),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = label_percent(scale = 1, accuracy = 0.01),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = "ClinGen HI genes only: gene-level rare PTV (AF<0.5%) carrier vs RNA underlier frequency",
    subtitle = str_wrap(
      sprintf("Restricted to ClinGen haploinsufficiency 'Sufficient Evidence' genes (%s of %s symbols plotted; rest not expression-tested or median TPM<=1). Each point = one such gene. Cohort n=%s. Shared x-axis, panel-specific y-axis; rho computed on genes with >=1 PTV carrier.",
              comma(length(matched_nv)), comma(length(clingen_genes)), comma(n_kept)),
      width = 110
    ),
    caption = str_wrap(
      "Dashed line = equal frequency (y = x). Compare rho/Pearson here against the all-gene 2x2 to see whether restricting to dominant HI disease genes strengthens the PTV vs underlier relationship.",
      width = 130
    ),
    x = "% of subjects carrying a rare PTV (AF<0.5%) in gene",
    y = "% of subjects who are RNA underliers for gene"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "grey25", lineheight = 1.15),
    plot.caption  = element_text(size = 7.5, color = "grey35", hjust = 0),
    strip.text    = element_text(face = "bold", size = 9),
    legend.position = "right",
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )

out_png <- "ptv_vs_underlier_gene_level_2x2_clingen_hi_haplo_z2_z3_z4.png"
out_pdf <- "ptv_vs_underlier_gene_level_2x2_clingen_hi_haplo_z2_z3_z4.pdf"
ggsave(file.path(OUT_DIR, out_png), p, width = 13.0, height = 9.8, dpi = 220)
ggsave(file.path(OUT_DIR, out_pdf), p, width = 13.0, height = 9.8, device = cairo_pdf)
cat("Wrote:", normalizePath(out_png), "\n")
cat("Wrote:", normalizePath(out_pdf), "\n")
