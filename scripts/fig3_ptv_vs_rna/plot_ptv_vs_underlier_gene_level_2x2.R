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
#
# Default colouring (FIG3B_COLOR=cnv, or unset): SegDup-pass genes with a
# primary gene-overlapping PASS muCNV SVTYPE (any AF) —
#   DEL / DUP / INV  = fig3b_segdup_svtype_genes.tsv (SegDup rule required;
#                      primary variant by |SVLEN| among gene-overlapping hits)
#   Other genes      = remainder (incl. curated-only w/o SegDup or w/o VCF hit)
# Point size (cnv mode): distinct gene-overlapping PASS muCNVs with >=1 alt
# carrier among RNA underliers for that gene (z<-3 panel; see SOURCE_fig3b_segdup_svtype).
#
# Alternate colouring: FIG3B_COLOR=nmd uses gene-level NMD escape status from
# Ensembl VEP NMD plugin annotations (see nmd_escape_vep/).
#
# Uniform style (no SV/CNV legend): NO_SV_COLOR=1 (default off) keeps all genes,
# single grey points, uniform size; output ..._no_sv_color (main panel unchanged).
#
# Test filter: FIG3B_OTHER_ONLY=1 (default off) drops DEL/DUP/INV points,
# plotting Other genes only. Output stem becomes
# ..._other_only_test so the main panel file is not overwritten.
#
# EXCLUDE_SV_CNV=1 (default off) removes genes in fig3b_segdup_svtype_genes.tsv
# (DEL/DUP/INV) and fig3b_cnv_multicopy_genes.tsv from the scatter, then replots
# with uniform grey points and recomputed Spearman. Output ..._excl_sv_cnv.
# (Do not combine EXCLUDE_SV_CNV with NO_SV_COLOR or FIG3B_OTHER_ONLY.)
#
# Outputs:
#   ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4.png / .pdf
#   (filename kept so Figure 3 assembly still finds this panel)
#   …_no_sv_color.png / .pdf when NO_SV_COLOR=1
#   …_other_only_test.png / .pdf when FIG3B_OTHER_ONLY=1
#   …_excl_sv_cnv.png / .pdf when EXCLUDE_SV_CNV=1

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

# Colour mode: "cnv" (default) or "nmd".
COLOR_MODE <- tolower(Sys.getenv("FIG3B_COLOR", unset = "cnv"))
if (!COLOR_MODE %in% c("cnv", "nmd")) {
  stop("FIG3B_COLOR must be 'cnv' or 'nmd' (got '", COLOR_MODE, "')")
}

# Test-only: keep Other genes and drop DEL/DUP/INV points.
# Default off so main Figure 3b colouring / outputs are unchanged.
OTHER_ONLY <- Sys.getenv("FIG3B_OTHER_ONLY", unset = "0") %in% c("1", "true", "TRUE", "yes", "YES")
if (OTHER_ONLY && !identical(COLOR_MODE, "cnv")) {
  stop("FIG3B_OTHER_ONLY=1 requires FIG3B_COLOR=cnv (got '", COLOR_MODE, "')")
}
NO_SV_COLOR <- Sys.getenv("NO_SV_COLOR", unset = "0") %in% c("1", "true", "TRUE", "yes", "YES")
EXCLUDE_SV_CNV <- Sys.getenv("EXCLUDE_SV_CNV", unset = "0") %in% c("1", "true", "TRUE", "yes", "YES")
if (NO_SV_COLOR && OTHER_ONLY) {
  stop("Set only one of NO_SV_COLOR=1 or FIG3B_OTHER_ONLY=1")
}
if (EXCLUDE_SV_CNV && (NO_SV_COLOR || OTHER_ONLY)) {
  stop("Set only one of EXCLUDE_SV_CNV=1, NO_SV_COLOR=1, or FIG3B_OTHER_ONLY=1")
}
if ((NO_SV_COLOR || EXCLUDE_SV_CNV) && !identical(COLOR_MODE, "cnv")) {
  stop("NO_SV_COLOR / EXCLUDE_SV_CNV require FIG3B_COLOR=cnv (got '", COLOR_MODE, "')")
}
UNIFORM_STYLE <- NO_SV_COLOR || EXCLUDE_SV_CNV
if (OTHER_ONLY) {
  cat("FIG3B_OTHER_ONLY=1: plotting Other genes only (SV/CNV points removed)\n")
}
if (NO_SV_COLOR) {
  cat("NO_SV_COLOR=1: uniform grey points, no SV/CNV colouring or size legend\n")
}
if (EXCLUDE_SV_CNV) {
  cat("EXCLUDE_SV_CNV=1: drop SegDup SV + CNV/multi-copy genes; uniform grey replot\n")
}

NO_SV_POINT_COLOR <- "grey55"
NO_SV_POINT_SIZE <- 1.3

# ---- Dosage colouring: SegDup-pass primary SVTYPE vs other (default) ------
CNV_CLASS_LEVELS <- c("DEL", "DUP", "INV", "Other genes")
CNV_CLASS_COLORS <- c(
  "DEL"          = unname(OI["blue"]),        # #0072B2
  "DUP"          = unname(OI["vermillion"]),  # #D55E00
  "INV"          = unname(OI["green"]),       # #009E73
  "Other genes"  = "grey70"
)
SEGDUP_SVTYPE_TSV <- file.path(ROOT, "data", "references", "fig3b_segdup_svtype_genes.tsv")
SEGDUP_HITS_TSV <- file.path(ROOT, "data", "references", "fig3b_segdup_sv_gene_hits.tsv")
SEGDUP_CARRIERS_TSV <- file.path(ROOT, "data", "references", "fig3b_segdup_sv_carriers.tsv")
CNV_MULTICOPY_TSV <- file.path(ROOT, "data", "references", "fig3b_cnv_multicopy_genes.tsv")

load_exclude_sv_cnv_genes <- function() {
  sv_nv <- character()
  if (file.exists(SEGDUP_SVTYPE_TSV)) {
    sv_nv <- read_tsv(SEGDUP_SVTYPE_TSV, show_col_types = FALSE) |>
      mutate(svtype_class = as.character(svtype_class)) |>
      filter(svtype_class %in% c("DEL", "DUP", "INV")) |>
      pull(gene_nv) |>
      unique()
  }
  cnv_nv <- character()
  if (file.exists(CNV_MULTICOPY_TSV)) {
    cnv_nv <- read_tsv(CNV_MULTICOPY_TSV, show_col_types = FALSE) |>
      pull(gene_nv) |>
      unique()
  } else {
    warning(
      paste(
        "Missing CNV/multi-copy list:", CNV_MULTICOPY_TSV,
        "\n  Run: Rscript scripts/fig3_ptv_vs_rna/build_cnv_multicopy_gene_set.R"
      ),
      call. = FALSE
    )
  }
  ex <- unique(c(sv_nv, cnv_nv))
  cat(sprintf(
    "  EXCLUDE_SV_CNV gene sets: DEL/DUP/INV=%d, CNV/multi-copy=%d, union=%d\n",
    length(sv_nv), length(cnv_nv), length(ex)
  ))
  ex
}

POINT_SIZE_MIN <- 0.9
POINT_SIZE_MAX <- 4.2

load_sv_underlier_rep_counts <- function(under_gene_pairs, kept_subjects) {
  if (!file.exists(SEGDUP_HITS_TSV) || !file.exists(SEGDUP_CARRIERS_TSV)) {
    stop(
      "Missing SV hit/carrier tables. Rebuild:\n",
      "  Rscript scripts/fig3_ptv_vs_rna/build_fig3b_segdup_svtype_genes.R"
    )
  }
  hits <- read_tsv(SEGDUP_HITS_TSV, show_col_types = FALSE) |>
    select(gene_nv, variant_id)
  carriers <- read_tsv(SEGDUP_CARRIERS_TSV, show_col_types = FALSE) |>
    filter(subjectid %in% kept_subjects)
  if (!nrow(under_gene_pairs)) {
    return(tibble(gene_nv = character(), n_sv_in_underliers = integer()))
  }
  vc <- hits |>
    inner_join(carriers, by = "variant_id", relationship = "many-to-many")
  under_sv <- under_gene_pairs |>
    inner_join(vc, by = c("gene_nv", "subjectid")) |>
    distinct(gene_nv, variant_id)
  under_sv |>
    count(gene_nv, name = "n_sv_in_underliers")
}

load_cnv_gene_map <- function() {
  if (!file.exists(SEGDUP_SVTYPE_TSV)) {
    stop(
      "Missing SegDup SVTYPE gene set: ", SEGDUP_SVTYPE_TSV,
      "\nBuild with: Rscript scripts/fig3_ptv_vs_rna/build_fig3b_segdup_svtype_genes.R"
    )
  }
  cat("Loading SegDup + primary SVTYPE gene set:", SEGDUP_SVTYPE_TSV, "\n")
  sv <- read_tsv(SEGDUP_SVTYPE_TSV, show_col_types = FALSE) |>
    distinct(gene_nv, .keep_all = TRUE) |>
    mutate(svtype_class = as.character(svtype_class)) |>
    filter(svtype_class %in% c("DEL", "DUP", "INV"))
  cat("  SegDup-pass genes with gene-overlap SV:", nrow(sv), "\n")
  print(as.data.frame(sort(table(sv$svtype_class), decreasing = TRUE)))

  sv |>
    transmute(
      gene_nv,
      cnv_symbol = as.character(symbol),
      cnv_reason = paste0(
        "segdup;",
        svtype_class, ":",
        as.character(variant_id)
      ),
      cnv_class = factor(svtype_class, levels = CNV_CLASS_LEVELS)
    )
}

# ---- NMD escape colouring (FIG3B_COLOR=nmd alternate) ---------------------
NMD_CLASS_LEVELS <- c("NMD-escaping", "NMD-triggering", "Unknown")
NMD_CLASS_COLORS <- c(
  "NMD-escaping"   = unname(OI["vermillion"]),  # #D55E00
  "NMD-triggering" = unname(OI["blue"]),        # #0072B2
  "Unknown"        = "grey55"
)
NMD_POINT_ALPHA <- 1
NMD_BY_VARIANT <- file.path(
  ROOT, "nmd_escape_vep", "fig3b_ptv_nmd_escape_by_variant.tsv"
)
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
cat("Loading expr_cpm_adj.rds (eligible gene universe) ...\n")
.expr_cpm_adj <- read_rds(EXPR_CPM_ADJ)
eligible_nv <- unique(str_remove(rownames(.expr_cpm_adj), "\\..*$"))
rm(.expr_cpm_adj); gc(verbose = FALSE)
cat("  eligible (expression-tested) genes:", length(eligible_nv), "\n")

TPM_MIN <- 1
tpm_med <- read_tsv(file.path(DATA_CORE, "gene_median_tpm_kept.tsv"), show_col_types = FALSE) |>
  distinct(gene_nv, .keep_all = TRUE)
tpm_ok_nv <- tpm_med |> filter(median_tpm > TPM_MIN) |> pull(gene_nv)
cat("  genes with median TPM >", TPM_MIN, ":", length(tpm_ok_nv), "\n")

ptv_gene <- ptv |>
  mutate(gene_nv = str_remove(gene, "\\..*$")) |>
  distinct(subjectid = sample, gene_nv) |>
  count(gene_nv, name = "n_ptv_carriers")

# Annotation maps for the active colour mode.
cnv_map <- NULL
nmd_map <- NULL
sv_size_map <- NULL
if (identical(COLOR_MODE, "cnv") && !UNIFORM_STYLE) {
  cnv_map <- load_cnv_gene_map()
  under_z3 <- read_rds(under_rds_for("z_-3")) |>
    filter(subjectid %in% kept_subjects) |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
    distinct(gene_nv, subjectid)
  sv_size_map <- load_sv_underlier_rep_counts(under_z3, kept_subjects)
  cat("  SV-in-underlier size map: ",
      nrow(sv_size_map), " genes with n>=1\n", sep = "")
} else if (identical(COLOR_MODE, "cnv") && UNIFORM_STYLE) {
  cat("Uniform-style plot: skipping SegDup SVTYPE / SV-in-underlier maps\n")
} else {
  nmd_map <- load_nmd_gene_map()
}

gene_table_for_mode <- function(m) {
  under <- read_rds(under_rds_for(m)) |>
    filter(subjectid %in% kept_subjects)
  under_gene <- under |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
    distinct(subjectid, gene_nv) |>
    count(gene_nv, name = "n_underlier_subjects")

  tab <- full_join(ptv_gene, under_gene, by = "gene_nv") |>
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

  if (identical(COLOR_MODE, "cnv") && !UNIFORM_STYLE) {
    tab <- tab |>
      left_join(cnv_map, by = "gene_nv") |>
      left_join(sv_size_map, by = "gene_nv") |>
      mutate(
        # Prefer LOEUF symbol; fall back to CNV-set / GENCODE symbol when LOEUF missing
        # (e.g. GATD3B, FP565260.1).
        symbol = coalesce(symbol, cnv_symbol),
        n_sv_in_underliers = replace_na(as.integer(n_sv_in_underliers), 0L),
        cnv_class = factor(
          ifelse(is.na(cnv_class), "Other genes", as.character(cnv_class)),
          levels = CNV_CLASS_LEVELS
        )
      )
  } else if (identical(COLOR_MODE, "cnv") && UNIFORM_STYLE) {
    tab
  } else {
    tab <- tab |>
      left_join(nmd_map, by = "gene_nv") |>
      mutate(
        nmd_class = factor(
          ifelse(is.na(nmd_class), "Unknown", as.character(nmd_class)),
          levels = NMD_CLASS_LEVELS
        )
      )
  }
  tab
}

all_df <- bind_rows(lapply(ORDER, function(m) {
  cat("Computing gene-level table for:", m, "\n")
  gene_table_for_mode(m)
})) |>
  mutate(mode = factor(mode, levels = ORDER, labels = LABELS[ORDER]))

if (identical(COLOR_MODE, "cnv") && !UNIFORM_STYLE) {
  all_df <- all_df |>
    mutate(
      n_sv_in_underliers = replace_na(as.integer(n_sv_in_underliers), 0L),
      n_sv_size = pmin(n_sv_in_underliers, 5L)
    ) |>
    arrange(
      cnv_class == "Other genes",
      cnv_class == "INV",
      cnv_class == "DUP",
      cnv_class == "DEL",
      cnv_class,
      n_sv_size
    )
} else if (identical(COLOR_MODE, "cnv") && UNIFORM_STYLE) {
  # uniform style: plotting order immaterial
} else {
  all_df <- all_df |>
    arrange(nmd_class == "Unknown", nmd_class == "NMD-triggering", nmd_class)
}

spearman_fit <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  r <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  tstat <- r * sqrt((n - 2) / (1 - r^2))
  logp  <- log(2) + stats::pt(-abs(tstat), df = n - 2, log.p = TRUE)
  tibble(rho = r, n_pairs = n, logp = logp, p_value = exp(logp))
}
fmt_logp <- function(logp) {
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
    .groups = "drop"
  ) |>
  mutate(p_label = fmt_logp(logp))

cat("\nPer-definition correlation (gene-level):\n")
print(summ |> select(mode, rho, p_label, pearson_all, n_ptv_genes, n_genes, n_under_only))

if (identical(COLOR_MODE, "cnv") && !UNIFORM_STYLE) {
  n_del <- sum(all_df$cnv_class == "DEL")
  n_dup <- sum(all_df$cnv_class == "DUP")
  n_inv <- sum(all_df$cnv_class == "INV")
  tl <- all_df$pct_under >= 0.5 & all_df$pct_ptv <= 0.25
  n_tl <- sum(tl)
  cat(sprintf(
    "\nDosage colouring: DEL=%d, DUP=%d, INV=%d / %d plotted; top-left: DEL=%d, DUP=%d, INV=%d / %d\n",
    n_del, n_dup, n_inv, nrow(all_df),
    sum(tl & all_df$cnv_class == "DEL"),
    sum(tl & all_df$cnv_class == "DUP"),
    sum(tl & all_df$cnv_class == "INV"),
    n_tl
  ))
  for (cl in c("DEL", "DUP", "INV")) {
    syms <- sort(unique(all_df$symbol[tl & all_df$cnv_class == cl]))
    if (length(syms)) {
      cat("  Top-left ", cl, ": ", paste(syms, collapse = ", "), "\n", sep = "")
    }
  }
}

if (OTHER_ONLY) {
  n_before <- nrow(all_df)
  all_df <- all_df |>
    filter(cnv_class == "Other genes") |>
    mutate(cnv_class = factor(as.character(cnv_class), levels = "Other genes"))
  cat(sprintf(
    "FIG3B_OTHER_ONLY: kept %d / %d genes (Other only)\n",
    nrow(all_df), n_before
  ))
  # Recompute correlation summary on the filtered set for the log.
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
  cat("\nPer-definition correlation (Other genes only):\n")
  print(summ |> select(mode, rho, p_label, pearson_all, n_ptv_genes, n_genes, n_under_only))
}


if (EXCLUDE_SV_CNV) {
  n_before <- nrow(all_df)
  n_ptv_before <- sum(all_df$n_ptv_carriers > 0, na.rm = TRUE)
  exclude_nv <- load_exclude_sv_cnv_genes()
  all_df <- all_df |> filter(!gene_nv %in% exclude_nv)
  n_removed <- n_before - nrow(all_df)
  cat(sprintf(
    "EXCLUDE_SV_CNV: removed %d genes; %d remain (%d with n_ptv>0 before -> %d after)\n",
    n_removed, nrow(all_df), n_ptv_before, sum(all_df$n_ptv_carriers > 0, na.rm = TRUE)
  ))
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
  cat("\nPer-definition correlation (SV/CNV genes excluded):\n")
  print(summ |> select(mode, rho, p_label, pearson_all, n_ptv_genes, n_genes, n_under_only))
}

# Log10 on both axes. Zeros cannot be logged — plot at half the smallest
# positive value on that axis (tick labeled "0%"). Correlations / top-left
# cuts still use raw pct_ptv / pct_under.
pos_ptv <- all_df$pct_ptv[is.finite(all_df$pct_ptv) & all_df$pct_ptv > 0]
pos_und <- all_df$pct_under[is.finite(all_df$pct_under) & all_df$pct_under > 0]
x_floor <- if (length(pos_ptv)) min(pos_ptv) / 2 else (50 / n_kept)
y_floor <- if (length(pos_und)) min(pos_und) / 2 else (50 / n_kept)
x_max <- max(all_df$pct_ptv, na.rm = TRUE)
y_max <- max(all_df$pct_under, na.rm = TRUE)
all_df <- all_df |>
  mutate(
    pct_ptv_plot = ifelse(pct_ptv > 0, pct_ptv, x_floor),
    pct_under_plot = ifelse(pct_under > 0, pct_under, y_floor)
  )
cat(sprintf(
  "X log10 floor for pct_ptv==0: %.5f%% (min positive=%.5f%%)\n",
  x_floor, if (length(pos_ptv)) min(pos_ptv) else NA_real_
))
cat(sprintf(
  "Y log10 floor for pct_under==0: %.5f%% (min positive=%.5f%%)\n",
  y_floor, if (length(pos_und)) min(pos_und) else NA_real_
))

# Log10 breaks on both axes: floor ("0%"), then decade-ish positives.
pct_breaks_cand <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
make_log_breaks <- function(floor_v, max_v) {
  cand <- pct_breaks_cand[
    pct_breaks_cand <= max_v * 1.05 &
      (log10(pct_breaks_cand) - log10(floor_v)) > 0.25
  ]
  br <- sort(unique(c(floor_v, cand)))
  if (length(br) < 2 || max_v > max(br) * 1.2) {
    br <- sort(unique(c(br, 10^ceiling(log10(max(max_v, floor_v * 10))))))
  }
  br[br > 0]
}
x_breaks <- make_log_breaks(x_floor, x_max)
y_breaks <- make_log_breaks(y_floor, y_max)
label_pct_log <- function(floor_v) {
  force(floor_v)
  function(x) {
    vapply(x, function(v) {
      if (!is.finite(v) || v <= 0) return("")
      if (abs(log10(v) - log10(floor_v)) < 1e-8) return("0%")
      acc <- if (v < 0.1) 0.01 else if (v < 1) 0.1 else 1
      label_percent(scale = 1, accuracy = acc)(v)
    }, character(1))
  }
}
label_pct_x <- label_pct_log(x_floor)
label_pct_y <- label_pct_log(y_floor)

if (UNIFORM_STYLE) {
  p <- ggplot(all_df, aes(x = pct_ptv_plot, y = pct_under_plot)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey55", linewidth = 0.35) +
    geom_point(
      color = NO_SV_POINT_COLOR, size = NO_SV_POINT_SIZE,
      alpha = 1, stroke = 0
    )
} else if (identical(COLOR_MODE, "cnv")) {
  p <- ggplot(all_df, aes(x = pct_ptv_plot, y = pct_under_plot)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey55", linewidth = 0.35) +
    geom_point(
      aes(color = cnv_class, size = n_sv_size),
      alpha = 1, stroke = 0
    ) +
    scale_color_manual(
      values = if (OTHER_ONLY) CNV_CLASS_COLORS["Other genes"] else CNV_CLASS_COLORS,
      name = NULL, drop = !OTHER_ONLY
    ) +
    scale_size_continuous(
      range = c(POINT_SIZE_MIN, POINT_SIZE_MAX),
      breaks = c(0, 1, 2, 5),
      labels = c("0", "1", "2", "5+"),
      limits = c(0, 5),
      name = "SVs in RNA underliers\n(distinct PASS, gene overlap)"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2.6)),
      size = guide_legend(order = 2, override.aes = list(color = "grey35", alpha = 1))
    )
} else {
  p <- ggplot(all_df, aes(x = pct_ptv_plot, y = pct_under_plot)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey55", linewidth = 0.35) +
    geom_point(aes(color = nmd_class), alpha = NMD_POINT_ALPHA, size = 1.3, stroke = 0) +
    scale_color_manual(values = NMD_CLASS_COLORS, name = "NMD status", drop = FALSE) +
    guides(color = guide_legend(override.aes = list(alpha = NMD_POINT_ALPHA, size = 2.4)))
}

# On-plot Spearman (genes with ≥1 rare PTV carrier); top-right.
# With log–log axes, geom_abline(slope=1) is again the true y=x diagonal.
rho_plot <- summ$rho[1]
p <- p +
  annotate(
    "text",
    x = max(all_df$pct_ptv_plot, na.rm = TRUE),
    y = max(all_df$pct_under_plot, na.rm = TRUE),
    label = sprintf("Spearman ρ = %.2f", rho_plot),
    hjust = 1.05, vjust = 1.25,
    size = 3.0, family = PAPER_FONT, colour = "grey15"
  ) +
  scale_x_log10(
    breaks = x_breaks,
    labels = label_pct_x,
    limits = c(x_floor, NA),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  scale_y_log10(
    breaks = y_breaks,
    labels = label_pct_y,
    limits = c(y_floor, NA),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "% of subjects carrying a rare PTV (AF<0.5%) in gene",
    y = "% of subjects with RNA\nunderliers in gene"
  ) +
  theme_paper() +
  theme(
    legend.position = if (UNIFORM_STYLE) "none" else "right",
    plot.margin = margin(t = 4, r = 6, b = 2, l = 16, unit = "mm"),
    axis.title.y = element_text(
      angle = 90, lineheight = 0.92,
      margin = margin(r = 3, unit = "mm")
    )
  )

.paper_plots <- list(a = p)
out_name <- if (EXCLUDE_SV_CNV) {
  "ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4_excl_sv_cnv"
} else if (NO_SV_COLOR) {
  "ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4_no_sv_color"
} else if (OTHER_ONLY) {
  "ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4_other_only_test"
} else {
  "ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4"
}
out_stem <- file.path(OUT_DIR, out_name)
paper_ggsave(p, out_stem, width = NATURE_DOUBLE_IN, height = 4.6)
if (!is_paper_assemble()) cat("Wrote:", out_stem, "\n")
