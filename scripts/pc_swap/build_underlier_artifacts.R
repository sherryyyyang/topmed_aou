#!/usr/bin/env Rscript
# Derive the PC-dependent underlier artifacts for whichever expression-PC
# setting PC_SETTING names, writing them into that setting's DATA_UNDERLIER.
#
# Everything downstream of the residualized z-score matrix is rebuilt here:
#
#   expr_z_join.rds                            gene_id, subjectid, z, log2CPM_adj
#   z_score_underlier_{haplo,z_-2,z_-3,z_-4}_clean.rds
#   z_score_tally_{...}_clean.rds              per-gene counts, rank, prevalence
#   z_score_tally_subj_{...}_clean.rds         per-subject counts, rank
#
# The recipe is lifted from scripts/analysis/run_strict_zcutoff_figures.R, which
# built the published 55-PC tables out of expr_z_join.rds. Running this with
# PC_SETTING=55 therefore has to reproduce those tables exactly; that round-trip
# is the correctness gate for the 2005-PC replicate and is asserted below.
#
# Only the z-score arm of the underlier definition moves with PC_SETTING. The
# log2CPM arm comes from expr_cpm_adj.rds, which is residualized on 55
# expression PCs and is the same matrix in both runs; see the README written
# into submission_package_nov_2026 for what that does and does not affect.
#
# Usage:
#   PC_SETTING=2005 Rscript scripts/pc_swap/build_underlier_artifacts.R
#   # round-trip gate (never touches data/underliers):
#   PC_SETTING=55 UNDERLIER_OUT_DIR=.../data/underliers_pc55_validation \
#     Rscript scripts/pc_swap/build_underlier_artifacts.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(data.table); library(tibble)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
cm <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")

PUBLISHED_UNDERLIER <- file.path(ROOT, "data", "underliers")
LOGCPM_DROP <- 1
Z_CUTOFFS   <- c(-2, -3, -4)

OUT_DIR <- Sys.getenv("UNDERLIER_OUT_DIR", unset = DATA_UNDERLIER)
# Hard stop: the published 55-PC artifacts are immutable.
if (normalizePath(OUT_DIR, mustWork = FALSE) ==
    normalizePath(PUBLISHED_UNDERLIER, mustWork = FALSE)) {
  stop("refusing to write into the published 55-PC directory ", PUBLISHED_UNDERLIER,
       "; set UNDERLIER_OUT_DIR to a scratch path")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

msg("PC_SETTING  = ", PC_SETTING)
msg("z long-form = ", Z_LONGFORM)
msg("output dir  = ", OUT_DIR)
stopifnot(file.exists(Z_LONGFORM))

# --------------------------------------------------------------------- axes
# The long-form table is one contiguous block of all subjects per gene. Rather
# than parse 45M key rows, read the subject block once and every ngene-th gene
# label, then check the recovered axes against expr_cpm_adj.
msg("loading expr_cpm_adj.rds ...")
cpm <- read_rds(EXPR_CPM_ADJ)
ngene <- nrow(cpm); nsub <- ncol(cpm)
msg("  grid: ", cm(ngene), " genes x ", cm(nsub), " subjects = ",
    cm(as.numeric(ngene) * nsub), " pairs")

AXES <- file.path(OUT_DIR, "longform_axes.rds")
axes <- if (file.exists(AXES)) {
  msg("reusing cached long-form axes"); read_rds(AXES)
} else {
  msg("recovering long-form gene/subject axes ...")
  a <- list(
    subjects = fread(cmd = sprintf("head -n %d %s | cut -f2 | tail -n +2", nsub + 1, Z_LONGFORM),
                     header = FALSE)[[1]],
    genes = fread(cmd = sprintf("awk -F'\\t' 'NR>1 && (NR-2)%%%d==0 {print $1}' %s", nsub, Z_LONGFORM),
                  header = FALSE)[[1]])
  write_rds(a, AXES); a
}
genes <- axes$genes; subjects <- axes$subjects
stopifnot(length(genes) == ngene, length(subjects) == nsub,
          !anyDuplicated(genes), !anyDuplicated(subjects),
          setequal(genes, rownames(cpm)), identical(subjects, colnames(cpm)))
cpm <- cpm[genes, , drop = FALSE]
msg("  axes verified; expr_cpm_adj aligned to them")

# ------------------------------------------------------------- expr_z_join
# Same schema as data/expression/expr_z_join.rds: one row per (gene, subject),
# gene-major with subject varying fastest.
EZJ <- file.path(OUT_DIR, "expr_z_join.rds")
if (file.exists(EZJ)) {
  msg("reusing cached expr_z_join: ", EZJ)
  expr_z_join <- read_rds(EZJ)
} else {
  msg("reading z column from the long-form table ...")
  t0 <- Sys.time()
  z <- fread(Z_LONGFORM, select = 3, showProgress = FALSE)[[1]]
  stopifnot(length(z) == as.numeric(ngene) * nsub, !anyNA(z))
  msg(sprintf("  %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  msg("assembling expr_z_join ...")
  expr_z_join <- tibble(
    gene_id           = rep(genes, each = nsub),
    subjectid         = rep(subjects, times = ngene),
    expression_zscore = z,
    # cpm is genes x subjects; t() makes subject vary fastest, matching the
    # long-form row order.
    expression_logcpm = as.vector(t(cpm)))
  rm(z); gc(verbose = FALSE)
  msg("  rows = ", cm(nrow(expr_z_join)))
  write_rds(expr_z_join, EZJ)          # uncompressed: this is a scratch cache
  msg("  wrote ", EZJ)
}

kept_subjects <- read_rds(KEPT_SUBJECTS)
N <- length(kept_subjects)
msg("kept subjects (extreme-subject filter): ", cm(N), " / ", cm(nsub))

# Per-gene mean log2CPM over every subject in the grid, exactly as
# run_strict_zcutoff_figures.R::derive_strict_underliers computed it.
gene_means <- expr_z_join |>
  group_by(gene_id) |>
  summarise(mean_logcpm = mean(expression_logcpm, na.rm = TRUE), .groups = "drop")
mean_logcpm <- gene_means$mean_logcpm[match(expr_z_join$gene_id, gene_means$gene_id)]
pass_logcpm <- expr_z_join$expression_logcpm < (mean_logcpm - LOGCPM_DROP)
rm(mean_logcpm); gc(verbose = FALSE)
msg("log2CPM-eligible pairs (< gene mean - ", LOGCPM_DROP, "): ", cm(sum(pass_logcpm)))

is_kept <- expr_z_join$subjectid %in% kept_subjects

# ------------------------------------------------------------------ derive
# Compare against the published 55-PC artifact when one exists. Under
# PC_SETTING=55 the comparison is an assertion, not a report.
check_against_published <- function(under, label) {
  ref <- file.path(PUBLISHED_UNDERLIER, sprintf("z_score_underlier_%s_clean.rds", label))
  if (!file.exists(ref)) { msg("    (no published counterpart for ", label, ")"); return(invisible()) }
  saved <- read_rds(ref)
  got  <- paste(under$gene_id, under$subjectid)
  want <- paste(saved$gene_id, saved$subjectid)
  msg("    vs published: ", cm(length(got)), " here, ", cm(length(want)),
      " saved, ", cm(length(intersect(got, want))), " shared")
  if (PC_IS_DEFAULT) {
    stopifnot(setequal(got, want))
    m <- match(want, got)
    stopifnot(max(abs(under$expression_zscore[m] - saved$expression_zscore)) < 1e-9,
              max(abs(under$expression_logcpm[m] - saved$expression_logcpm)) < 1e-9)
    msg("    ROUND-TRIP OK: reproduces the published ", label, " table exactly")
  }
  invisible()
}

summary_rows <- list()

emit <- function(label, keep) {
  sel <- which(keep & is_kept)
  under <- tibble(
    gene_id           = expr_z_join$gene_id[sel],
    subjectid         = expr_z_join$subjectid[sel],
    expression_zscore = expr_z_join$expression_zscore[sel],
    expression_logcpm = expr_z_join$expression_logcpm[sel])
  msg("  ", label, ": ", cm(nrow(under)), " underlier pairs after the extreme-subject filter")
  check_against_published(under, label)

  tally_gene <- under |>
    group_by(gene_id) |> tally(name = "n") |> arrange(desc(n)) |>
    mutate(pct_underlier = n / N, rank = row_number())
  tally_subj <- under |>
    group_by(subjectid) |> tally(name = "n") |> arrange(desc(n)) |>
    mutate(rank = row_number())

  saveRDS(under,      file.path(OUT_DIR, sprintf("z_score_underlier_%s_clean.rds", label)))
  saveRDS(tally_gene, file.path(OUT_DIR, sprintf("z_score_tally_%s_clean.rds", label)))
  saveRDS(tally_subj, file.path(OUT_DIR, sprintf("z_score_tally_subj_%s_clean.rds", label)))

  summary_rows[[label]] <<- tibble(
    pc_setting = PC_SETTING, definition = label,
    n_underlier_pairs = nrow(under),
    n_genes_with_underlier = nrow(tally_gene),
    n_subjects_with_underlier = nrow(tally_subj),
    median_per_subject = if (nrow(tally_subj)) median(tally_subj$n) else NA_real_,
    max_per_gene = if (nrow(tally_gene)) max(tally_gene$n) else NA_integer_)
  invisible()
}

msg("deriving underlier sets ...")
emit("haplo", pass_logcpm)                       # log2CPM drop alone, no z threshold
for (zc in Z_CUTOFFS) emit(sprintf("z_%s", zc), pass_logcpm & (expr_z_join$expression_zscore < zc))

# The published tree also carries an unsuffixed "strict" alias for z < -4.
for (stem in c("z_score_underlier", "z_score_tally", "z_score_tally_subj")) {
  file.copy(file.path(OUT_DIR, sprintf("%s_z_-4_clean.rds", stem)),
            file.path(OUT_DIR, sprintf("%s_clean.rds", stem)), overwrite = TRUE)
}

summ <- bind_rows(summary_rows)
write_tsv(summ, file.path(OUT_DIR, "underlier_counts_summary.tsv"))
cat("\n"); print(as.data.frame(summ)); cat("\n")
msg("wrote ", OUT_DIR)
msg("DONE")
