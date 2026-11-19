#!/usr/bin/env Rscript
# Re-run the published HVG -> Freeman-Tukey -> connectivity -> re-HVG chain
# WITHOUT the GENCODE protein-coding filter. Stop before PCA / z-scores.
#
# Input counts: gene_counts.csv, same 5,560 NWD samples as expr_cpm_adj
# (the sample set the original HVG was run on). All genes kept.
#
# HVG = scran Poisson HVGs intersect CV2 HVGs with ratio > 1
# (same filter_expression() as outlier_processing.R / run_underlier_pipeline.R).
#
# Writes under pc2005_work/fig3_ptv_vs_rna/allgenes_hvg/ and a composition figure
# to submission_package_nov_2026/. Does not touch published expression matrices.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(scran)
  library(scuttle)
  library(corral)
  library(WGCNA)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
Sys.setenv(PC_SETTING = "2005")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

OUT <- file.path(PC_WORK, "fig3_ptv_vs_rna", "allgenes_hvg")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
COUNTS_CSV <- file.path(SHEYANG, "counts", "gene_counts.csv")
GFF <- file.path(SHEYANG, "annotation", "gencode.v30.annotation.gff3")
CONN_Z <- -3
allowWGCNAThreads()
options(warn = 1)
LNC_G30 <- c(
  "lincRNA", "antisense", "sense_intronic", "sense_overlapping",
  "processed_transcript", "3prime_overlapping_ncRNA", "macro_lncRNA",
  "bidirectional_promoter_lncRNA", "non_coding"
)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}
cache_rds <- function(path, expr) {
  if (file.exists(path)) {
    msg("cache load ", basename(path))
    return(read_rds(path))
  }
  val <- force(expr)
  write_rds(val, path)
  msg("wrote ", path)
  val
}

filter_expression <- function(counts) {
  gv_pois <- modelGeneVarByPoisson(counts)
  top_pois <- getTopHVGs(gv_pois)
  gv_cv2 <- modelGeneCV2(counts)
  top_cv2 <- getTopHVGs(gv_cv2, var.field = "ratio", var.threshold = 1.0)
  keep <- intersect(top_pois, top_cv2)
  msg("  HVG Poisson=", length(top_pois), " CV2=", length(top_cv2),
      " intersect=", length(keep), " / ", nrow(counts))
  counts[keep, , drop = FALSE]
}

g30_class_of <- function(types) {
  case_when(
    types == "protein_coding" ~ "Protein-coding",
    types %in% LNC_G30 ~ "lncRNA",
    grepl("pseudo", types) ~ "Pseudogene",
    is.na(types) ~ "Unmapped",
    TRUE ~ "Other"
  )
}

PAL <- c(
  "Protein-coding" = "grey70",
  "lncRNA"         = unname(OI["purple"]),
  "Pseudogene"     = unname(OI["green"]),
  "Other"          = unname(OI["sky"]),
  "Unmapped"       = unname(OI["vermillion"])
)
CLASS_LVLS <- names(PAL)

# --------------------------------------------------------------------------
msg("loading published sample/gene axes and GENCODE types ...")
cpm_pub <- read_rds(EXPR_CPM_ADJ)
samples <- colnames(cpm_pub)
pub_genes <- rownames(cpm_pub)
pub_nv <- unique(str_remove(pub_genes, "\\..*$"))
rm(cpm_pub)
gc(verbose = FALSE)
msg("  published HVG genes=", length(pub_genes), " samples=", length(samples))

g30_path <- file.path(OUT, "gencode_v30_gene_types.tsv")
g30 <- cache_rds(sub("\\.tsv$", ".rds", g30_path), {
  if (file.exists(g30_path)) {
    read_tsv(g30_path, show_col_types = FALSE)
  } else {
    cmd <- sprintf(
      "awk -F'\\t' '$3==\"gene\"{id=\"\";ty=\"\";n=split($9,a,\";\");for(i=1;i<=n;i++){if(a[i]~/^gene_id=/){split(a[i],b,\"=\");id=b[2]}if(a[i]~/^gene_type=/){split(a[i],b,\"=\");ty=b[2]}}gsub(/\\..*$/,\"\",id);print id\"\\t\"ty}' %s",
      shQuote(GFF)
    )
    d <- fread(cmd = cmd, header = FALSE, col.names = c("gene_nv", "g30_type"))
    d <- as_tibble(d) |> distinct(gene_nv, .keep_all = TRUE)
    write_tsv(d, g30_path)
    d
  }
})

meta <- read_tsv(METADATA, show_col_types = FALSE)

# --------------------------------------------------------------------------
counts <- cache_rds(file.path(OUT, "counts_allgenes_5560samples.rds"), {
  msg("reading gene_counts.csv ...")
  dt <- fread(COUNTS_CSV, sep = ",", header = TRUE, data.table = FALSE)
  gene_ids <- dt[[1]]
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- gene_ids
  rm(dt)
  gc(verbose = FALSE)
  nwd <- meta$WGS[match(colnames(mat), meta$RNA)]
  keep_col <- !is.na(nwd) & nwd %in% samples & !duplicated(nwd)
  mat <- mat[, keep_col, drop = FALSE]
  colnames(mat) <- nwd[keep_col]
  stopifnot(setequal(colnames(mat), samples))
  mat <- mat[, samples, drop = FALSE]
  msg("  counts: ", nrow(mat), " genes x ", ncol(mat), " samples")
  mat
})

annotate_genes <- function(ids, step) {
  tibble(gene_id = ids, gene_nv = str_remove(ids, "\\..*$"), step = step) |>
    left_join(g30, by = "gene_nv") |>
    mutate(
      g30_class = factor(g30_class_of(g30_type), levels = CLASS_LVLS),
      in_published_hvg = gene_nv %in% pub_nv
    )
}

input_ann <- annotate_genes(rownames(counts), "Input (all genes)")

hvg1 <- cache_rds(file.path(OUT, "hvg_pass1_counts.rds"), {
  msg("HVG pass 1 (all samples, all genes) ...")
  filter_expression(counts)
})
hvg1_ann <- annotate_genes(rownames(hvg1), "HVG pass 1")

norm1 <- cache_rds(file.path(OUT, "normalized_pass1.rds"), {
  msg("Freeman-Tukey pass 1 ...")
  n <- corral_preproc(hvg1, rtype = "freemantukey")
  colnames(n) <- colnames(hvg1)
  n
})
keep_samples <- cache_rds(file.path(OUT, "qc_kept_samples.rds"), {
  msg("connectivity QC ...")
  adj <- 0.5 + 0.5 * cor(norm1)
  conn <- as.numeric(fundamentalNetworkConcepts(adj)$Connectivity)
  sid <- colnames(norm1)
  z <- as.numeric((conn - mean(conn)) / sd(conn))
  keep <- z > CONN_Z
  msg(sprintf("  connectivity: removed %d of %d samples (Z < %s)",
              sum(!keep), length(keep), CONN_Z))
  tibble(subjectid = sid, connectivity = conn, z = z, keep = keep)
})
kept_ids <- keep_samples$subjectid[keep_samples$keep]
msg("  samples kept: ", length(kept_ids), " / ", nrow(keep_samples),
    "  overlap published 5560: ", sum(kept_ids %in% samples))

hvg2 <- cache_rds(file.path(OUT, "hvg_pass2_counts.rds"), {
  msg("HVG pass 2 + Freeman-Tukey on QC samples ...")
  qc_counts <- counts[, kept_ids, drop = FALSE]
  h <- filter_expression(qc_counts)
  norm2 <- corral_preproc(h, rtype = "freemantukey")
  colnames(norm2) <- colnames(h)
  write_rds(norm2, file.path(OUT, "qc_normalized.rds"))
  h
})
hvg2_ann <- annotate_genes(rownames(hvg2), "HVG pass 2 (final)")

pub_ann <- annotate_genes(pub_genes, "Published HVG (coding only)")

all_ann <- bind_rows(input_ann, hvg1_ann, hvg2_ann, pub_ann)
write_tsv(all_ann, file.path(OUT, "genes_by_step.tsv"))

comp <- all_ann |>
  dplyr::count(step, g30_class, name = "n") |>
  group_by(step) |>
  mutate(n_tot = sum(n), pct = 100 * n / n_tot) |>
  ungroup() |>
  mutate(
    step = factor(
      step,
      levels = c(
        "Input (all genes)", "HVG pass 1", "HVG pass 2 (final)",
        "Published HVG (coding only)"
      )
    ),
    g30_class = factor(g30_class, levels = CLASS_LVLS)
  )
write_tsv(comp, file.path(OUT, "composition_by_step.tsv"))

input_n <- input_ann |> dplyr::count(g30_class, name = "n_input")
retain <- hvg2_ann |>
  dplyr::count(g30_class, name = "n_final") |>
  full_join(input_n, by = "g30_class") |>
  mutate(
    n_input = replace(n_input, is.na(n_input), 0L),
    n_final = replace(n_final, is.na(n_final), 0L),
    pct_kept = 100 * n_final / n_input,
    g30_class = factor(g30_class, levels = CLASS_LVLS)
  )
write_tsv(retain, file.path(OUT, "retention_by_class.tsv"))

final_nv <- unique(hvg2_ann$gene_nv)
recov <- tibble(
  published_hvg = length(pub_nv),
  new_final = length(final_nv),
  shared = length(intersect(pub_nv, final_nv)),
  published_only = length(setdiff(pub_nv, final_nv)),
  new_only = length(setdiff(final_nv, pub_nv))
)
write_tsv(recov, file.path(OUT, "published_hvg_recovery.tsv"))

# --------------------------------------------------------------------------
msg("plotting composition sanity check ...")
pa <- ggplot(comp, aes(x = step, y = n, fill = g30_class)) +
  geom_col(width = 0.72) +
  geom_text(
    data = comp |> group_by(step) |> summarise(n = sum(n), .groups = "drop"),
    aes(x = step, y = n, label = comma(n)),
    inherit.aes = FALSE, vjust = -0.25, size = 2.4, family = PAPER_FONT
  ) +
  scale_fill_manual(values = PAL, name = NULL, drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.14))) +
  labs(x = NULL, y = "Genes") +
  theme_paper() +
  theme(axis.text.x = element_text(size = 7.2, angle = 22, hjust = 1))

pb <- ggplot(retain, aes(x = g30_class, y = pct_kept, fill = g30_class)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%s / %s\n%.1f%%", comma(n_final), comma(n_input), pct_kept)),
    vjust = -0.12, size = 2.2, family = PAPER_FONT, lineheight = 0.95
  ) +
  scale_fill_manual(values = PAL, guide = "none", drop = FALSE) +
  scale_y_continuous(
    limits = c(0, max(120, max(retain$pct_kept, na.rm = TRUE) * 1.25)),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(v) paste0(v, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "% of input class kept in final HVG") +
  theme_paper() +
  theme(axis.text.x = element_text(size = 7.5, angle = 25, hjust = 1))

recov_lab <- sprintf(
  "Published coding HVG recovered: %s / %s (%.1f%%)\nNew final only: %s   published dropped: %s",
  comma(recov$shared), comma(recov$published_hvg),
  100 * recov$shared / recov$published_hvg,
  comma(recov$new_only), comma(recov$published_only)
)
pc <- ggplot(hvg2_ann, aes(x = "HVG pass 2 (final)", fill = g30_class)) +
  geom_bar(width = 0.55) +
  scale_fill_manual(values = PAL, guide = "none", drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = NULL, y = "Final HVG genes") +
  theme_paper() +
  annotate(
    "text", x = 1, y = Inf, label = recov_lab,
    vjust = 1.15, hjust = 0.5, size = 2.4, family = PAPER_FONT, lineheight = 0.95
  )

fig <- (pa / (pb | pc)) +
  plot_layout(heights = c(1.15, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = PAPER_FONT, face = "bold", size = 14, colour = "black"
    ),
    legend.position = "bottom",
    plot.margin = margin(6, 10, 8, 6)
  )

draft <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
src <- file.path(draft, "source_panels")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
stem <- "Fig_hvg_allgenes_composition"
paper_ggsave(fig, file.path(OUT, stem), width = NATURE_DOUBLE_IN, height = 6.2)
paper_ggsave(fig, file.path(draft, stem), width = NATURE_DOUBLE_IN, height = 6.2)
paper_ggsave(fig, file.path(src, stem), width = NATURE_DOUBLE_IN, height = 6.2)

cat("\nComposition by step:\n")
print(comp, n = 30)
cat("\nRetention of each input class into final HVG:\n")
print(retain, n = 20)
cat("\nPublished HVG recovery:\n")
print(recov)
cat("\nSample QC: kept", length(kept_ids), "of", nrow(keep_samples),
    "; overlap with published 5560 =", sum(kept_ids %in% samples), "\n")
msg("DONE  figure ", file.path(draft, paste0(stem, ".png")))
