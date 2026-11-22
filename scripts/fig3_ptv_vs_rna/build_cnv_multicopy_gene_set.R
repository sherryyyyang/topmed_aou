#!/usr/bin/env Rscript
# Build the Fig3b CNV / multi-copy gene highlight set.
#
# A gene is flagged if EITHER:
#   (1) >=50% of its GENCODE v30 gene body overlaps a UCSC hg38 genomicSuperDups
#       interval with fracMatch >= 0.95, OR
#   (2) its symbol matches a curated multi-copy / polymorphic-CNV family prefix
#       (SIGLEC, DEFA/DEFB, TBC1D3, GOLGA8/6, AGAP4-12, NOMO, GATD3, NBPF, …).
#
# See data/references/SOURCE_fig3b_cnv_multicopy.txt for the full rule and
# provenance of the SegDup track.
#
# Outputs (under data/references/):
#   genomicSuperDups_hg38.txt.gz / .bed   (downloaded if missing)
#   gencode.v30.genes.bed                 (from annotation/gencode.v30.annotation.gff3)
#   gencode.v30_genes_segdup_overlap.tsv
#   fig3b_cnv_multicopy_genes.tsv
#   SOURCE_fig3b_cnv_multicopy.txt

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))

ANN <- file.path(ROOT, "data", "references")
dir.create(ANN, recursive = TRUE, showWarnings = FALSE)
SHEYANG_ANN <- file.path(SHEYANG, "annotation")
GFF <- file.path(SHEYANG_ANN, "gencode.v30.annotation.gff3")
SEGDUP_GZ <- file.path(ANN, "genomicSuperDups_hg38.txt.gz")
SEGDUP_URL <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"
GENE_BED <- file.path(ANN, "gencode.v30.genes.bed")
SEGDUP_BED <- file.path(ANN, "genomicSuperDups_hg38.bed")
OVERLAP_TSV <- file.path(ANN, "gencode.v30_genes_segdup_overlap.tsv")
OUT_TSV <- file.path(ANN, "fig3b_cnv_multicopy_genes.tsv")
SOURCE_TXT <- file.path(ANN, "SOURCE_fig3b_cnv_multicopy.txt")

CURATED_PREFIX_RE <- paste0(
  "^(",
  paste(
    c(
      "SIGLEC", "DEFA", "DEFB", "TBC1D3", "GOLGA8", "GOLGA6",
      "AGAP([4-9]|1[0-2])$", "NOMO", "GATD3", "NBPF", "WASH[0-9]",
      "SPDYE", "FAM90A", "RGPD", "GTF2I$", "GTF2IRD1?$",
      "SMN[12]$", "FCGR3", "CCL3L", "CCL4L", "UGT2B17$",
      "GSTM1$", "GSTT1$", "AMY1", "OPN1[LM]W", "RHD$", "RHCE$",
      "DUX4", "CT47A", "TSPY([0-9]+|P[0-9]*)?$", "PRAMEF",
      "SSX[0-9]+P?$", "NPIP", "FRG[12]", "ROCK1P", "SERF1"
    ),
    collapse = "|"
  ),
  ")"
)

if (!file.exists(SEGDUP_GZ)) {
  cat("Downloading UCSC genomicSuperDups ...\n")
  download.file(SEGDUP_URL, SEGDUP_GZ, mode = "wb", quiet = FALSE)
}
if (!file.exists(GENE_BED)) {
  stopifnot(file.exists(GFF))
  cat("Building gene-body BED from", GFF, "\n")
  cmd <- sprintf(
    "awk -F'\\t' '$3==\"gene\"{match($9,/gene_id=([^;]+)/,a); match($9,/gene_name=([^;]+)/,b); match($9,/gene_type=([^;]+)/,c); gid=a[1]; gsub(/\\..*/,\"\",gid); print $1\"\\t\"($4-1)\"\\t\"$5\"\\t\"gid\"\\t\"b[1]\"\\t\"c[1]}' %s > %s",
    shQuote(GFF), shQuote(GENE_BED)
  )
  stopifnot(system(cmd) == 0)
}
if (!file.exists(SEGDUP_BED)) {
  cat("Writing SegDup BED (both pair loci) ...\n")
  cmd <- sprintf(
    "zcat %s | awk -F'\\t' 'BEGIN{OFS=\"\\t\"} {print $2,$3,$4,$5,$27; print $8,$9,$10,$5,$27}' | awk -F'\\t' 'BEGIN{OFS=\"\\t\"} $1!=\"\" && $3>$2' | sort -k1,1 -k2,2n | uniq > %s",
    shQuote(SEGDUP_GZ), shQuote(SEGDUP_BED)
  )
  stopifnot(system(cmd) == 0)
}

if (!file.exists(OVERLAP_TSV) || file.info(OVERLAP_TSV)$mtime < file.info(SEGDUP_BED)$mtime) {
  cat("Intersecting gene bodies with SegDups (bedtools) ...\n")
  tmp <- tempfile(fileext = ".tsv")
  cmd <- sprintf(
    "bedtools intersect -a %s -b %s -wa -wb | awk -F'\\t' 'BEGIN{OFS=\"\\t\"; print \"gene_nv\",\"symbol\",\"gene_type\",\"gene_chr\",\"gene_start\",\"gene_end\",\"segdup_chr\",\"segdup_start\",\"segdup_end\",\"segdup_name\",\"fracMatch\",\"overlap_bp\"} {gs=$2; ge=$3; ss=$8; se=$9; o_start=(gs>ss?gs:ss); o_end=(ge<se?ge:se); ov=o_end-o_start; if(ov<0) ov=0; print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,ov}' > %s",
    shQuote(GENE_BED), shQuote(SEGDUP_BED), shQuote(tmp)
  )
  stopifnot(system(cmd) == 0)
  ov <- read_tsv(tmp, show_col_types = FALSE) |>
    mutate(gene_len = gene_end - gene_start) |>
    group_by(gene_nv, symbol, gene_type, gene_chr, gene_start, gene_end, gene_len) |>
    summarise(
      n_segdup_intervals = n_distinct(segdup_name),
      max_fracMatch = max(fracMatch, na.rm = TRUE),
      max_overlap_bp = max(overlap_bp, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(max_overlap_frac = max_overlap_bp / gene_len)
  write_tsv(ov, OVERLAP_TSV)
  unlink(tmp)
  cat("  genes overlapping any SegDup:", nrow(ov), "\n")
}

genes <- read_tsv(
  GENE_BED,
  col_names = c("chr", "start", "end", "gene_nv", "symbol", "gene_type"),
  show_col_types = FALSE
)
ov <- read_tsv(OVERLAP_TSV, show_col_types = FALSE) |>
  mutate(
    segdup_strict = !is.na(max_overlap_frac) &
      max_overlap_frac >= 0.5 &
      max_fracMatch >= 0.95
  ) |>
  select(gene_nv, max_overlap_frac, max_fracMatch, segdup_strict)

out <- genes |>
  left_join(ov, by = "gene_nv") |>
  mutate(
    curated = str_detect(coalesce(symbol, ""), CURATED_PREFIX_RE),
    segdup_strict = coalesce(segdup_strict, FALSE),
    reason = case_when(
      curated & segdup_strict ~ "curated_family;segdup_ge50pct_fracMatch_ge0.95",
      curated ~ "curated_family",
      segdup_strict ~ "segdup_ge50pct_fracMatch_ge0.95",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(reason)) |>
  transmute(
    gene_nv, symbol, gene_type, chr, start, end,
    cnv_class = "CNV / multi-copy",
    reason,
    max_overlap_frac, max_fracMatch
  ) |>
  arrange(symbol, gene_nv)

write_tsv(out, OUT_TSV)
cat("Wrote", OUT_TSV, "n=", nrow(out),
    " protein_coding=", sum(out$gene_type == "protein_coding"), "\n")
print(as.data.frame(sort(table(out$reason), decreasing = TRUE)))

if (!file.exists(SOURCE_TXT)) {
  writeLines(c(
    "Fig3b CNV / multi-copy gene highlight set",
    "See build_cnv_multicopy_gene_set.R header and comments for the rule.",
    paste("SegDup URL:", SEGDUP_URL),
    paste("Built:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ), SOURCE_TXT)
}
