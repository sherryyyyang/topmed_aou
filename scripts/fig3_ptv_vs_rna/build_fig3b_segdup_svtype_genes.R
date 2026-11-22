#!/usr/bin/env Rscript
# Build Fig3b dosage colouring: SegDup-passing genes with a primary gene-overlapping
# PASS muCNV SVTYPE (DEL / DUP / INV). Any AF; primary variant = largest |SVLEN|.
#
# Outputs (data/references/):
#   fig3b_segdup_svtype_genes.tsv
#   SOURCE_fig3b_segdup_svtype.txt

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))

ANN <- file.path(ROOT, "data", "references")
CNV_GENES_TSV <- file.path(ANN, "fig3b_cnv_multicopy_genes.tsv")
OUT_TSV <- file.path(ANN, "fig3b_segdup_svtype_genes.tsv")
HITS_TSV <- file.path(ANN, "fig3b_segdup_sv_gene_hits.tsv")
CARRIERS_TSV <- file.path(ANN, "fig3b_segdup_sv_carriers.tsv")
SOURCE_TXT <- file.path(ANN, "SOURCE_fig3b_segdup_svtype.txt")
VDIR <- Sys.getenv("SV_VCF_DIR", unset = "")
if (!nzchar(VDIR)) VDIR <- file.path(ROOT, "data", "vcfs-pass-variants-only")

path_bin <- function(cmd) {
  p <- Sys.which(cmd)
  if (!nzchar(p)) stop("Required executable not found on PATH: ", cmd)
  p
}
BCF <- path_bin("bcftools")
BT <- path_bin("bedtools")

stopifnot(file.exists(CNV_GENES_TSV))
segdup_genes <- read_tsv(CNV_GENES_TSV, show_col_types = FALSE) |>
  filter(str_detect(coalesce(reason, ""), "segdup_ge50pct_fracMatch_ge0.95")) |>
  transmute(
    gene_nv, symbol, gene_type, chr, start, end,
    segdup_reason = reason,
    max_overlap_frac, max_fracMatch
  )
cat("SegDup-pass genes:", nrow(segdup_genes), "\n")
if (!nrow(segdup_genes)) stop("No SegDup-pass genes in ", CNV_GENES_TSV)

tmpdir <- tempfile("fig3b_segdup_sv_")
dir.create(tmpdir)

hits_list <- list()
for (ch in sort(unique(segdup_genes$chr))) {
  vcf <- file.path(VDIR, paste0(ch, ".vcf.gz"))
  if (!file.exists(vcf)) {
    warning("Missing VCF for ", ch, " — skipping")
    next
  }
  sv_bed <- file.path(tmpdir, paste0(ch, "_sv.bed"))
  cmd_q <- sprintf(
    "%s query -f '%%CHROM\\t%%POS\\t%%INFO/END\\t%%ID\\t%%INFO/SVTYPE\\t%%INFO/SVLEN\\t%%INFO/AF\\n' %s",
    BCF, shQuote(vcf)
  )
  cmd_awk <- paste(
    "awk -F'\\t' 'BEGIN{OFS=\"\\t\"}",
    "$5==\"DEL\"||$5==\"DUP\"||$5==\"INV\" {",
    "end=$3; if(end==\".\"||end==\"\") end=$2;",
    "n=$6+0; if(n==0) n=end-$2; if(n<0) n=-n;",
    "print $1,$2-1,end,$4,$5,n,$7",
    "}'"
  )
  full <- paste(cmd_q, "|", cmd_awk, ">", shQuote(sv_bed))
  if (system(full) != 0) stop("bcftools query failed for ", ch)

  gsub <- segdup_genes |> filter(chr == ch)
  if (!nrow(gsub)) next
  gchr <- file.path(tmpdir, paste0(ch, "_genes.bed"))
  write_tsv(
    gsub |> mutate(start0 = pmax(0L, as.integer(start) - 1L)) |>
      select(chr, start0, end, gene_nv, symbol),
    gchr,
    col_names = FALSE
  )
  out_inter <- file.path(tmpdir, paste0(ch, "_inter.tsv"))
  cmd_i <- sprintf(
    "%s intersect -a %s -b %s -wa -wb > %s",
    BT, shQuote(gchr), shQuote(sv_bed), shQuote(out_inter)
  )
  if (system(cmd_i) != 0) stop("bedtools intersect failed for ", ch)
  if (!file.exists(out_inter) || !file.info(out_inter)$size) next

  h <- read_tsv(
    out_inter,
    col_names = c(
      "chr", "g_start", "g_end", "gene_nv", "symbol",
      "sv_chr", "sv_start", "sv_end", "sv_id", "svtype", "svlen", "af"
    ),
    show_col_types = FALSE
  )
  if (nrow(h)) hits_list[[ch]] <- h
}

unlink(tmpdir, recursive = TRUE)

if (!length(hits_list)) {
  stop("No gene-overlapping PASS DEL/DUP/INV hits for SegDup-pass genes")
}

hits <- bind_rows(hits_list) |>
  mutate(
    svlen = abs(as.numeric(svlen)),
    af = as.numeric(af),
    svtype = as.character(svtype)
  ) |>
  filter(is.finite(svlen))

gene_hits <- hits |>
  transmute(
    gene_nv, symbol, chr, variant_id = sv_id, svtype, svlen, af,
    sv_start = as.integer(sv_start) + 1L,
    sv_end = as.integer(sv_end)
  ) |>
  distinct()
write_tsv(gene_hits, HITS_TSV)
cat("Wrote", HITS_TSV, "n=", nrow(gene_hits), "\n")

# Alt carriers per variant (kept-subject sample names as in VCF columns).
kept <- readRDS(KEPT_SUBJECTS)
carrier_list <- list()
tmpdir2 <- tempfile("fig3b_sv_car_")
dir.create(tmpdir2)
for (ch in sort(unique(gene_hits$chr))) {
  vcf <- file.path(VDIR, paste0(ch, ".vcf.gz"))
  if (!file.exists(vcf)) next
  ids <- gene_hits |> filter(chr == ch) |> pull(variant_id) |> unique()
  idf <- file.path(tmpdir2, paste0(ch, "_ids.txt"))
  writeLines(ids, idf)
  outc <- file.path(tmpdir2, paste0(ch, "_car.tsv"))
  cmd <- paste(
    sprintf(
      "%s query -i 'GT=\"alt\"' -f '[%%ID\\t%%SAMPLE\\n]' %s",
      BCF, shQuote(vcf)
    ),
    "| awk -F'\\t' 'NR==FNR{v[$1]=1; next} $1 in v {print}'",
    shQuote(idf), "- >", shQuote(outc)
  )
  if (system(cmd) != 0) stop("bcftools carrier query failed for ", ch)
  if (file.exists(outc) && file.info(outc)$size > 0) {
    carrier_list[[ch]] <- read_tsv(
      outc, col_names = c("variant_id", "subjectid"), show_col_types = FALSE
    )
  }
}
unlink(tmpdir2, recursive = TRUE)
carriers <- bind_rows(carrier_list) |>
  filter(subjectid %in% kept) |>
  distinct()
write_tsv(carriers, CARRIERS_TSV)
cat("Wrote", CARRIERS_TSV, "n=", nrow(carriers),
    " variants=", n_distinct(carriers$variant_id), "\n")

primary <- hits |>
  arrange(gene_nv, desc(svlen), desc(af)) |>
  group_by(gene_nv) |>
  slice(1) |>
  ungroup() |>
  left_join(
    segdup_genes |> select(gene_nv, gene_type, segdup_reason, max_overlap_frac, max_fracMatch),
    by = "gene_nv"
  ) |>
  transmute(
    gene_nv, symbol, gene_type, chr,
    svtype_class = svtype,
    variant_id = sv_id,
    svlen, af,
    sv_start = sv_start + 1L,
    sv_end,
    segdup_reason, max_overlap_frac, max_fracMatch,
    note = "SegDup-pass; gene-overlap PASS muCNV; primary by |SVLEN|; any AF"
  ) |>
  arrange(symbol, gene_nv)

write_tsv(primary, OUT_TSV)
cat("Wrote", OUT_TSV, "n=", nrow(primary), "\n")
print(as.data.frame(sort(table(primary$svtype_class), decreasing = TRUE)))

writeLines(c(
  "Fig3b SegDup + primary SVTYPE gene highlight set",
  "===============================================",
  "",
  "A gene is coloured DEL / DUP / INV on Fig3b only if:",
  "  (1) reason contains segdup_ge50pct_fracMatch_ge0.95 (see build_cnv_multicopy_gene_set.R), AND",
  "  (2) it has at least one PASS muCNV with SVTYPE in {DEL,DUP,INV} overlapping the gene body,",
  "      with no AF filter.",
  "",
  "Primary overlapping variant per gene: largest |SVLEN| among gene-overlapping hits.",
  "Genes that are curated-only or SegDup-pass without such a VCF hit -> Other genes.",
  "",
  "Also: fig3b_segdup_sv_gene_hits.tsv (all gene-overlap hits),",
  "      fig3b_segdup_sv_carriers.tsv (variant x kept-subject alt carriers).",
  "Fig3b point size = distinct gene-overlap PASS SVs with >=1 carrier among",
  "RNA underliers for that gene (strict z<-3 panel; built at plot time).",
  "",
  paste("VCF directory:", VDIR),
  paste("Built:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
), SOURCE_TXT)