#!/usr/bin/env Rscript
# Per-gene summed gnomAD AF_exomes for TOPMed-overlapping PTVs.
#
# Pipeline (Gencode gene body -> VCF subset -> overlap join -> TOPMed PTV filter):
#   1. GENCODE v30 GFF3: keep type=="gene" for the 624 shared rare-outlier genes,
#      write a BED.
#   2. bcftools view -R BED on the local gnomAD v4.1 joint-exomes VCF (INFO/AF_exomes).
#   3. Variant GRanges x gene GRanges, plyranges::join_overlap_left, so a variant
#      is assigned by genomic overlap, not by the VEP gene column.
#   4. Keep only sites that are also a PTV in ptv_all.tsv (TOPMed).
#   5. Collapse: sum(AF_exomes) per gene. gnomAD has no individual-level genotypes,
#      so this slightly undercounts if one person carries more than one PTV in the gene.
#
# VariantAnnotation is not installed in this R library; the VCF is read with
# bcftools query and coerced to GRanges, which is the same overlap object
# readVcf() would have given for this stripped (AF_exomes/AC/AN only) VCF.

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicRanges)
  library(plyranges)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG4D

HERE      <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
ANN       <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang"
GFF       <- file.path(ANN, "annotation/gencode.v30.annotation.gff3.gz")
VCF       <- file.path(ANN, "annotation/gnomad/gnomad.joint.v4.1.sites.exomes.vcf.bgz")
PTV_TSV   <- file.path(ANN, "ptv_all.tsv")
GENES_TSV <- file.path(HERE, "shared_rare_outlier_genes_4dot.tsv")
BED       <- file.path(HERE, "gencode_v30_shared_rare_outlier_genes.bed")
VCF_SUB   <- file.path(HERE, "gnomad_af_exomes_shared_rare_outlier_genes.vcf.gz")
OUT       <- file.path(HERE, file.path(FIG4D, "sum_af_exome_topmed_overlap.tsv"))
SITES_OUT <- file.path(HERE, "topmed_ptv_gnomad_af_exome_by_gene_overlap.tsv")

pairs <- read_tsv(GENES_TSV, show_col_types = FALSE)
gene_set <- unique(pairs$gene_nv)
cat(sprintf("target genes: %d\n", length(gene_set)))

# --- 1. Gencode gene features -> BED ------------------------------------------------
cat("loading GENCODE v30 GFF3 (gene features)...\n")
gff <- rtracklayer::import(GFF)
stopifnot("type" %in% names(mcols(gff)))
genes_all <- gff[gff$type == "gene"]
gid <- as.character(genes_all$gene_id)
gene_nv <- sub("\\..*$", "", gid)
keep <- gene_nv %in% gene_set
genes <- genes_all[keep]
genes$gene_nv <- gene_nv[keep]
genes$symbol  <- as.character(genes$gene_name)
genes <- genes[!duplicated(genes$gene_nv)]
cat(sprintf("Gencode gene bodies matching target set: %d / %d\n",
            length(genes), length(gene_set)))
miss <- setdiff(gene_set, genes$gene_nv)
if (length(miss)) cat("missing from GFF gene features:", paste(miss, collapse = ", "), "\n")

bed_tbl <- tibble(
  chrom = as.character(seqnames(genes)),
  start = start(genes) - 1L,   # BED is 0-based
  end   = end(genes),
  name  = genes$gene_nv,
  score = 0,
  strand = as.character(strand(genes))
) |> arrange(chrom, start, end)
write_tsv(bed_tbl, BED, col_names = FALSE)
cat(sprintf("wrote %s (%d intervals)\n", BED, nrow(bed_tbl)))

# --- 2. subset VCF to those gene bodies --------------------------------------------
if (file.exists(VCF_SUB) && file.exists(paste0(VCF_SUB, ".tbi"))) {
  cat("reusing existing subset VCF\n")
} else {
  cat("bcftools view -R BED ...\n")
  st <- system2("bcftools", c("view", "-O", "z", "-R", BED, "-o", VCF_SUB, VCF))
  if (!identical(st, 0L)) stop("bcftools view failed")
  st <- system2("bcftools", c("index", "-f", "-t", VCF_SUB))
  if (!identical(st, 0L)) stop("bcftools index of subset failed")
}
n_var <- as.integer(system2("bcftools", c("index", "-n", VCF_SUB), stdout = TRUE))
cat(sprintf("subset VCF: %s variants overlapping the %d gene bodies\n",
            format(n_var, big.mark = ","), length(genes)))

# --- 3. VCF -> GRanges, overlap-join onto genes ------------------------------------
qfile <- tempfile(fileext = ".tsv")
# Pass \t/\n as two-character escapes (shQuote so the shell does not eat the backslash).
fmt <- "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF_exomes\\n"
st <- system2("bcftools",
              c("query", "-f", shQuote(fmt, type = "sh"), "-o", qfile, VCF_SUB))
if (!identical(st, 0L)) stop("bcftools query failed")
vcf_tbl <- read_tsv(qfile, col_names = c("chrom", "pos", "ref", "alt", "AF_exomes"),
                    col_types = "ciccc", show_col_types = FALSE)
unlink(qfile)
vcf_tbl <- vcf_tbl |>
  mutate(alt = str_split(alt, ","), AF_exomes = str_split(AF_exomes, ",")) |>
  unnest(c(alt, AF_exomes)) |>
  mutate(AF_exomes = suppressWarnings(as.numeric(AF_exomes)),
         pos = as.integer(pos)) |>
  filter(!is.na(AF_exomes))
cat(sprintf("variant alleles after splitting: %s\n", format(nrow(vcf_tbl), big.mark = ",")))

vcf_gr <- GRanges(
  seqnames = vcf_tbl$chrom,
  ranges   = IRanges(start = vcf_tbl$pos,
                     end   = vcf_tbl$pos + pmax(nchar(vcf_tbl$ref), 1L) - 1L),
  ref = vcf_tbl$ref,
  alt = vcf_tbl$alt,
  AF_exomes = vcf_tbl$AF_exomes,
  vcf_pos = vcf_tbl$pos
)
seqlevelsStyle(vcf_gr) <- seqlevelsStyle(genes)[1]
ol <- join_overlap_left(vcf_gr, genes)
ol_tbl <- as_tibble(ol) |>
  filter(!is.na(gene_nv)) |>
  transmute(chrom = as.character(seqnames), pos = vcf_pos, ref, alt,
            AF_exomes, gene_nv, symbol)
cat(sprintf("variant-gene overlap rows: %s  (unique genes: %d)\n",
            format(nrow(ol_tbl), big.mark = ","), n_distinct(ol_tbl$gene_nv)))

# --- 4. restrict to TOPMed PTVs ----------------------------------------------------
ptv <- read_tsv(PTV_TSV, show_col_types = FALSE,
                col_select = c(chrom, pos, ref, alt, gene, symbol)) |>
  mutate(gene_nv_vep = str_remove(gene, "\\..*$"),
         chrom = as.character(chrom), pos = as.integer(pos),
         ref = as.character(ref), alt = as.character(alt)) |>
  distinct(chrom, pos, ref, alt, .keep_all = TRUE)
cat(sprintf("unique TOPMed PTV sites (genome-wide): %s\n",
            format(nrow(ptv), big.mark = ",")))

matched <- ol_tbl |>
  inner_join(ptv |> transmute(chrom, pos, ref, alt, gene_nv_vep, symbol_vep = symbol),
             by = c("chrom", "pos", "ref", "alt"))
cat(sprintf("overlap variants also in TOPMed PTVs: %s sites, %d genes\n",
            format(nrow(distinct(matched, chrom, pos, ref, alt)), big.mark = ","),
            n_distinct(matched$gene_nv)))

# --- 5. sum AF_exomes per gene -----------------------------------------------------
per_gene <- matched |>
  distinct(gene_nv, chrom, pos, ref, alt, .keep_all = TRUE) |>
  group_by(gene_nv) |>
  summarise(n_ptv_overlap = dplyr::n(),
            sum_AF_exomes = sum(AF_exomes),
            max_AF_exomes = max(AF_exomes),
            .groups = "drop")

out <- pairs |>
  dplyr::select(gene_nv, symbol) |>
  left_join(per_gene, by = "gene_nv") |>
  mutate(n_ptv_overlap = replace_na(n_ptv_overlap, 0L))

write_tsv(matched, SITES_OUT)
write_tsv(out, OUT)
cat(sprintf("wrote %s\n", OUT))
cat(sprintf("wrote %s\n", SITES_OUT))
have <- !is.na(out$sum_AF_exomes)
cat(sprintf("genes with summed AF_exomes: %d / %d\n", sum(have), nrow(out)))
if (any(have)) {
  s <- out$sum_AF_exomes[have]
  cat(sprintf("sum AF_exomes: median=%.4g  min=%.4g  max=%.4g\n",
              median(s), min(s), max(s)))
}
