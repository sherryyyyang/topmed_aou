#!/usr/bin/env Rscript
# Gene-count matrix used by run_underlier_pipeline.R.
# Read the cohort gene-count parquet (index column Name = Ensembl gene id,
# remaining columns = RNA sample ids) and write genes x samples CSV.
#
#   export GENE_COUNTS_PARQUET=/path/to/gene-counts.parquet
#   export OUT_CSV=$DATA_ROOT/counts/gene_counts.csv
#   Rscript scripts/counts/parquet_to_counts.R

suppressPackageStartupMessages(library(arrow))

parquet_path <- Sys.getenv("GENE_COUNTS_PARQUET", unset = "")
out_csv <- Sys.getenv("OUT_CSV", unset = "")
if (!nzchar(parquet_path) || !nzchar(out_csv)) {
  stop("Set GENE_COUNTS_PARQUET and OUT_CSV (see header).")
}
if (!file.exists(parquet_path)) stop("Parquet not found: ", parquet_path)

counts <- read_parquet(parquet_path)
if (!"Name" %in% names(counts)) stop("Parquet has no Name column: ", parquet_path)

counts <- counts[, c("Name", setdiff(names(counts), "Name"))]
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(counts, out_csv, row.names = FALSE, quote = FALSE)
message("Wrote ", nrow(counts), " genes x ", ncol(counts) - 1L, " samples -> ", out_csv)
