# Adapted from the user's outlier_processing.R (provided 2026-11-09).
# Bulk RNA-seq POC: count filtering, Freeman-Tukey normalization, connectivity QC,
# expression/genotype PC removal, gene-wise z scores and a logCPM underlier gate.
# No phenotype or rare-variant association is used to choose PCs.

read_table <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  } else {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  }
}

need_columns <- function(x, columns, label) {
  if (!all(columns %in% names(x))) stop(label, ": required columns ", paste(columns, collapse = ", "))
}

unique_key <- function(x, column) {
  if (anyNA(x[[column]]) || any(!nzchar(x[[column]])) || anyDuplicated(x[[column]])) {
    stop("Missing or duplicate ", column)
  }
}

strip_gene_version <- function(x) sub("\\.[0-9]+$", "", x)

pc_count <- function(x, label) {
  if (is.null(x) || length(x) != 1 || !is.numeric(x) || !is.finite(x) || x < 0 || x != as.integer(x)) {
    stop("Set ", label, " explicitly to a nonnegative integer (0 is allowed)")
  }
  as.integer(x)
}

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")

residualize_and_scale <- function(expression, covariates) {
  # Both matrices have exactly the same ordered samples. Fit with an intercept,
  # retain the gene mean (like limma::removeBatchEffect), then scale within gene.
  if (!identical(colnames(expression), rownames(covariates))) stop("Covariate/sample order mismatch")
  design <- cbind(Intercept = 1, covariates)
  if (any(!is.finite(design)) || qr(design)$rank != ncol(design)) stop("Covariates are missing, constant, or linearly dependent")
  if (nrow(design) - ncol(design) < 3) stop("Too few residual degrees of freedom")
  fit <- lm.fit(design, t(expression))
  residuals <- t(fit$residuals)
  adjusted <- sweep(residuals, 1, rowMeans(expression), "+")
  spread <- apply(residuals, 1, sd)
  valid <- is.finite(spread) & spread > sqrt(.Machine$double.eps)
  z <- matrix(NA_real_, nrow(expression), ncol(expression), dimnames = dimnames(expression))
  if (any(valid)) z[valid, ] <- sweep(sweep(adjusted[valid, , drop = FALSE], 1,
                                         rowMeans(adjusted[valid, , drop = FALSE]), "-"),
                                    1, spread[valid], "/")
  list(adjusted = adjusted, z = z, callable = valid)
}

underlier_gate <- function(z, adjusted_logcpm, z_cutoff, logcpm_drop) {
  if (length(z_cutoff) != 1 || !is.finite(z_cutoff) || z_cutoff >= 0) stop("underlier_z must be negative")
  if (length(logcpm_drop) != 1 || !is.finite(logcpm_drop) || logcpm_drop < 0) stop("logcpm_drop must be nonnegative")
  z < z_cutoff & adjusted_logcpm < (rowMeans(adjusted_logcpm) - logcpm_drop)
}

prepare_expression <- function(config) {
  for (p in c("edgeR", "corral", "scran")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Install Bioconductor package: ", p)
  }
  if (!identical(config$release, "v9")) stop("This POC config must explicitly specify release v9")
  if (is.null(config$genes_of_interest) || !length(config$genes_of_interest) || anyDuplicated(config$genes_of_interest)) stop("Specify unique genes_of_interest")
  output <- config$output_dir
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  metadata <- read_table(config$sample_metadata_file)
  need_columns(metadata, c("sample_id", "subject_id"), "metadata")
  unique_key(metadata, "sample_id")
  unique_key(metadata, "subject_id") # one selected RNA sample per participant
  annotations <- read_table(config$gene_annotation_file)
  need_columns(annotations, c("gene_id", "gene", "gene_type"), "annotations")
  annotations$gene_id <- strip_gene_version(annotations$gene_id)
  unique_key(annotations, "gene_id")
  counts_df <- read_table(config$counts_file)
  need_columns(counts_df, "gene_id", "counts")
  counts_df$gene_id <- strip_gene_version(counts_df$gene_id)
  unique_key(counts_df, "gene_id")
  if (anyDuplicated(names(counts_df))) stop("Duplicated sample columns")
  sample_ids <- setdiff(names(counts_df), "gene_id")
  if (!setequal(sample_ids, metadata$sample_id)) stop("Counts and metadata sample sets must match exactly; subset them explicitly first")
  samples <- metadata$sample_id
  counts <- as.matrix(counts_df[, samples, drop = FALSE])
  suppressWarnings(storage.mode(counts) <- "double")
  rownames(counts) <- counts_df$gene_id
  colnames(counts) <- metadata$subject_id
  if (any(!is.finite(counts)) || any(counts < 0) || any(abs(counts - round(counts)) > 1e-6)) stop("Input must be finite, nonnegative raw integer counts; not TPM/logcounts")
  if (ncol(counts) < 8 || any(colSums(counts) == 0)) stop("Need >=8 participants and nonzero libraries")
  annotated <- intersect(rownames(counts), annotations$gene_id[annotations$gene_type == "protein_coding"])
  counts <- counts[annotated, , drop = FALSE]
  if (nrow(counts) < 3) stop("Too few annotated protein-coding genes")
  adequately_expressed <- function(x) {
    if (config$min_expression_fraction <= 0 || config$min_expression_fraction > 1 || config$min_cpm <= 0) stop("Invalid expression filter")
    rowSums(edgeR::cpm(x) >= config$min_cpm) >= max(3, ceiling(ncol(x) * config$min_expression_fraction))
  }
  select_hvg <- function(x) {
    poisson <- scran::getTopHVGs(scran::modelGeneVarByPoisson(x))
    cv2 <- scran::getTopHVGs(scran::modelGeneCV2(x), var.field = "ratio", var.threshold = 1)
    genes <- intersect(poisson, cv2)
    if (length(genes) < 3) stop("Too few variable genes for the original scran filter; revise QC/filter settings explicitly")
    genes
  }
  keep <- adequately_expressed(counts)
  initial <- counts[keep, , drop = FALSE]
  hvg_initial <- select_hvg(initial)
  normalized_initial <- corral::corral_preproc(initial, rtype = "freemantukey")
  network <- 0.5 + 0.5 * cor(normalized_initial[hvg_initial, , drop = FALSE])
  diag(network) <- 0
  if (any(!is.finite(network))) stop("Undefined sample correlations")
  connectivity <- colSums(network)
  connectivity_z <- if (sd(connectivity) > 0) as.numeric(scale(connectivity)) else rep(0, length(connectivity))
  retained <- connectivity_z > config$connectivity_z_min
  write_csv(data.frame(subject_id = colnames(counts), connectivity_z = connectivity_z, retained = retained), file.path(output, "sample_qc.csv"))
  counts <- counts[, retained, drop = FALSE]
  if (ncol(counts) < 8) stop("Too few samples after connectivity QC")
  # Recompute filtering and normalization after sample exclusion; never overwrite
  # retained column names with the original unfiltered sample names.
  keep <- adequately_expressed(counts)
  counts <- counts[keep, , drop = FALSE]
  hvg <- select_hvg(counts)
  normalized <- corral::corral_preproc(counts, rtype = "freemantukey")
  dimnames(normalized) <- dimnames(counts)
  pca_input <- t(normalized[hvg, , drop = FALSE])
  max_pcs <- pc_count(config$max_expression_pcs, "max_expression_pcs")
  max_pcs <- min(max_pcs, ncol(pca_input), nrow(pca_input) - 1L)
  if (max_pcs < 1) stop("Prepare at least one expression PC")
  fit <- prcomp(pca_input, center = TRUE, scale. = TRUE, rank. = max_pcs)
  pcs <- fit$x
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  write_csv(data.frame(subject_id = rownames(pcs), pcs, check.names = FALSE), file.path(output, "expression_pcs.csv"))
  write_csv(data.frame(component = seq_along(fit$sdev), variance_fraction = fit$sdev^2 / sum(fit$sdev^2)), file.path(output, "expression_pc_variance.csv"))
  # TMM is an explicit adaptation for this bulk-RNA POC.
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  # All callable requested genes remain eligible; HVGs are for PCA/QC only.
  annotations <- annotations[match(rownames(normalized), annotations$gene_id), , drop = FALSE]
  targets <- annotations$gene %in% config$genes_of_interest
  duplicate_symbols <- unique(annotations$gene[targets & (duplicated(annotations$gene) | duplicated(annotations$gene, fromLast = TRUE))])
  if (length(duplicate_symbols)) stop("Resolve multiple gene IDs for target symbols: ", paste(duplicate_symbols, collapse = ", "))
  target_annotations <- annotations[targets, , drop = FALSE]
  target_normalized <- normalized[targets, , drop = FALSE]
  target_logcpm <- logcpm[targets, , drop = FALSE]
  rownames(target_normalized) <- target_annotations$gene
  rownames(target_logcpm) <- target_annotations$gene
  coverage <- data.frame(gene = config$genes_of_interest,
                         expression_status = ifelse(config$genes_of_interest %in% target_annotations$gene, "callable_pending_residualization", "not_detected_or_filtered"))
  write_csv(coverage, file.path(output, "expression_gene_coverage.csv"))
  saveRDS(list(normalized = target_normalized, logcpm = target_logcpm,
               expression_pcs = pcs, genes = target_annotations, config = config),
          file.path(output, "prepared_expression.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output, "R_session_prepare.txt"))
  jsonlite::write_json(config, file.path(output, "prepare_config.json"), pretty = TRUE, auto_unbox = TRUE)
  message("Preparation complete. Inspect QC and PC variance; set both PC counts explicitly before calling underliers.")
}

call_underliers <- function(config) {
  n_expr <- pc_count(config$n_expression_pcs, "n_expression_pcs")
  n_geno <- pc_count(config$n_genotype_pcs, "n_genotype_pcs")
  output <- config$output_dir
  prepared <- readRDS(file.path(output, "prepared_expression.rds"))
  stable_fields <- setdiff(names(prepared$config), c("n_expression_pcs", "n_genotype_pcs", "genotype_pcs_file", "additional_covariates_file", "underlier_z", "logcpm_drop"))
  if (!identical(prepared$config[stable_fields], config[stable_fields])) stop("Preparation inputs/config changed; rerun prepare")
  expression <- prepared$normalized
  if (!nrow(expression)) stop("None of the requested genes passed expression filtering")
  ids <- colnames(expression)
  if (n_expr > ncol(prepared$expression_pcs)) stop("Requested more expression PCs than prepared; increase max_expression_pcs and rerun prepare")
  covariates <- prepared$expression_pcs[ids, seq_len(n_expr), drop = FALSE]
  if (n_geno > 0) {
    if (is.null(config$genotype_pcs_file)) stop("genotype_pcs_file required when n_genotype_pcs > 0")
    genotype <- read_table(config$genotype_pcs_file)
    required <- paste0("Genotype_PC", seq_len(n_geno))
    need_columns(genotype, c("subject_id", required), "genotype PCs")
    unique_key(genotype, "subject_id")
    if (!all(ids %in% genotype$subject_id)) stop("Missing retained participants in genotype PCs")
    matrix <- as.matrix(genotype[match(ids, genotype$subject_id), required, drop = FALSE])
    suppressWarnings(storage.mode(matrix) <- "double")
    covariates <- cbind(covariates, matrix)
  }
  if (!is.null(config$additional_covariates_file)) {
    extra <- read_table(config$additional_covariates_file)
    need_columns(extra, "subject_id", "additional covariates")
    unique_key(extra, "subject_id")
    if (!all(ids %in% extra$subject_id)) stop("Missing retained participants in additional covariates")
    matrix <- as.matrix(extra[match(ids, extra$subject_id), setdiff(names(extra), "subject_id"), drop = FALSE])
    suppressWarnings(storage.mode(matrix) <- "double")
    covariates <- cbind(covariates, matrix)
  }
  rownames(covariates) <- ids
  scaled <- residualize_and_scale(expression, covariates)
  cpm_adjusted <- residualize_and_scale(prepared$logcpm, covariates)$adjusted
  under <- underlier_gate(scaled$z, cpm_adjusted, config$underlier_z, config$logcpm_drop)
  # Samples x genes, vectorized in column-major order.
  long <- expand.grid(subject_id = ids, gene = rownames(expression), stringsAsFactors = FALSE)
  long$expression_zscore <- as.vector(t(scaled$z))
  long$expression_logcpm <- as.vector(t(cpm_adjusted))
  long$underlier <- as.vector(t(under))
  long <- long[is.finite(long$expression_zscore), , drop = FALSE]
  if (!nrow(long)) stop("No genes have residual variance after adjustment")
  write_csv(long, file.path(output, "expression_calls.csv"))
  write_csv(long[, c("subject_id", "gene")], file.path(output, "eligible.csv"))
  write_csv(long[long$underlier, c("subject_id", "gene"), drop = FALSE], file.path(output, "underliers.csv"))
  write_csv(data.frame(gene = rownames(expression), callable = scaled$callable), file.path(output, "residual_gene_coverage.csv"))
  jsonlite::write_json(config, file.path(output, "call_config.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(output, "R_session_call.txt"))
  message("Underliers and gene-specific eligible populations written to ", output)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || !args[1] %in% c("prepare", "call")) stop("Usage: Rscript R/outlier_processing.R prepare|call config.json")
  config <- jsonlite::read_json(args[2], simplifyVector = TRUE)
  if (args[1] == "prepare") prepare_expression(config) else call_underliers(config)
}
