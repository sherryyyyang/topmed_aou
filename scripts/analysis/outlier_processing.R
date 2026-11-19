ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
library(edgeR)
library(biomaRt)
library(limma)
library(RNOmni)
library(uwot)
library(gdsfmt)
library(SNPRelate)
# library(GWASTools)
library(corral)
# library(SingleCellExperiment)
# library(tidySingleCellExperiment)
library(scran)
library(scuttle)
library(vroom)
library(PCAtools)
library(rtracklayer)
library(plyranges)
library(WGCNA)
library(broom)
library(tidyverse)
library(magrittr)
library(readr)
library(RMTstat)
library(entropy)

pseudobulk_celltype <- function(cell_type, sce_object) {
    sce_subset <- filter(sce_object, celltype == !!cell_type)
    subsetted_counts <- sumCountsAcrossCells(sce_subset, id = sce_subset$PIDN)
    subsetted_counts
}

filter_expression <- function(counts) {
    genevar_poisson <- modelGeneVarByPoisson(counts)
    top_genes_poisson <- getTopHVGs(genevar_poisson)
    genevar_cv2 <- modelGeneCV2(counts)
    top_genes_cv2 <- getTopHVGs(genevar_cv2, var.field="ratio", var.threshold=1.0)
    top_var_intersect <- intersect(top_genes_poisson, top_genes_cv2)
    counts_subset <- counts[top_var_intersect,]
}

choose_highest <- function(grouped_rows, grouping_df) {
    arrange(grouped_rows, desc(ncells)) |> slice(1)
}

residualize_phenotypes <- function(n_expression_pcs, n_genotype_pcs, normalized_expression, original_metadata = NULL, expression_pcs, genotype_pcs) {
    expression_pc_names <- str_c("PC", 1:n_expression_pcs)
    expression_pc_subset <- dplyr::select(expression_pcs, ID, one_of(expression_pc_names))
    genotype_pc_names <- str_c("Genotype_PC", 1:n_genotype_pcs)
    genotype_pc_subset <- dplyr::select(genotype_pcs, ID, one_of(genotype_pc_names))
    # combined_phenodata <- left_join(original_metadata, expression_pc_subset) |>
    #     left_join(genotype_pc_subset) |> dplyr::select(-genotype_id)
    combined_phenodata <- left_join(expression_pc_subset, genotype_pc_subset)

    combined_phenotype_matrix <- column_to_rownames(combined_phenodata, "ID") |> as.matrix()
    # combined_phenotype_matrix <- select(original_metadata, -genotype_id) |> column_to_rownames("PIDN") |> as.matrix()
    # combined_phenotype_qr <- qr(combined_phenotype_matrix)
    # combined_phenotype_predictors <- combined_phenotype_qr$qr
    normalized_expression_ordered <- as.matrix(normalized_expression[,rownames(combined_phenotype_matrix)])

    expr_residuals <- removeBatchEffect(normalized_expression_ordered, covariates = combined_phenotype_matrix) 
    ## row(sd) / row mean --> CV 
    expr_scaled <- t(expr_residuals) |> scale() |> t()
    expr_scaled
}

map_expression_npcs <- function(n_genotype_pcs, n_pcs_list, ...) {
    scaled_list <- map(n_pcs_list, residualize_phenotypes, n_genotype_pcs = n_genotype_pcs, ...)
    names(scaled_list) <- str_c("expression_pc", n_pcs_list, sep = "_")
    scaled_list
}

subset_genotypes <- function(genes, svs_gr, flanking_window, af_threshold) {
    start(genes) <- start(genes) - flanking_window
    end(genes) <- end(genes) + flanking_window

    genes_subset <- select(genes, type, gene_id, gene_type, gene_name)

    svs_gr_join <- join_overlap_left(svs_gr, genes_subset)
    svs_gr_df <- as_tibble(svs_gr_join) |>
        filter(gene_type == "protein_coding") |>
        filter(AF <= af_threshold)
    # svs_gr_df$gene_id <- str_replace_all(svs_gr_df$gene_id.x, "\\..*$", "")
    # svs_gr_df$ID <- as.character(svs_gr_df$ID)
    svs_gr_df

    # sv_genes <- na.omit(unique(svs_gr_df$gene_id))
    # sv_samples <- na.omit(unique(svs_gr_df$ID))
    # sv_sample_gene_mat <- matrix(nrow = length(sv_genes), ncol=length(sv_samples), 0)
    # colnames(sv_sample_gene_mat) <- sv_samples
    # rownames(sv_sample_gene_mat) <- sv_genes
    # sv_sample_gene_mat
}

collapse_variants <- function(grouped_rows, grouped_df) {
    if (nrow(grouped_rows) > 1) {
        collapsed_df <- tibble(variant_id = str_c(grouped_rows$variant_id, collapse = ";"))
    } else {
        collapsed_df <- grouped_rows
    }
    collapsed_df
}

outlier_enrichment <- function(scaled_matrix, sv_subset, z_score) {
    scaled_matrix_longer <- as.data.frame(scaled_matrix) |>
        rownames_to_column("gene") |> pivot_longer(-gene, names_to = "sample", values_to = "expression_zscore")
    scaled_matrix_longer$sample <- str_remove_all(scaled_matrix_longer$sample, "_.*$")

    longer_filter <- filter(scaled_matrix_longer, is_in(sample, unique(sv_subset$sample)))
    colnames(longer_filter)[1] <- "gene_id"

    sv_join <- left_join(longer_filter, sv_subset)
    sv_join$underlier <- sv_join$expression_zscore < -z_score
    sv_join$overlier <- sv_join$expression_zscore > z_score
    sv_join$rv <- !is.na(sv_join$variant_id)

    risk_ratio_over <- fisher.test(sv_join$overlier, sv_join$rv)
    risk_ratio_over_df <- tibble(odds_ratio = risk_ratio_over$estimate, p_value = risk_ratio_over$p.value, ci_95_lo = risk_ratio_over$conf.int[1], ci_95_hi = risk_ratio_over$conf.int[2], direction = "overlier")
    risk_ratio_under <- fisher.test(sv_join$underlier, sv_join$rv)
    risk_ratio_under_df <- tibble(odds_ratio = risk_ratio_under$estimate, p_value = risk_ratio_under$p.value, ci_95_lo = risk_ratio_under$conf.int[1], ci_95_hi = risk_ratio_under$conf.int[2], direction = "underlier")

    bind_rows(risk_ratio_over_df, risk_ratio_under_df)
}

meta_data <- read_tsv("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/metadata.filtered.txt")

## gene counts
expr_data <- vroom("../counts/gene_counts.csv")
mapping <- setNames(meta_data$WGS, meta_data$RNA)
colnames(expr_counts_matrix) <- mapping[colnames(expr_counts_matrix)]
# colnames(expr_counts_matrix)<- gsub('"', '', colnames(expr_data))


# Expression data
write_rds(expr_data, "../counts/gene_counts.rds")

expr_names <- colnames(expr_data)
# Check that sample names match expression data
setdiff(expr_names, meta_data$RNA)
setdiff(meta_data$RNA, expr_names)

wgs_names <- read_tsv("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/genotypes/plink/chr1.fam", col_names = FALSE)
setdiff(wgs_names$X2, meta_data$WGS)
setdiff(meta_data$WGS, wgs_names$X2)

#sv_calls
sv_names <- read_lines("../sv/sv_subjects")
setdiff(sv_names, meta_data$WGS)
setdiff(meta_data$WGS, sv_names)

meta_data$SV_subjects <- ifelse(meta_data$WGS %in% sv_names, meta_data$WGS, NA)
write_tsv(meta_data, "metadata.filtered_with_sv.txt")

expression_sv_overlap <- filter(meta_data, !is.na(SV_subjects))
map_vec <- setNames(meta_data$WGS, meta_data$RNA)
library(UpSetR)

# Example vectors (replace with yours)
set1 <- unique(sv_names)
set2 <- unique(meta_data$RNA)
set3 <- unique(meta_data$WGS)

# Put into a list
input_list <- list(
  SV_vcfs = set1,
  TPM = set2,
  WGS = set3
)

# Convert list → UpSetR input
upset_data <- fromList(input_list)

png("upset_plot.png", width = 2000, height = 1500, res = 200)

upset(
    upset_data,
    nsets = 3,
    nintersects = NA,
    order.by = "freq",
    sets.x.label = "Set size"
)

dev.off()

# Genotype PCs

convert_bed_to_gds <- function(chr_num) {
    chr_name <- str_c("chr", chr_num)
    bed_name <- str_c(chr_name, ".bed")
    bim_name <- str_c(chr_name, ".bim")
    fam_name <- str_c(chr_name, ".fam")
    gds_name <- str_c(chr_name, ".gds")
    snpgdsBED2GDS(bed.fn = str_c("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/genotypes/plink/", bed_name),
        bim.fn = str_c("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/genotypes/plink/", bim_name),
        fam.fn = str_c("/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/genotypes/plink/", fam_name),
        out.gdsfn = str_c("../data/", gds_name)
    )
}

chrs <- c(1:21, "X")
catch <- map(chrs, convert_bed_to_gds)
gds_files <- list.files("../data", pattern = "*.gds", full.names = TRUE)
gds_output <- "../data/snv_genotypes.gds"
# Combine GDS files across chromosomes - significantly faster than merging with plink
snpgdsCombineGeno(gds_files, gds_output, method = "exact", compress.annotation = "", compress.geno = "")

genotype_gds <- snpgdsOpen(gds_output)
genotype_sample_ids <- snpgdsSummary(genotype_gds)
write_rds(genotype_sample_ids, "genotype_summary.rds")
shared_ids <- expression_sv_overlap$WGS

# structural_variant_ids <- read_lines("/oak/stanford/groups/smontgom/dnachun/data/1kg/vcf_structural/structural_variant_vcf_ids")
# shared_ids2 <- intersect(shared_ids, structural_variant_ids)
# write_lines(shared_ids2, "expr_ids")
# shared_ids2 <- read_lines("./expr_ids")

# Protein coding
ensembl_96 <- readGFF("../annotation/gencode.v30.annotation.gff3")
# ensembl <- useEnsembl(biomart = "genes", version = "96", dataset = "hsapiens_gene_ensembl")
# annotations <- getBM(
#       attributes = c("ensembl_gene_id_version", "hgnc_symbol", "gene_biotype"),
#       filters = "ensembl_gene_id_version",
#       values = expr_data$Name,
#       mart = ensembl
#     )
ensembl_96_df <- as_tibble(ensembl_96)
ensembl_96_df <- read_rds("./ensembl_96_df.rds")
ensembl_96_genes <- filter(ensembl_96_df, type == "gene")
annotations_coding <- filter(ensembl_96_genes, gene_type == "protein_coding")
expr_coding <- filter(expr_data, is_in(Name, annotations_coding$gene_id))
expr_coding_geno <- select(expr_coding, Name, one_of(expression_sv_overlap$RNA))
expr_coding_counts <- column_to_rownames(expr_coding_geno, "Name") |> as.matrix()
colnames(expr_coding_counts) <- map_vec[colnames(expr_coding_counts)]


# Identify highly variable genes using a mean variance trend analysis
# Also normalize with Freeman-Tukey residuals using corral

# IMPORTANT!!!! These steps and all subsequent analyses depend upon which samples you have chosen
# If you need to exclude samples later, you must eventually return to this point and do the exclusion here!
expr_filtered_counts <- filter_expression(expr_coding_counts)
expr_normalized <- corral_preproc(expr_filtered_counts, rtype = "freemantukey")
write_rds(expr_normalized, "expr_normalized.rds")
expr_normalized <- read_rds("expr_normalized.rds")

# Connectivity-based outlier identification
# Compute correlation between samples and identify ones
# whose correlation with all other samples is suspiciously low
# Typically remove samples with connectivity Z < -3
expr_norm_adj <- (0.5 + 0.5 * cor(expr_normalized))
colnames(expr_norm_adj) <- colnames(expr_normalized)
rownames(expr_norm_adj) <- colnames(expr_normalized)
expr_net_summary <- fundamentalNetworkConcepts(expr_norm_adj)
expr_net_connectivity <- expr_net_summary$Connectivity
expr_connectivity_zscore <- (expr_net_connectivity - mean(expr_net_connectivity)) / sd(expr_net_connectivity)
names(expr_connectivity_zscore) <- colnames(expr_normalized)
expr_connectivity_z3 <- names(expr_connectivity_zscore[expr_connectivity_zscore > -3])

n_remove_samples <- length(expr_connectivity_zscore) - length(expr_connectivity_z3)
print(str_c("Removed ", n_remove_samples, " with connectivity Z < -3"))

expr_qc_counts <- expr_coding_counts[,is_in(colnames(expr_coding_counts), expr_connectivity_z3)]
expr_qc_filtered_counts <- filter_expression(expr_qc_counts)
expr_qc_normalized <- corral_preproc(expr_qc_filtered_counts, rtype = "freemantukey")
colnames(expr_qc_normalized) <- colnames(expr_coding_counts)
write_rds(expr_qc_normalized, "expr_qc_normalized.rds")
expr_qc_normalized <-  read_rds("expr_qc_normalized.rds")

expr_pca <- as.matrix(expr_qc_normalized) |> pca(scale = TRUE)

write_rds(expr_pca, "expr_pca.rds")
expr_pca <- read_rds("expr_pca.rds")
sum(expr_pca$sdev^2)
sum(expr_pca$sdev[1:2000]^2)/sum(expr_pca$sdev^2)



eigenvalue_entropy <- entropy((expr_pca$sdev^2)/sum(expr_pca$sdev^2), method = "ML")

p <- screeplot(expr_pca, 
# Generate the scree plot
          axisLabSize = 5, 
          titleLabSize = 12,
          components = 1:100,
          colBar = "dodgerblue",
          drawCumulativeSumLine = TRUE,
          colCumulativeSumLine = "red2")

ggsave(file.path(OUT_DIR, "topmed_screeplot.png"), plot = p)

expr_eigenvectors <- expr_pca$rotated
rownames(expr_eigenvectors) <- colnames(expr_qc_normalized)

expr_eigenvectors_subset <- expr_eigenvectors[,1:55]

expr_counts_matrix <- column_to_rownames(expr_data, "Name") 
expr_counts_matrix_align <- as.matrix(expr_counts_matrix[rownames(expr_qc_normalized),rownames(expr_eigenvectors_subset)])
expr_dgelist <- DGEList(expr_counts_matrix_align)
expr_cpm <- edgeR::cpm(expr_dgelist, log = TRUE)
expr_cpm_adj <- removeBatchEffect(expr_cpm, covariates = expr_eigenvectors_subset)
write_rds(expr_cpm_adj, "expr_cpm_adj.rds")
expr_cpm_df <- as_tibble(expr_cpm_adj, rownames = "gene_id") |> pivot_longer(contains("NWD"), names_to = "subjectid", values_to = "expression_id")

expr_norm_adj <- removeBatchEffect(expr_qc_normalized, covariates = expr_eigenvectors_subset)
expr_norm_adj_pca <- pca(expr_norm_adj, scale = TRUE, center = TRUE)
write_rds(expr_norm_adj_pca, "expr_norm_adj_pca.rds")
expr_min_k_norm_adj <- chooseGavishDonoho(expr_norm_adj, var.explained = expr_norm_adj_pca$sdev^2, noise = 1)

expr_umap <- umap(expr_norm_adj_pca$rotated[,1:expr_min_k_norm_adj])
colnames(expr_umap) <- str_c("UMAP_", 1:ncol(expr_umap))
expr_umap_df <- as.data.frame(expr_umap) |> rownames_to_column(var = "WGS") |> left_join(meta_data)

umap_plot <- ggplot(expr_umap_df, aes(UMAP_1, UMAP_2, color = COHORT)) +
    geom_point() +
    theme_classic()

ggsave(file.path(OUT_DIR, "topmed_sample_umap_cohort.pdf"), umap_plot, device = cairo_pdf(), height = 4, width = 4)

get_underlier <- function(grouped_data, grouping_df, z_cutoff = -4) {
    mean_expr <- mean(grouped_data$expression_logcpm)
    filtered_data <- filter(grouped_data, expression_zscore < z_cutoff & expression_logcpm < (mean_expr - 1))
    filtered_data
}

expr_z <- vroom("./expr_scaled_genotype_pc_5_expression_pc_55_longer_form.tsv", num_threads = 1)
expr_z_join <- left_join(expr_z, expr_cpm_df)
colnames(expr_z_join)[4] <- "expression_logcpm"
write_rds(expr_z_join, "./expr_z_join.rds")

z_score_underlier <- group_by(expr_z_join, gene_id) |> group_modify(get_underlier)
write_rds(z_score_underlier, "z_score_underlier.rds")
z_score_tally <- group_by(z_score_underlier, gene_id) |> tally() |> arrange(desc(n))
z_score_tally$pct_underlier <- z_score_tally$n / ncol(expr_cpm)

z_score_tally_subj <- group_by(z_score_counts, subjectid) |> tally() |> arrange(desc(n))

expr_eigenvectors_df <- as.data.frame(expr_eigenvectors) |> rownames_to_column(var = "ID")
expr_eigenvalues <- expr_pca$sdev^2
expr_min_k <- chooseGavishDonoho(expr_normalized, var.explained = expr_eigenvalues, noise = 1)
expr_min_k_mp <- chooseMarchenkoPastur(expr_normalized, var.explained = expr_eigenvalues, noise = 1)

# Demonstration - do more formally later!
expr_pc_min <- expr_eigenvectors[,1:expr_min_k]
expr_pc_min_norm <- apply(expr_pc_min, 2, RankNorm)
expr_umap <- umap(expr_eigenvectors_subset)
colnames(expr_umap) <- str_c("UMAP_", 1:ncol(expr_umap))
expr_umap_df <- as.data.frame(expr_umap) |> rownames_to_column(var = "WGS") |> left_join(meta_data)

cluster_labels <- kmeans(expr_umap, centers = 2)$cluster

umap_df <- data.frame(
  umap_coords,
  Cluster = factor(cluster_labels)
)


umap_plot <- ggplot(expr_umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point() +
    theme_classic()

ggsave(file.path(OUT_DIR, "topmed_sample_umap.pdf"), umap_plot, device = cairo_pdf(), height = 4, width = 4)

umap_plot <- ggplot(expr_umap_df, aes(UMAP_1, UMAP_2, color = COHORT)) +
    geom_point() +
    theme_classic()

ggsave(file.path(OUT_DIR, "topmed_sample_umap_cohort.pdf"), umap_plot, device = cairo_pdf(), height = 4, width = 4)

# cd8t_genotype <- readRDS("./cd8t_genotype.rds")
# Covariates - explore if ncells as covariate max sense?
# Diagnosis is messy - revisit
# The "Batch" variable in celltype_metadata has one duplicate - maybe drop sample?
# cd8t_genotype_colData <- colData(cd8t_coding) |> as_tibble()
# cd8t_genotype_colData$Sample_Batch %<>% factor()
# cd8t_genotype_colData$Processing_Batch %<>% factor()
# cd8t_genotype_colData$Sex %<>% factor()
# colnames(cd8t_genotype_colData)[1] <- "PIDN"
# cd8t_design_matrix <- model.matrix( ~ Sample_Batch + Processing_Batch + Sex + Age, cd8t_genotype_colData)
# cd8t_covariates_df <- as.matrix(cd8t_design_matrix[,-1]) |> as_tibble()
# cd8t_covariates_df$PIDN <- cd8t_genotype_colData$PIDN

# Genotype PCA
# Note - converting between GDS and BED probably not needed
# First, subset genotypes to subjects with expression data
expr_snps_ldprune <- snpgdsLDpruning(genotype_gds, method = "corr", slide.max.bp = 10e6, maf = 0.05, missing.rate = 0.10, ld.threshold = sqrt(0.1), verbose = TRUE, num.thread = 16L, sample.id = shared_ids)
expr_snp_ids <- unlist(expr_snps_ldprune, use.names = TRUE)
write(expr_snp_ids, "expr_snp_ids.rda")
expr_snp_ids <- read_rds("expr_snp_ids.rda")
expr_genotype_pca <- snpgdsPCA(genotype_gds, snp.id = expr_snp_ids, sample.id = shared_ids, algorithm = "exact", eigen.method = "DSPEV", num.thread = 16L, eigen.cnt = -1)
expr_genotype_eigenvectors <- expr_genotype_pca$eigenvect
rownames(expr_genotype_eigenvectors) <- expr_genotype_pca$sample.id
colnames(expr_genotype_eigenvectors) <- str_c("Genotype_PC", 1:ncol(expr_genotype_eigenvectors))
expr_genotype_eigenvectors_df <- as.data.frame(expr_genotype_eigenvectors) |> rownames_to_column(var = "ID")
colnames(expr_genotype_eigenvectors_df)[1] <- "ID"

expr_genotype_min_k <- chooseGavishDonoho(.dim = c(length(expr_genotype_pca$eigenval), length(expr_genotype_pca$snp.id)), var.explained = expr_genotype_pca$eigenval, noise = 1)

# expr_combined_meta <- dplyr::select(expr_genotype_colData, -Age) |> left_join(expr_covariates_df) |>
#     dplyr::select(genotype_id, one_of(colnames(expr_covariates_df)))
# write_rds(expr_combined_meta, "expr_combined_meta.rds")
colnames(expr_eigenvectors_df)[1] <- "ID"
read_rds(expr_eigenvectors_df, "expr_eigenvectors_df.rds")
write_rds(expr_genotype_eigenvectors_df, "expr_genotype_eigenvectors_df.rds")
# expr_combined_meta <- read_rds("./expr_combined_meta.rds")
expr_eigenvectors_df <- read_rds("./expr_eigenvectors_df.rds")
expr_genotype_eigenvectors_df <- read_rds("./expr_genotype_eigenvectors_df.rds")

# Start here if you want to change the number of expression/genotype PCs to remove
n_expression_pcs <- seq(55, 4000, 200)
# n_expression_pcs <- seq(100, 500, 100)
# n_expression_pcs <- c(20)
# n_expression_pcs <- c(20, 25, 30, 35)
# n_expression_pcs <- 1
# n_expression_pcs <- seq(1, 10, 1)
n_genotype_pcs <- c(5)
# original_metadata = cd8t_combined_meta,
expr_scaled_lists <- map(n_genotype_pcs,
    map_expression_npcs,
    n_expression_pcs,
    normalized_expression = expr_qc_normalized,
    expression_pcs = expr_eigenvectors_df,
    genotype_pcs = expr_genotype_eigenvectors_df
)
names(expr_scaled_lists) <- str_c("genotype_pc", n_genotype_pcs, sep = "_")
expr_scaled_flat_list <- list_flatten(expr_scaled_lists)
write_rds(expr_scaled_flat_list, "expr_scaled_flat_list.rds")

chr_names <- c(1:22, "X", "Y")
gff <- readGFFAsGRanges("/oak/stanford/groups/smontgom/dnachun/data/AFGR/rnaseqc/Homo_sapiens.GRCh38.110.collapse.gtf")
seqlevels(gff)[1:24] <- c(str_c("chr", 1:22), "chrX", "chrY")
gff_genes <- gff[gff$type == "gene" & gff$gene_type == "protein_coding",]

write_rds(gff_genes, "gff_genes.rds")

# Also load SV data here and remove samples that are missing long reads
# svs <- vroom("/oak/stanford/groups/smontgom/syang/ADRC/unprocessed_data/pre_zscore_filtered_dna.txt", num_threads = 1)
# svs_nonzero <- filter(svs, gt > 0)
# svs_nonzero$length <- svs_nonzero$end - svs_nonzero$start
# write_rds(svs_nonzero, "svs_nonzero.rds")
# svs_nonzero <- read_rds("./svs_nonzero.rds")
# svs_nonzero$SCMD <- str_remove_all(svs_nonzero$ID, "^.*_")

sv_names <- list.files("/oak/stanford/groups/smontgom/dnachun/data/1kg/vcf_structural/individual_svs") |> str_subset("X", negate = TRUE)
svs_list <- list.files("/oak/stanford/groups/smontgom/dnachun/data/1kg/vcf_structural/individual_svs", full.names = TRUE) |>
    str_subset("X", negate = TRUE) |>
    map(vroom, col_names = FALSE)
names(svs_list) <- sv_names
svs <- bind_rows(svs_list, .id = "sample")
colnames(svs)[2:8] <- c("chr", "start", "end", "variant_id", "alt", "svlen", "AF")
write_rds(svs, "svs.rds")

svs_1mb <- filter(svs, abs(svlen) < 10000000)
svs_gr <- makeGRangesFromDataFrame(svs_1mb, keep.extra.columns = TRUE, start.field = "start", end.field = "end")
write_rds(svs_gr, "svs_gr.rds")

# Start here to tweak flanking window and AF for SVs
gff_genes <- read_rds("gff_genes.rds")
svs_gr <- read_rds("svs_gr.rds")
flanking_window <- 100000
af_threshold <- 0.01

subset_genotypes_100000_0.01 <- subset_genotypes(gff_genes, svs_gr, flanking_window, af_threshold)
subset_genotypes_100000_0.01_select <- select(subset_genotypes_100000_0.01, sample, gene_id, variant_id) |>
    group_by(sample, gene_id) |> group_modify(collapse_variants)
group_by(subset_genotypes_100000_0.01_select, sample) |> tally()  
group_by(subset_genotypes_100000_0.01_select, sample) |>
  tally() |>
  summarise(q25 = quantile(n, 0.25),
    median = quantile(n, 0.5),
    q75 = quantile(n, 0.75))

group_by(subset_genotypes_100000_0.01_select, sample) |>
  tally() |>
  summarise(mean_n = mean(n))


subset_genotypes_100000_0.01_select |> group_by(variant_id) |> summarise(n_samples = n_distinct(sample)) |> summarise(average_samples = mean(n_samples)) 
subset_genotypes_100000_0.01_select |> group_by(gene_id) |> summarise(n_samples = n_distinct(sample)) |> summarise(average_samples = mean(n_samples))
# subset_genotypes_100000_0.01_select$SCMD <- str_remove_all(subset_genotypes_100000_0.01_select$ID, "^.*_")
write_rds(subset_genotypes_100000_0.01_select, "./subset_genotypes_100000_0.01_select.rds")
subset_genotypes_100000_0.01_select <- read_rds("./subset_genotypes_100000_0.01_select.rds")
scaled_matrix_longer <- as.data.frame(expr_scaled_flat_list[[1]]) |>
    rownames_to_column("gene") |> pivot_longer(-gene, names_to = "sample", values_to = "expression_zscore")
group_by(scaled_matrix_longer, sample) |> tally()

# RV-outlier enrichment analysis begins here
subset_genotypes_100000_0.01_select <- read_rds("./subset_genotypes_100000_0.01_select.rds")
expr_scaled_flat_list <- read_rds("./expr_scaled_flat_list.rds")
# expr_enrichments <- map(expr_scaled_flat_list, outlier_enrichment, subset_genotypes_100000_0.01_select, 4)
# expr_enrichments_df <- bind_rows(expr_enrichments, .id = "npcs")
# write_rds(expr_enrichments_df, "expr_enrichments_npc_10_130.rds")
# filter(expr_enrichments_df, direction == "underlier") |> arrange(desc(odds_ratio))
# filter(expr_enrichments_df, direction == "underlier") 
# filter(expr_enrichments_df, direction == "overlier") |> arrange(desc(odds_ratio))
# expr_enrichments_df$npcs_format <- as.integer(str_remove_all(expr_enrichments_df$npcs, "^.*_expression_pc_"))
# expr_enrichments_underliers <- filter(expr_enrichments_df, direction == "underlier")

# enrichment_plot <- ggplot(expr_enrichments_underliers, aes(npcs_format, odds_ratio, ymin = ci_95_lo, ymax = ci_95_hi)) +
#     geom_point() +
#     geom_line(group = 1) +
#     geom_errorbar() +
#     geom_hline(yintercept = 1.5, color = "red") +
#     theme_classic() +
#     xlab("# Expression PCs") +
#     ylab("RV/outlier odds ratio")
# ggsave(file.path(OUT_DIR, "enrichment_npc_35_295_3.5.pdf"), enrichment_plot, device = cairo_pdf(), width = 4, height = 3)

cutoffs <- c(3.5, 4, 5)

# run enrichment for each cutoff
expr_enrichments_all <- map_dfr(cutoffs, function(cutoff) {
  expr_enrichments <- map(expr_scaled_flat_list, outlier_enrichment,
                          subset_genotypes_100000_0.01_select, cutoff)
  expr_enrichments_df <- bind_rows(expr_enrichments, .id = "npcs")
  expr_enrichments_df$npcs_format <- as.integer(str_remove_all(expr_enrichments_df$npcs, "^.*_expression_pc_"))
  expr_enrichments_df$cutoff <- cutoff
  return(expr_enrichments_df)
})

# save combined data
write_rds(expr_enrichments_all, "./expr_enrichments_all_cutoffs.rds")
expr_enrichments_all <- read_rds("./expr_enrichments_all_cutoffs.rds")
# prepare underliers for plotting
expr_enrichments_underliers <- filter(expr_enrichments_all, direction == "underlier")
expr_enrichments_underliers |> filter(cutoff == 4) |> arrange(desc(odds_ratio))

# make combined plot
enrichment_plot <- ggplot(expr_enrichments_underliers,
                          aes(x = npcs_format, y = odds_ratio, color = as.factor(cutoff))) +
  geom_point() +
  geom_line(aes(group = cutoff)) +
  geom_errorbar(aes(ymin = ci_95_lo, ymax = ci_95_hi), width = 0.2) +
  geom_hline(yintercept = 1.5, color = "red", linetype = "dashed") +
  theme_classic() +
  xlab("# Expression PCs") +
  ylab("RV/outlier odds ratio") +
  labs(color = "Z-score cutoff") +
  theme(legend.position = "right")

# save figure
ggsave(file.path(OUT_DIR, "enrichment_npc_all_cutoffs.pdf"), enrichment_plot, device = cairo_pdf(), width = 5, height = 4)

## selected based on previous function
expr_scaled_flat_list <- read_rds("./expr_scaled_flat_list.rds")
selected_names <- str_subset(names(expr_scaled_flat_list), "genotype_pc_5_expression_pc_2005$")  # optimal pc selected using gavish donoho method
selected_data <- expr_scaled_flat_list[selected_names]
selected_data_longer <- as.data.frame(selected_data[[1]]) |>
    rownames_to_column("gene_id") |>
    pivot_longer(-gene_id, names_to = "subjectid", values_to = "expression_zscore")
#selected_data_longer %<>% left_join(annotations_minimal)
selected_data_longer$gene_id %<>% str_remove_all("_.*$")
#colnames(selected_data_longer) <- c("gene_id", "subjectid","expression_zscore")

filename <- str_c("expr_scaled_", selected_names, "_longer_form.tsv")
write_tsv(selected_data_longer, filename)
# expr_scaled_geno1_expr10_longer_form <- vroom("expr_scaled_geno1_expr10_longer_form.tsv")

gff_genes_df <- as.data.frame(mcols(gff_genes)) %>% select(gene_name, gene_id) 
gff_genes_df$gene_id <- sub("\\..*$", "", gff_genes_df$gene_id)
gff_genes_df %>% distinct()

selected_data_longer <- left_join(
    selected_data_longer,
    gff_genes_df,
    by = c("gene" = "gene_name")
)

# annotations_duplicated <- filter(annotations_coding, duplicated(hgnc_symbol) | duplicated(hgnc_symbol, fromLast = TRUE)) |> 
#     pull("hgnc_symbol") |>
#     unique()

## outliers from esther, limit to the sv and subjects matched to those from the RNA outliers list
splicing_outliers_fraser <- vroom("/oak/stanford/groups/smontgom/erobb/data/annotations/outliers/AFGR_MAGE.fraser2_filtered.tsv")
splicing_outliers_fraser_filtered <- splicing_outliers_fraser %>% filter (SubjectID %in% sv_names)
colnames(splicing_outliers_fraser_filtered) <- c("SubjectID","gene_id","zscore", "deltaPsi")
write_tsv(splicing_outliers_fraser_filtered, "/oak/stanford/groups/smontgom/syang/AFGR/AFGR_MAGE.fraser2_splicing_filtered_less_subjects.tsv")
