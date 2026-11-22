# TOPMed / AoU underlier figures

Scripts that build main Figures 1–4. Inputs are not in this repository. Set `DATA_ROOT` to the analysis directory (it holds `data/` and these `scripts/`). Set `PC_SETTING` to the number of expression PCs removed before z-scores. That count is an input; the selection rule is not in this repository.

An underlier is a gene–subject pair with expression z < −3 AND log2CPM_adj < gene_mean − 1, where gene_mean is the per-gene mean of the same PC-adjusted log2CPM.

Export `TOPMED_OUTLIERS_ROOT=$DATA_ROOT` and `PC_SETTING`. Panel scripts that assign `ROOT` before sourcing `scripts/helpers/_paths.R` must point `ROOT` at `$DATA_ROOT`. Only the z-score table changes with `PC_SETTING`. log2CPM_adj is read from `$DATA_ROOT/data/expression/expr_cpm_adj.rds`.

```bash
export DATA_ROOT=/path/to/analysis
export TOPMED_OUTLIERS_ROOT="$DATA_ROOT"
export PC_SETTING=<expression PC count>
cd "$DATA_ROOT"
./scripts/run_pipeline_pc.sh
```

## 0. Source VCFs to PLINK

`scripts/genotypes/vcf_to_plink.sh` is the standalone genotype fetch. It reads `VCF_DIR/chr{1..22,X}.vcf.gz` and a sample list. For each chromosome, bcftools keeps those samples, SNPs or indels with `|ILEN| < 50`, sites with minor allele count at least 1, and `FILTER` of `PASS` or `.`. `plink2 --vcf` then writes `$OUT_DIR/plink/chr*.bed`, `.bim`, and `.fam`. `run_underlier_pipeline.R` reads that PLINK directory.

`scripts/genotypes/filter_cohort_by_fmissing_pass_snp_indel_30bp_array.sbatch` is a stricter bcftools pass on chromosome VCFs: `FILTER=PASS`, then `F_MISSING < 0.1`, then SNP or indel with `|ILEN| <= 30`. Set `COHORT_VCF_DIR` and `OUTDIR`. It does not read PLINK and it does not call underliers.

## 0b. Structural variants

`scripts/genotypes/process_svs.sh` concatenates the per-chromosome PASS SV VCFs in `SV_VCF_DIR` (`chr*.vcf.gz`, muCNV PASS callset) into `$OUT_DIR/svs_pass.vcf.gz`. For each id in `SV_SUBJECTS` it writes `$OUT_DIR/individual_svs/<sample>`: sites with `AF < 0.01`, `F_MISSING < 0.1`, and an alternate genotype (`CHROM`, `POS`, `END`, `ID`, `ALT`, `SVLEN`, `AF`).

Figure 3b color is a separate step. `scripts/fig3_ptv_vs_rna/build_fig3b_segdup_svtype_genes.R` reads those same per-chromosome PASS SV VCFs, keeps SegDup-pass genes, and labels each gene by its gene-overlapping DEL, DUP, or INV with the largest `|SVLEN|` (any allele frequency).

## 0c. Gene-count parquet to the counts matrix

`scripts/counts/parquet_to_counts.R` is the expression matrix the underlier run starts from. It reads `GENE_COUNTS_PARQUET` (column `Name` = Ensembl gene id, remaining columns = RNA sample ids) and writes `OUT_CSV`, genes × samples. `run_underlier_pipeline.R` reads that CSV. It does not read the parquet.

```bash
export GENE_COUNTS_PARQUET=/path/to/gene-counts.parquet
export OUT_CSV=$DATA_ROOT/counts/gene_counts.csv
Rscript scripts/counts/parquet_to_counts.R
```

## 1. Counts and genotypes to z and adjusted log2CPM

`scripts/analysis/run_underlier_pipeline.R` reads that CSV (`CFG$counts_csv`), a GENCODE gene model, and the PLINK files from section 0. Those PLINK files are the genotypes made from the source VCFs. Set `n_expr_pcs` to `PC_SETTING`. It writes PC-adjusted log2CPM and z-scores in its output directory. Copy `expr_cpm_adj.rds` to `$DATA_ROOT/data/expression/expr_cpm_adj.rds`.

If `$DATA_ROOT/data/expression/expr_scaled_flat_list.rds` already exists, `scripts/pc_swap/write_z_longform_from_flat.R` writes `$DATA_ROOT/data/expression/expr_scaled_genotype_pc_5_expression_pc_${PC_SETTING}_longer_form.tsv`.

## 2. Underlier calls

`scripts/pc_swap/build_underlier_artifacts.R` reads that long-form z table and `expr_cpm_adj.rds`. It keeps pairs with z < −3 AND log2CPM_adj < gene_mean − 1 and writes `$DATA_ROOT/data/underliers_pc${PC_SETTING}/` (`z_score_underlier_z_-3_clean.rds` and the per-gene and per-subject tallies). Figures 1–3 read those objects. Also required: `$DATA_ROOT/data/core/kept_subjects.rds`.

## 3. TOPMed vs AoU prevalence

`scripts/pc_swap/build_topmed_vs_aou_prevalence.R` joins the TOPMed tallies to the AoU prevalence release under `$DATA_ROOT/call-RunRnaUnderlier/rna_outputs`. It writes `topmed_vs_new_z_-3_prevalence_per_gene.tsv`. Figures 1d, 2d, and 4 read that table through `scripts/helpers/load_topmed_aou_prevalence_scatter.R`.

## 4. Figure 3 dosage tables

`scripts/fig3_ptv_vs_rna/build_cnv_multicopy_gene_set.R` overlaps GENCODE gene bodies with the UCSC segmental-duplication track. `scripts/fig3_ptv_vs_rna/build_fig3b_segdup_svtype_genes.R` labels SegDup-pass genes by the primary gene-overlapping PASS SV (DEL, DUP, or INV, largest absolute SVLEN) in the PASS SV VCFs. Figure 3b reads `fig3b_segdup_svtype_genes.tsv`, `fig3b_segdup_sv_gene_hits.tsv`, and `fig3b_segdup_sv_carriers.tsv` under `$DATA_ROOT/data/references/`. The scatter also reads the PTV table (`PTV_TSV`, AF already < 0.5%), `loeuf_lookup.rds`, and `expr_cpm_adj.rds`.

## 5. Figure 4 gnomAD column

`scripts/fig4_acmg/sum_af_exome_topmed_overlap.R` sums gnomAD exome allele frequency for TOPMed PTVs on the shared-gene list. `scripts/fig4_acmg/plot_shared_rare_outlier_genes_4dot.R` reads `sum_af_exome_topmed_overlap.tsv` plus the prevalence table and the PTV table. `scripts/fig4_acmg/plot_acmg_4dot_no_ci.R` reads that table and the ACMG SF v3.2 CSV (`TOPMED_ACMG_CSV`).

## 6. Draw Figures 1–4

`scripts/assemble_nature_figures.R` stacks the panels (`ONLY_FIG=1` through `4`). `scripts/run_pipeline_pc.sh` runs the stages above, with Figure 2 at `DEFINITION_LABEL=z_-3`. Shared code: `scripts/helpers/_paths.R`, `paper_theme.R`, `_underlier_def.R`.

Figure 1 reads the underlier tallies and the prevalence table: `fig1_overview/make_fig1ab_cohort_overviews.R` (calls `recreate_publication_figure.py`), `plot_underlier_definition_demo.R`, `plot_underlier_tallies_topmed_aou.R`, `plot_topmed_vs_aou_prevalence_uncolored.R`.

Figure 2 reads underlier calls, `loeuf_lookup.rds`, and the prevalence table: `fig2_constraint/plot_loeuf_calibration.R`, `plot_underlier_tallies_loeuf_constrained_tolerant.R`, `plot_topmed_vs_aou_prevalence_loeuf.R`.

Figure 3 reads underlier calls, the PTV table, and the dosage tables. The histogram hardcodes eight gene symbols and does not open the ClinGen HI CSV. Scripts: `fig3_ptv_vs_rna/plot_ptv_vs_nonptv_zscore.R`, `plot_ptv_vs_underlier_gene_level_2x2.R`, `plot_candidate_gene_zscore_dist_clingen_z3.R`.

Figure 4 reads the prevalence table, the PTV table, the gnomAD AF sum, and the ACMG list: `fig4_acmg/plot_shared_rare_outlier_genes_4dot.R`, `plot_acmg_4dot_no_ci.R`.

The assembler writes Figure 1–4 png/pdf files into the delivery directory for this run. Those images are not part of this repository.

## External annotations

Download URLs, versions, and the filename under `$DATA_ROOT` for each annotation Figures 1–4 actually read are in [docs/ANNOTATIONS.md](docs/ANNOTATIONS.md). ClinVar is not one of those inputs.
