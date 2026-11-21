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

## 1. Counts and genotypes to z and adjusted log2CPM

`scripts/analysis/run_underlier_pipeline.R` reads gene counts, a GENCODE gene model, and PLINK genotypes from the paths in its CONFIG block. Set `n_expr_pcs` to `PC_SETTING`. It writes PC-adjusted log2CPM and z-scores in its output directory. Copy `expr_cpm_adj.rds` to `$DATA_ROOT/data/expression/expr_cpm_adj.rds`.

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

Figure 3 reads underlier calls, the PTV table, the dosage tables, and the ClinGen HI list: `fig3_ptv_vs_rna/plot_ptv_vs_nonptv_zscore.R`, `plot_ptv_vs_underlier_gene_level_2x2.R`, `plot_candidate_gene_zscore_dist_clingen_z3.R`.

Figure 4 reads the prevalence table, the PTV table, the gnomAD AF sum, and the ACMG list: `fig4_acmg/plot_shared_rare_outlier_genes_4dot.R`, `plot_acmg_4dot_no_ci.R`.

The assembler writes Figure 1–4 png/pdf files into the delivery directory for this run. Those images are not part of this repository.
