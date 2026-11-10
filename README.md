# All of Us expression-underlier phenotype enrichment

Proof of concept for All of Us v9. The R script follows the user's outlier-processing approach: raw bulk RNA-seq counts, count filtering, Freeman–Tukey normalization, connectivity QC, expression-PC adjustment, optional genotype-PC adjustment, gene-wise standardized expression, and a low-expression gate. It writes one eligibility row per participant and selected gene plus the called underliers.

## Important input requirements

No v9 multiomics files were supplied with the request, so this repo contains no fabricated v9 schema or participant data. Before running, adapt the example input paths/column names to the actual v9 files in the Workbench. Counts must be raw integer gene counts with one gene per row and one sample per column. Metadata must map each sample to exactly one participant. Annotation must map Ensembl gene IDs to symbols and `protein_coding`. The script deliberately errors when identifiers, sample sets, genes, or PC counts do not match.

Set `genes_of_interest` in `config/poc.example.json` to only the genes in the planned analysis. The example names RBM20 and LMNA only to show the format. After inspecting sample QC and the expression-PC variance output, set `n_expression_pcs` and `n_genotype_pcs` explicitly. Zero genotype PCs is valid. Do not use disease labels to pick PCs.

```sh
Rscript R/outlier_processing.R prepare config/poc.example.json
# Review results/poc_v9/sample_qc.csv and expression_pc_variance.csv first.
# Set both PC counts in a private config copy, then:
Rscript R/outlier_processing.R call private_inputs/poc.v9.private.json
```

Outputs are written under `results/` and ignored by Git. Keep subject-level data and Workbench outputs inside the Workbench. Do not upload them to GitHub.

## Phenotype enrichment

`data/gene_disease_associations.csv` is the ACMG SF v3.3 gene-disease source (84 genes, 100 rows). `data/gene_phenotype_map.csv` is an explicit, versioned PhecodeX v1.0 crosswalk. Disease-level matches and proxy manifestations are labeled separately; primary-only is the default. Five gene-disease pairs currently have no selected crosswalk. Review `data/mapping_coverage.csv` and the mapping file before analysis. Do not interpret proxy outcomes as diagnosis of the inherited syndrome.

In the Workbench, use the PhecodeX ICD-CM crosswalk to make subject-level calls from `condition_occurrence` and the source vocabulary fields. Preserve the PhecodeX version, source ICD vocabulary and observed EHR window. Distinguish participants with no matching code but observed EHR from participants without observable EHR; the latter must have missing phenotype status, not controls. For repeated diagnoses, collapse to one case status per participant and phenotype using a prespecified rule. Verify catalog-specific ICD-10-only and sex rules before defining controls.

Export three temporary CSVs inside the Workbench: `eligible.csv` (`subject_id,gene`, all retained callable samples per selected gene), `underliers.csv` (`subject_id,gene`, positive calls only), and `phenotypes.csv` (`subject_id,category,phenoname,status`, status 1/0/blank). Then run:

```sh
python -m pip install -e .
aou-enrichment --eligible results/eligible.csv \
  --underliers results/underliers.csv \
  --phenotypes results/phenotypes.csv \
  --selected-genes RBM20 LMNA \
  --output results/enrichment.csv
```

Replace the example genes with exactly the selected genes. The output contains per-gene/per-phenotype 2×2 counts, case fractions, odds ratios, exact Fisher p-values and Benjamini–Hochberg q-values. The default one-sided test asks whether phenotype cases are enriched among underliers. This is an unadjusted association screen; interpret it with the clinical specificity of each mapped phenotype and the blood/tissue relevance of the expression measurement in mind.

## Provenance

- ACMG SF v3.3 reference: ClinGen secondary findings resource, `data/sources/clingen_acmg_sf_v3_3.csv`.
- All of Us phenotype catalog: user-provided v7 workbook. It describes phenotype names and analysis availability; it is not participant-level phenotype data and does not determine an individual's case status.
- ICD-to-PhecodeX definitions: PheWAS PhecodeX v1.0 source files in `data/sources/`.

The reference builder is `scripts/build_references.py`. Public references can be regenerated from the ClinGen and All of Us workbook exports. Inspect its mapping rules before regenerating.
