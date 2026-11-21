# External annotations for Figures 1–4

An underlier is a gene–subject pair with z < −3 AND log2CPM_adj < gene_mean − 1. `PC_SETTING` (how many expression PCs were removed) stays an input from a future selection rule. It is not one of the annotations below.

`$DATA_ROOT` is the analysis directory (`TOPMED_OUTLIERS_ROOT`). Files the scripts place beside that directory are written as `dirname($DATA_ROOT)/...`. Only annotations that a Figure 1–4 script, or a script that builds a table those figures read, actually opens are listed. Download URLs are copied from those scripts or from a note stored next to the file. Where no command was recorded, that is stated and the filename is quoted.

## GENCODE v30

**Figures.** Figure 1 rank labels read `$DATA_ROOT/data/core/gencode_v30_gene_symbols.tsv` (no header: Ensembl gene id without a version, then symbol) in `scripts/fig1_overview/plot_underlier_tallies_topmed_aou.R`. Figure 3 SegDup overlap and the Figure 4 gnomAD overlap use the GFF3 gene models: `dirname($DATA_ROOT)/annotation/gencode.v30.annotation.gff3` (`build_cnv_multicopy_gene_set.R`) and `dirname($DATA_ROOT)/annotation/gencode.v30.annotation.gff3.gz` (`sum_af_exome_topmed_overlap.R`).

**How it was obtained.** Source not recorded in the scripts. The GFF3 header on the file those builders read is `#description: evidence-based annotation of the human genome (GRCh38), version 30 (Ensembl 96)`, `#provider: GENCODE`, `#date: 2019-03-28`. `plot_ptv_vs_underlier_gene_level.R` says the symbol table was built once with `zcat gencode.v30.annotation.gff3.gz | awk ... > gencode_v30_gene_symbols.tsv`; the `awk` body is not stored. `build_cnv_multicopy_gene_set.R` builds `gencode.v30.genes.bed` from the GFF3 with an `awk` that keeps `type=="gene"` and prints `gene_id` (version stripped), `gene_name`, and `gene_type`.

## gnomAD LOEUF

**Figures.** Figure 2. `plot_loeuf_calibration.R` and `plot_underlier_tallies_loeuf_constrained_tolerant.R` read `$DATA_ROOT/data/core/loeuf_lookup.rds` and bin genes as Constrained (LOEUF < 0.35), Mid-low (0.35 ≤ LOEUF < 0.6), Mid-high (0.6 ≤ LOEUF < 1.0), Tolerant (LOEUF ≥ 1.0), or No LOEUF. `build_topmed_vs_aou_prevalence.R` joins the same RDS onto `topmed_vs_new_z_-3_prevalence_per_gene.tsv`, which Figure 2d colors by `LOEUF_tier`. The Figure 3 gene-level scatter also joins the RDS, for the symbol column.

**How it was obtained.** Source not recorded in the scripts. The plot comments name `build_loeuf_lookup.R`, and that file is not in this repository. The figure reads `loeuf_lookup.rds` (columns `gene_id`, `gene_id_no_version`, `symbol`, `gene_type`, `LOEUF`, `LOEUF_tier`).

## PTV table (LOFTEE `lof`, allele frequency)

**Figures.** Figure 3 and Figure 4. Default path in `scripts/helpers/_paths.R` is `dirname($DATA_ROOT)/ptv_all.tsv` (`TOPMED_PTV_TSV`). Columns the panels use are `gene`, `sample`, `lof`, and `af`.

- `plot_ptv_vs_nonptv_zscore.R` keeps `lof == "HC"`. The drawn row is `af < 0.001` (the script also caches `af < 0.005`).
- `plot_candidate_gene_zscore_dist_clingen_z3.R` keeps `af < 0.001`. Its comment calls `af` the gnomAD allele fraction in `ptv_all.tsv`.
- `plot_ptv_vs_underlier_gene_level_2x2.R` (default `FIG3B_COLOR=cnv`) keeps every row whose `sample` is in `kept_subjects.rds`. It does not apply a second `af` cut. The SegDup / DEL / DUP / INV color is a separate table (below).
- `plot_shared_rare_outlier_genes_4dot.R` counts carriers from the same file.

The VEP `consequence` string is not a filter in Figures 1–4.

**How it was obtained.** Source not recorded in the scripts. The figures read `ptv_all.tsv`.

## UCSC genomicSuperDups (hg38)

**Figures.** Figure 3 gene-level scatter, as the SegDup arm of the DEL / DUP / INV color. `build_cnv_multicopy_gene_set.R` flags a gene when at least 50% of its GENCODE v30 body overlaps a SegDup interval with `fracMatch >= 0.95`. `build_fig3b_segdup_svtype_genes.R` keeps those genes (reason contains `segdup_ge50pct_fracMatch_ge0.95`) and joins them to PASS SVs. The scatter reads `$DATA_ROOT/data/references/fig3b_segdup_svtype_genes.tsv`, `fig3b_segdup_sv_gene_hits.tsv`, and `fig3b_segdup_sv_carriers.tsv`.

**How it was obtained.** If the gzip is missing, `build_cnv_multicopy_gene_set.R` runs:

```r
SEGDUP_URL <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"
download.file(SEGDUP_URL, SEGDUP_GZ, mode = "wb", quiet = FALSE)
```

Expected file: `$DATA_ROOT/data/references/genomicSuperDups_hg38.txt.gz`. The script does not pin a UCSC release date.

## muCNV PASS SV VCFs (DEL / DUP / INV)

**Figures.** Same Figure 3 color and point size. `build_fig3b_segdup_svtype_genes.R` queries PASS records with `SVTYPE` of DEL, DUP, or INV (`bcftools query` on `INFO/END`, `ID`, `SVTYPE`, `SVLEN`, `AF`), intersects them with SegDup-pass gene bodies, and picks the gene-overlapping variant with the largest `|SVLEN|`. Point size is the number of those SVs with an alternate carrier among that gene’s underliers.

**How it was obtained.** Source not recorded in the scripts. The builder reads `{chr}.vcf.gz` from `file.path(dirname(SHEYANG), "data", "vcfs-pass-variants-only")`, with `SHEYANG = dirname($DATA_ROOT)`, so the directory is `dirname(dirname($DATA_ROOT))/data/vcfs-pass-variants-only/`. The VCF header on `chr1.vcf.gz` is `##source=muCNV_pipeline_v0.9.6`.

## gnomAD v4.1 joint exome AF (`AF_exomes`)

**Figures.** Figure 4. `sum_af_exome_topmed_overlap.R` subsets `dirname($DATA_ROOT)/annotation/gnomad/gnomad.joint.v4.1.sites.exomes.vcf.bgz` to GENCODE v30 bodies of the shared genes (`bcftools view -R`), overlap-joins sites to genes, keeps sites that are also in `ptv_all.tsv`, and writes `sum(AF_exomes)` per gene to `sum_af_exome_topmed_overlap.tsv`. `plot_shared_rare_outlier_genes_4dot.R` copies that column into `shared_rare_outlier_genes_4dot.tsv` as the fourth measure, `gnomAD sum AF_exomes (TOPMed-overlapping PTVs)`. `plot_acmg_4dot_no_ci.R` reads that table. It sets `drop_gnomad <- TRUE`, so the drawn RBM20 and MYH11 strips omit the gnomAD point.

**How it was obtained.** Source not recorded in the scripts. The figure pipeline reads `gnomad.joint.v4.1.sites.exomes.vcf.bgz`. The VCF header has `##hailversion=0.2.130-bea04d9c79b5` and contig assembly `gnomAD_GRCh38`, and `INFO/AF_exomes` is the exome alternate allele frequency. No download URL is in `sum_af_exome_topmed_overlap.R`.

## ACMG SF v3.2 gene list

**Figures.** Figure 4 opens `acmg_sf_v3.2.csv` and sets `is_acmg` from the `Gene` column. That flag is not used to drop genes or to choose the panels. The two panels are the hardcoded symbols RBM20 and MYH11, and the ranks keep genes with `n_ptv > 0`. `_paths.R` default is `dirname($DATA_ROOT)/annotation/variant_annotations/acmg_sf_v3.2.csv` (`TOPMED_ACMG_CSV`). `plot_acmg_4dot_no_ci.R` then assigns `ACMG_CSV` itself.

**How it was obtained.** A `commands` note stored next to `acmg_sf_v3.2.csv` records the spreadsheet as supplementary Table 1 of <https://www.sciencedirect.com/science/article/pii/S1098360023008791?via%3Dihub#appsec1> (`acmg_sf_v3.2.xlsx`). The command that wrote the CSV is not in the scripts. The figure reads `acmg_sf_v3.2.csv`.

## ClinGen haploinsufficiency list

**Figures.** The assembled histogram, `plot_candidate_gene_zscore_dist_clingen_z3.R`, does not open a ClinGen file. It hardcodes eight genes: NOG (`ENSG00000183691`), CDKN1C (`ENSG00000129757`), FGFR1 (`ENSG00000077782`), SLC6A8 (`ENSG00000130821`), DDX3X (`ENSG00000215301`), RB1 (`ENSG00000139687`), PUF60 (`ENSG00000179950`), CHD7 (`ENSG00000171316`). `$DATA_ROOT/data/references/ClinGen_HI_Sufficient_evidence.csv` is the `_paths.R` object `CLINGEN_HI`. `assemble_nature_figures.R` does not call the script that reads it (`plot_ptv_vs_underlier_gene_level_2x2_clingen_hi.R`).

**How it was obtained.** Source not recorded in the scripts. The assembled panel does not read `ClinGen_HI_Sufficient_evidence.csv`.

## Checked and not used by Figures 1–4

**ClinVar.** No Figure 1–4 script reads a ClinVar VCF or uses a ClinVar label as a filter. The chromosome-wise overlap job is a separate analysis (`slurm_clinvar_x_cohort_full_array.sbatch`); it is not called by `assemble_nature_figures.R`. That job’s default VCF is `annotation/variant_annotations/clinvar.vcf.gz`. Its header is `##fileDate=2024-02-15`, `##source=ClinVar`, `##reference=GRCh38`. The download command is not in the scripts.

**Ancestry and cohort labels.** The assembled Figure 1 cohort panel copies a provided overview image (`make_fig1ab_cohort_overviews.R`) or, if that image is absent, draws hardcoded cohort counts in `recreate_publication_figure.py`. It does not read an ancestry table. TOPMed vs AoU on the prevalence panels is a column already present on `topmed_vs_new_z_-3_prevalence_per_gene.tsv`.
