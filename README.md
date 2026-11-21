# TOPMed / AoU RNA underlier figures

Scripts that go from raw counts and genotypes to the figure files in `figures/`. The scripts are the copies from `submission_package_nov_2026/scripts` that assembled the panels now in `figures/`. They still point at cluster paths. Large inputs (counts, VCFs, RDS, `ptv_all.tsv`) are not in this repository.

Cluster roots used by the scripts:

| Name | Path |
|------|------|
| Analysis root (`ROOT`) | `/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts` |
| Submission package | `/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026` |
| Parent (`SHEYANG`) | `/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang` |

`scripts/helpers/_paths.R` is the path switch. Several panel scripts assign `ROOT` to the analysis-root path above *before* sourcing `_paths.R`, so `TOPMED_OUTLIERS_ROOT` is ignored in those scripts. To run off this cluster, edit those `ROOT` assignments and the defaults in `_paths.R`, `scripts/run_pipeline_pc.sh`, `scripts/assemble_nature_figures.R`, and the two `export_fig3b_*.R` wrappers.

## PC setting for the figures in this repo

**The files in `figures/` were produced with `PC_SETTING=1055`** (5 genotype PCs + 1055 expression PCs on the *z* arm).

The script default is different. `scripts/helpers/_paths.R` and `scripts/run_pipeline_pc.sh` default to **`PC_SETTING=55`**. A run that does not set the variable will not regenerate these figures. Set it explicitly:

```bash
export PC_SETTING=1055
export TOPMED_SUBMISSION_PACKAGE=/path/to/this/clone   # see "Where scripts are read from"
```

What on disk shows 1055:

- `pipeline_pc1055_20261119_021923.log` starts with `PC_SETTING=1055`, derives `data/underliers_pc1055/`, and records `wrote Figure1`, `wrote Figure2`, and `wrote Figure3`. The follow-up log `pipeline_pc1055_20261119_023053.log` (also `PC_SETTING=1055`) rewrote Figure 4 and the supplementary stage. Both logs live in the submission package `logs/` directory and are not part of this repo.
- Figure 1 and Figure 2 panels were refreshed on 2026-11-20 under `pre_Watershed_scripts/pc1055_work/` and copied into `figures/main` the same morning.
- Figure 3 was reassembled on 2026-11-20 from panels written that evening under `pc1055_work/fig3_ptv_vs_rna/` (`panel_a` 19:37, `panel_b` 19:39, `figures/main/Figure3.png` 19:40).
- `export_fig3b_no_sv_color.R` and `export_fig3b_excl_sv_cnv.R` default `PC_SETTING` to **1055** and write `figures/main/Figure3b_no_sv_color.*` and `Figure3b_excl_sv_cnv.*`.
- Supplementary Figure 1 marks a dashed guide at 1055 (`plot_enrichment_npc_sv_vs_ptv.R`). Supplementary Figure 6 keeps the enrichment-OR row where `npcs == PC_SETTING` (`make_supplementary_s2_s5.R`).

The log2CPM arm of the underlier call does **not** move with `PC_SETTING`. `scripts/pc_swap/build_underlier_artifacts.R` reads `data/expression/expr_cpm_adj.rds`, which was residualized on 55 expression PCs. Only the *z* arm uses the 1055 long-form table.

## Underlier definition

From `scripts/pc_swap/build_underlier_artifacts.R` and `scripts/helpers/_underlier_def.R`:

- Primary underlier: expression *z* < −3 **and** log2CPM_adj < gene mean(log2CPM_adj) − 1.
- Haplo: log2CPM arm only (unchanged across PC settings).

## Raw inputs (on the cluster, not in git)

| Input | Where the scripts read it |
|-------|---------------------------|
| Gene counts | `scripts/analysis/outlier_processing.R` reads `../counts/gene_counts.csv` relative to the working directory. `scripts/analysis/run_underlier_pipeline.R` CONFIG sets `counts_csv = "../counts/gene_counts.csv"` and `gencode_gff = "../annotation/gencode.v30.annotation.gff3"`. |
| PLINK genotypes | `run_underlier_pipeline.R`: `/net/topmed11/working/porchard/for-others/julie-and-esther/2025-06-11/data/genotypes/plink` (`chr{1..22,X}.bed/bim/fam`). |
| Multi-NPC residualized matrices | `data/expression/expr_scaled_flat_list.rds` under `ROOT` (~28 GB). `scripts/pc_swap/write_z_longform_from_flat.R` writes `data/expression/expr_scaled_genotype_pc_5_expression_pc_<PC_SETTING>_longer_form.tsv`. |
| Adjusted log2CPM | `ROOT/data/expression/expr_cpm_adj.rds` |
| Kept subjects, LOEUF, metadata | `ROOT/data/core/kept_subjects.rds`, `loeuf_lookup.rds`, `metadata.filtered_with_sv.txt`, `gene_median_tpm_kept.tsv` |
| PTV table | `SHEYANG/ptv_all.tsv` (override `TOPMED_PTV_TSV` or `PTV_TSV`). Every row already has `af < 0.005` (maximum AF in the file is 0.00498). |
| PASS SV VCFs | `sv_processing.sh` lists `/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/data/vcfs-pass-variants-only/`. `build_fig3b_segdup_svtype_genes.R` reads that same directory (`dirname(SHEYANG)/data/vcfs-pass-variants-only`). |
| SegDup track | `build_cnv_multicopy_gene_set.R` downloads UCSC `genomicSuperDups.txt.gz` if missing and overlaps GENCODE v30 (`SHEYANG/annotation/gencode.v30.annotation.gff3`). |
| AoU prevalence | External release under `ROOT/call-RunRnaUnderlier/rna_outputs`, joined by `scripts/pc_swap/build_topmed_vs_aou_prevalence.R`. |
| ACMG SF v3.2 | `SHEYANG/annotation/variant_annotations/acmg_sf_v3.2.csv` |
| ClinGen HI | `ROOT/data/references/ClinGen_HI_Sufficient_evidence.csv` |

`scripts/analysis/processing.sh` is an earlier 1KGP SV filter (csh `>!` redirects). It is not the path that feeds the figures in `figures/`.

## From raw files to `figures/`

Run from the analysis root with `PC_SETTING=1055`. `scripts/run_pipeline_pc.sh` changes directory to `ROOT` and calls `$ROOT/scripts/...`. On the cluster those are the `pre_Watershed_scripts` copies. **Figure 3’s gene-level panel in this repo is newer than that tree:** `assemble_nature_figures.R` prefers `TOPMED_SUBMISSION_PACKAGE/scripts` when that directory exists, and that is the copy exported here. To redraw Figure 3 as in `figures/main/Figure3.png`, point `TOPMED_SUBMISSION_PACKAGE` at a checkout of these scripts (on the cluster, the submission package).

### 1. Residualize expression

`scripts/analysis/run_underlier_pipeline.R` is the one-shot counts → *z* pipeline. Its CONFIG default is `n_expr_pcs = 55` and `n_geno_pcs = 5`. The 1055-PC *z* matrix used for these figures is taken from the existing flat list, not by re-running that CONFIG:

```bash
PC_SETTING=1055 Rscript scripts/pc_swap/write_z_longform_from_flat.R
```

`scripts/analysis/outlier_processing.R` is the original notebook that function was distilled from (counts via `vroom("../counts/gene_counts.csv")`, edgeR CPM, PCAtools, SNPRelate genotype PCs, `limma::removeBatchEffect`).

### 2. Call underliers and cross-cohort prevalence

```bash
PC_SETTING=1055 Rscript scripts/pc_swap/build_underlier_artifacts.R
PC_SETTING=1055 Rscript scripts/pc_swap/build_topmed_vs_aou_prevalence.R
```

Writes `ROOT/data/underliers_pc1055/` (`expr_z_join.rds`, `z_score_underlier_*_clean.rds`, tallies). At `PC_SETTING=55` the same script must reproduce `data/underliers/` exactly; `run_pipeline_pc.sh` treats that as a gate and will not overwrite the 55-PC tree unless `ALLOW_PC55_WRITE=1`.

### 3. Fig. 3b dosage tables (before the scatter)

```bash
Rscript scripts/fig3_ptv_vs_rna/build_cnv_multicopy_gene_set.R
Rscript scripts/fig3_ptv_vs_rna/build_fig3b_segdup_svtype_genes.R
```

Outputs under `ROOT/data/references/`: `fig3b_cnv_multicopy_genes.tsv`, `fig3b_segdup_svtype_genes.tsv`. A gene is coloured only if it passes the SegDup rule (`segdup_ge50pct_fracMatch_ge0.95`) **and** has a gene-overlapping PASS muCNV. The primary call is the overlapping DEL, DUP, or INV with the largest `|SVLEN|` (any AF). Curated-only genes and SegDup genes with no such VCF hit stay in “Other genes”.

### 4. Draw panels and assemble

`scripts/run_pipeline_pc.sh` stages, in order: `preflight`, `gate`, `derive`, `seed`, `carryover`, `prevalence`, `fig1`, `fig2`, `fig3`, `fig4`, `supplement`, `tables`, `sourcedata`, `manuscript`, `prosedocs`, `readme`, `verify`.

```bash
PC_SETTING=1055 ./scripts/run_pipeline_pc.sh
# or one figure:
PC_SETTING=1055 ONLY_FIG=3 Rscript scripts/assemble_nature_figures.R
```

`assemble_nature_figures.R` writes composites to `figures/main/` and `figures/supplementary/` (and also to `paper/composites/` and `source_panels/` on the package). Panel scripts:

| Figure | Script | What it draws |
|--------|--------|----------------|
| 1 | `fig1_overview/make_fig1ab_cohort_overviews.R`, `plot_underlier_definition_demo.R`, `plot_underlier_tallies_topmed_aou.R`, `plot_topmed_vs_aou_prevalence_uncolored.R` | Cohort schematic, underlier definition, gene/subject ranks, TOPMed vs AoU prevalence |
| 2 | `fig2_constraint/plot_loeuf_calibration.R`, `plot_underlier_tallies_loeuf_constrained_tolerant.R`, `plot_topmed_vs_aou_prevalence_loeuf.R` | LOEUF forest and deciles, constrained vs tolerant ranks, prevalence × LOEUF |
| 3a | `fig3_ptv_vs_rna/plot_ptv_vs_nonptv_zscore.R` | PTV-carrier vs non-carrier *z* (top row of Figure 3) |
| 3b | `fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level_2x2.R` | Gene-level PTV % vs underlier % |
| 3c | `fig3_ptv_vs_rna/plot_candidate_gene_zscore_dist_clingen_z3.R` | ClinGen example histograms |
| 4 | `fig4_acmg/plot_shared_rare_outlier_genes_4dot.R`, `plot_acmg_4dot_no_ci.R` | RBM20 and MYH11 rank curves and 4-dot strips. Focus genes are force-included. |

Figure 3 row order in `assemble_nature_figures.R` is *z*-score (`f3b`), then gene-level scatter (`f3a`), then ClinGen (`f3c`).

### 5. Supplementary figures

| Output in `figures/supplementary/` | Script | PC |
|------------------------------------|--------|----|
| Supplementary Figure 1 (SV vs PTV enrichment across NPC; guide at 1055), `enrichment_npc_sv_vs_ptv.*`, no-overlay copies | `supplement/plot_enrichment_npc_sv_vs_ptv.R` | Sweep across NPC; not a single-PC recompute |
| Supplementary Figures 4, 5, 6 and `FigS4`–`FigS6` | `supplement/make_supplementary_s2_s5.R` | S4–S6 use `PC_SETTING`. S6 filters `npcs == PC_SETTING`. |
| Supplementary Figures 2 and 3 | same script | Skipped when `Supplementary_Figure2.png` and `Supplementary_Figure3.png` already exist (`SKIP_S23`). The copies here are the 2026-09-18 QC / UMAP figures (connectivity *z*, PCA scree, PERMANOVA, sequencing centre, cohort). They do not use the 1055 *z* table. |
| Supplementary Figure 7, `FigS7*`, `FigS_celltype_pc_eval` | `supplement/plot_s1_celltype_overlay_s7.R` | Cell-type vs PC composition (PC1–80), not the 1055 underlier calls |
| Supplementary Figure 8, `FigS8*`, and the S1 a\|b reassembly plus Supplementary Figure 10 | `supplement/plot_s8_gene_z_vs_celltype.R` | Gene *z* vs cell type uses the PC setting |
| Supplementary Figure 9, `FigS9*` | `supplement/plot_s9_participant_characteristics.R` | Demographics; no TOPMed *z* |
| Supplementary Figures 11–16 | `supplement/organize_canonical_supp_figs.py` | Copies earlier renders (HVG composition, UMAP, coverage tracks, PC55 vs PC2005 overlap, clinical-class prevalence, all-genes constraint sanity). Not recomputed at 1055. |

`supplement/plot_ptv_consequence_summary.R` draws the PTV-consequence panel used as Supplementary Figure 10 / `FigS10a`.

## Figure 3 details (as coded)

**3a** (`plot_ptv_vs_nonptv_zscore.R`). The script caches both AF < 0.5% (`af < 0.005`) and AF < 0.1% (`af < 0.001`) inside PTV genes with LOFTEE `lof == "HC"`. The paper row is AF < 0.1% only:

```r
fig <- (rows[["AF<0.1%"]]$pA | rows[["AF<0.1%"]]$pB | rows[["AF<0.1%"]]$pC)
```

The right-hand panel (`pC`, title “Outlier rate (log % axis)”) plots percent of pairs on `scale_y_log10`. Values are drawn with `pmax(pct, 0.005)` so *z* < −3 and *z* < −4 bars stay visible. Fold annotations still use the unfloored percents.

**3b** (`plot_ptv_vs_underlier_gene_level_2x2.R`). One point per gene. Both axes are `scale_*_log10`. Zeros are placed at half the smallest positive percent and labelled “0%”. Spearman ρ is computed on genes with `n_ptv_carriers > 0` (`stats::cor(..., method = "spearman")`).

- *x* = percent of kept subjects carrying a PTV in the gene. The loader keeps `ptv_all.tsv` rows for kept subjects and does not apply a second `af <` cut. That file is already AF < 0.5%, which matches the axis title.
- *y* = percent of kept subjects with an RNA underlier in the gene (primary definition, *z* < −3 and the log2CPM drop).
- Colour (`FIG3B_COLOR` unset or `cnv`): SegDup-pass primary PASS SVTYPE — DEL / DUP / INV — versus other genes.
- Point size: number of distinct gene-overlapping PASS SVs with at least one alt carrier among that gene’s *z* < −3 RNA underliers. The scale caps the aesthetic at 5 (`n_sv_size = pmin(n_sv_in_underliers, 5)`, breaks 0, 1, 2, 5+).

**3c** (`plot_candidate_gene_zscore_dist_clingen_z3.R`). Candidate genes are filtered with `af < 0.001` (`PTV_AF_MAX <- 0.001`).

### Alternate Figure 3b exports

These do not replace `figures/main/Figure3.png`. Flags are read by `plot_ptv_vs_underlier_gene_level_2x2.R`:

| Flag | Wrapper | Output |
|------|---------|--------|
| `NO_SV_COLOR=1` | `scripts/fig3_ptv_vs_rna/export_fig3b_no_sv_color.R` | `figures/main/Figure3b_no_sv_color.png` and `.pdf` — all genes, grey, one size, log–log, Spearman |
| `EXCLUDE_SV_CNV=1` | `scripts/fig3_ptv_vs_rna/export_fig3b_excl_sv_cnv.R` | `figures/main/Figure3b_excl_sv_cnv.png` and `.pdf` — drops SegDup DEL/DUP/INV genes and `fig3b_cnv_multicopy_genes.tsv`, then replots in grey and recomputes Spearman |
| `FIG3B_COLOR=nmd` | no wrapper; set the variable when sourcing the panel script | NMD-escape colouring from `nmd_escape_vep/`. Not the colouring in `Figure3.png`. |
| `FIG3B_OTHER_ONLY=1` | `plot_fig3b_no_overlap_test.R` / the panel script | Test output only. Those PNGs are not in `figures/`. |

`nmd_escape_vep/VEP_NMD_COMMAND.sh` runs Ensembl VEP offline with the NMD plugin. The input VCF and VEP text are on the cluster and are not in this repo. `nmd_escape_vep/README_NMD_escape.md` records that annotation; its ECDF note is a PC 2005 summary, not Figure 3.

## What is in `figures/`

- `figures/main/Figure{1,2,3,4}.{png,pdf}` — assembled main figures.
- `figures/main/Figure3b_no_sv_color.{png,pdf}` and `Figure3b_excl_sv_cnv.{png,pdf}` — the two alternate 3b exports.
- `figures/supplementary/` — supplementary composites (`Supplementary_Figure1`–`16`, including `Supplementary_Figure1_no_overlay`) and the panel files (`FigS*`, `enrichment_npc_sv_vs_ptv`).

Test renders (`Figure3b_*_test`, `figures/pc55_scratch/`) are not included. Manuscript docx/pptx files are not included.

## Scripts in the tree that do not write `figures/`

Exploratory or earlier-draft scripts are still here so the directory matches the package, including all-genes / lncRNA Fig. 3a checks, `test_high_underlier_low_ptv_sv*.R`, `plot_fig3b_prevalence_boxplot_test.R`, and `plot_loeuf_calibration.original_from_transcript.R` / `.pre_legend_change.R`. They are not the commands that produced `figures/main/Figure{1–4}`.
