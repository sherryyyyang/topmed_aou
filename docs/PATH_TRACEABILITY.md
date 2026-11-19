# Path traceability — submission package (Nov 2026)

This folder is the **canonical delivery set** for the 2005-expression-PC paper figures.
Scripts to rebuild live in `scripts/` here. **Large raw / derived data stay** under
`pre_Watershed_scripts/` and `sheyang/` (not duplicated). Absolute paths below are the
authoritative inputs.

**Package root**

`/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026`

**Data / pipeline root (ROOT)**

`/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts`

**Sheyang root**

`/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang`

Rebuild (2005-PC figures into this package):

```bash
cd /net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026
PC_SETTING=2005 Rscript scripts/assemble_nature_figures.R
PC_SETTING=2005 python3 scripts/supplement/organize_canonical_supp_figs.py
```

---

## Layer 0 — raw / primary inputs

| Asset | Absolute path | Role |
|---|---|---|
| Residualized expression z (long form, 2005 PCs) | `ROOT/data/expression/expr_scaled_genotype_pc_5_expression_pc_2005_longer_form.tsv` | Gene×subject z after 5 genotype + 2005 expression PCs |
| Adjusted CPM | `ROOT/data/expression/expr_cpm_adj.rds` | Expression level filter (underlier: log2CPM &lt; gene_mean−1) |
| Kept subjects | `ROOT/data/core/kept_subjects.rds` | Analysis cohort (n≈5489) |
| Metadata | `ROOT/data/core/metadata.filtered_with_sv.txt` | Cohort / sequencing covariates |
| Gene symbols (GENCODE v30) | `ROOT/data/core/gencode_v30_gene_symbols.tsv` | ID ↔ symbol |
| LOEUF lookup | `ROOT/data/core/loeuf_lookup.rds` | Constraint |
| ClinGen HI | `ROOT/data/references/ClinGen_HI_Sufficient_evidence.csv` | Fig3c examples |
| TOPMed PTVs | `SHEYANG/ptv_all.tsv` | Fig3 / enrichment |
| ACMG SF v3.2 | `SHEYANG/annotation/variant_annotations/acmg_sf_v3.2.csv` | Fig4 gene set |
| Genotype GDS / SNV | `SHEYANG/gds/variants/` | Variant genotypes (live tree; large) |
| AoU RNA underlier tables | `ROOT/call-RunRnaUnderlier/rna_outputs/` | Fig1 AoU overlays |
| Cell-type / PC refs | `ROOT/pc_optimization/` | Supp S7–S8 |

Underlier definition used throughout: **z &lt; −3 AND log2CPM &lt; gene_mean − 1**.

---

## Layer 1 — derived underlier / enrichment artifacts (2005-PC)

| Asset | Absolute path | Built by |
|---|---|---|
| Underlier calls (z−3 clean) | `ROOT/data/underliers_pc2005/z_score_underlier_z_-3_clean.rds` | `scripts/pc_swap/build_underlier_artifacts.R` / underlier pipeline |
| Underlier tallies / pipeline_out | `ROOT/data/underliers_pc2005/pipeline_out/` | same |
| expr_z_join (2005) | `ROOT/data/underliers_pc2005/expr_z_join.rds` | PC-swap join |
| Enrichment RDS | `ROOT/data/enrichment/` | `scripts/supplement/sweep_enrichment_*.R` |
| Panel workdir outputs | `ROOT/pc2005_work/fig{1–4}_*/` and `ROOT/pc2005_work/supplement/` | Individual panel scripts |

Path constants: `scripts/helpers/_paths.R` and `_paths.py` (`PC_SETTING=2005`).

---

## Layer 2 — panel scripts → intermediate PNGs/PDFs

### Main Figure 1

| Panel | Script | Key inputs | Intermediate / panel export |
|---|---|---|---|
| a cohort overview | `scripts/fig1_overview/make_fig1ab_cohort_overviews.R` | metadata, kept subjects, AoU counts | `source_panels/Fig1a_cohort_overview_combined.png` |
| b underlier definition demo | `scripts/fig1_overview/plot_underlier_definition_demo.R` | illustrative z / CPM | `source_panels/Fig1b_underlier_definition_demo.{png,pdf}` |
| c–d gene/subject ranks | `scripts/fig1_overview/plot_underlier_tallies_topmed_aou.R` | underlier tallies; AoU gene prevalence | `Fig1c_*`, `Fig1d_*`; also `pc2005_work/fig1_overview/panel_b_*` |
| e TOPMed vs AoU prevalence | `scripts/fig1_overview/plot_topmed_vs_aou_prevalence_uncolored.R` | per-gene prevalence TSVs | `Fig1e_topmed_vs_aou_uncolored.*` |

### Main Figure 2

| Panel | Script | Key inputs |
|---|---|---|
| a–b LOEUF calibration | `scripts/fig2_constraint/plot_loeuf_calibration.R` | underliers + LOEUF (`DEFINITION_LABEL=z_-3`) |
| c constrained vs tolerant | `scripts/fig2_constraint/plot_underlier_tallies_loeuf_constrained_tolerant.R` | same |
| d TOPMed vs AoU × LOEUF | `scripts/fig2_constraint/plot_topmed_vs_aou_prevalence_loeuf.py` | prevalence + LOEUF |

### Main Figure 3

| Panel | Script | Key inputs |
|---|---|---|
| a PTV vs underlier (gene) | `scripts/fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level_2x2.R` | `ptv_all.tsv` + underliers |
| b PTV vs non-PTV z | `scripts/fig3_ptv_vs_rna/plot_ptv_vs_nonptv_zscore.R` | genotypes / PTV flags + z |
| c ClinGen histograms | `scripts/fig3_ptv_vs_rna/plot_candidate_gene_zscore_dist_clingen_z3.R` | ClinGen HI + z |

### Main Figure 4

| Panel | Script | Key inputs |
|---|---|---|
| a–d RBM20 / MYH11 | `scripts/fig4_acmg/plot_acmg_4dot_no_ci.R` | ACMG set, underliers, gnomAD AF helpers, rank curves |

### Supplementary (assembled or dedicated)

| Fig | Script(s) | Notes |
|---|---|---|
| S1 a\|b | `scripts/supplement/plot_enrichment_npc_sv_vs_ptv.R` | SV \| PTV enrichment vs n PCs |
| S2–S5, S6–S9 | `scripts/supplement/make_supplementary_s2_s5.R`, `plot_s1_celltype_overlay_s7.R`, `plot_s8_gene_z_vs_celltype.R`, `plot_s9_participant_characteristics.R`, enrichment / LOEUF helpers | See `supp_figs/Supplementary_Figures_provenance.md` |
| S10 | `scripts/supplement/plot_ptv_consequence_summary.R` | PTV VEP consequences (also sourced in assemble) |
| S11–S12 | `scripts/fig3_ptv_vs_rna/run_hvg_allgenes.R` (+ related) | HVG composition / UMAP → promoted via `organize_canonical_supp_figs.py` |
| S13 | coverage exploratory under `ROOT/supplement/exploratory/` | ClinGen coverage tracks |
| S14 | `scripts/supplement/plot_pc55_vs_pc2005_outlier_overlap.R` | PC55 vs PC2005 overlap |
| S15 | `scripts/fig2_constraint/plot_topmed_vs_aou_prevalence_clinical_class.py` | Clinical dosage class |
| S16 | `scripts/fig3_ptv_vs_rna/plot_fig3a_constraint_sanity_allgenes_pc2005.R` | All-genes constraint sanity |

---

## Layer 3 — composites written into this package

| Output | Built by | Location in package |
|---|---|---|
| Figure1–4 | `scripts/assemble_nature_figures.R` | `figures/main/Figure{1–4}.{png,pdf}` |
| S1, S10 | same assemble | `figures/supplementary/Supplementary_Figure{1,10}.*` |
| S11–S16 sync | `scripts/supplement/organize_canonical_supp_figs.py` | `figures/supplementary/Supplementary_Figure{11–16}.*` |
| Per-panel exports | assemble + panel scripts | `./source_panels/` |
| LaTeX-linked copies | symlinks | `./latex/figures/` |
| Index | | `./supp_figs/CANONICAL_LATEST.md` |

Old draft directory names were **not** kept as symlinks; use this package path only.

---

## What is *not* in this package (by design)

| Location | Why |
|---|---|
| `ROOT/data/expression/` (~tens of GB) | Raw/long-form matrices |
| `SHEYANG/gds/variants/` | Genotype GDS |
| `SHEYANG/archive/` | Archived intermediates / exploratory |
| `SHEYANG/writing_and_documentations/` | Older Aug-2026 prose + script snapshot (not this delivery) |
| `ROOT/final_paper_draft_pc55/` | Published 55-PC draft (separate) |

See `CHANGELOG_reorganization_nov_2026.md` for the move/archive history.
