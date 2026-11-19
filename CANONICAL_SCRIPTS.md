# Canonical scripts to regenerate figures & tables

## Orchestrators

| Script | Role |
|--------|------|
| `scripts/run_pipeline_pc.sh` | Full PC_SETTING rebuild (derive → figures → tables → verify) |
| `scripts/assemble_nature_figures.R` | Nature-tagged composites (Figs 1–4, S1) |
| `scripts/helpers/_paths.R` / `_paths.py` | Shared path + `PC_SETTING` resolution |
| `scripts/helpers/paper_theme.R` / `paper_mpl.py` | Shared typography / Okabe–Ito |

## Main figures (pipeline `fig1`–`fig4`)

| Stage | Scripts |
|-------|---------|
| fig1 | `fig1_overview/plot_underlier_tallies_topmed_aou.R` + assemble `ONLY_FIG=1` |
| fig2 | `fig2_constraint/plot_loeuf_calibration.R` (defs haplo/z_-2/z_-3/z_-4), `plot_underlier_tallies_loeuf_constrained_tolerant.R`, `analysis/assess_wilson_vs_subject_variance.R` + assemble |
| fig3 | `plot_ptv_vs_underlier_gene_level_2x2.R`, `plot_ptv_vs_underlier_gene_level.R`, `plot_ptv_vs_nonptv_zscore.R`, `plot_candidate_gene_zscore_dist_clingen_z3.R` + assemble |
| fig4 | `plot_shared_rare_outlier_genes_4dot.R`, `plot_acmg_4dot_no_ci.R` + assemble |

## PC swap / derive

| Script | Role |
|--------|------|
| `pc_swap/build_underlier_artifacts.R` | z-longform → underlier RDS / tallies |
| `pc_swap/build_topmed_vs_aou_prevalence.R` | TOPMed vs AoU prevalence table |
| `pc_swap/write_z_longform_from_flat.R` | Optional flat→long helper |
| `pc_swap/copy_pc_independent.sh` | Carry PC-invariant panels |
| `pc_swap/build_source_data.R` | Nature source-data tables |
| `pc_swap/verify_outputs.R` / `pc55_integrity.sh` | Checks |

## Supplement & tables

| Script | Role |
|--------|------|
| `supplement/plot_enrichment_npc_sv_vs_ptv.R` | S1a/b |
| `supplement/plot_ptv_consequence_summary.R` | S1c / S10 |
| `supplement/make_supplementary_s2_s5.R` | S2–S6 family |
| `supplement/plot_s1_celltype_overlay_s7.R` | S7 |
| `supplement/plot_s8_gene_z_vs_celltype.R` | S8 |
| `supplement/plot_s9_participant_characteristics.R` | S9 |
| `supplement/organize_canonical_supp_figs.py` | Canonical S11–S16 placement |
| `supplement/pack_supplementary_tables.py` | Supplementary tables pack |
| `supplement/sweep_enrichment_*.R` | PC-sweep enrichment (inputs to S1) |

## NMD (Fig3b-related)

| Path | Role |
|------|------|
| `nmd_escape_vep/VEP_NMD_COMMAND.sh` | Offline VEP + NMD plugin |
| `nmd_escape_vep/fig3b_ptv_nmd_escape_by_variant.tsv` | Per-variant status (external or `--with-nmd-tables`) |
| `fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level_2x2.R` | Consumes NMD map for gene colouring |

## Included but not always on the critical path

Older experiments / alternate cutoffs live under the same `fig*` folders
(e.g. `*.pre_legend_change.R`, `*_pc2005.R`, DeGoede comparisons). Keep for
provenance; the pipeline stages above define the paper path.

## Analysis precursors

`scripts/analysis/run_underlier_pipeline*.R`, `outlier_processing.R`,
`export_underlier_constrained_tolerant_tables.R` — earlier / all-genes
pipelines that feed some tables; prefer `pc_swap/build_underlier_artifacts.R`
for PC-swapped deliveries.
