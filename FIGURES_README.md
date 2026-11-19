# TOPMed RNA-outlier figure & table scripts

Code to regenerate main Figures 1–4, Supplementary Figures, and supporting
tables for the TOPMed whole-blood RNA-outlier / underlier analyses.

This repository is **scripts + docs only**. Expression matrices, underlier
RDS objects, VEP caches, and rendered figure PNGs/PDFs stay on analysis
storage (see [PATHS.md](PATHS.md) and [data/README.md](data/README.md)).

## Layout

```
scripts/
  helpers/           # _paths.R/_paths.py, paper_theme, underlier defs
  analysis/          # underlier pipeline, prevalence helpers
  fig1_overview/ … fig4_acmg/
  supplement/        # S1–S16 panel + table packers
  pc_swap/           # PC_SETTING rebuild (derive, verify, source data)
  assemble_nature_figures.R
  run_pipeline_pc.sh
nmd_escape_vep/      # VEP+NMD docs/command (tables optional; see DATA_EXTERNAL.txt)
config/paths.example.sh
docs/                # figure sources, PC plateau, NMD notes
```

## Quick start (off-cluster)

```bash
cp config/paths.example.sh config/paths.sh
# edit TOPMED_OUTLIERS_ROOT, TOPMED_SHEYANG, PAPER_PYTHON, …
source config/paths.sh

# Assemble Nature composites from existing panel outputs:
Rscript scripts/assemble_nature_figures.R

# Or rebuild a non-55 PC delivery end-to-end (needs long-form z + data tree):
PC_SETTING=505 ./scripts/run_pipeline_pc.sh
PC_SETTING=505 STAGES=fig3,fig4 ./scripts/run_pipeline_pc.sh   # subset
./scripts/run_pipeline_pc.sh --list
```

Published analysis is **PC_SETTING=55**. The Nov 2026 submission delivery
uses the AoU-style z&lt;−3 enrichment plateau at **PC_SETTING=505**
(see `docs/PC_PLATEAU_z3.md`).

## Pipeline stages

`run_pipeline_pc.sh` stages (in order):

`preflight` → `gate` → `derive` → `seed` → `carryover` → `prevalence` →
`fig1` → `fig2` → `fig3` → `fig4` → `supplement` → `tables` →
`sourcedata` → `manuscript` → `prosedocs` → `readme` → `verify`

- **gate**: 55-PC round-trip must reproduce 48,632 z&lt;−3 underlier pairs.
- **fig\*** stages call panel scripts then `assemble_nature_figures.R` with `ONLY_FIG=N`.
- At PC_SETTING=55, render stages are refused unless `ALLOW_PC55_WRITE=1`
  (protects the checksummed published tree).

## Canonical panel scripts (main figures)

| Figure | Scripts sourced by `assemble_nature_figures.R` |
|--------|--------------------------------------------------|
| 1 | `fig1_overview/make_fig1ab_cohort_overviews.R`, `plot_underlier_definition_demo.R`, `plot_underlier_tallies_topmed_aou.R`, `plot_topmed_vs_aou_prevalence_uncolored.R` (+ Python study-design / scatter helpers) |
| 2 | `fig2_constraint/plot_loeuf_calibration.R`, `plot_underlier_tallies_loeuf_constrained_tolerant.R`, `plot_topmed_vs_aou_prevalence_loeuf.R` |
| 3 | `fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level_2x2.R`, `plot_ptv_vs_nonptv_zscore.R`, `plot_candidate_gene_zscore_dist_clingen_z3.R` |
| 4 | `fig4_acmg/plot_acmg_4dot_no_ci.R` (data via `plot_shared_rare_outlier_genes_4dot.R` in the pipeline) |
| S1 | `supplement/plot_enrichment_npc_sv_vs_ptv.R`, `plot_ptv_consequence_summary.R` |

Full inventory: [CANONICAL_SCRIPTS.md](CANONICAL_SCRIPTS.md).

## Fig. 3b / NMD escape

Gene-level NMD colouring (when enabled) reads
`nmd_escape_vep/fig3b_ptv_nmd_escape_by_variant.tsv` produced with Ensembl
VEP + NMD plugin. Re-export that table with
`--with-nmd-tables`, or keep it external. See `nmd_escape_vep/README_NMD_escape.md`
and `docs/README_Fig3b_NMD_escape.md`.

## Environment

R packages and Python deps: [ENVIRONMENT.md](ENVIRONMENT.md).

## License

See [LICENSE](LICENSE). Data access remains under TOPMed / dbGaP terms;
this repo does not redistribute individual-level data.
