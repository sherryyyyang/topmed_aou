# Path configuration

All R/Python entrypoints resolve paths through
`scripts/helpers/_paths.R` and `scripts/helpers/_paths.py`.

Cluster defaults remain the original absolute paths when env vars are
**unset**, so existing runs keep working without changes.

## Environment variables

| Variable | Purpose |
|----------|---------|
| `TOPMED_OUTLIERS_ROOT` / `ANALYSIS_ROOT` | Analysis tree root (`pre_Watershed_scripts` layout: `data/`, `fig*/`, `scripts/`, …) |
| `TOPMED_SHEYANG` | Parent directory holding `ptv_all.tsv` and `annotation/` |
| `TOPMED_SUBMISSION_PACKAGE` | Delivery folder for non-55-PC outputs (figures, manuscript swaps) |
| `TOPMED_PTV_TSV` | Absolute path to PTV table (default `$TOPMED_SHEYANG/ptv_all.tsv`) |
| `TOPMED_ACMG_CSV` | ACMG SF list CSV |
| `PC_SETTING` | Expression PCs residualized (default `55`; paper delivery `505`) |
| `PAPER_PYTHON` | Python 3 for matplotlib / docx / pptx helpers |
| `RSCRIPT` | Rscript binary (pipeline) |
| `DEFINITION_LABEL` | Underlier definition for LOEUF panels (`z_-3`, `z_-2`, …) |
| `ONLY_FIG` / `SKIP_PANELS` | Restrict / skip panel regen in `assemble_nature_figures.R` |
| `VEP_ENV` / `VEP_CACHE` / `NMD_OUT` | Offline VEP + NMD plugin (Fig3b annotations) |

## Recommended setup

```bash
cp config/paths.example.sh config/paths.sh
$EDITOR config/paths.sh
source config/paths.sh
```

`config/paths.sh` is gitignored (see `.gitignore`).

## Expected on-disk layout under `TOPMED_OUTLIERS_ROOT`

```
data/
  core/                 # kept_subjects.rds, loeuf_lookup.rds, metadata, …
  expression/           # long-form z TSV(s), expr_z_join.rds, …  (LARGE)
  underliers/           # PC=55 calls
  underliers_pc<N>/     # non-55 calls
  enrichment/ references/
call-RunRnaUnderlier/rna_outputs/   # AoU RNA tables
fig1_overview/ … fig4_acmg/ supplement/ paper/   # panel outputs (or pc<N>_work/)
nmd_escape_vep/         # NMD annotations (optional clone of scripts + TSV)
```

## What must stay external (do not git-commit)

| Asset | Typical size | Why |
|-------|-------------:|-----|
| `data/expression/*` long-form z / matrices | tens of GB | individual-level |
| `data/underliers*/` RDS | GB-scale | derived individual-level |
| `nmd_escape_vep/vep_cache/` | multi-GB | Ensembl cache |
| `nmd_escape_vep/fig3b_ptv_vep_nmd.txt` | ~70 MB | raw VEP dump |
| Pipeline logs, intermediate PNG/PDF panels | varies | regenerable |
| Manuscript `.docx` / large PPTX with embeds | MB–GB | not code |

Optional small exports (`--with-nmd-tables`, `--with-fig1-assets`) are
documented in the export script; review before pushing.
