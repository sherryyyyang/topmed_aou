> **PROVENANCE — 505-EXPRESSION-PC DELIVERY (z<−3 AoU-style 95% plateau).**
> Package figures are regenerated at **PC_SETTING=505**. Rule: earliest NPC with OR ≥ 95% of max
> on each Supp Fig1 z<−3 curve (SV→255, PTV→505); paper PC = max = **505**. See `docs/PC_PLATEAU_z3.md`.
>
> Regenerate with: `PC_SETTING=505 ./scripts/run_pipeline_pc.sh`


# What is actually in the latest Word file

Live Word file (this folder):
`Outliers_in_TOPMed_paper_draft.docx`

Rebuild everything (Arial, Nature lowercase tags a/b/c, headlines in
captions not on the graphics, Okabe–Ito palette):

```
Rscript ../scripts/assemble_nature_figures.R
python3 ../scripts/helpers/update_docx_figures.py
```

`assemble_nature_figures.R` uses **patchwork** to tag ggplot panels and
**cowplot** to wrap Python PNG panels and to stack Supplementary Figure 1.
Shared type is `scripts/helpers/paper_theme.R` (R) and `paper_mpl.py`
(Python). Double-column width is 183 mm. `DEFINITION_LABEL=z_-3` is set
for Figure 2.

Composites (`Figure1.png` … `Supplementary_Figure1.png`) are written here,
to `source_panels/`, and to `source_panels/as_in_docx/`. Individual panels
use sequential lowercase names (`Fig1a_…`, `Fig2d_…`).

Nature source tables (Figs 1–4, S1, robustness, provenance):
`supplemental_data/` (`MANIFEST.tsv`, `PROVENANCE.tsv`).

Figure 1a is a short two-line banner (not a tall hero strip). Raster rows
are sized to the PNG aspect so they are not stretched or padded.

## Palette (Okabe–Ito)

| Role | Hex |
|---|---|
| TOPMed | `#0072B2` |
| AoU | `#CC79A7` |
| PTV / underlier highlight | `#D55E00` |
| Constrained → Tolerant | `#D55E00`, `#E69F00`, `#56B4E9`, `#0072B2` |
| No LOEUF | `#000000` |
| Sequential | viridis |

## Panel letters (sequential lowercase Nature tags)

| Figure | Panels |
|---|---|
| 1 | **a** study design; **b** gene ranks; **c** subject ranks; **d** TOPMed vs AoU (uncolored) |
| 2 | **a** LOEUF-tier forest; **b** LOEUF-decile trend; **c** constrained vs tolerant; **d** TOPMed vs AoU (LOEUF color) |
| 3 | **a** gene-level PTV vs underlier; **b** 2×3 carrier z-score grid (one tag); **c** ClinGen histograms |
| 4 | **a** RBM20 ranks; **b** RBM20 4-dot; **c** MYH11 ranks; **d** MYH11 4-dot |
| S1 | **a** SV enrichment vs n PCs; **b** PTV enrichment vs n PCs; **c** PTV consequences |

Python panels: `scripts/fig1_overview/plot_study_design_topmed_aou.py`,
`plot_topmed_vs_aou_prevalence_uncolored.py`,
`scripts/fig2_constraint/plot_topmed_vs_aou_prevalence_loeuf.py`.
