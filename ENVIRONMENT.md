# Software environment

## R

Suggested packages (as used by panel / assemble scripts):

- tidyverse (`dplyr`, `readr`, `tidyr`, `ggplot2`, `stringr`, `purrr`)
- `patchwork`, `cowplot`, `grid`, `png`
- `data.table` (some burden / enrichment scripts)
- `scales`

Run with a recent R (≥ 4.2 recommended).

```r
install.packages(c(
  "tidyverse", "patchwork", "cowplot", "png", "data.table", "scales"
))
```

## Python

Used for study-design / scatter panels, docx/pptx swaps, and table packing:

- Python ≥ 3.9
- `matplotlib`, `numpy`, `pandas`
- `python-docx`, `python-pptx` (manuscript / slide helpers)

```bash
pip install matplotlib numpy pandas python-docx python-pptx
# or use your cluster pixi/conda env and set PAPER_PYTHON
```

## Optional: Ensembl VEP + NMD plugin

Only needed to **regenerate** Fig3b NMD annotations (not needed if you ship
`fig3b_ptv_nmd_escape_by_variant.tsv`). See `nmd_escape_vep/README_NMD_escape.md`.
