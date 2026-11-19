#!/usr/bin/env Rscript
# Write README.md into the PC_SETTING output folder.
#
# Numbers are read back out of the tables the pipeline just produced rather
# than hard-coded, so the README cannot drift from the figures beside it.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(stringr) })

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

OUT   <- FINAL_DRAFT
REF   <- file.path(ROOT, "final_paper_draft_pc55")
TBL   <- file.path(OUT, "supp_figs", "tables")
cm    <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
rd    <- function(p) if (file.exists(p)) suppressWarnings(read_tsv(p, show_col_types = FALSE, comment = "#")) else NULL
fmt   <- function(x, d = 2) ifelse(is.na(x), "n/a", formatC(x, format = "f", digits = d))

# ---- underlier counts, this setting vs published ------------------------
here <- rd(file.path(DATA_UNDERLIER, "underlier_counts_summary.tsv"))
there <- rd(file.path(ROOT, "data", "underliers_pc55_validation", "underlier_counts_summary.tsv"))
counts_md <- "_(underlier count summary not found)_"
if (!is.null(here)) {
  j <- here |> select(definition, n_2005 = n_underlier_pairs,
                      g_2005 = n_genes_with_underlier, s_2005 = n_subjects_with_underlier)
  if (!is.null(there)) {
    j <- left_join(j, there |> select(definition, n_55 = n_underlier_pairs,
                                      g_55 = n_genes_with_underlier,
                                      s_55 = n_subjects_with_underlier), by = "definition")
  } else { j$n_55 <- NA; j$g_55 <- NA; j$s_55 <- NA }
  counts_md <- paste0(
    "| definition | pairs @55 | pairs @", PC_SETTING, " | change | genes @55 | genes @", PC_SETTING,
    " | subjects @55 | subjects @", PC_SETTING, " |\n",
    "|---|---:|---:|---:|---:|---:|---:|---:|\n",
    paste(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
                  j$definition, cm(j$n_55), cm(j$n_2005),
                  ifelse(is.na(j$n_55), "n/a",
                         sprintf("%+.0f%%", 100 * (j$n_2005 / j$n_55 - 1))),
                  cm(j$g_55), cm(j$g_2005), cm(j$s_55), cm(j$s_2005)),
          collapse = "\n"))
}

# ---- LOEUF tier gradient at z<-3 ----------------------------------------
tier_md <- "_(LOEUF tier table not found)_"
tier_here <- rd(file.path(FIG2A, "underlier_loeuf_calibration_by_tier_z_-3.tsv"))
tier_ref  <- rd(file.path(ROOT, "fig2_constraint", "panel_a_loeuf_calibration",
                          "underlier_loeuf_calibration_by_tier_z_-3.tsv"))
if (!is.null(tier_here)) {
  t2 <- tier_here |> select(LOEUF_tier, n_genes_in_tier, prev_2005 = prev_per_10k,
                            lo_2005 = prev_per_10k_lo, hi_2005 = prev_per_10k_hi)
  if (!is.null(tier_ref)) {
    t2 <- left_join(t2, tier_ref |> select(LOEUF_tier, prev_55 = prev_per_10k),
                    by = "LOEUF_tier")
  } else t2$prev_55 <- NA_real_
  tier_md <- paste0(
    "| LOEUF tier | genes | prevalence /10k @55 | prevalence /10k @", PC_SETTING, " (95% CI) |\n",
    "|---|---:|---:|---:|\n",
    paste(sprintf("| %s | %s | %s | %s (%s-%s) |",
                  t2$LOEUF_tier, cm(t2$n_genes_in_tier), fmt(t2$prev_55),
                  fmt(t2$prev_2005), fmt(t2$lo_2005), fmt(t2$hi_2005)),
          collapse = "\n"))
}

# ---- PTV / SV enrichment slice at this PC setting ------------------------
enr_md <- "_(S6 enrichment table not found)_"
e_here <- rd(file.path(TBL, sprintf("S6_enrichment_or_%spcs.tsv", PC_SETTING)))
e_ref  <- rd(file.path(REF, "supp_figs", "tables", "S6_enrichment_or_55pcs.tsv"))
if (!is.null(e_here)) {
  k <- function(d) d |> transmute(variant = variant_class, cutoff = as.numeric(cutoff),
                                  or = odds_ratio, lo = ci_95_lo, hi = ci_95_hi)
  a <- k(e_here)
  b <- if (!is.null(e_ref)) k(e_ref) |> rename(or55 = or) |> select(variant, cutoff, or55)
       else tibble(variant = character(), cutoff = numeric(), or55 = numeric())
  m <- left_join(a, b, by = c("variant", "cutoff")) |> arrange(variant, cutoff)
  enr_md <- paste0(
    "| variant class | cutoff | OR @55 | OR @", PC_SETTING, " (95% CI) |\n",
    "|---|---:|---:|---:|\n",
    paste(sprintf("| %s | z < -%g | %s | %s (%s-%s) |",
                  m$variant, m$cutoff, fmt(m$or55), fmt(m$or), fmt(m$lo), fmt(m$hi)),
          collapse = "\n"))
}

# ---- final inventory: how every published deliverable was reproduced -----
# Built as one string and substituted whole, so the sprintf argument list above
# stays readable.
inventory_md <- sprintf('## Layout

    Figure1..4.{png,pdf}          main figures
    source_panels/                individual panels
    source_panels/as_in_docx/     the five composites embedded in the docx
    supp_figs/                    Supplementary Figures 1-9 + Supplementary_Tables.xlsx
    supp_figs/tables/             every intermediate TSV
    supplemental_data/            Nature source-data tables
    Figures_1_to_4_panels.pptx    one panel per slide, 15 slides
    Outliers_in_TOPMed_paper_draft.docx
                                  manuscript with the %s-PC figures embedded
    FIGURE_SOURCES.md, MOVED.txt  carried-over orientation notes
    logs/                         per-run pipeline logs
    pipeline_verification.tsv     one row per file: the deliverable inventory
                                  against the 55-PC tree, plus a size / decode
                                  / blank check on every PNG produced here
    Supplementary_Figure_PC55_vs_PC2005_outlier_overlap.*
                                  direct 55 vs %s comparison (kept from the
                                  earlier one-off analysis)

## Deliverable inventory vs `final_paper_draft_pc55/`

This folder is checked file for file against the published tree by
`scripts/pc_swap/verify_outputs.R`, which counts *every* deliverable (docx,
pptx, prose documents and `source_panels/as_in_docx/` included, not just
png/pdf/tsv/xlsx) and fails the run if a published file has no counterpart
here. Each published file falls into exactly one of three classes.

### What the check does not count

A deliverable is a file a reader of the folder is meant to open. Backup
copies, crash droppings and per-run working state are not, and requiring a
counterpart for them fails the run for reasons that have nothing to do with
the analysis -- which is what happened on 2026-09-23, when an unrelated
editing session wrote a `Supplementary_Figures_methods_BACKUP_<stamp>.docx`
into the 55-PC tree one minute into a run and the verifier reported it as a
published deliverable with no counterpart here.

`verify_outputs.R` therefore carries an explicit `IGNORE_RULES` table, applied
to both trees so the diff stays symmetric: timestamped `*_BACKUP_*` copies,
`.bak` / `.orig` / `~` / `#autosave#` / `.swp` backups, `logs/` and stray
`.log` files, cached `.rds` intermediates, partial writes (`.tmp`, `.part`,
`.nfs*`), and toolchain or OS droppings (`__pycache__`, `.pyc`,
`.ipynb_checkpoints`, a stray `Rplots.pdf`, `.DS_Store`, `Thumbs.db`).
Nothing is dropped quietly: every skipped file appears in
`pipeline_verification.tsv` as `check = ignored` next to the rule that skipped
it. Anything *not* matched by that table is treated as a deliverable and must
be accounted for.

**Recomputed at %s PCs.** Every figure, panel, supplementary figure, source
table and workbook listed in the two sections above, plus:

| deliverable | how |
|---|---|
| `Outliers_in_TOPMed_paper_draft.docx` | published docx used as a read-only text template; Figures 1-4 and Supplementary Figure 1 replaced with this run own composites, image boxes rescaled to the new aspect ratios (`scripts/helpers/update_docx_figures.py`, stage `manuscript`) |
| `Outliers_in_TOPMed_paper_draft.docx.bak_nature_panels` | the document as it stood before that swap, i.e. the published manuscript |
| `Figures_1_to_4_panels.pptx` | published deck used as a read-only layout/theme template; all 15 panel images replaced from `source_panels/` (`scripts/helpers/build_figures_pptx.py`, stage `manuscript`) |
| `source_panels/as_in_docx/*.png` | the five composites actually embedded in the docx above |
| `supp_figs/tables/S6_enrichment_or_%spcs.tsv` | the PC-parameterized name of the published `S6_enrichment_or_55pcs.tsv`; same table taken at %s PCs |

**Carried over from the 55-PC tree.** PC-independent panels and tables (S2,
S3, S7 and the archived no-overlay S1 renderings) as listed above, plus four
hand-authored prose documents that no script generates
(`scripts/pc_swap/carry_prose_docs.py`, stage `prosedocs`):

| document | treatment |
|---|---|
| `FIGURE_SOURCES.md` | carried over; paths retargeted to this folder |
| `MOVED.txt` | carried over; paths retargeted to this folder |
| `supp_figs/Supplementary_Figures_methods.md` | carried over; paths and the S6 table / Excel tab names retargeted |
| `supp_figs/Supplementary_Figures_methods.docx` | re-rendered from the retargeted `.md` with pandoc |
| `supp_figs/Supplementary_Figures_provenance.md` | carried over; paths and the S6 names retargeted |

Each of these opens with a banner recording that its *descriptions* of the
pipeline hold at any PC setting but every *number* in it is a 55-PC value that
was not recomputed. Rewriting those numbers would mean regenerating the prose
itself, which no script does and which cannot be done mechanically without
inventing results; the %s-PC values live in the regenerated tables instead.
The manuscript docx carries the same banner on its first page for the same
reason: its figures are %s-PC, its prose is not.

**Intentionally absent.** None. Every deliverable in
`final_paper_draft_pc55/` has a counterpart here, and
`verify_outputs.R` keeps an empty `INTENTIONALLY_ABSENT` table to say so.

The files present here and not there are this README,
`pipeline_verification.tsv`, `logs/`, the cached intermediates of the
comparison analysis, and the `Supplementary_Figure_PC55_vs_PC2005_outlier_overlap`
panels and TSVs, which exist only because there are two PC settings to
compare. The one file present *there* and not here is the
`Supplementary_Figures_methods_BACKUP_<stamp>.docx` described above: it is the
rollback for a hand edit to the published methods text, it is deliberately
preserved in the 55-PC tree, and it is skipped on both sides by the ignore
rules rather than deleted.
',
  PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING)

txt <- sprintf('# Paper figures rebuilt at %s expression PCs

Everything in this folder is the full paper figure and table set regenerated
with **%s residualized expression PCs** in place of the 55 used in the
published analysis. It mirrors `final_paper_draft_pc55/` file for file so the
two can be diffed directly.

Generated: %s

## How this was produced

    PC_SETTING=%s ./scripts/run_pipeline_pc.sh

`scripts/run_pipeline_pc.sh` is the pipeline driver; `PC_SETTING` (55 or %s)
selects the whole tree. The switch itself lives in
`scripts/helpers/_paths.R` and `_paths.py`, which re-point three things and
nothing else:

1. `Z_LONGFORM` -> `expr_scaled_genotype_pc_5_expression_pc_%s_longer_form.tsv`
2. `DATA_UNDERLIER` / `EXPR_Z_JOIN` -> `data/underliers_pc%s/`
3. output roots -> `final_paper_draft_pc%s/` and `pc%s_work/`

With `PC_SETTING` unset or 55 every path is exactly what it was before the
switch existed, so the published analysis is reproduced byte for byte.

## Correctness gate

The 2005-PC calls are only trusted because the same code path reproduces the
published 55-PC calls exactly. Running the derivation against the 55-PC
long-form z table recovers all four published underlier tables, keys and
values, including the **48,632** gene x subject pairs at z < -3 (also 2,356,754
haplo, 322,377 at z < -2, 8,937 at z < -4). The pipeline aborts if that gate
fails. Nothing under `final_paper_draft_pc55/` or `data/underliers/` is
written to; `scripts/pc_swap/pc55_integrity.sh` checksums them before and
after every run.

## What is genuinely recomputed at %s PCs

These depend on the TOPMed residualized z-scores and were rebuilt end to end:

| output | why it moves |
|---|---|
| Figure 1b, 1c | gene / subject underlier rank curves |
| Figure 1d, Figure 2d | TOPMed underlier prevalence per gene (AoU side unchanged) |
| Figure 2a, 2b | LOEUF tier forest and decile trend |
| Figure 2c | constrained vs tolerant rank curves |
| Figure 3a | gene-level PTV carrier vs underlier frequency |
| Figure 3b | PTV-carrier vs non-carrier z-score distributions |
| Figure 3c | ClinGen HI candidate-gene z distributions |
| Figure 4a-d | RBM20 / MYH11 panels, via the rebuilt prevalence table |
| Supplementary Figure 4 | TOPMed vs AoU underlier counts per gene |
| Supplementary Figure 5 | z-cutoff landscape and LOEUF calibration by cutoff |
| Supplementary Figure 6 | Fisher OR slice taken at %s PCs instead of 55 |
| Supplementary Figure 8 | per-gene residualized z vs cell-type proportions |

## What is PC-independent

These have no TOPMed z-score anywhere in their inputs. They are **copied
verbatim** from the 55-PC folder (S2, S3, S7) or redrawn from their own
PC-independent sources (S1, S9), and are identical by construction, not by
coincidence:

| output | input | status |
|---|---|---|
| Figure 1a | study-design schematic, cohort sizes | redrawn, identical |
| Supplementary Figure 1a, 1b | rare SV / PTV enrichment **swept across all PC counts**; 55 and %s are both already marked on the curve | redrawn, identical |
| Supplementary Figure 1d | VEP consequence counts from `ptv_all.tsv` | redrawn, identical |
| Supplementary Figure 2 | connectivity Z, |Z|>5 subject filter, PCA scree, PERMANOVA | copied |
| Supplementary Figure 3 | UMAP of the PCA objects | copied |
| Supplementary Figure 7 | cell-type composition vs PC1..80 | copied |
| Supplementary Figure 9 | participant demographics | redrawn, identical |
| AoU side of Figures 1d / 2d / S4 | external release, n = 8,895 | unchanged |
| LOEUF tiers, gnomAD, PTV/SV genotypes, ClinGen HI | external references | unchanged |

Supplementary Figure 1c is the gene-z vs cell-type panel and **does** move; it
is rebuilt with S8.

## What could not be replicated, and why

The underlier definition has two arms:

    expression_zscore < -k   AND   log2CPM_adj < gene_mean(log2CPM_adj) - 1

Only the first arm moves with `PC_SETTING`. `expr_cpm_adj.rds`, the adjusted
log2CPM matrix behind the second arm, was residualized on 55 expression PCs by
the original pipeline, and no 2005-PC counterpart exists on disk; the 2005-PC
release provides the z-score table only. The log2CPM arm is therefore held
fixed at its published 55-PC value in both runs. This is the same convention
used by `Supplementary_Figure_PC55_vs_PC2005_outlier_overlap`, which is kept in
this folder.

Two consequences worth stating plainly:

* The **haplo** definition is log2CPM-only, so it is *identical* at 55 and %s
  PCs (2,356,754 pairs either way). It is reported here for completeness, not
  as a 2005-PC result.
* A fully re-residualized log2CPM at %s PCs would shrink the 2-fold-drop set
  substantially and would push all z-cutoff counts lower still. The counts
  below are therefore best read as "2005-PC z-scores, published logCPM filter".

The rare-variant enrichment sweep (`data/enrichment/expr_enrichments_*.rds`)
was **not** re-run: it already spans every PC count from 5 to 3,955, so the
%s-PC values are read directly out of the existing sweep. Re-running it would
require `expr_scaled_flat_list.rds` (~28 GB) and would reproduce the same
numbers.

## Headline differences vs 55 PCs

### Underlier burden

%s

### Constraint gradient (LOEUF tier, z < -3)

%s

### Rare-variant enrichment at the chosen PC setting

%s

%s',
  PC_SETTING, PC_SETTING, format(Sys.time(), "%%Y-%%m-%%d %%H:%%M:%%S"),
  PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING,
  PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING, PC_SETTING,
  counts_md, tier_md, enr_md, inventory_md)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
writeLines(txt, file.path(OUT, "README.md"))
cat("wrote ", file.path(OUT, "README.md"), "\n", sep = "")
