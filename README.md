# TOPMed / AoU RNA underlier analysis — standardized pipeline

Scripts and documentation to regenerate main Figures 1–4 and Supplementary Figure 1 for whole-blood RNA **underlier** analyses.

**Repository:** https://github.com/sherryyyyang/topmed_aou

This tree is primarily **scripts + docs** (plus rendered figures in the local submission package). Large expression matrices, genotype callsets, and individual-level data remain on analysis storage under TOPMed / dbGaP terms and are **not** redistributed here.

Use this README as the **harmonization checklist** when applying the same procedure to an independent WGS + bulk RNA-seq cohort.

---

## 0. Analysis cohort (this paper)

| Cohort | Role | N (kept) |
|--------|------|----------|
| **TOPMed** | Matched WGS + whole-blood bulk RNA-seq (5 studies) | 5,489 |
| **AoU** | External RNA-outlier prevalence release (cross-cohort panels only) | 8,895 in Source_Data tables |

Paper delivery setting: **`PC_SETTING=1055`** (see §4).

---

## 1. Key definitions (lock these for a new cohort)

| Term | Definition |
|------|------------|
| **Underlier** (primary) | Expression *z* &lt; **−3** **and** log2CPM_adj &lt; gene_mean(log2CPM_adj) − **1** |
| **Haplo** (secondary) | log2CPM arm only (PC-independent given fixed CPM matrix) |
| **Rare PTV** (Fig. 3b *x*) | High-confidence protein-truncating variant; carrier AF **&lt; 0.5%** (AF &lt; 0.005) |
| **Rare SV** (enrichment / Supp. Fig. 1) | Structural variant AF **≤ 1%**, gene ± **100 kb** |
| **Fig. 3b SV (DEL)** | Gene-overlapping PASS muCNV **DEL** whose **primary** (largest \|SVLEN\|) DEL has **AF &lt; 0.5%** |

Do **not** put common DELs (AF ≥ 0.5%) in the Fig. 3b **SV (DEL)** class.

---

## 2. End-to-end procedure (raw → figures)

Apply these stages in order for a new dataset. Paths below use:

- `ROOT` = analysis root (`pre_Watershed_scripts` locally)
- `PKG` = this submission package
- `SHEYANG` = sibling data/annotation tree

### Stage A — Genotypes (WGS)

1. **Align / call** WGS to GRCh38 (cohort-standard pipeline).
2. **SNV/indel PTVs:** annotate LoF / HC PTVs (VEP + LOFTEE or equivalent); export gene–subject carrier table analogous to `SHEYANG/ptv_all.tsv`.
3. **SVs/CNVs:** PASS callset (here: muCNV) with `ALT=<DEL|DUP|INV|…>`, `INFO/SVTYPE`, `AF`, `END`.
   - Local PASS VCFs: `…/2025-11-17/data/vcfs-pass-variants-only/chr*.vcf.gz`
4. **Rare SV catalog for enrichment:** filter AF ≤ 1%, genotype missingness, map to genes ± 100 kb → object like `data/expression/subset_genotypes_100000_0.01_select.rds`.

### Stage B — RNA-seq

1. Quantify gene expression (GENCODE **v30** gene model used here).
2. Build adjusted CPM / log2CPM matrix (`expr_cpm_adj.rds`).
3. Keep analysis subjects (`kept_subjects.rds`; metadata with cohort covariates).

### Stage C — Residualization (expression PCs + genotype PCs)

1. Compute expression PCs on the analysis matrix; genotype PCs (here: **5**).
2. Residualize expression; write **long-form z** table:
   - `data/expression/expr_scaled_genotype_pc_5_expression_pc_<NPC>_longer_form.tsv`
3. Point the pipeline at that table via `PC_SETTING=<NPC>` (`scripts/helpers/_paths.R`).

### Stage D — Call underliers

```bash
PC_SETTING=<NPC> Rscript scripts/pc_swap/build_underlier_artifacts.R
```

Writes under `data/underliers_pc<NPC>/` (or `data/underliers/` when NPC=55):

| Artifact | Content |
|----------|---------|
| `expr_z_join.rds` | gene × subject: *z*, log2CPM_adj |
| `z_score_underlier_z_-3_clean.rds` | primary underlier pairs |
| `z_score_tally_*_clean.rds` | per-gene / per-subject ranks & prevalences |

**Note:** Only the *z* arm moves with `PC_SETTING`. The log2CPM arm uses the fixed `expr_cpm_adj.rds` matrix in this package’s PC-swap design.

### Stage E — Choose NPC (`PC_SETTING`)

Sweep rare-SV and rare-PTV enrichment vs number of expression PCs (Supp. Fig. 1):

```bash
# enrichment sweeps → plot
Rscript scripts/supplement/plot_enrichment_npc_sv_vs_ptv.R
```

| Variant class | Peak OR (this study) | Rule used for paper |
|---------------|----------------------|---------------------|
| Rare SV (AF ≤ 1%, ±100 kb) | ~805 | — |
| Rare PTV / LoF (AF ≤ 0.5%) | **1055** | **Chosen NPC = PTV enrichment peak = 1055** |

Documented in `docs/PC_PLATEAU_z3.md`. Mark a single dashed guide at **1055** on Supp. Fig. 1a/b.

For a new cohort: re-run the enrichment sweep and pick the **PTV-curve peak** (or pre-register an AoU-style 95%-of-max plateau if that is the study design).

### Stage F — Cross-cohort prevalence (optional)

AoU side is a **fixed external release** in this paper (not re-residualized with TOPMed PCs). Rebuild TOPMed vs AoU tables with:

```bash
PC_SETTING=<NPC> Rscript scripts/pc_swap/build_topmed_vs_aou_prevalence.R
```

Shared loader for Fig. 1d / 2d scatter geometry: `scripts/helpers/load_topmed_aou_prevalence_scatter.R`.

### Stage G — PTV vs RNA underlier comparison (Fig. 3)

1. **Carrier *z* distributions (Fig. 3a):** rare PTV carriers vs non-carriers (`plot_ptv_vs_nonptv_zscore.R`).
2. **Gene-level scatter (Fig. 3b):** *x* = % subjects with rare PTV (AF &lt; 0.5%); *y* = % subjects with underliers.
3. **ClinGen HI examples (Fig. 3c):** candidate histograms (`plot_candidate_gene_zscore_dist_clingen_z3.R`).

#### Fig. 3b dosage colouring (default)

| Legend | Rule |
|--------|------|
| **SV (DEL)** | Primary gene-overlapping PASS **DEL** has **AF &lt; 0.5%** (`fig3b_sv_del_genes.tsv`) |
| **CNV / multi-copy** | SegDup (≥50% gene body, fracMatch ≥ 0.95) **or** curated multi-copy family prefixes; not already SV (DEL) |
| **Other genes** | Remainder |

**Labels:** all **SV (DEL)** genes + all **top-left** CNV / multi-copy genes (`pct_under ≥ 0.5%` and `pct_ptv ≤ 0.25%`).

**AF gate:** any gene whose primary gene-overlapping DEL has **AF ≥ 0.5%** is **excluded** from **SV (DEL)** (it may still appear as **CNV / multi-copy** via SegDup/family annotation).

Reference files (`ROOT/data/references/` and `PKG/supporting_data/`):

- `fig3b_sv_del_genes.tsv` / `Source_Data_Fig3b_sv_del_genes.tsv`
- `fig3b_sv_annotation_af0.005.tsv` / `Source_Data_Fig3b_sv_annotation_af0.005.tsv` (saved AF &lt; 0.5% SV table)
- `fig3b_cnv_multicopy_genes.tsv`
- `SOURCE_fig3b_sv_del.txt`, `SOURCE_fig3b_cnv_multicopy.txt`

Optional: `FIG3B_COLOR=nmd` for NMD-escape colouring (not default).

### Stage H — Constraint calibration (Fig. 2)

LOEUF-tier forest, decile trend, constrained vs tolerant ranks, TOPMed vs AoU coloured by LOEUF (`scripts/fig2_constraint/`).

### Stage I — ACMG exemplars (Fig. 4)

1. Build shared rare-outlier gene table (`plot_shared_rare_outlier_genes_4dot.R`):
   - Genes with **−log10** RNA-outlier prevalence **≥ 3** in **both** TOPMed and AoU (CUT = 3).
   - Annotate ACMG SF v3.2 membership (`acmg_sf_v3.2.csv`).
   - **Force-include** paper focus genes: **RBM20**, **MYH11** (even if they fall outside the cut).
2. Plot rank curves + 4-dot summaries (`plot_acmg_4dot_no_ci.R`) for **RBM20** and **MYH11**.

**Why these two:** ACMG secondary-findings genes that illustrate contrasting RNA vs PTV prevalence patterns in the shared high-prevalence set; force-included so low TOPMed underlier counts (e.g. MYH11 at PC1055) do not drop them from the figure.

### Stage J — Assemble Nature composites

```bash
PC_SETTING=1055 ./scripts/run_pipeline_pc.sh          # full rebuild
PC_SETTING=1055 Rscript scripts/assemble_nature_figures.R
PC_SETTING=1055 ONLY_FIG=3 Rscript scripts/assemble_nature_figures.R
```

Outputs: `figures/main/Figure{1–4}.{png,pdf}`, `figures/supplementary/`, `source_panels/`.

---

## 3. Path map (this deployment)

| Role | Path |
|------|------|
| Package (`PKG`) | `…/sheyang/submission_package_nov_2026` |
| Analysis root (`ROOT`) | `…/sheyang/pre_Watershed_scripts` |
| Path switch | `scripts/helpers/_paths.R` / `_paths.py` |
| Detailed absolute inputs | `docs/PATH_TRACEABILITY.md` |
| PC choice notes | `docs/PC_PLATEAU_z3.md` |

`PC_SETTING` re-points:

1. `Z_LONGFORM` → `expr_scaled_genotype_pc_5_expression_pc_<NPC>_longer_form.tsv`
2. `DATA_UNDERLIER` → `data/underliers_pc<NPC>/` (or published `data/underliers/` for 55)
3. Work / delivery dirs → `pc<NPC>_work/`, submission package figures

---

## 4. Checklist for an independent cohort

- [ ] WGS: PTV table + PASS SV VCF with `SVTYPE` / `AF`
- [ ] RNA: GENCODE gene quant + adjusted CPM + kept-subject list
- [ ] Residualize (genotype PCs + expression PCs); write long-form *z*
- [ ] Call underliers: *z* &lt; −3 **and** log2CPM &lt; mean − 1
- [ ] Enrichment vs NPC → choose **PC_SETTING** (document peak / plateau rule)
- [ ] Gene-level PTV % vs underlier %; colour **SV (DEL, AF&lt;0.5%)** vs **CNV/multi-copy** vs other
- [ ] Optional: external cohort prevalence scatter (lock external release)
- [ ] ACMG (or disease-panel) exemplars with pre-registered focus genes + force-include rule
- [ ] Assemble figures; deposit Source_Data tables

---

## 5. Figure overview (this package)

| Figure | Content |
|--------|---------|
| **Fig. 1** | Design; gene/subject underlier ranks; TOPMed vs AoU prevalence |
| **Fig. 2** | LOEUF calibration; constrained vs tolerant; prevalence × LOEUF |
| **Fig. 3** | PTV-carrier *z*; gene-level PTV vs underlier (dosage colours); ClinGen HI |
| **Fig. 4** | RBM20 / MYH11 ACMG exemplars |
| **Supp. Fig. 1** | Rare SV / PTV enrichment vs NPC (guide at 1055) |

---

## 6. Rebuild commands (quick)

```bash
cd PKG   # or ROOT — assemble resolves scripts via _paths.R

PC_SETTING=1055 ./scripts/run_pipeline_pc.sh
PC_SETTING=1055 ONLY_FIG=1 Rscript scripts/assemble_nature_figures.R
PC_SETTING=1055 ONLY_FIG=2 Rscript scripts/assemble_nature_figures.R
PC_SETTING=1055 ONLY_FIG=3 Rscript scripts/assemble_nature_figures.R
PC_SETTING=1055 ONLY_FIG=4 Rscript scripts/assemble_nature_figures.R
```

Published checksummed tree uses `PC_SETTING=55`. Non-55 rebuilds must not overwrite it (pipeline gate).

---

## 7. Further docs

| Doc | Topic |
|-----|--------|
| `docs/PATH_TRACEABILITY.md` | Absolute raw/derived paths |
| `docs/PC_PLATEAU_z3.md` | Why PC1055 |
| `docs/FIGURE_SOURCES.md` | Panel script map |
| `docs/README_Fig3b_NMD_escape.md` | Optional NMD colouring |
| `docs/Supplementary_Figures_methods.md` | Supp methods |
| `supporting_data/MANIFEST.tsv` | Source_Data inventory |
