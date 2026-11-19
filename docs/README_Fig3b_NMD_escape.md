# Fig3b PTV NMD-escape annotation (VEP NMD plugin)

## Variant set
Unique LOFTEE high-confidence PTVs from `sheyang/ptv_all.tsv` (chrom/pos/ref/alt),
the same table underlying the Fig3 PTV vs non-PTV / carrier z-score panel
(`plot_ptv_vs_nonptv_zscore.R`; AF<0.1% for the paper panel).

## Method
Ensembl VEP **v113.0** offline with cache **109_GRCh38** and plugin **NMD**:

```bash
# see VEP_NMD_COMMAND.sh
vep -i fig3b_ptv_unique.vcf \
  --offline --cache --dir_cache $CACHE --cache_version 109 \
  --assembly GRCh38 --species homo_sapiens \
  --plugin NMD --dir_plugins $PLUGIN_DIR \
  --symbol --canonical --mane_select --fork 8 \
  -o fig3b_ptv_vep_nmd.txt
```

- VEP binary / plugins: `~/.pixi/envs/ensembl-vep` (`NMD.pm` in `lib/perl5/site_perl`)
- Cache: `sheyang/annotation/vep_cache/homo_sapiens/109_GRCh38`
- Plugin flag: `--plugin NMD` (outputs Extra `NMD=NMD_escaping_variant` when rules match)
- Per-variant call prefers **MANE Select** transcript with a PTV SO term, else canonical, else first PTV transcript

NMD-escape rules (plugin): last exon; ≤50 bp from end of penultimate exon; first ≤100 coding bases; intronless transcript.

## Key counts (unique variants)
| Set | NMD_escaping | NMD_triggering | inconclusive |
|-----|-------------:|---------------:|-------------:|
| all HC unique | 8156 (22.8%) | 27493 (77.0%) | 67 (0.2%) |
| AF<0.1% (Fig3b) | 7895 (22.7%) | 26885 (77.2%) | 64 (0.2%) |

Note: `ptv_all.tsv` is already LOFTEE HC (END_TRUNC filtered), so residual VEP-NMD escapes
largely reflect first-100bp / intronless / transcript-choice differences vs LOFTEE.

## Expression pattern (PC2005, AF<0.1%, kept subjects)
NMD-triggering carriers show the low-expression shift; NMD-escaping carriers resemble non-carriers.

See `fig3b_zscore_by_nmd_status_AF0.1_pc2005.tsv` and supplemental ECDF.

## Outputs
- `fig3b_ptv_nmd_escape_by_variant.tsv` — per-variant NMD status
- `fig3b_ptv_vep_nmd.txt` — raw VEP
- `Fig3b_supp_nmd_escape_ecdf_pc2005.{pdf,png}` — supplemental panel
- Copied into submission package `supplemental_data/` and `supporting_data/` as `Source_Data_Fig3b_*`

Main Nature Fig3b was **not** recolored (GO-class work left intact). NMD status can be layered later as shape/facet.
