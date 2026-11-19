# Expression-PC plateau at z < −3 (AoU-style)

**Chosen PC_SETTING: 505**

## Rule

AoU RareVariantEnrichment–style 95%-of-max plateau, applied to Supplementary
Figure 1 underlier enrichment curves at **cutoff z < −3**:

1. For each variant class (rare SV; rare PTV/LoF), take the earliest NPC where
   Fisher OR ≥ **0.95 × max(OR)** on that curve.
2. Paper PC count = **max** of those earliest-plateau NPCs (minimal count at
   which *both* definitions are on plateau).

## Result on TOPMed sweep (`enrichment_npc_sv_vs_ptv_underlier.tsv`)

| Variant class | Peak OR at | 95% plateau entry |
|---|---|---|
| Rare SV within 100 kb (AF ≤ 1%) | 805 | **255** |
| LoF / PTV in gene (AF ≤ 0.5%) | 1055 | **505** |

→ **PC_SETTING = max(255, 505) = 505**

Source sweep: `pre_Watershed_scripts/data/enrichment/expr_enrichments_{sv_full,ptv}.rds`
via `scripts/supplement/plot_enrichment_npc_sv_vs_ptv.R`.
