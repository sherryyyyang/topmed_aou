#!/usr/bin/env bash
# Rebuild main Figures 1–4.
#
#   DATA_ROOT=/path/to/analysis PC_SETTING=<expression PC count> ./scripts/run_pipeline_pc.sh
#
# PC_SETTING is an input. The rule that chooses it is not in this repository.
# DATA_ROOT is the analysis directory that contains data/ and these scripts.
# Optional: STAGES=derive,prevalence,fig3sv,fig1,fig2,fig3,fig4
set -euo pipefail

ROOT="${DATA_ROOT:-${TOPMED_OUTLIERS_ROOT:-${ANALYSIS_ROOT:-}}}"
if [[ -z "$ROOT" ]]; then
  echo "Set DATA_ROOT to the analysis directory." >&2
  exit 2
fi
if [[ -z "${PC_SETTING:-}" ]]; then
  echo "Set PC_SETTING to the expression-PC count." >&2
  exit 2
fi
if ! [[ "$PC_SETTING" =~ ^[1-9][0-9]*$ ]]; then
  echo "PC_SETTING must be a positive integer (got '$PC_SETTING')" >&2
  exit 2
fi

cd "$ROOT"
export DATA_ROOT="$ROOT"
export TOPMED_OUTLIERS_ROOT="$ROOT"
export PC_SETTING
export DEFINITION_LABEL="${DEFINITION_LABEL:-z_-3}"
RSCRIPT="${RSCRIPT:-Rscript}"

ALL_STAGES=(longform derive prevalence fig3sv fig4af fig1 fig2 fig3 fig4)
if [[ "${1:-}" == "--list" ]]; then
  printf '%s\n' "${ALL_STAGES[@]}"
  exit 0
fi

STAGES="${STAGES:-}"
read -r -a WANTED <<< "${STAGES//,/ }"
[[ ${#WANTED[@]} -eq 0 ]] && WANTED=("${ALL_STAGES[@]}")
want() { local s; for s in "${WANTED[@]}"; do [[ "$s" == "$1" ]] && return 0; done; return 1; }
for s in "${WANTED[@]}"; do
  printf '%s\n' "${ALL_STAGES[@]}" | grep -qx "$s" || { echo "unknown stage '$s'" >&2; exit 2; }
done

run_stage() {
  local name="$1"; shift
  if ! want "$name"; then
    echo "-- skip $name"
    return 0
  fi
  echo
  echo "== $name =="
  "$@"
}

stage_longform() {
  "$RSCRIPT" "$ROOT/scripts/pc_swap/write_z_longform_from_flat.R"
}
stage_derive() {
  "$RSCRIPT" "$ROOT/scripts/pc_swap/build_underlier_artifacts.R"
}
stage_prevalence() {
  "$RSCRIPT" "$ROOT/scripts/pc_swap/build_topmed_vs_aou_prevalence.R"
}
stage_fig3sv() {
  "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/build_cnv_multicopy_gene_set.R"
  "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/build_fig3b_segdup_svtype_genes.R"
}
stage_fig4af() {
  "$RSCRIPT" "$ROOT/scripts/fig4_acmg/sum_af_exome_topmed_overlap.R"
}
stage_fig1() {
  "$RSCRIPT" "$ROOT/scripts/fig1_overview/plot_underlier_tallies_topmed_aou.R"
  ONLY_FIG=1 "$RSCRIPT" "$ROOT/scripts/assemble_nature_figures.R"
}
stage_fig2() {
  DEFINITION_LABEL=z_-3 "$RSCRIPT" "$ROOT/scripts/fig2_constraint/plot_loeuf_calibration.R"
  DEFINITION_LABEL=z_-3 "$RSCRIPT" "$ROOT/scripts/fig2_constraint/plot_underlier_tallies_loeuf_constrained_tolerant.R"
  ONLY_FIG=2 DEFINITION_LABEL=z_-3 "$RSCRIPT" "$ROOT/scripts/assemble_nature_figures.R"
}
stage_fig3() {
  "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level_2x2.R"
  "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/plot_ptv_vs_nonptv_zscore.R"
  "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/plot_candidate_gene_zscore_dist_clingen_z3.R"
  ONLY_FIG=3 "$RSCRIPT" "$ROOT/scripts/assemble_nature_figures.R"
}
stage_fig4() {
  "$RSCRIPT" "$ROOT/scripts/fig4_acmg/plot_shared_rare_outlier_genes_4dot.R"
  "$RSCRIPT" "$ROOT/scripts/fig4_acmg/plot_acmg_4dot_no_ci.R"
  ONLY_FIG=4 "$RSCRIPT" "$ROOT/scripts/assemble_nature_figures.R"
}

run_stage longform    stage_longform
run_stage derive      stage_derive
run_stage prevalence  stage_prevalence
run_stage fig3sv      stage_fig3sv
run_stage fig4af      stage_fig4af
run_stage fig1        stage_fig1
run_stage fig2        stage_fig2
run_stage fig3        stage_fig3
run_stage fig4        stage_fig4
echo
echo "Figures 1–4 stages finished."
