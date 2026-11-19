#!/usr/bin/env bash
#
# run_pipeline_pc.sh -- rebuild the whole paper figure/table set at a chosen
# number of residualized expression PCs.
#
# The published analysis uses 55 expression PCs. Any other PC_SETTING (e.g.
# 505 = z<-3 AoU-style 95% plateau; 2005 = Gavish-Donoho) reruns every
# z-score-dependent step against that PC's long-form z table and writes the
# result to submission_package_nov_2026/, leaving the 55-PC tree untouched.
#
#   PC_SETTING=505 ./scripts/run_pipeline_pc.sh
#   PC_SETTING=505 STAGES=fig3,fig4 ./scripts/run_pipeline_pc.sh    # subset
#   PC_SETTING=505 FORCE=1 ./scripts/run_pipeline_pc.sh             # ignore caches
#   ./scripts/run_pipeline_pc.sh --list                              # show stages
#
# Environment:
#   PC_SETTING         55 (default) or any positive integer (505, 2005, ...)
#   STAGES             comma-separated subset of stages; default = all
#   FORCE              1 = delete derived caches first and recompute everything
#   SKIP_GATE          1 = skip the 55-PC round-trip gate (not recommended)
#   ALLOW_PC55_WRITE   1 = permit render stages at PC_SETTING=55, which would
#                      overwrite the published tree. Off by default on purpose.
#
# Every stage prints a banner with start time, end time, elapsed seconds and
# exit status. The whole run is teed to
#   final_paper_draft_pc<PC_SETTING>/logs/pipeline_pc<PC_SETTING>_<stamp>.log
#
# set -euo pipefail: a failing stage aborts the run rather than leaving a
# half-built figure set that looks complete.
set -euo pipefail

DEFAULT_ROOT="/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
DEFAULT_PKG="/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
ROOT="${TOPMED_OUTLIERS_ROOT:-${ANALYSIS_ROOT:-$DEFAULT_ROOT}}"
cd "$ROOT"

PC_SETTING="${PC_SETTING:-55}"
FORCE="${FORCE:-0}"
SKIP_GATE="${SKIP_GATE:-0}"
ALLOW_PC55_WRITE="${ALLOW_PC55_WRITE:-0}"
export PC_SETTING
export TOPMED_OUTLIERS_ROOT="$ROOT"

RSCRIPT="${RSCRIPT:-Rscript}"
if [[ -z "${PAPER_PYTHON:-}" ]]; then
  if [[ -x /net/wonderland/home/sheyang/.pixi/envs/python/bin/python3 ]]; then
    PAPER_PYTHON=/net/wonderland/home/sheyang/.pixi/envs/python/bin/python3
  else
    PAPER_PYTHON=python3
  fi
fi
export PAPER_PYTHON

if ! [[ "$PC_SETTING" =~ ^[1-9][0-9]*$ ]]; then
  echo "PC_SETTING must be a positive integer (got '$PC_SETTING')" >&2; exit 2
fi
if [[ "$PC_SETTING" == "55" ]]; then
  OUT_ROOT="$ROOT/final_paper_draft_pc55"
  WORK="$ROOT"
  UNDERLIER_DIR="$ROOT/data/underliers"
else
  OUT_ROOT="${TOPMED_SUBMISSION_PACKAGE:-$DEFAULT_PKG}"
  WORK="$ROOT/pc${PC_SETTING}_work"
  UNDERLIER_DIR="$ROOT/data/underliers_pc${PC_SETTING}"
fi
export TOPMED_SUBMISSION_PACKAGE="$OUT_ROOT"
VALIDATION_DIR="$ROOT/data/underliers_pc55_validation"

ALL_STAGES=(preflight gate derive seed carryover prevalence fig1 fig2 fig3 fig4 supplement tables sourcedata manuscript prosedocs readme verify)

if [[ "${1:-}" == "--list" ]]; then
  printf '%s\n' "${ALL_STAGES[@]}"; exit 0
fi

STAGES="${STAGES:-}"
read -r -a WANTED <<< "${STAGES//,/ }"
[[ ${#WANTED[@]} -eq 0 ]] && WANTED=("${ALL_STAGES[@]}")
want() { local s; for s in "${WANTED[@]}"; do [[ "$s" == "$1" ]] && return 0; done; return 1; }
for s in "${WANTED[@]}"; do
  printf '%s\n' "${ALL_STAGES[@]}" | grep -qx "$s" || { echo "unknown stage '$s'" >&2; exit 2; }
done

# Logs for the published setting go outside final_paper_draft_pc55/, which is
# checksum-protected and must not gain files.
if [[ "$PC_SETTING" == "55" ]]; then LOGDIR="$ROOT/logs"; else LOGDIR="$OUT_ROOT/logs"; fi
mkdir -p "$LOGDIR"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="$LOGDIR/pipeline_pc${PC_SETTING}_${STAMP}.log"

# Everything below is teed to $LOG as well as the terminal.
exec > >(tee -a "$LOG") 2>&1

STAGE_STATUS=()
t_run_start=$SECONDS

banner() { printf '\n%s\n== %-68s ==\n%s\n' "$(printf '=%.0s' {1..76})" "$1" "$(printf '=%.0s' {1..76})"; }

run_stage() {
  local name="$1"; shift
  if ! want "$name"; then
    echo "-- skip stage '$name' (not in STAGES)"
    STAGE_STATUS+=("$name:SKIPPED")
    return 0
  fi
  banner "STAGE $name  |  start $(date '+%F %T')"
  local t0=$SECONDS rc=0
  # Run in a subshell so a stage cannot leak cd / env into the next one.
  #
  # The subshell must NOT be written as `( ... ) || rc=$?`: bash ignores `set -e`
  # inside any command that is an operand of && or ||, so a stage function would
  # carry on past its first failing command and report the status of its *last*
  # one. That is not hypothetical -- it hid a pc55_integrity.sh failure behind a
  # later verify_outputs.R failure in the 2026-09-23 run, and would have reported
  # the stage OK had verify_outputs.R passed. Disabling `set -e` around a plain
  # subshell keeps it honoured inside, so a stage stops at its first failure.
  set +e
  ( set -euo pipefail; "$@" )
  rc=$?
  set -e
  local dt=$((SECONDS - t0))
  if [[ $rc -eq 0 ]]; then
    printf -- '-- STAGE %s OK  (%ds)  end %s\n' "$name" "$dt" "$(date '+%F %T')"
    STAGE_STATUS+=("$name:OK:${dt}s")
  else
    printf -- '-- STAGE %s FAILED rc=%d  (%ds)  end %s\n' "$name" "$rc" "$dt" "$(date '+%F %T')"
    STAGE_STATUS+=("$name:FAILED(rc=$rc):${dt}s")
    summary; exit "$rc"
  fi
}

summary() {
  banner "PIPELINE SUMMARY  PC_SETTING=$PC_SETTING  total $((SECONDS - t_run_start))s"
  local s
  for s in "${STAGE_STATUS[@]}"; do printf '   %s\n' "$s"; done
  printf '   log: %s\n' "$LOG"
}

# --------------------------------------------------------------------------
# stages
# --------------------------------------------------------------------------

stage_preflight() {
  echo "ROOT          = $ROOT"
  echo "PC_SETTING    = $PC_SETTING"
  echo "output root   = $OUT_ROOT"
  echo "work tree     = $WORK"
  echo "underlier dir = $UNDERLIER_DIR"
  echo "Rscript       = $($RSCRIPT --version 2>&1 | head -1)"
  echo "python        = $PAPER_PYTHON"
  echo "FORCE=$FORCE SKIP_GATE=$SKIP_GATE"
  echo "stages requested: ${WANTED[*]}"

  if [[ "$PC_SETTING" == "55" && "$ALLOW_PC55_WRITE" != "1" ]]; then
    echo
    echo "REFUSING to run render stages at PC_SETTING=55: that would overwrite the"
    echo "published tree. Re-run with ALLOW_PC55_WRITE=1 if this is really intended,"
    echo "or use STAGES=gate to only exercise the round-trip validation."
    for s in "${WANTED[@]}"; do
      case "$s" in preflight|gate|verify) ;; *) exit 3 ;; esac
    done
  fi

  mkdir -p "$OUT_ROOT" "$LOGDIR"
  # Baseline of every file the run must not disturb.
  if [[ ! -f "$ROOT/data/pc55_integrity_baseline.sha256" ]]; then
    "$ROOT/scripts/pc_swap/pc55_integrity.sh" snapshot
  fi
  "$ROOT/scripts/pc_swap/pc55_integrity.sh" verify

  if [[ "$FORCE" == "1" ]]; then
    echo "FORCE=1: clearing derived caches for PC_SETTING=$PC_SETTING"
    [[ "$PC_SETTING" != "55" ]] && rm -rf "$UNDERLIER_DIR"
    rm -rf "$VALIDATION_DIR"
    rm -f "$WORK/fig3_ptv_vs_rna/panel_b_carrier_zscore/ptv_vs_nonptv_zscore_cache_v2.rds"
  fi
}

# The correctness gate. Rebuilds the 55-PC underlier calls straight from the
# 55-PC long-form z table using the same code path the 2005 run uses, and
# requires them to reproduce the published tables (48,632 pairs at z < -3).
stage_gate() {
  if [[ "$SKIP_GATE" == "1" ]]; then echo "SKIP_GATE=1 -- gate not run"; return 0; fi
  PC_SETTING=55 UNDERLIER_OUT_DIR="$VALIDATION_DIR" \
    "$RSCRIPT" "$ROOT/scripts/pc_swap/build_underlier_artifacts.R"

  local n
  n=$(V="$VALIDATION_DIR" "$RSCRIPT" -e \
        'cat(nrow(readRDS(file.path(Sys.getenv("V"), "z_score_underlier_z_-3_clean.rds"))))' \
      | tail -1)
  echo "round-trip z<-3 underlier pairs: $n (expected 48632)"
  if [[ "$n" != "48632" ]]; then
    echo "GATE FAILED: expected 48632 underlier pairs at z<-3, got '$n'" >&2
    exit 1
  fi
  echo "GATE PASSED: 55-PC round-trip reproduces the published underlier calls"

  PC_SETTING=55 FIG1C_OUT="$VALIDATION_DIR/fig1c" \
    "$RSCRIPT" "$ROOT/scripts/pc_swap/build_topmed_vs_aou_prevalence.R"
}

# z-scores -> underlier calls -> per-gene / per-subject tallies at PC_SETTING.
stage_derive() {
  "$RSCRIPT" "$ROOT/scripts/pc_swap/build_underlier_artifacts.R"
  cat "$UNDERLIER_DIR/underlier_counts_summary.tsv"
}

# Populate the PC-specific work tree with the inputs the figure scripts expect
# (reference tables, sweep outputs, static panels). Regenerated outputs
# overwrite these in later stages; anything not regenerated is PC-independent.
stage_seed() {
  [[ "$PC_SETTING" == "55" ]] && { echo "nothing to seed at 55 PCs"; return 0; }
  mkdir -p "$WORK"
  local d
  for d in fig1_overview fig2_constraint fig3_ptv_vs_rna fig4_acmg supplement paper; do
    if [[ ! -d "$WORK/$d" ]]; then
      echo "seeding $d"; cp -a "$ROOT/$d" "$WORK/"
    else
      echo "$d already present"
    fi
  done
  # This cache is keyed to the 55-PC z-scores; it must not be reused.
  rm -f "$WORK/fig3_ptv_vs_rna/panel_b_carrier_zscore/ptv_vs_nonptv_zscore_cache_v2.rds"
}

stage_prevalence() {
  "$RSCRIPT" "$ROOT/scripts/pc_swap/build_topmed_vs_aou_prevalence.R"
}

# assemble_nature_figures.R sources the panel scripts with paper.assemble=TRUE,
# which suppresses their write_tsv() calls so the composite build cannot clobber
# a panel table. That means the tables only get refreshed if the panel script is
# also run on its own, which is what the standalone calls below are for. Without
# them the figures would be redrawn at 2005 PCs while the tables beside them
# stayed at their seeded 55-PC values.
stage_fig1() {
  "$RSCRIPT" "$ROOT/scripts/fig1_overview/plot_underlier_tallies_topmed_aou.R"
  ONLY_FIG=1 "$RSCRIPT" "$ROOT/scripts/assemble_nature_figures.R"
}

# The LOEUF calibration tables are needed at every cutoff, because
# Supplementary Figure 5b dodges the tiers across haplo / z<-2 / -3 / -4.
stage_fig2() {
  local def
  for def in haplo z_-2 z_-3 z_-4; do
    echo "--- LOEUF calibration, DEFINITION_LABEL=$def"
    DEFINITION_LABEL="$def" "$RSCRIPT" "$ROOT/scripts/fig2_constraint/plot_loeuf_calibration.R"
  done
  "$RSCRIPT" "$ROOT/scripts/fig2_constraint/plot_underlier_tallies_loeuf_constrained_tolerant.R"
  "$RSCRIPT" "$ROOT/scripts/analysis/assess_wilson_vs_subject_variance.R"
  ONLY_FIG=2 "$RSCRIPT" "$ROOT/scripts/assemble_nature_figures.R"
}

stage_fig3() {
  "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level_2x2.R"
  # Writes fig3_ptv_vs_rna/extra/ptv_vs_underlier_gene_level_<def>.tsv, which is
  # the source table behind Figure 3a; the 2x2 variant above writes elsewhere.
  UNDERLIER_LABEL=z_-3 "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level.R"
  "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/plot_ptv_vs_nonptv_zscore.R"
  "$RSCRIPT" "$ROOT/scripts/fig3_ptv_vs_rna/plot_candidate_gene_zscore_dist_clingen_z3.R"
  ONLY_FIG=3 "$RSCRIPT" "$ROOT/scripts/assemble_nature_figures.R"
}

# Figure 4's dot data is rebuilt first: it reads the PC-dependent TOPMed-vs-AoU
# prevalence table written by stage_prevalence.
stage_fig4() {
  "$RSCRIPT" "$ROOT/scripts/fig4_acmg/plot_shared_rare_outlier_genes_4dot.R"
  "$RSCRIPT" "$ROOT/scripts/fig4_acmg/plot_acmg_4dot_no_ci.R"
  ONLY_FIG=4 "$RSCRIPT" "$ROOT/scripts/assemble_nature_figures.R"
}

# S1a/S1b (PC sweep) and S1c/S1d (VEP) are PC-independent but cheap to redraw,
# so they are rebuilt in place rather than copied. S4/S5/S6 and S8 genuinely
# move with PC_SETTING.
stage_supplement() {
  "$RSCRIPT" "$ROOT/scripts/supplement/plot_enrichment_npc_sv_vs_ptv.R"
  "$RSCRIPT" "$ROOT/scripts/supplement/plot_ptv_consequence_summary.R"
  "$RSCRIPT" "$ROOT/scripts/supplement/make_supplementary_s2_s5.R"
  "$RSCRIPT" "$ROOT/scripts/supplement/plot_s1_celltype_overlay_s7.R"
  "$RSCRIPT" "$ROOT/scripts/supplement/plot_s8_gene_z_vs_celltype.R"
  "$RSCRIPT" "$ROOT/scripts/supplement/plot_s9_participant_characteristics.R"
}

stage_tables() {
  "$PAPER_PYTHON" "$ROOT/scripts/supplement/pack_supplementary_tables.py"
}

# Panels whose inputs contain no TOPMed z-scores are carried over verbatim from
# the 55-PC folder rather than being redrawn, and are listed in the README.
# This runs *before* the supplement stage so that make_supplementary_s2_s5.R
# takes its documented S2/S3 skip path instead of recomputing WGCNA and UMAP.
stage_carryover() {
  [[ "$PC_SETTING" == "55" ]] && { echo "n/a at 55 PCs"; return 0; }
  "$ROOT/scripts/pc_swap/copy_pc_independent.sh"
}

stage_sourcedata() {
  "$RSCRIPT" "$ROOT/scripts/pc_swap/build_source_data.R"
}

# The two bundled manuscript artifacts. Both take the published 55-PC file as a
# read-only text/layout template and swap in the panels this run just drew, so
# the figures inside them are never stale 55-PC renders. The docx also gets a
# banner saying its prose is inherited and not recomputed.
# This also (re)populates source_panels/as_in_docx/, which is defined as "the
# composites that are inside the docx" and so can only be written here:
# assemble_nature_figures.R is always invoked with ONLY_FIG=<n> by the figure
# stages above, which returns before its Supplementary Figure 1 block.
stage_manuscript() {
  [[ "$PC_SETTING" == "55" ]] && { echo "n/a at 55 PCs (published manuscript is not rebuilt)"; return 0; }
  "$PAPER_PYTHON" "$ROOT/scripts/helpers/update_docx_figures.py"
  "$PAPER_PYTHON" "$ROOT/scripts/helpers/build_figures_pptx.py"
}

# Hand-authored prose (FIGURE_SOURCES.md, MOVED.txt, the supplementary methods
# and provenance documents). No script generates these; they are carried over
# with path/table-name retargeting plus a banner marking their numbers as
# 55-PC values. Runs after tables so the retargeted S6 name can be checked.
stage_prosedocs() {
  [[ "$PC_SETTING" == "55" ]] && { echo "n/a at 55 PCs (prose documents are the originals)"; return 0; }
  "$PAPER_PYTHON" "$ROOT/scripts/pc_swap/carry_prose_docs.py"
}

stage_readme() {
  PC_SETTING="$PC_SETTING" "$RSCRIPT" "$ROOT/scripts/pc_swap/write_readme.R"
}

stage_verify() {
  "$ROOT/scripts/pc_swap/pc55_integrity.sh" verify
  [[ "$PC_SETTING" == "55" ]] && return 0
  "$RSCRIPT" "$ROOT/scripts/pc_swap/verify_outputs.R"
}

# --------------------------------------------------------------------------
banner "PIPELINE START  PC_SETTING=$PC_SETTING  $(date '+%F %T')"
run_stage preflight        stage_preflight
run_stage gate             stage_gate
run_stage derive           stage_derive
run_stage seed             stage_seed
run_stage carryover        stage_carryover
run_stage prevalence       stage_prevalence
run_stage fig1             stage_fig1
run_stage fig2             stage_fig2
run_stage fig3             stage_fig3
run_stage fig4             stage_fig4
run_stage supplement       stage_supplement
run_stage tables           stage_tables
run_stage sourcedata       stage_sourcedata
run_stage manuscript       stage_manuscript
run_stage prosedocs        stage_prosedocs
run_stage readme           stage_readme
run_stage verify           stage_verify
summary
echo "PIPELINE COMPLETE"
