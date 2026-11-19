#!/usr/bin/env bash
# Checksum guard for everything the 2005-PC run must leave alone.
#
#   pc55_integrity.sh snapshot   write the baseline manifest (once)
#   pc55_integrity.sh verify     re-hash and diff against the baseline
#
# Covers the published figure folder plus every 55-PC derived data artifact and
# panel directory that a figure script might write into by accident.
set -euo pipefail

ROOT="/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
MANIFEST="${ROOT}/data/pc55_integrity_baseline.sha256"

TARGETS=(
  "final_paper_draft_pc55"
  "data/underliers"
  "fig1_overview"
  "fig2_constraint"
  "fig3_ptv_vs_rna"
  "fig4_acmg"
  "supplement"
  "paper"
)

hash_all() {
  cd "$ROOT"
  for t in "${TARGETS[@]}"; do
    [[ -e "$t" ]] || continue
    find "$t" -type f -print0 | sort -z | xargs -0 -r sha256sum
  done
}

case "${1:-verify}" in
  snapshot)
    hash_all > "$MANIFEST"
    echo "baseline written: $MANIFEST ($(wc -l < "$MANIFEST") files)"
    ;;
  verify)
    if [[ ! -f "$MANIFEST" ]]; then
      echo "NO BASELINE at $MANIFEST -- run 'pc55_integrity.sh snapshot' first" >&2
      exit 2
    fi
    tmp="$(mktemp)"; trap 'rm -f "$tmp"' EXIT
    hash_all > "$tmp"
    if diff -q "$MANIFEST" "$tmp" >/dev/null; then
      echo "pc55 integrity OK: all $(wc -l < "$MANIFEST") protected files unchanged"
    else
      echo "pc55 INTEGRITY FAILURE -- protected files changed:" >&2
      diff "$MANIFEST" "$tmp" | head -40 >&2
      exit 1
    fi
    ;;
  *)
    echo "usage: $0 {snapshot|verify}" >&2; exit 2 ;;
esac
