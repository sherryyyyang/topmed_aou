#!/usr/bin/env bash
# Standalone structural-variant processing.
# Concatenate per-chromosome PASS SV VCFs, then write one table per subject
# of alternate-carrying SVs with AF < 0.01 and F_MISSING < 0.1.
#
#   export SV_VCF_DIR=/path/to/vcfs-pass-variants-only   # chr*.vcf.gz
#   export SV_SUBJECTS=/path/to/sv_subjects             # one sample id per line
#   export OUT_DIR=$DATA_ROOT/sv/sv_genotypes
#   bash scripts/genotypes/process_svs.sh
#
# Writes:
#   $OUT_DIR/svs_pass.vcf.gz
#   $OUT_DIR/individual_svs/<sample>
#     columns: CHROM POS END ID ALT SVLEN AF
#
# Figure 3b DEL/DUP/INV labels are not built here. They are built by
# scripts/fig3_ptv_vs_rna/build_fig3b_segdup_svtype_genes.R from the same
# per-chromosome PASS SV VCFs (any AF; primary call = largest |SVLEN|).

set -euo pipefail

if [ -z "${SV_VCF_DIR:-}" ] || [ -z "${SV_SUBJECTS:-}" ] || [ -z "${OUT_DIR:-}" ]; then
  echo "Set SV_VCF_DIR, SV_SUBJECTS, and OUT_DIR (see header)." >&2
  exit 1
fi
if [ ! -s "$SV_SUBJECTS" ]; then
  echo "Subject list not found: $SV_SUBJECTS" >&2
  exit 1
fi
command -v bcftools >/dev/null

mkdir -p "$OUT_DIR/individual_svs"
list="${OUT_DIR}/vcf_list"
find "$SV_VCF_DIR" -maxdepth 1 -type f -name 'chr*.vcf.gz' | sort -V > "$list"
if [ ! -s "$list" ]; then
  echo "No chr*.vcf.gz under $SV_VCF_DIR" >&2
  exit 1
fi

concat="${OUT_DIR}/svs_pass.vcf.gz"
echo "Concatenating PASS SV VCFs -> $concat"
bcftools concat --file-list "$list" --output "$concat"
bcftools index -f -t "$concat"

echo "Per-subject SV tables (AF < 0.01, F_MISSING < 0.1, GT = alt)"
while IFS= read -r sample || [ -n "$sample" ]; do
  [ -z "$sample" ] && continue
  out="${OUT_DIR}/individual_svs/${sample}"
  bcftools view --include 'AF<0.01 && F_MISSING<0.1' -Ou -s "$sample" "$concat" \
    | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%ALT\t%SVLEN\t%AF\n' --include "GT=\"alt\"" \
    > "$out"
done < "$SV_SUBJECTS"

echo "Done. Per-subject tables: ${OUT_DIR}/individual_svs/"
