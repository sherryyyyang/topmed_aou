#!/usr/bin/env bash
# Standalone replacement for the genotype fetch workflow.
# Source per-chromosome VCFs go through bcftools, then plink2 writes bed/bim/fam.
#
#   export VCF_DIR=/path/to/source_vcfs          # chr{1..22,X}.vcf.gz
#   export SAMPLES=/path/to/samples.txt          # one sample id per line
#   export OUT_DIR=$DATA_ROOT/data/genotypes
#   bash scripts/genotypes/vcf_to_plink.sh
#
# Writes:
#   $OUT_DIR/vcf/chr*.vcf.gz          PASS or missing FILTER, tabix index
#   $OUT_DIR/plink/chr*.{bed,bim,fam} from those VCFs
#
# bcftools rules, in order:
#   samples in $SAMPLES
#   SNP, or indel with |ILEN| < 50
#   minor allele count >= 1
#   FILTER is PASS or '.'

set -euo pipefail

if [ -z "${VCF_DIR:-}" ] || [ -z "${SAMPLES:-}" ] || [ -z "${OUT_DIR:-}" ]; then
  echo "Set VCF_DIR, SAMPLES, and OUT_DIR (see header)." >&2
  exit 1
fi
if [ ! -s "$SAMPLES" ]; then
  echo "Sample list not found: $SAMPLES" >&2
  exit 1
fi

command -v bcftools >/dev/null
command -v tabix >/dev/null
command -v plink2 >/dev/null

vcf_out="${OUT_DIR}/vcf"
plink_out="${OUT_DIR}/plink"
mkdir -p "$vcf_out" "$plink_out"

chroms=( $(seq 1 22) X )
for CHR in "${chroms[@]}"; do
  in_vcf="${VCF_DIR}/chr${CHR}.vcf.gz"
  if [ ! -f "$in_vcf" ]; then
    echo "Missing input VCF: $in_vcf" >&2
    exit 1
  fi

  subset="${vcf_out}/chr${CHR}.subset.vcf.gz"
  pass="${vcf_out}/chr${CHR}.vcf.gz"

  echo "[chr${CHR}] subset samples, SNPs/indels, minor AC >= 1"
  bcftools view --types snps,indels --samples-file "$SAMPLES" -Ou \
    -i 'TYPE="snp" || (ILEN>-50 && ILEN<50)' "$in_vcf" \
    | bcftools view --min-ac 1:minor -Oz -o "$subset"

  echo "[chr${CHR}] keep FILTER PASS or missing"
  bcftools view -f 'PASS,.' -Oz -o "$pass" "$subset"
  tabix -f -p vcf "$pass"
  rm -f "$subset"

  echo "[chr${CHR}] plink2 --make-bed"
  plink2 --vcf "$pass" --make-bed --out "${plink_out}/chr${CHR}"
done

echo "Done. PLINK prefixes: ${plink_out}/chr{1..22,X}"
