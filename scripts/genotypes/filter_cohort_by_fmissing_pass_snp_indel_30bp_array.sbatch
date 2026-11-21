#!/usr/bin/env bash
# SLURM array: filter per-chromosome cohort VCFs.
#   1. FILTER=PASS
#   2. F_MISSING < 0.1
#   3. SNP, or indel with |ILEN| <= 30
#
# Inputs:  $COHORT_VCF_DIR/chr{1..22,X}.vcf.gz
# Outputs: $OUTDIR/chr{1..22,X}.vcf.gz (+ tabix index)
#
#   export DATA_ROOT=/path/to/analysis
#   export COHORT_VCF_DIR=$DATA_ROOT/data/genotypes/vcf
#   export OUTDIR=$DATA_ROOT/data/genotypes/vcf/filtered_by_f_missing_01_pass_snp_indel_30bp
#   sbatch scripts/genotypes/filter_cohort_by_fmissing_pass_snp_indel_30bp_array.sbatch

#SBATCH --job-name=cohort-fmissing-filter
#SBATCH --array=1-23%23
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-64}"

if [ -z "${COHORT_VCF_DIR:-}" ] || [ -z "${OUTDIR:-}" ]; then
  echo "Set COHORT_VCF_DIR and OUTDIR (see header)." >&2
  exit 1
fi

mkdir -p "$OUTDIR" logs

if [ "${SLURM_ARRAY_TASK_ID}" -le 22 ]; then
  CHR="${SLURM_ARRAY_TASK_ID}"
else
  CHR="X"
fi

in_vcf="${COHORT_VCF_DIR}/chr${CHR}.vcf.gz"
base="$(basename "$in_vcf")"
out_vcf="${OUTDIR}/${base}"

if [ ! -f "$in_vcf" ]; then
  echo "Missing input VCF: $in_vcf" >&2
  exit 1
fi

if [ -s "$out_vcf" ]; then
  echo "Already exists (skip): $out_vcf"
  exit 0
fi

echo "[chr${CHR}] Filtering cohort VCF -> $out_vcf"

bcftools view -f PASS "$in_vcf" -r "chr${CHR}" -Ou --threads "$THREADS" \
  | bcftools view -i 'F_MISSING<0.1' -Ou --threads "$THREADS" \
  | bcftools view -i 'TYPE=="snp" || (ILEN <= 30 && ILEN >= -30)' -Oz -o "$out_vcf" --threads "$THREADS"

tabix -f -p vcf "$out_vcf"

echo "[chr${CHR}] Done"
