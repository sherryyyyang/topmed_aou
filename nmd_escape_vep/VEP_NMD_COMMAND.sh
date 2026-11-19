#!/usr/bin/env bash
# Ensembl VEP + NMD plugin for Fig3b PTV unique variants (offline cache 109_GRCh38)
set -euo pipefail
# Override with env vars for off-cluster use; see PATHS.md / config/paths.example.sh
ENV="${VEP_ENV:-/net/wonderland/home/sheyang/.pixi/envs/ensembl-vep}"
CACHE="${VEP_CACHE:-/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/annotation/vep_cache}"
OUT="${NMD_OUT:-${TOPMED_OUTLIERS_ROOT:-/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts}/nmd_escape_vep}"
PLUGIN_DIR="${VEP_PLUGIN_DIR:-$ENV/lib/perl5/site_perl}"
export PATH="$ENV/bin:$PATH"
export PERL5LIB="$ENV/lib/perl5/site_perl:$ENV/lib/perl5/vendor_perl:$ENV/lib/perl5/core_perl:$ENV/lib/perl5/5.32/site_perl:$ENV/lib/perl5/5.32/vendor_perl:$ENV/lib/perl5/5.32/core_perl"
# Input VCF uses non-chr contigs to match cache; ptv_all.tsv is chr-prefixed and stripped when building VCF.
vep -i "$OUT/fig3b_ptv_unique.vcf" \
  --offline --cache --dir_cache "$CACHE" --cache_version 109 \
  --assembly GRCh38 --species homo_sapiens \
  --plugin NMD --dir_plugins "$PLUGIN_DIR" \
  --symbol --canonical --mane_select \
  --fork 8 \
  --force_overwrite -o "$OUT/fig3b_ptv_vep_nmd.txt"
