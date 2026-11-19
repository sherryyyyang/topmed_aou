# Copy to config/paths.sh and edit. Then: source config/paths.sh
# config/paths.sh is gitignored.

# Analysis tree (contains data/, fig*/, scripts/ when running in place)
export TOPMED_OUTLIERS_ROOT="${TOPMED_OUTLIERS_ROOT:-/path/to/pre_Watershed_scripts}"

# Parent folder with ptv_all.tsv and annotation/
export TOPMED_SHEYANG="${TOPMED_SHEYANG:-$(dirname "$TOPMED_OUTLIERS_ROOT")}"

# Where non-55-PC figures / manuscript swaps are written
export TOPMED_SUBMISSION_PACKAGE="${TOPMED_SUBMISSION_PACKAGE:-/path/to/submission_package_nov_2026}"

# Optional explicit inputs
# export TOPMED_PTV_TSV="$TOPMED_SHEYANG/ptv_all.tsv"
# export TOPMED_ACMG_CSV="$TOPMED_SHEYANG/annotation/variant_annotations/acmg_sf_v3.2.csv"

# Paper PC setting: 55 = published; 505 = Nov 2026 plateau delivery
export PC_SETTING="${PC_SETTING:-505}"

# Interpreters
export PAPER_PYTHON="${PAPER_PYTHON:-python3}"
export RSCRIPT="${RSCRIPT:-Rscript}"

# VEP / NMD (only if regenerating Fig3b annotations)
# export VEP_ENV="$HOME/.pixi/envs/ensembl-vep"
# export VEP_CACHE="/path/to/vep_cache"
# export NMD_OUT="$TOPMED_OUTLIERS_ROOT/nmd_escape_vep"
