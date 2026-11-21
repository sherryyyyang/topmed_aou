#!/usr/bin/env Rscript
# Fig 3b with SegDup SV (DEL/DUP/INV) and CNV/multi-copy genes removed;
# uniform grey points, log–log axes, Spearman on remaining genes (n_ptv>0).
# Does not replace main Figure3b.
#
# Usage (default PC_SETTING=1055):
#   Rscript scripts/fig3_ptv_vs_rna/export_fig3b_excl_sv_cnv.R

.PKG <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/submission_package_nov_2026"
.ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
panel_script <- file.path(.PKG, "scripts/fig3_ptv_vs_rna/plot_ptv_vs_underlier_gene_level_2x2.R")
stopifnot(file.exists(panel_script))

pc <- Sys.getenv("PC_SETTING", unset = "1055")
Sys.setenv(
  PC_SETTING = pc,
  TOPMED_OUTLIERS_ROOT = .ROOT,
  EXCLUDE_SV_CNV = "1"
)

cat("Running Fig3b panel (EXCLUDE_SV_CNV=1, PC_SETTING=", pc, ") ...\n", sep = "")
status <- system2("Rscript", panel_script, stdout = "", stderr = "")
if (!identical(status, 0L)) {
  stop("plot_ptv_vs_underlier_gene_level_2x2.R failed with status ", status)
}

panel_stem <- file.path(
  .ROOT, paste0("pc", pc, "_work"),
  "fig3_ptv_vs_rna", "panel_a_gene_level",
  "ptv_vs_underlier_gene_level_2x2_haplo_z2_z3_z4_excl_sv_cnv"
)
for (ext in c("png", "pdf")) {
  src <- paste0(panel_stem, ".", ext)
  if (!file.exists(src)) stop("Missing panel output: ", src)
}

out_dir <- file.path(.PKG, "figures/main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dest_stem <- file.path(out_dir, "Figure3b_excl_sv_cnv")
file.copy(paste0(panel_stem, ".png"), paste0(dest_stem, ".png"), overwrite = TRUE)
file.copy(paste0(panel_stem, ".pdf"), paste0(dest_stem, ".pdf"), overwrite = TRUE)

sp <- file.path(.PKG, "source_panels")
dir.create(sp, recursive = TRUE, showWarnings = FALSE)
file.copy(paste0(dest_stem, ".png"), file.path(sp, "Fig3b_excl_sv_cnv.png"), overwrite = TRUE)
file.copy(paste0(dest_stem, ".pdf"), file.path(sp, "Fig3b_excl_sv_cnv.pdf"), overwrite = TRUE)

cat("Wrote:", paste0(dest_stem, ".png"), "\n")
cat("Wrote:", paste0(dest_stem, ".pdf"), "\n")
cat("Panel source:", paste0(panel_stem, ".png"), "\n")
