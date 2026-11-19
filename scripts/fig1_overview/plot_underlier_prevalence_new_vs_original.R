#!/usr/bin/env Rscript
# Per-gene RNA-outlier (underlier) prevalence, strict z<-3: NEW dataset
# (call-RunRnaUnderlier, n=8,895) vs the ORIGINAL TOPMed dataset (n=5,489).
# Joins on the UNVERSIONED Ensembl id because the two runs use different
# GENCODE versions (e.g. ENSG00000168903.9 vs .8).
suppressPackageStartupMessages({
  library(tidyverse); library(scales); library(patchwork); library(ggrepel)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG1C

NEW  <- "call-RunRnaUnderlier/glob-b354ee3aa408dbeefbcd9f3f19b78ae6/rna_outlier_prevalence_per_gene_z_-3.tsv"
ORIG <- "underlier_pipeline_out/rna_outlier_prevalence_per_gene_z_-3.tsv"
OUT  <- "underlier_prevalence_new_vs_original_z_-3"

rd <- function(f, tag) read_tsv(f, show_col_types = FALSE) |>
  transmute(gene_nv, symbol,
            !!paste0("n_", tag)    := n_underlier_subjects,
            !!paste0("N_", tag)    := n_subjects_tested,
            !!paste0("prev_", tag) := rna_outlier_prevalence)

new <- rd(NEW, "new"); orig <- rd(ORIG, "orig")
n_new_tested  <- max(new$N_new);  n_orig_tested <- max(orig$N_orig)

# NOTE: these tables list only genes with >=1 underlier, so a gene missing from
# one table has prevalence 0 there (or was not tested). A log-log scatter can
# only show genes non-zero in BOTH; the one-sided genes are counted separately.
both <- inner_join(dplyr::select(new,  gene_nv, symbol, n_new,  prev_new),
                   dplyr::select(orig, gene_nv, symbol, n_orig, prev_orig),
                   by = "gene_nv", suffix = c("_new", "_orig")) |>
  mutate(symbol = coalesce(symbol_new, symbol_orig),
         ratio  = prev_new / prev_orig)

only_new  <- anti_join(new,  orig, by = "gene_nv")
only_orig <- anti_join(orig, new,  by = "gene_nv")

cat(sprintf("new: %s genes with >=1 underlier (N=%s subjects)\n", comma(nrow(new)),  comma(n_new_tested)))
cat(sprintf("orig: %s genes with >=1 underlier (N=%s subjects)\n", comma(nrow(orig)), comma(n_orig_tested)))
cat(sprintf("shared (non-zero in both): %s | only new: %s | only original: %s\n",
            comma(nrow(both)), comma(nrow(only_new)), comma(nrow(only_orig))))

sp <- cor(both$prev_new, both$prev_orig, method = "spearman")
pe <- cor(log10(both$prev_new), log10(both$prev_orig), method = "pearson")
pv <- cor.test(both$prev_new, both$prev_orig, method = "spearman")$p.value
cat(sprintf("Spearman rho=%.3f (p=%.2g) | Pearson on log10=%.3f\n", sp, pv, pe))
cat(sprintf("median prevalence: new=%.2e  orig=%.2e  (new/orig ratio median=%.2f)\n",
            median(both$prev_new), median(both$prev_orig), median(both$ratio)))
write_tsv(both, paste0(OUT, "_merged.tsv"))

# label the extremes: highest in each dataset + biggest off-diagonal both ways
lab <- bind_rows(
  slice_max(both, prev_new,  n = 6),
  slice_max(both, prev_orig, n = 6),
  slice_max(both, ratio,     n = 4),
  slice_min(both, ratio,     n = 4)) |> distinct(gene_nv, .keep_all = TRUE)

lims <- range(c(both$prev_new, both$prev_orig))

pA <- ggplot(both, aes(prev_orig, prev_new)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_point(alpha = 0.45, size = 1.5, color = "#2166AC") +
  geom_text_repel(data = lab, aes(label = symbol), size = 2.7, color = "grey15",
                  min.segment.length = 0, segment.size = 0.2, segment.color = "grey55",
                  box.padding = 0.4, max.overlaps = Inf, seed = 5) +
  scale_x_log10(labels = label_percent(accuracy = 0.001), limits = lims) +
  scale_y_log10(labels = label_percent(accuracy = 0.001), limits = lims) +
  annotation_logticks(sides = "bl", size = 0.25) +
  coord_fixed() +
  labs(x = sprintf("Original TOPMed prevalence  (n=%s subjects)", comma(n_orig_tested)),
       y = sprintf("New dataset prevalence  (n=%s subjects)", comma(n_new_tested)),
       title = "A. Per-gene underlier prevalence, new vs original",
       subtitle = sprintf("%s genes non-zero in both | Spearman rho=%.2f, Pearson(log10)=%.2f | dashed = y=x",
                          comma(nrow(both)), sp, pe)) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.6, color = "grey30"))

# panel B: where the genes live (shared vs dataset-specific)
ov <- tibble(cls = factor(c("non-zero in both", "only new", "only original"),
                          levels = c("non-zero in both", "only new", "only original")),
             n = c(nrow(both), nrow(only_new), nrow(only_orig)))
pB <- ggplot(ov, aes(cls, n, fill = cls)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = comma(n)), vjust = -0.4, size = 3.1) +
  scale_fill_manual(values = c("#2166AC", "#D95F02", "grey60")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "genes",
       title = "B. Gene overlap",
       subtitle = "genes with >=1 underlier in each dataset") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.6, color = "grey30"),
        axis.text.x = element_text(angle = 20, hjust = 1))

p <- (pA | pB) + plot_layout(widths = c(2.1, 1)) +
  plot_annotation(
    title = "RNA-outlier prevalence per gene (strict z<-3 + log2CPM): new dataset vs original TOPMed",
    subtitle = str_wrap(sprintf(
      "Both axes are the fraction of tested subjects called an underlier for that gene. Joined on unversioned Ensembl id. The new dataset yields far fewer underlier genes (%s vs %s) and a lower ceiling (max %.2f%% vs %.2f%%), so most genes sit below y=x; median new/original prevalence ratio = %.2f. Genes absent from a table have zero underliers there and cannot appear on log axes - they are counted in panel B.",
      comma(nrow(new)), comma(nrow(orig)),
      100*max(new$prev_new), 100*max(orig$prev_orig), median(both$ratio)), width = 150),
    theme = theme(plot.title = element_text(face = "bold", size = 12.5),
                  plot.subtitle = element_text(size = 8.8, color = "grey25")))

ggsave(paste0(OUT, ".png"), p, width = 12.2, height = 6.2, dpi = 200)
ggsave(paste0(OUT, ".pdf"), p, width = 12.2, height = 6.2)
message("wrote ", OUT, ".{png,pdf} + _merged.tsv")
