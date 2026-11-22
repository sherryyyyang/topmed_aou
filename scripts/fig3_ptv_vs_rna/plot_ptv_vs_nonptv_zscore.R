#!/usr/bin/env Rscript
# Sanity check: do PTV carriers have lower expression z-scores than non-carriers?
# Pool over all (subject x gene) pairs, within PTV-carrying genes. z-scores in
# expr_z_join are per-gene standardized (mean 0, sd 1), so non-carriers center on
# 0 and PTV carriers should shift negative / be enriched in the low tail.
# Cache both AF<0.5% and AF<0.1%; the paper panel is AF<0.1% only.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(ggplot2); library(scales); library(patchwork)
})

ROOT <- Sys.getenv("TOPMED_OUTLIERS_ROOT", unset = "")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("DATA_ROOT", unset = "")
if (!nzchar(ROOT)) stop("Set TOPMED_OUTLIERS_ROOT or DATA_ROOT to the analysis directory.")
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))
OUT_DIR <- FIG3B
PW <- ROOT
CACHE <- file.path(FIG3B, "ptv_vs_nonptv_zscore_cache_v2.rds")
CUTOFFS <- c(`AF<0.5%`=0.005, `AF<0.1%`=0.001)

kept <- read_rds(KEPT_SUBJECTS); n_kept <- length(kept)

if (file.exists(CACHE)) {
  cat("Loading cache:", CACHE, "\n"); res <- read_rds(CACHE)
} else {
  ptv_all <- read_tsv(PTV_TSV, show_col_types=FALSE) |>
    filter(sample %in% kept, lof=="HC") |>
    mutate(gene_nv = str_remove(gene, "\\..*$"))
  cat("Loading expr_z_join.rds ...\n")
  z <- read_rds(EXPR_Z_JOIN) |>
    filter(subjectid %in% kept) |>
    mutate(gene_nv = str_remove(gene_id, "\\..*$")) |>
    select(gene_nv, subjectid, expression_zscore)

  build <- function(af_max) {
    ptv <- ptv_all |> filter(af < af_max) |>
      distinct(gene_nv, subjectid = sample) |> mutate(is_ptv=TRUE)
    ptv_genes <- unique(ptv$gene_nv)
    zg <- z |> filter(gene_nv %in% ptv_genes) |>
      left_join(ptv, by=c("gene_nv","subjectid")) |>
      mutate(grp = ifelse(coalesce(is_ptv, FALSE), "PTV carrier", "non-carrier"))
    summ <- zg |> group_by(grp) |>
      summarise(n=n(), mean_z=mean(expression_zscore), median_z=median(expression_zscore),
                pct_z_lt_m2=100*mean(expression_zscore < -2),
                pct_z_lt_m3=100*mean(expression_zscore < -3),
                pct_z_lt_m4=100*mean(expression_zscore < -4), .groups="drop")
    set.seed(1)
    pcz <- zg |> filter(grp=="PTV carrier") |> pull(expression_zscore)
    ncz <- zg |> filter(grp=="non-carrier") |> slice_sample(n=1e6) |> pull(expression_zscore)
    wt_p <- wilcox.test(c(pcz,ncz) ~ c(rep(TRUE,length(pcz)), rep(FALSE,length(ncz))))$p.value

    # ---- gene-level (independent-unit) normalization ----
    # one row per gene: mean z of its PTV carriers vs its non-carriers.
    gene_stats <- zg |> group_by(gene_nv) |>
      summarise(n_ptv = sum(grp=="PTV carrier"),
                n_non = sum(grp=="non-carrier"),
                mz_ptv = mean(expression_zscore[grp=="PTV carrier"]),
                mz_non = mean(expression_zscore[grp=="non-carrier"]),
                .groups="drop") |>
      filter(n_ptv >= 1, n_non >= 5) |>
      mutate(diff = mz_ptv - mz_non)
    # paired tests across genes (H0: PTV carriers not lower than non-carriers)
    gl_wilcox <- wilcox.test(gene_stats$diff, alternative="less")$p.value      # signed-rank
    gl_t      <- t.test(gene_stats$diff, alternative="less")$p.value
    gl_sign_p <- binom.test(sum(gene_stats$diff < 0), nrow(gene_stats), 0.5,
                            alternative="greater")$p.value
    gene_level <- list(n_genes=nrow(gene_stats),
                       median_diff=median(gene_stats$diff),
                       mean_diff=mean(gene_stats$diff),
                       pct_genes_lower=100*mean(gene_stats$diff < 0),
                       wilcox_p=gl_wilcox, t_p=gl_t, sign_p=gl_sign_p)

    list(ptv_carrier_z=pcz, non_sample_z=ncz, summ=summ, wt_p=wt_p,
         n_genes=length(ptv_genes), n_ptv=sum(zg$grp=="PTV carrier"),
         n_non=sum(zg$grp=="non-carrier"),
         gene_level=gene_level, gene_stats=gene_stats)
  }
  res <- lapply(CUTOFFS, build)
  write_rds(res, CACHE, compress="gz")
}

COL <- c("PTV carrier"=COL_PTV,"non-carrier"=COL_NON)

make_row <- function(r, aflab) {
  summ <- r$summ
  med_ptv <- summ$median_z[summ$grp=="PTV carrier"]; med_non <- summ$median_z[summ$grp=="non-carrier"]
  zp <- bind_rows(tibble(grp="PTV carrier", z=r$ptv_carrier_z),
                  tibble(grp="non-carrier", z=r$non_sample_z))
  fold <- tibble(threshold=factor(c("z<-2","z<-3","z<-4"), levels=c("z<-2","z<-3","z<-4")),
                 ptv=c(summ$pct_z_lt_m2[summ$grp=="PTV carrier"],summ$pct_z_lt_m3[summ$grp=="PTV carrier"],summ$pct_z_lt_m4[summ$grp=="PTV carrier"]),
                 non=c(summ$pct_z_lt_m2[summ$grp=="non-carrier"],summ$pct_z_lt_m3[summ$grp=="non-carrier"],summ$pct_z_lt_m4[summ$grp=="non-carrier"])) |>
    mutate(fold=ptv/non)
  fl <- fold |> pivot_longer(c(ptv,non), names_to="g", values_to="pct") |>
    mutate(grp=ifelse(g=="ptv","PTV carrier","non-carrier"))

  pA <- ggplot(zp, aes(z, fill=grp, color=grp)) +
    geom_density(alpha=0.35, linewidth=0.6, adjust=1.2) +
    geom_vline(xintercept=c(-2,-3,-4), linetype="dotted", color="grey60", linewidth=0.3) +
    scale_fill_manual(values=COL, name=NULL, guide="none") +
    scale_color_manual(values=COL, name=NULL, guide="none") +
    coord_cartesian(xlim=c(-6,4)) +
    labs(title=sprintf("%s  density", aflab),
         x="expression z-score", y="density") +
    theme_paper() +
    theme_paper_inner_title() +
    theme(legend.position = "none")

  # Explicit ECDF: stat_ecdf() emits y = 0 at the left edge, which is -Inf
  # once the cumulative axis is logged. A log axis gives each decade of the
  # left tail equal height, so the carrier excess stays visible all the way
  # down to z < -4 while both curves still run up to 100%.
  cdf <- zp |>
    arrange(grp, z) |>
    group_by(grp) |>
    mutate(cum = row_number() / n()) |>
    group_by(grp, z) |>
    slice_tail(n = 1) |>
    ungroup() |>
    filter(cum > 0)
  cdf_brks <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)

  pB <- ggplot(cdf, aes(z, cum, color=grp)) +
    geom_step(linewidth=0.8, direction="hv") +
    geom_vline(xintercept=c(-2,-3,-4), linetype="dotted", color="grey60", linewidth=0.3) +
    scale_color_manual(values=COL, name=NULL, guide="none") +
    coord_cartesian(xlim=c(-6,4)) +
    scale_y_log10(breaks=cdf_brks,
                  labels=c("0.001%","0.01%","0.1%","1%","10%","100%"),
                  expand=expansion(mult=c(0.03, 0.05))) +
    annotation_logticks(sides="l", linewidth=0.25, short=unit(0.04,"cm"),
                        mid=unit(0.07,"cm"), long=unit(0.11,"cm")) +
    labs(title="Cumulative distribution (log axis)",
         x="expression z-score", y="cumulative \u2264 z") +
    theme_paper() +
    theme_paper_inner_title() +
    theme(legend.position = "none")

  pC <- ggplot(fl, aes(threshold, pct, fill=grp)) +
    geom_col(position=position_dodge(0.7), width=0.62) +
    geom_text(
      data = fold,
      aes(threshold, ptv, label = sprintf("%.0f\u00d7", fold)),
      inherit.aes = FALSE, vjust = -0.40, fontface = "bold", size = 2.8,
      color = COL_PTV, family = PAPER_FONT, nudge_x = 0.18, show.legend = FALSE
    ) +
    scale_fill_manual(values=COL, name=NULL) +
    scale_y_continuous(
      labels = label_percent(scale = 1),
      expand = expansion(mult = c(0, 0.22))
    ) +
    coord_cartesian(clip = "off") +
    labs(title="Outlier rate (fold-enrichment)",
         x=NULL, y="% of pairs") +
    theme_paper() +
    theme_paper_inner_title() +
    theme(legend.position = "none",
          plot.title = element_text(margin = margin(b = 6)),
          plot.margin = margin(8, 6, 4, 8))
  list(pA=pA,pB=pB,pC=pC, fold=fold, wt_p=r$wt_p, n_ptv=r$n_ptv, n_non=r$n_non, n_genes=r$n_genes)
}

rows <- Map(make_row, res, names(res))
for (nm in names(res)) {
  gl <- res[[nm]]$gene_level
  cat(sprintf("\n== %s ==  PTV pairs=%s non=%s genes=%s  pooled-MW p=%.1e\n",
      nm, comma(rows[[nm]]$n_ptv), comma(rows[[nm]]$n_non), comma(rows[[nm]]$n_genes), rows[[nm]]$wt_p))
  print(as.data.frame(rows[[nm]]$fold), digits=4)
  cat(sprintf("  GENE-LEVEL (n=%s independent genes): median(meanZ_PTV - meanZ_non)=%.3f | %.1f%% of genes have PTV carriers lower | signed-rank p=%.2e | paired-t p=%.2e | sign-test p=%.2e\n",
      comma(gl$n_genes), gl$median_diff, gl$pct_genes_lower, gl$wilcox_p, gl$t_p, gl$sign_p))
}
fig <- (rows[["AF<0.1%"]]$pA | rows[["AF<0.1%"]]$pB | rows[["AF<0.1%"]]$pC) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig <- fig + theme(plot.margin = margin(t = 8, r = 4, b = 1, l = 10))
.paper_plots <- list(b = fig)
paper_ggsave(fig, file.path(OUT_DIR, "ptv_vs_nonptv_zscore"),
             width = NATURE_DOUBLE_IN, height = 3.6)
if (!is_paper_assemble()) {
  bind_rows(lapply(names(res), function(nm) mutate(res[[nm]]$summ, cutoff=nm))) |>
    write_tsv(file.path(OUT_DIR, "ptv_vs_nonptv_zscore_summary.tsv"))
  bind_rows(lapply(names(res), function(nm) as_tibble(res[[nm]]$gene_level) |> mutate(cutoff=nm))) |>
    write_tsv(file.path(OUT_DIR, "ptv_vs_nonptv_zscore_genelevel.tsv"))
  bind_rows(lapply(names(res), function(nm) mutate(res[[nm]]$gene_stats, cutoff=nm))) |>
    write_tsv(file.path(OUT_DIR, "ptv_vs_nonptv_zscore_per_gene.tsv"))
  cat("\nWrote ptv_vs_nonptv_zscore\n")
}
