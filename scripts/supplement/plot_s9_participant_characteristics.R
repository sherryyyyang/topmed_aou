#!/usr/bin/env Rscript
# Supplementary Figure 9: participant characteristics of the TOPMed analysis
# cohort (the 5,489 QC-kept subjects used for every underlier call).
#
#   a  age distribution (subjects with an age on record)
#   b  sex
#   c  contributing TOPMed study
#   d  sequencing centre
#   e  predominant genetic ancestry
#   f  mean genetic ancestry proportions
#
# Panels e/f are GENETIC ANCESTRY proportions from ancestry.9b.txt, not
# self-reported race or ethnicity; no self-reported field exists in either
# metadata file on disk.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(scales); library(patchwork)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
source(file.path(ROOT, "scripts/helpers/paper_theme.R"))

ANCESTRY_TSV <- "/net/topmed11/working/porchard/for-others/sherry/2026-03-03/ancestry.9b.txt"
OUT_DIR <- file.path(FINAL_DRAFT, "supp_figs")
TAB_DIR <- file.path(OUT_DIR, "tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

kept <- read_rds(KEPT_SUBJECTS)
n_kept <- length(kept)
cat("kept subjects:", n_kept, "\n")

meta <- read_tsv(METADATA, show_col_types = FALSE) |>
  filter(WGS %in% kept) |>
  distinct(WGS, .keep_all = TRUE) |>
  mutate(AGE = suppressWarnings(as.numeric(AGE)))
cat("metadata rows matched to kept subjects:", nrow(meta), "\n")

# ---------------------------------------------------------------- a. age ----
age_ok <- meta |> filter(!is.na(AGE))
n_age <- nrow(age_ok)
n_age_missing <- n_kept - n_age
age_med <- median(age_ok$AGE)
missing_studies <- meta |>
  filter(is.na(AGE)) |>
  count(COHORT, sort = TRUE)
cat("age available:", n_age, " missing:", n_age_missing,
    " median:", age_med, "\n")
print(missing_studies)

p_age <- ggplot(age_ok, aes(AGE)) +
  geom_histogram(binwidth = 5, boundary = 0, fill = COL_TOPMED,
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = age_med, linetype = "dashed",
             color = COL_PTV, linewidth = 0.5) +
  annotate("text", x = age_med, y = Inf,
           label = sprintf(" median %.0f y", age_med),
           hjust = 0, vjust = 1.6, size = 2.6, color = COL_PTV,
           family = PAPER_FONT) +
  scale_x_continuous(breaks = pretty_breaks(6)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Age (years)", y = "Subjects", title = "Age") +
  theme_paper() + theme_paper_inner_title()

# ---------------------------------------------------------------- b. sex ----
sex_tbl <- meta |>
  count(SEX, name = "n") |>
  mutate(pct = 100 * n / sum(n)) |>
  arrange(desc(n))

# Counts ride on the category axis; only the percent sits beside the bar, so
# long labels cannot overflow these one-third-width panels.
p_sex <- ggplot(sex_tbl, aes(reorder(sprintf("%s (n=%s)", SEX, comma(n)), n),
                             pct, fill = SEX)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            hjust = -0.15, size = 2.4, family = PAPER_FONT, color = "grey15") +
  scale_fill_manual(values = c(Female = unname(OI["purple"]),
                               Male = unname(OI["sky"]))) +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.26))) +
  coord_flip() +
  labs(x = NULL, y = "% of subjects", title = "Sex") +
  theme_paper() + theme_paper_inner_title()

# -------------------------------------------------------------- c. study ----
study_tbl <- meta |>
  count(COHORT, name = "n") |>
  mutate(pct = 100 * n / sum(n)) |>
  arrange(desc(n))

p_study <- ggplot(study_tbl, aes(reorder(sprintf("%s (n=%s)", COHORT, comma(n)), n),
                                 pct)) +
  geom_col(width = 0.7, fill = unname(OI["green"])) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            hjust = -0.15, size = 2.4, family = PAPER_FONT, color = "grey15") +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.26))) +
  coord_flip() +
  labs(x = NULL, y = "% of subjects", title = "Contributing study") +
  theme_paper() + theme_paper_inner_title()

# ------------------------------------------------------- d. seq. centre ----
ctr_tbl <- meta |>
  count(SEQUENCING_CENTER, name = "n") |>
  mutate(pct = 100 * n / sum(n)) |>
  arrange(desc(n))

p_ctr <- ggplot(ctr_tbl, aes(
  reorder(sprintf("%s (n=%s)", SEQUENCING_CENTER, comma(n)), n), pct)) +
  geom_col(width = 0.6, fill = unname(OI["orange"])) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            hjust = -0.15, size = 2.4, family = PAPER_FONT, color = "grey15") +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.26))) +
  coord_flip() +
  labs(x = NULL, y = "% of subjects", title = "Sequencing centre") +
  theme_paper() + theme_paper_inner_title()

# ----------------------------------------------------- e/f. ancestry ----
ANC <- c("AFR", "AMR", "EAS", "EUR", "MES", "SAS", "OCN")
ANC_COLORS <- c(
  AFR = unname(OI["orange"]),
  AMR = unname(OI["vermillion"]),
  EAS = unname(OI["green"]),
  EUR = unname(OI["blue"]),
  MES = unname(OI["purple"]),
  SAS = unname(OI["sky"]),
  OCN = unname(OI["yellow"])
)

anc <- read_tsv(ANCESTRY_TSV, show_col_types = FALSE) |>
  filter(NWD_ID %in% kept) |>
  distinct(NWD_ID, .keep_all = TRUE)
n_anc <- nrow(anc)
cat("kept subjects with genetic ancestry:", n_anc,
    " missing:", n_kept - n_anc, "\n")

anc_long <- anc |>
  pivot_longer(all_of(ANC), names_to = "ancestry", values_to = "prop")

# Predominant component; below 50% there is no majority ancestry.
anc_major <- anc_long |>
  group_by(NWD_ID) |>
  slice_max(prop, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(group = ifelse(prop >= 0.5, ancestry, "No majority"))

major_tbl <- anc_major |>
  count(group, name = "n") |>
  mutate(pct = 100 * n / sum(n)) |>
  arrange(desc(n))

p_major <- ggplot(major_tbl, aes(
  reorder(sprintf("%s (n=%s)", group, comma(n)), n), pct, fill = group)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            hjust = -0.15, size = 2.4, family = PAPER_FONT, color = "grey15") +
  scale_fill_manual(values = c(ANC_COLORS, `No majority` = "grey60")) +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.26))) +
  coord_flip() +
  labs(x = NULL, y = "% of subjects",
       title = "Predominant genetic ancestry (>50%)") +
  theme_paper() + theme_paper_inner_title()

mean_tbl <- anc_long |>
  group_by(ancestry) |>
  summarise(mean_prop = mean(prop), .groups = "drop") |>
  mutate(pct = 100 * mean_prop) |>
  arrange(desc(pct))

p_mean <- ggplot(mean_tbl, aes(reorder(ancestry, pct), pct, fill = ancestry)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            hjust = -0.15, size = 2.3, family = PAPER_FONT, color = "grey15") +
  scale_fill_manual(values = ANC_COLORS) +
  scale_y_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.34))) +
  coord_flip() +
  labs(x = NULL, y = "Mean ancestry proportion",
       title = "Mean genetic ancestry") +
  theme_paper() + theme_paper_inner_title()

# ------------------------------------------------------------ assemble ----
fig <- (p_age | p_sex) / (p_study | p_ctr) / (p_major | p_mean) +
  plot_annotation(tag_levels = "a") & theme_paper_tag()

paper_ggsave(fig, file.path(OUT_DIR, "Supplementary_Figure9"),
             width = NATURE_DOUBLE_IN, height = 7.6)

panels <- list(
  FigS9a_age = p_age, FigS9b_sex = p_sex, FigS9c_study = p_study,
  FigS9d_sequencing_centre = p_ctr, FigS9e_predominant_ancestry = p_major,
  FigS9f_mean_ancestry = p_mean
)
for (nm in names(panels)) {
  paper_ggsave(panels[[nm]], file.path(OUT_DIR, nm),
               width = NATURE_DOUBLE_IN / 2, height = 2.6)
}

# -------------------------------------------------------- source table ----
src <- bind_rows(
  tibble(section = "Cohort", level = "TOPMed QC-kept subjects",
         n = n_kept, percent = 100, extra = "kept_subjects.rds"),
  sex_tbl |> transmute(section = "Sex", level = SEX, n, percent = pct,
                       extra = NA_character_),
  tibble(section = "Age (years)", level = "Median (IQR)", n = n_age,
         percent = NA_real_,
         extra = sprintf("%.1f (%.1f-%.1f); n_with_age=%d; n_missing=%d",
                         age_med, quantile(age_ok$AGE, 0.25),
                         quantile(age_ok$AGE, 0.75), n_age, n_age_missing)),
  tibble(section = "Age (years)", level = "Mean (SD)", n = n_age,
         percent = NA_real_,
         extra = sprintf("%.2f (%.2f)", mean(age_ok$AGE), sd(age_ok$AGE))),
  study_tbl |> transmute(section = "Study", level = COHORT, n, percent = pct,
                         extra = NA_character_),
  ctr_tbl |> transmute(section = "Sequencing centre", level = SEQUENCING_CENTER,
                       n, percent = pct, extra = NA_character_),
  major_tbl |> transmute(section = "Predominant genetic ancestry", level = group,
                         n, percent = pct,
                         extra = "max component; genetic ancestry, not self-reported"),
  mean_tbl |> transmute(section = "Mean genetic ancestry proportion",
                        level = ancestry, n = n_anc, percent = pct,
                        extra = "mean over subjects with ancestry data"),
  tibble(section = "Self-reported race / ethnicity", level = "Not available",
         n = NA_integer_, percent = NA_real_,
         extra = "absent from metadata.filtered_with_sv.txt and 2025-06-11 metadata.filtered.txt")
)
out_tsv <- file.path(TAB_DIR, "S9_participant_characteristics.tsv")
write_tsv(src, out_tsv)
cat("\nWrote:", out_tsv, "\n")
print(as.data.frame(src), max = 400)
