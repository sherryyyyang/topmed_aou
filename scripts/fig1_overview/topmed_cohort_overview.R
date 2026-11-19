ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))
OUT_DIR <- FIG1A
library(vcfR)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel) 
library(vroom)
library(ggpattern)
library(tidyverse)
library(magrittr)

ancestry <- vroom("/net/topmed11/working/porchard/for-others/sherry/2026-03-03/ancestry.9b.txt") |> filter(is_in(NWD_ID, rawdata$SV_subjects)) |> left_join(counts_df)
ancestry$pan_EAS <- ancestry$EAS + ancestry$AMR
ancestry$pan_WES <- ancestry$EUR + ancestry$SAS + ancestry$MES
ancestry_mat <- select(ancestry, -SV_Type, -Count) |> distinct() |> column_to_rownames("NWD_ID") |> as.matrix()
ancestry_max <- apply(ancestry_mat[rawdata$SV_subjects,], 1, max)
ancestry
test()

ancestry_eur <- filter(ancestry, EUR > 0.9) |> group_by(SV_Type) |> summarise(median_n_sv = median(Count), mad_n_sv = mad(Count)) 
ancestry_noneur <- filter(ancestry, EUR <= 0.9) |> group_by(SV_Type) |> summarise(median_n_sv = median(Count), mad_n_sv = mad(Count)) 

rawdata<- vroom("./metadata.filtered_with_sv.txt") 
rawdata <- rawdata[!is.na(rawdata$SV_subjects),]
rawdata$age_bins <- cut(as.numeric(rawdata$AGE), breaks = seq(0, 90, by = 10), right = FALSE, labels = c("0-9", "10-19", "20-29","30-39","40-49","50-59","60-69","70-79","80-89"))
plot_data <- rawdata %>%
  group_by(age_bins, COHORT, SEX) %>%
  summarise(Count = n(), .groups = "drop") 

# 1. LOAD DATA
# Replace 'your_variants.vcf' with your actual VCF file path
vcf <- read.vcfR("~/topmed11/sv_genotypes/svs_af0.01_missing0.1.vcf.gz")

# 2. EXTRACT INFORMATION
# Extract the fixed info (to get SV Type) and genotypes (to identify carriers)
fix_data <- as.data.frame(getFIX(vcf))
gt_data <- extract.gt(vcf) # Matrix of samples (columns) x variants (rows)

# Create a mapping of Variant Index to SV Type
# Assuming SVTYPE is in the INFO column
sv_types <- fix_data$ALT %>% 
  str_replace_all("<|>", "")
gt_matrix <- extract.gt(vcf)

counts_df <- as.data.frame(gt_matrix) %>%
  mutate(SV_Type = sv_types) %>%
  pivot_longer(cols = -SV_Type, names_to = "Subject", values_to = "GT") %>%
  # Filter for variants where the subject has at least one ALT allele
  # This handles both "0/1" and "0|1" styles shown in your data
  filter(!is.na(GT) & GT != "0/0" & GT != "0|0" & GT != "./." & GT != "0") %>%
  group_by(Subject, SV_Type) %>%
  summarise(Count = n(), .groups = 'drop')
colnames(counts_df)[1] <- "NWD_ID"


p <- ggplot(plot_data, aes(x = age_bins, y = Count, fill = COHORT, pattern = SEX)) +
  geom_bar_pattern(
    stat = "identity",
    position = "stack",
    pattern_fill = "black",
    pattern_density = 0.1,
    pattern_spacing = 0.02
  ) +
  scale_pattern_manual(values = c(Male = "stripe", Female = "none")) +
  labs(x = "Age Group", y = "Count", fill = "Cohort", pattern = "Sex") +
  ggtitle("TOPMed Cohort Summary") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.line = element_line(size = 1.5, color = "black")  # thicker axis lines
  )

ggsave(filename = "TOPMed_cohort_summary_plot.png", plot = p, width = 10, height = 6, dpi = 300) 


# Load necessary libraries
library(ggplot2)
library(patchwork) # For side-by-side juxtaposition

# 1. DATA
all_counts  <- c(Deletion = 47686, Duplication = 8072, Inversion = 8186)
rare_counts <- c(Deletion = 39533, Duplication = 7975, Inversion = 7795)

df_all <- data.frame(
  Variant = factor(names(all_counts), levels = c("Deletion", "Duplication", "Inversion")),
  Count   = all_counts
)

df_rare <- data.frame(
  Variant = factor(names(rare_counts), levels = c("Deletion", "Duplication", "Inversion")),
  Count   = rare_counts
)

# 2. SHARED VISUALS
max_y_limit <- 55000 
palette <- c("Deletion" = "#A65A5A", "Duplication" = "#E9E5B5", "Inversion" = "#6E8B6E")

# Updated theme to specifically make axis lines visible
fixed_axis_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    # This makes the X and Y axis lines solid black
    axis.line = element_line(color = "black", linewidth = 0.8),
    # This adds small ticks on the axes
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# 3. GENERATE PLOTS
p1 <- ggplot(df_all, aes(x = Variant, y = Count, fill = Variant)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = paste0("N=", Count)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = palette) +
  scale_y_continuous(limits = c(0, max_y_limit), expand = c(0, 0)) + # expand=c(0,0) puts bars flush against X-axis
  labs(title = "All Structural Variants", x = "Variant Type", y = "Number of Counts (N)") +
  fixed_axis_theme

p2 <- ggplot(df_rare, aes(x = Variant, y = Count, fill = Variant)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = paste0("N=", Count)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = palette) +
  scale_y_continuous(limits = c(0, max_y_limit), expand = c(0, 0)) +
  labs(title = "Rare Variants (AF < 0.01)", x = "Variant Type", y = "") + 
  fixed_axis_theme

# 4. COMBINE
final_plot <- p1 + p2 + 
  plot_annotation(title = 'Structural Variant Frequency Comparison')

ggsave(filename = "SV_rare_sv_distribution.png", plot = final_plot, width = 10, height = 6, dpi = 300) 



# Load necessary libraries
library(vcfR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

# 1. LOAD DATA
# Replace 'your_variants.vcf' with your actual VCF file path
vcf <- read.vcfR("~/topmed11/sv_genotypes/svs_af0.01_missing0.1.vcf.gz")

# 2. EXTRACT INFORMATION
# Extract the fixed info (to get SV Type) and genotypes (to identify carriers)
fix_data <- as.data.frame(getFIX(vcf))
gt_data <- extract.gt(vcf) # Matrix of samples (columns) x variants (rows)

# Create a mapping of Variant Index to SV Type
# Assuming SVTYPE is in the INFO column
sv_types <- fix_data$ALT %>% 
  str_replace_all("<|>", "")
gt_matrix <- extract.gt(vcf)

counts_df <- as.data.frame(gt_matrix) %>%
  mutate(SV_Type = sv_types) %>%
  pivot_longer(cols = -SV_Type, names_to = "Subject", values_to = "GT") %>%
  # Filter for variants where the subject has at least one ALT allele
  # This handles both "0/1" and "0|1" styles shown in your data
  filter(!is.na(GT) & GT != "0/0" & GT != "0|0" & GT != "./." & GT != "0") %>%
  group_by(Subject, SV_Type) %>%
  summarise(Count = n(), .groups = 'drop')

library(ggplot2)

# 1. SETUP THE PLOT
p <- ggplot(counts_df, aes(x = SV_Type, y = Count, fill = SV_Type)) +
  geom_boxplot(outlier.shape = 16, outlier.alpha = 0.5, linewidth = 1.2) + # Thicker lines for the box
  scale_y_log10(breaks = c(1, 10, 100, 1000), labels = scales::comma) +
  scale_fill_manual(values = c("DEL"="#A65A5A", "INS"="#E9E5B5", "DUP"="#6E8B6E", "INV"="#4682B4")) +
  
  # 2. LABELS
  labs(
    title = "Rare SV Count per Individual",
    x = "SV Type",
    y = "Count per Individual (Log Scale)"
  ) +
  
  # 3. MEGA THEME
  theme_classic(base_size = 20) + # Increases all base text to 20pt
  theme(
    # Plot Title: Huge and Bold
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    
    # Axis Titles: Large and Bold
    axis.title.x = element_text(size = 26, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 26, face = "bold", margin = margin(r = 15)),
    
    # Axis Tick Labels: Large
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    
    # Legend: (Optional, if you decide to show it)
    legend.position = "none",
    
    # Make axis lines thicker to match the big text
    axis.line = element_line(linewidth = 1.5),
    axis.ticks = element_line(linewidth = 1.5, size = 1)
  )

ggsave(filename = "rare_sv_types_per_subject.png", plot = p, width = 14, height = 10, dpi = 300) 

