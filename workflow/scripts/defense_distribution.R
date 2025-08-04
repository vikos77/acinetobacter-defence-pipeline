# Defence System Distribution Analysis: DefenseFinder and PADLOC

# Analysis of defence system distribution and co-occurrence patterns

# Load required libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
library(patchwork)
library(circlize)
library(png)

# Get input and output paths from snakemake
defense_file <- snakemake@input[["defense"]]

# Output files
report_file <- snakemake@output[["report"]]
figures_dir <- snakemake@output[["figures_dir"]]

# Create figures directory
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Set theme for consistent visualization
theme_set(theme_bw(base_size = 12))
custom_theme <- theme(
  axis.title = element_text(face = "bold"),
  axis.text = element_text(size = 10),
  plot.title = element_text(size = 12, face = "bold"),
  legend.title = element_text(size = 10, face = "bold"),
  legend.text = element_text(size = 9),
  panel.grid.minor = element_blank()
)



#-----------------------------------------------------------------
# Data Loading and Preparation
#-----------------------------------------------------------------

# Load defense systems data
defense_df <- read_tsv(defense_file, show_col_types = FALSE)


# Clean DefenseFinder data to ensure one system per genome-type combination
defense_systems_data <- defense_df %>% 
  filter(!activity == "Antidefense") %>% 
  group_by(Genome_ID) %>% 
  distinct(type) %>%
  ungroup()

#-----------------------------------------------------------------
# FIGURE 1: Histogram
#-----------------------------------------------------------------

# Count systems per genome for DefenseFinder
defense_systems_per_genome <- defense_systems_data %>%
  group_by(Genome_ID) %>%
  summarize(num_systems = n()) %>%
  ungroup()


# Create histogram for DefenseFinder
hist_defensefinder <- ggplot(defense_systems_per_genome, aes(x = num_systems)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#377EB8", alpha = 0.8) +
  labs(
    title = "Distribution of Defence System Counts per Genome",
    x = "Number of Defense Systems per Genome",
    y = "Number of Genomes"
  ) +
  scale_x_continuous(breaks = seq(0, max(defense_systems_per_genome$num_systems), by = 1)) +
  custom_theme


# Save Figure 1
ggsave(file.path(figures_dir, "figure1_defense_count_distribution.png"), hist_defensefinder, 
       width = 12, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "figure1_defense_count_distribution.pdf"), hist_defensefinder, 
       width = 12, height = 6)

#-----------------------------------------------------------------
# Bar Charts (DefenseFinder and PADLOC side by side)
#-----------------------------------------------------------------

# Count occurrences for DefenseFinder
defense_type_counts <- defense_systems_data %>%
  count(type, name = "count") %>%
  arrange(desc(count))


# Create bar chart for DefenseFinder (top 20)
bar_defensefinder <- defense_type_counts %>% 
  head(20) %>%
  ggplot(aes(x = count, y = reorder(type, count), fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), hjust = -0.2, size = 3) +
  scale_fill_viridis_d(option = "turbo", begin = 0, end = 0.8, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Prevalence of Defence System Types",
    x = "Number of Genomes",
    y = NULL
  ) +
  custom_theme +
  theme(axis.text.y = element_text(face = "bold"))


# Save Figure 2
ggsave(file.path(figures_dir, "figure2_defense_type_prevalence.png"), bar_defensefinder, 
       width = 14, height = 8, dpi = 300)
ggsave(file.path(figures_dir, "figure2_defense_type_prevalence.pdf"), bar_defensefinder, 
       width = 14, height = 8)


# Combine Figure 1 and 2 side by side
figure12 <- hist_defensefinder + 
  bar_defensefinder + 
  plot_layout(ncol = 2, widths = c(1, 1)) + 
  plot_annotation(
    title = "A. Defence System Count Distribution    |    B. Defense System Type Prevalence",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 10, b = 10)
    )
  )

# Save the combined figure
ggsave(
  filename = file.path(figures_dir, "figure12_combined.png"),
  plot     = figure12,
  width    = 12,    # adjust to taste
  height   = 10,
  dpi      = 300
)
ggsave(
  filename = file.path(figures_dir, "figure12_combined.pdf"),
  plot     = figure12,
  width    = 12,
  height   = 10
)

#-----------------------------------------------------------------
# Circos Plot and Fisher's Exact Test Matrix
#-----------------------------------------------------------------

# Use DefenseFinder data for co-occurrence analysis (can be changed to PADLOC if preferred)
system_order <- defense_type_counts %>% 
  head(20) %>% 
  pull(type)

# Create presence-absence matrix
presence_matrix <- defense_systems_data %>%
  filter(type %in% system_order) %>%
  distinct(Genome_ID, type) %>% 
  mutate(
    type = factor(type, levels = system_order),
    present = 1
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = present,
    values_fill = 0
  ) %>%
  column_to_rownames("Genome_ID")

# Panel A: Circos Plot
set.seed(777) # For reproducibility

# Generate colors using a colorblind-friendly palette
system_colors_dynamic <- setNames(
  viridis::viridis_pal(option = "turbo")(length(system_order)),
  system_order
)
# Calculate co-occurrence counts
cooccur_counts <- crossprod(as.matrix(presence_matrix))

# Prepare data for circos plot
threshold <- max(cooccur_counts) * 0.05
filtered_counts <- cooccur_counts
filtered_counts[filtered_counts < threshold] <- 0

# Create the circos plot as a PNG file
png(file.path(figures_dir, "temp_circos_plot.png"), width = 2000, height = 1600, res = 200)
circos.clear()
circos.par(gap.after = 2, start.degree = 90)
par(mar = c(1, 1, 2, 1))

# Plot the chords diagram
chordDiagram(
  filtered_counts,
  grid.col = system_colors_dynamic,
  transparency = 0.6,
  directional = 0,
  annotationTrack = c("grid", "axis"),
  preAllocateTracks = list(track.height = 0.15)
)

# Add labels to the sectors
circos.trackPlotRegion(
  track.index = 1, 
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    # Add system labels
    circos.text(
      mean(xlim), ylim[1] + 0.2, 
      sector.name, 
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.9,
      font = 2,  # Bold font
      family = "Times New Roman"
    )
  }, 
  bg.border = NA
)

# Add title with Times New Roman
title(paste("A.", "Defence System Co-occurrence Network"), 
      line = 0.5,  
      cex.main = 1.4, 
      font.main = 2,
      family = "Times New Roman")
# Add legend directly
legend(x = 1.02, y = -0.1, 
       legend = names(system_colors_dynamic),
       fill = system_colors_dynamic,
       title = "Defence Systems",
       cex = 0.8,
       xpd = TRUE
)
dev.off()

# Convert PNG to a grob for inclusion in the patchwork
circos_grob <- rasterGrob(readPNG(file.path(figures_dir, "temp_circos_plot.png")), 
                         interpolate = TRUE)

# Panel B: Fisher's Exact Test Matrix
# Initialize matrices for odds ratios and p-values
odds_ratio_matrix <- matrix(NA, nrow = length(system_order), ncol = length(system_order))
pval_matrix <- matrix(NA, nrow = length(system_order), ncol = length(system_order))
rownames(odds_ratio_matrix) <- colnames(odds_ratio_matrix) <- system_order
rownames(pval_matrix) <- colnames(pval_matrix) <- system_order

# Perform Fisher's exact test for each pair
for (i in 1:length(system_order)) {
  for (j in 1:length(system_order)) {
    if (i == j) {
      odds_ratio_matrix[i, j] <- 100
      pval_matrix[i, j] <- 0.0001
    } else {
      sys1 <- system_order[i]
      sys2 <- system_order[j]
      
      n11 <- sum(presence_matrix[, sys1] == 1 & presence_matrix[, sys2] == 1)
      n10 <- sum(presence_matrix[, sys1] == 1 & presence_matrix[, sys2] == 0)
      n01 <- sum(presence_matrix[, sys1] == 0 & presence_matrix[, sys2] == 1)
      n00 <- sum(presence_matrix[, sys1] == 0 & presence_matrix[, sys2] == 0)
      
      contingency_table <- matrix(c(n11, n01, n10, n00), nrow = 2, byrow = TRUE)
      test_result <- fisher.test(contingency_table)
      
      odds_ratio_matrix[i, j] <- test_result$estimate
      pval_matrix[i, j] <- test_result$p.value
    }
  }
}

# Apply FDR correction
pval_adjusted_matrix <- matrix(p.adjust(pval_matrix, method = "BH"), 
                               nrow = length(system_order), 
                               ncol = length(system_order))
rownames(pval_adjusted_matrix) <- colnames(pval_adjusted_matrix) <- system_order

# Calculate log2 odds ratio
log2_odds_ratio_matrix <- log2(odds_ratio_matrix)
log2_odds_ratio_matrix[is.infinite(log2_odds_ratio_matrix) & log2_odds_ratio_matrix > 0] <- 16
log2_odds_ratio_matrix[is.infinite(log2_odds_ratio_matrix) & log2_odds_ratio_matrix < 0] <- -16
log2_odds_ratio_matrix[is.na(log2_odds_ratio_matrix)] <- 0

# Convert to data frame for ggplot
fisher_df <- data.frame(
  System1 = rep(system_order, each = length(system_order)),
  System2 = rep(system_order, times = length(system_order)),
  OddsRatio = as.vector(odds_ratio_matrix),
  Log2OddsRatio = as.vector(log2_odds_ratio_matrix),
  Pvalue = as.vector(pval_matrix),
  Padj = as.vector(pval_adjusted_matrix)
)

# Create Fisher's test plot
panel_fisher <- ggplot(fisher_df, aes(x = factor(System1, levels = system_order), 
                                      y = factor(System2, levels = rev(system_order)))) +
  geom_tile(
    aes(fill = Padj < 0.05),
    alpha = 0.2,
    color = "grey90"
  ) +
  scale_fill_manual(values = c("white", "yellow"), guide = "none") +
  geom_point(
    aes(size = -log10(Padj), color = Log2OddsRatio),
    shape = 16
  ) +
  scale_color_distiller(
    palette = "RdBu",
    direction = 1,
    limits = c(-16, 16),
    name = expression(log[2](OR))
  ) +
  scale_size_continuous(
    range = c(1, 10),
    name = expression(-log[10](p)),
    limits = c(0, 5),
    breaks = c(1, 2, 3, 4, 5)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title = element_blank(),
    panel.grid = element_line(color = "gray95"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "B. Defence System Co-occurrence Matrix",
    subtitle = "Yellow background indicates statistical significance (p < 0.05, FDR-corrected)"
  )

# Combine circos and Fisher's test plots
figure3 <- wrap_elements(full = circos_grob) / panel_fisher +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = 'Defense System Co-occurrence Analysis',
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Save Figure 3
ggsave(file.path(figures_dir, "defence_cooccurrence_analysis.png"), figure3, 
       width = 10, height = 14, dpi = 300)
ggsave(file.path(figures_dir, "defence_cooccurrence_analysis.pdf"), figure3, 
       width = 10, height = 14)



# Save summary data
write_tsv(defense_type_counts, 
          file.path(figures_dir, "defensefinder_type_counts.tsv"))
write_tsv(defense_systems_per_genome, 
          file.path(figures_dir, "defensefinder_systems_per_genome.tsv"))
write_tsv(fisher_df, 
          file.path(figures_dir, "fisher_test_results.tsv"))

# Clean up temporary files
if (file.exists(file.path(figures_dir, "temp_circos_plot.png"))) {
  file.remove(file.path(figures_dir, "temp_circos_plot.png"))
}

cat("Defense distribution analysis completed successfully!\n")
cat("Report saved to:", report_file, "\n")
cat("Figures saved to:", figures_dir, "\n")
