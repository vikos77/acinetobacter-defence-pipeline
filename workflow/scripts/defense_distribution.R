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

# Get input and output paths - check if running via snakemake or standalone
if (exists("snakemake")) {
  # Running via Snakemake
  defense_file <- snakemake@input[["defense"]]
  padloc_file <- snakemake@input[["padloc"]]
  resfinder_file <- snakemake@input[["resfinder"]]
  figures_dir <- snakemake@output[["figures_dir"]]
  report_file <- if ("report" %in% names(snakemake@output)) snakemake@output[["report"]] else file.path(figures_dir, "defense_analysis_report.html")
} else {
  # Running standalone - use default paths
  defense_file <- "/app/data/consolidated/defense_finder_consolidated.tsv"
  padloc_file <- "/app/data/consolidated/padloc_consolidated.tsv"
  resfinder_file <- "/app/data/consolidated/resfinder_consolidated.tsv"
  figures_dir <- "/app/results/defense_figures"
  report_file <- file.path(figures_dir, "defense_analysis_report.html")
  cat("Running in standalone mode\n")
  cat("Defense file:", defense_file, "\n")
  cat("PADLOC file:", padloc_file, "\n")
  cat("ResFinder file:", resfinder_file, "\n")
  cat("Figures dir:", figures_dir, "\n")
}

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
padloc_df <- read_tsv(padloc_file, show_col_types = FALSE)
resfinder_df <- read_tsv(resfinder_file, show_col_types = FALSE)

# Clean DefenseFinder data to ensure one system per genome-type combination
defense_systems_data <- defense_df %>% 
  filter(!activity == "Antidefense") %>% 
  group_by(Genome_ID) %>% 
  distinct(type) %>%
  ungroup() %>%
  mutate(tool = "DefenseFinder")

# Clean PADLOC data to ensure one system per genome-type combination  
padloc_systems_data <- padloc_df %>%
  group_by(seqid) %>%
  distinct(system) %>%
  ungroup() %>%
  rename(Genome_ID = seqid, type = system) %>%
  mutate(tool = "PADLOC")

#-----------------------------------------------------------------
# FIGURE 1: Dual Histograms (DefenseFinder + PADLOC)
#-----------------------------------------------------------------

# Count systems per genome for DefenseFinder
defense_systems_per_genome <- defense_systems_data %>%
  group_by(Genome_ID) %>%
  summarize(num_systems = n(), tool = "DefenseFinder") %>%
  ungroup()

# Count systems per genome for PADLOC
padloc_systems_per_genome <- padloc_systems_data %>%
  group_by(Genome_ID) %>%
  summarize(num_systems = n(), tool = "PADLOC") %>%
  ungroup()

# Combine both datasets
combined_systems_per_genome <- bind_rows(defense_systems_per_genome, padloc_systems_per_genome)

# Create dual histogram
hist_combined <- ggplot(combined_systems_per_genome, aes(x = num_systems, fill = tool)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.8, position = "dodge") +
  scale_fill_manual(values = c("DefenseFinder" = "#377EB8", "PADLOC" = "#E41A1C")) +
  labs(
    title = "Distribution of Defence System Counts per Genome",
    x = "Number of Defense Systems per Genome",
    y = "Number of Genomes",
    fill = "Tool"
  ) +
  scale_x_continuous(breaks = seq(0, max(combined_systems_per_genome$num_systems), by = 1)) +
  custom_theme +
  theme(legend.position = "top")

# Create individual histograms for DefenseFinder and PADLOC
hist_defensefinder <- ggplot(defense_systems_per_genome, aes(x = num_systems)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#377EB8", alpha = 0.8) +
  labs(
    title = "A. DefenseFinder",
    x = "Number of Defense Systems per Genome",
    y = "Number of Genomes"
  ) +
  scale_x_continuous(breaks = seq(0, max(defense_systems_per_genome$num_systems), by = 1)) +
  custom_theme

hist_padloc <- ggplot(padloc_systems_per_genome, aes(x = num_systems)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#E41A1C", alpha = 0.8) +
  labs(
    title = "B. PADLOC", 
    x = "Number of Defense Systems per Genome",
    y = "Number of Genomes"
  ) +
  scale_x_continuous(breaks = seq(0, max(padloc_systems_per_genome$num_systems), by = 1)) +
  custom_theme

# Save individual and combined histograms
ggsave(file.path(figures_dir, "figure1_defense_count_distribution_combined.png"), hist_combined, 
       width = 12, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "figure1_defense_count_distribution_combined.pdf"), hist_combined, 
       width = 12, height = 6)

ggsave(file.path(figures_dir, "figure1a_defensefinder_distribution.png"), hist_defensefinder, 
       width = 10, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "figure1b_padloc_distribution.png"), hist_padloc, 
       width = 10, height = 6, dpi = 300)

#-----------------------------------------------------------------
# Bar Charts (DefenseFinder and PADLOC side by side)
#-----------------------------------------------------------------

# Count occurrences for DefenseFinder
defense_type_counts <- defense_systems_data %>%
  count(type, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(tool = "DefenseFinder")

# Count occurrences for PADLOC
padloc_type_counts <- padloc_systems_data %>%
  count(type, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(tool = "PADLOC")

# Create bar chart for DefenseFinder (top 20)
bar_defensefinder <- defense_type_counts %>% 
  head(20) %>%
  ggplot(aes(x = count, y = reorder(type, count), fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), hjust = -0.2, size = 3) +
  scale_fill_viridis_d(option = "turbo", begin = 0, end = 0.8, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "A. DefenseFinder",
    x = "Number of Genomes",
    y = NULL
  ) +
  custom_theme +
  theme(axis.text.y = element_text(face = "bold"))

# Create bar chart for PADLOC (top 20)  
bar_padloc <- padloc_type_counts %>% 
  head(20) %>%
  ggplot(aes(x = count, y = reorder(type, count), fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), hjust = -0.2, size = 3) +
  scale_fill_viridis_d(option = "plasma", begin = 0, end = 0.8, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "B. PADLOC",
    x = "Number of Genomes", 
    y = NULL
  ) +
  custom_theme +
  theme(axis.text.y = element_text(face = "bold"))

# Save individual bar charts
ggsave(file.path(figures_dir, "figure2a_defensefinder_prevalence.png"), bar_defensefinder, 
       width = 14, height = 8, dpi = 300)
ggsave(file.path(figures_dir, "figure2b_padloc_prevalence.png"), bar_padloc, 
       width = 14, height = 8, dpi = 300)


# Create Figure 1: Side-by-side histograms comparison
figure1_histograms <- hist_defensefinder + hist_padloc + 
  plot_layout(ncol = 2) + 
  plot_annotation(
    title = "Figure 1: Defence System Count Distribution Comparison",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 15, b = 10)
    )
  )

# Create Figure 2: Side-by-side bar charts comparison  
figure2_barplots <- bar_defensefinder + bar_padloc + 
  plot_layout(ncol = 2) + 
  plot_annotation(
    title = "Figure 2: Defence System Type Prevalence Comparison",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 15, b = 10)
    )
  )

# Save the comparison figures
ggsave(
  filename = file.path(figures_dir, "figure1_histograms_comparison.png"),
  plot     = figure1_histograms,
  width    = 20,
  height   = 8,
  dpi      = 300
)
ggsave(
  filename = file.path(figures_dir, "figure1_histograms_comparison.pdf"),
  plot     = figure1_histograms,
  width    = 20,
  height   = 8
)
ggsave(
  filename = file.path(figures_dir, "figure2_barplots_comparison.png"), 
  plot     = figure2_barplots,
  width    = 24,
  height   = 10,
  dpi      = 300
)
ggsave(
  filename = file.path(figures_dir, "figure2_barplots_comparison.pdf"), 
  plot     = figure2_barplots,
  width    = 24,
  height   = 10
)

#-----------------------------------------------------------------
# Circos Plot and Fisher's Exact Test Matrix (Both Tools)
#-----------------------------------------------------------------

# DEFENSEFINDER Co-occurrence Analysis
defense_system_order <- defense_type_counts %>% 
  head(15) %>%  # Reduced to 15 for better visualization
  pull(type)

# Create DefenseFinder presence-absence matrix
defense_presence_matrix <- defense_systems_data %>%
  filter(type %in% defense_system_order) %>%
  distinct(Genome_ID, type) %>% 
  mutate(
    type = factor(type, levels = defense_system_order),
    present = 1
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = present,
    values_fill = 0
  ) %>%
  column_to_rownames("Genome_ID")

# PADLOC Co-occurrence Analysis  
padloc_system_order <- padloc_type_counts %>% 
  head(15) %>%  # Reduced to 15 for better visualization
  pull(type)

# Create PADLOC presence-absence matrix
padloc_presence_matrix <- padloc_systems_data %>%
  filter(type %in% padloc_system_order) %>%
  distinct(Genome_ID, type) %>% 
  mutate(
    type = factor(type, levels = padloc_system_order),
    present = 1
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = present,
    values_fill = 0
  ) %>%
  column_to_rownames("Genome_ID")

# Generate Circos Plots for Both Tools
set.seed(777) # For reproducibility

# DEFENSEFINDER Circos Plot
defense_colors <- setNames(
  viridis::viridis_pal(option = "turbo")(length(defense_system_order)),
  defense_system_order
)

# Calculate DefenseFinder co-occurrence counts
defense_cooccur_counts <- crossprod(as.matrix(defense_presence_matrix))
defense_threshold <- max(defense_cooccur_counts) * 0.05
defense_filtered_counts <- defense_cooccur_counts
defense_filtered_counts[defense_filtered_counts < defense_threshold] <- 0

# Create DefenseFinder circos plot
png(file.path(figures_dir, "temp_defensefinder_circos.png"), width = 2000, height = 1600, res = 200)
circos.clear()
circos.par(gap.after = 2, start.degree = 90)
par(mar = c(1, 1, 2, 1))

chordDiagram(
  defense_filtered_counts,
  grid.col = defense_colors,
  transparency = 0.6,
  directional = 0,
  annotationTrack = c("grid", "axis"),
  preAllocateTracks = list(track.height = 0.15)
)

circos.trackPlotRegion(
  track.index = 1, 
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    circos.text(
      mean(xlim), ylim[1] + 0.2, 
      sector.name, 
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.8,
      font = 2
    )
  }, 
  bg.border = NA
)

title("A. DefenseFinder Co-occurrence Network", 
      line = 0.5, cex.main = 1.4, font.main = 2)
dev.off()

# PADLOC Circos Plot
padloc_colors <- setNames(
  viridis::viridis_pal(option = "plasma")(length(padloc_system_order)),
  padloc_system_order
)

# Calculate PADLOC co-occurrence counts
padloc_cooccur_counts <- crossprod(as.matrix(padloc_presence_matrix))
padloc_threshold <- max(padloc_cooccur_counts) * 0.05
padloc_filtered_counts <- padloc_cooccur_counts
padloc_filtered_counts[padloc_filtered_counts < padloc_threshold] <- 0

# Create PADLOC circos plot
png(file.path(figures_dir, "temp_padloc_circos.png"), width = 2000, height = 1600, res = 200)
circos.clear()
circos.par(gap.after = 2, start.degree = 90)
par(mar = c(1, 1, 2, 1))

chordDiagram(
  padloc_filtered_counts,
  grid.col = padloc_colors,
  transparency = 0.6,
  directional = 0,
  annotationTrack = c("grid", "axis"),
  preAllocateTracks = list(track.height = 0.15)
)

circos.trackPlotRegion(
  track.index = 1, 
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    circos.text(
      mean(xlim), ylim[1] + 0.2, 
      sector.name, 
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.8,
      font = 2
    )
  }, 
  bg.border = NA
)

title("B. PADLOC Co-occurrence Network", 
      line = 0.5, cex.main = 1.4, font.main = 2)
dev.off()

# Convert PNGs to grobs for patchwork
defensefinder_circos_grob <- rasterGrob(readPNG(file.path(figures_dir, "temp_defensefinder_circos.png")), 
                                       interpolate = TRUE)
padloc_circos_grob <- rasterGrob(readPNG(file.path(figures_dir, "temp_padloc_circos.png")), 
                                interpolate = TRUE)

# Fisher's Exact Test Matrices for Both Tools

# DEFENSEFINDER Fisher's Test
defense_odds_ratio_matrix <- matrix(NA, nrow = length(defense_system_order), ncol = length(defense_system_order))
defense_pval_matrix <- matrix(NA, nrow = length(defense_system_order), ncol = length(defense_system_order))
rownames(defense_odds_ratio_matrix) <- colnames(defense_odds_ratio_matrix) <- defense_system_order
rownames(defense_pval_matrix) <- colnames(defense_pval_matrix) <- defense_system_order

# Perform DefenseFinder Fisher's exact test
for (i in 1:length(defense_system_order)) {
  for (j in 1:length(defense_system_order)) {
    if (i == j) {
      defense_odds_ratio_matrix[i, j] <- 100
      defense_pval_matrix[i, j] <- 0.0001
    } else {
      sys1 <- defense_system_order[i]
      sys2 <- defense_system_order[j]
      
      n11 <- sum(defense_presence_matrix[, sys1] == 1 & defense_presence_matrix[, sys2] == 1)
      n10 <- sum(defense_presence_matrix[, sys1] == 1 & defense_presence_matrix[, sys2] == 0)
      n01 <- sum(defense_presence_matrix[, sys1] == 0 & defense_presence_matrix[, sys2] == 1)
      n00 <- sum(defense_presence_matrix[, sys1] == 0 & defense_presence_matrix[, sys2] == 0)
      
      contingency_table <- matrix(c(n11, n01, n10, n00), nrow = 2, byrow = TRUE)
      test_result <- fisher.test(contingency_table)
      
      defense_odds_ratio_matrix[i, j] <- test_result$estimate
      defense_pval_matrix[i, j] <- test_result$p.value
    }
  }
}

# DefenseFinder FDR correction and log2 transformation
defense_pval_adjusted_matrix <- matrix(p.adjust(defense_pval_matrix, method = "BH"), 
                                      nrow = length(defense_system_order), 
                                      ncol = length(defense_system_order))
rownames(defense_pval_adjusted_matrix) <- colnames(defense_pval_adjusted_matrix) <- defense_system_order

defense_log2_odds_ratio_matrix <- log2(defense_odds_ratio_matrix)
defense_log2_odds_ratio_matrix[is.infinite(defense_log2_odds_ratio_matrix) & defense_log2_odds_ratio_matrix > 0] <- 16
defense_log2_odds_ratio_matrix[is.infinite(defense_log2_odds_ratio_matrix) & defense_log2_odds_ratio_matrix < 0] <- -16
defense_log2_odds_ratio_matrix[is.na(defense_log2_odds_ratio_matrix)] <- 0

# DefenseFinder data frame for ggplot
defense_fisher_df <- data.frame(
  System1 = rep(defense_system_order, each = length(defense_system_order)),
  System2 = rep(defense_system_order, times = length(defense_system_order)),
  OddsRatio = as.vector(defense_odds_ratio_matrix),
  Log2OddsRatio = as.vector(defense_log2_odds_ratio_matrix),
  Pvalue = as.vector(defense_pval_matrix),
  Padj = as.vector(defense_pval_adjusted_matrix)
)

# PADLOC Fisher's Test
padloc_odds_ratio_matrix <- matrix(NA, nrow = length(padloc_system_order), ncol = length(padloc_system_order))
padloc_pval_matrix <- matrix(NA, nrow = length(padloc_system_order), ncol = length(padloc_system_order))
rownames(padloc_odds_ratio_matrix) <- colnames(padloc_odds_ratio_matrix) <- padloc_system_order
rownames(padloc_pval_matrix) <- colnames(padloc_pval_matrix) <- padloc_system_order

# Perform PADLOC Fisher's exact test
for (i in 1:length(padloc_system_order)) {
  for (j in 1:length(padloc_system_order)) {
    if (i == j) {
      padloc_odds_ratio_matrix[i, j] <- 100
      padloc_pval_matrix[i, j] <- 0.0001
    } else {
      sys1 <- padloc_system_order[i]
      sys2 <- padloc_system_order[j]
      
      n11 <- sum(padloc_presence_matrix[, sys1] == 1 & padloc_presence_matrix[, sys2] == 1)
      n10 <- sum(padloc_presence_matrix[, sys1] == 1 & padloc_presence_matrix[, sys2] == 0)
      n01 <- sum(padloc_presence_matrix[, sys1] == 0 & padloc_presence_matrix[, sys2] == 1)
      n00 <- sum(padloc_presence_matrix[, sys1] == 0 & padloc_presence_matrix[, sys2] == 0)
      
      contingency_table <- matrix(c(n11, n01, n10, n00), nrow = 2, byrow = TRUE)
      test_result <- fisher.test(contingency_table)
      
      padloc_odds_ratio_matrix[i, j] <- test_result$estimate
      padloc_pval_matrix[i, j] <- test_result$p.value
    }
  }
}

# PADLOC FDR correction and log2 transformation
padloc_pval_adjusted_matrix <- matrix(p.adjust(padloc_pval_matrix, method = "BH"), 
                                     nrow = length(padloc_system_order), 
                                     ncol = length(padloc_system_order))
rownames(padloc_pval_adjusted_matrix) <- colnames(padloc_pval_adjusted_matrix) <- padloc_system_order

padloc_log2_odds_ratio_matrix <- log2(padloc_odds_ratio_matrix)
padloc_log2_odds_ratio_matrix[is.infinite(padloc_log2_odds_ratio_matrix) & padloc_log2_odds_ratio_matrix > 0] <- 16
padloc_log2_odds_ratio_matrix[is.infinite(padloc_log2_odds_ratio_matrix) & padloc_log2_odds_ratio_matrix < 0] <- -16
padloc_log2_odds_ratio_matrix[is.na(padloc_log2_odds_ratio_matrix)] <- 0

# PADLOC data frame for ggplot
padloc_fisher_df <- data.frame(
  System1 = rep(padloc_system_order, each = length(padloc_system_order)),
  System2 = rep(padloc_system_order, times = length(padloc_system_order)),
  OddsRatio = as.vector(padloc_odds_ratio_matrix),
  Log2OddsRatio = as.vector(padloc_log2_odds_ratio_matrix),
  Pvalue = as.vector(padloc_pval_matrix),
  Padj = as.vector(padloc_pval_adjusted_matrix)
)

# Create DefenseFinder Fisher's test plot
defense_fisher_plot <- ggplot(defense_fisher_df, aes(x = factor(System1, levels = defense_system_order), 
                                                    y = factor(System2, levels = rev(defense_system_order)))) +
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
    range = c(1, 8),
    name = expression(-log[10](p)),
    limits = c(0, 5),
    breaks = c(1, 2, 3, 4, 5)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title = element_blank(),
    panel.grid = element_line(color = "gray95"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "A. DefenseFinder Co-occurrence Matrix")

# Create PADLOC Fisher's test plot
padloc_fisher_plot <- ggplot(padloc_fisher_df, aes(x = factor(System1, levels = padloc_system_order), 
                                                  y = factor(System2, levels = rev(padloc_system_order)))) +
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
    range = c(1, 8),
    name = expression(-log[10](p)),
    limits = c(0, 5),
    breaks = c(1, 2, 3, 4, 5)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title = element_blank(),
    panel.grid = element_line(color = "gray95"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  labs(title = "B. PADLOC Co-occurrence Matrix")

# Create Figure 3: Side-by-side Circos plots
figure3_circos <- wrap_elements(full = defensefinder_circos_grob) + wrap_elements(full = padloc_circos_grob) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = 'Figure 3: Defence System Co-occurrence Networks Comparison',
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 15, b = 10)
    )
  )

# Create Figure 4: Side-by-side Fisher's exact test matrices
figure4_fisher <- defense_fisher_plot + padloc_fisher_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = 'Figure 4: Defence System Co-occurrence Statistical Analysis',
    subtitle = "Yellow background indicates statistical significance (p < 0.05, FDR-corrected)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.margin = margin(t = 15, b = 10)
    )
  )

# Save the co-occurrence comparison figures
ggsave(
  filename = file.path(figures_dir, "figure3_circos_networks_comparison.png"),
  plot     = figure3_circos,
  width    = 20,
  height   = 10,
  dpi      = 300
)
ggsave(
  filename = file.path(figures_dir, "figure3_circos_networks_comparison.pdf"),
  plot     = figure3_circos,
  width    = 20,
  height   = 10
)
ggsave(
  filename = file.path(figures_dir, "figure4_fisher_matrices_comparison.png"),
  plot     = figure4_fisher,
  width    = 20,
  height   = 12,
  dpi      = 300
)
ggsave(
  filename = file.path(figures_dir, "figure4_fisher_matrices_comparison.pdf"),
  plot     = figure4_fisher,
  width    = 20,
  height   = 12
)



#-----------------------------------------------------------------
# ARG-Defense Correlation Analysis (Both Tools)
#-----------------------------------------------------------------

# Prepare ResFinder data
resfinder_clean <- resfinder_df %>%
  select(Genome_ID, `Resistance gene`) %>%
  distinct(Genome_ID, `Resistance gene`) %>%
  filter(!is.na(`Resistance gene`), `Resistance gene` != "")

# Get top ARGs for analysis (top 10 most prevalent)
top_args_counts <- resfinder_clean %>%
  count(`Resistance gene`, name = "count") %>%
  arrange(desc(count)) %>%
  head(10)

top_args_list <- top_args_counts %>% pull(`Resistance gene`)

cat("Top ARGs for correlation analysis:", paste(head(top_args_list, 5), collapse = ", "), "...\n")

# Function to perform Fisher's exact test between defense systems and ARGs
perform_arg_defense_correlation <- function(defense_data, defense_systems, tool_name) {
  
  # Get top defense systems (top 10 for cleaner visualization)
  top_defense_systems <- defense_data %>%
    count(type, name = "count") %>%
    arrange(desc(count)) %>%
    head(10) %>%
    pull(type)
  
  # Create presence-absence matrices
  defense_matrix <- defense_data %>%
    filter(type %in% top_defense_systems) %>%
    select(Genome_ID, type) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(
      id_cols = Genome_ID,
      names_from = type,
      values_from = present,
      values_fill = 0
    )
  
  arg_matrix <- resfinder_clean %>%
    filter(`Resistance gene` %in% top_args_list) %>%
    select(Genome_ID, `Resistance gene`) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(
      id_cols = Genome_ID,
      names_from = `Resistance gene`,
      values_from = present,
      values_fill = 0
    )
  
  # Join matrices
  all_genomes <- unique(c(defense_matrix$Genome_ID, arg_matrix$Genome_ID))
  joined_matrix <- data.frame(Genome_ID = all_genomes) %>%
    left_join(defense_matrix, by = "Genome_ID") %>%
    left_join(arg_matrix, by = "Genome_ID") %>%
    mutate(across(everything(), ~replace_na(., 0)))
  
  # Extract relevant columns
  defense_cols <- top_defense_systems[top_defense_systems %in% colnames(joined_matrix)]
  arg_cols <- top_args_list[top_args_list %in% colnames(joined_matrix)]
  
  # Initialize result matrix
  fisher_results <- data.frame(
    Defense_System = character(),
    ARG = character(),
    OddsRatio = numeric(),
    Pvalue = numeric(),
    Tool = character(),
    stringsAsFactors = FALSE
  )
  
  # Perform Fisher's exact test for each pair
  for (defense in defense_cols) {
    for (arg in arg_cols) {
      defense_data_vec <- joined_matrix[[defense]]
      arg_data_vec <- joined_matrix[[arg]]
      
      contingency_table <- table(
        factor(defense_data_vec, levels = c(0, 1)),
        factor(arg_data_vec, levels = c(0, 1))
      )
      
      # Add small constant if table has zeros
      if (any(contingency_table == 0)) {
        contingency_table <- contingency_table + 0.5
      }
      
      # Perform Fisher's exact test
      test_result <- fisher.test(contingency_table)
      
      # Store results
      fisher_results <- rbind(fisher_results, data.frame(
        Defense_System = defense,
        ARG = arg,
        OddsRatio = test_result$estimate,
        Pvalue = test_result$p.value,
        Tool = tool_name,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Apply FDR correction and process results
  fisher_results <- fisher_results %>%
    mutate(
      Padj = p.adjust(Pvalue, method = "BH"),
      Significant = Padj < 0.05,
      LogOddsRatio = log2(OddsRatio),
      # Handle infinite values
      LogOddsRatio = case_when(
        is.infinite(LogOddsRatio) & LogOddsRatio > 0 ~ 10,
        is.infinite(LogOddsRatio) & LogOddsRatio < 0 ~ -10,
        TRUE ~ LogOddsRatio
      )
    )
  
  return(list(
    results = fisher_results,
    defense_systems = defense_cols,
    args = arg_cols
  ))
}

# Perform correlation analysis for DefenseFinder
cat("Performing ARG correlation analysis for DefenseFinder...\n")
defensefinder_arg_results <- perform_arg_defense_correlation(defense_systems_data, defense_system_order, "DefenseFinder")

# Perform correlation analysis for PADLOC
cat("Performing ARG correlation analysis for PADLOC...\n")
padloc_arg_results <- perform_arg_defense_correlation(padloc_systems_data, padloc_system_order, "PADLOC")

# Create DefenseFinder ARG correlation plot
defensefinder_fisher_long <- defensefinder_arg_results$results %>%
  mutate(
    Defense_System = factor(Defense_System, levels = rev(defensefinder_arg_results$defense_systems)),
    ARG = factor(ARG, levels = defensefinder_arg_results$args)
  )

defensefinder_arg_plot <- ggplot(defensefinder_fisher_long, aes(x = ARG, y = Defense_System, fill = LogOddsRatio)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(
    low = "navy", 
    mid = "white", 
    high = "firebrick",
    midpoint = 0,
    name = expression(log[2](OR)),
    limits = c(-10, 10)
  ) +
  geom_text(aes(label = ifelse(Significant, 
                               ifelse(Padj < 0.001, "***", 
                                      ifelse(Padj < 0.01, "**", "*")), 
                               "")),
            color = "black", size = 3) +
  labs(
    title = "A. DefenseFinder - ARG Correlations",
    x = "Antibiotic Resistance Genes",
    y = "Defence Systems"
  ) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )

# Create PADLOC ARG correlation plot
padloc_fisher_long <- padloc_arg_results$results %>%
  mutate(
    Defense_System = factor(Defense_System, levels = rev(padloc_arg_results$defense_systems)),
    ARG = factor(ARG, levels = padloc_arg_results$args)
  )

padloc_arg_plot <- ggplot(padloc_fisher_long, aes(x = ARG, y = Defense_System, fill = LogOddsRatio)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(
    low = "navy", 
    mid = "white", 
    high = "firebrick",
    midpoint = 0,
    name = expression(log[2](OR)),
    limits = c(-10, 10)
  ) +
  geom_text(aes(label = ifelse(Significant, 
                               ifelse(Padj < 0.001, "***", 
                                      ifelse(Padj < 0.01, "**", "*")), 
                               "")),
            color = "black", size = 3) +
  labs(
    title = "B. PADLOC - ARG Correlations",
    x = "Antibiotic Resistance Genes",
    y = "Defence Systems"
  ) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )

# Create Figure 5: Side-by-side ARG correlation matrices
figure5_arg_correlations <- defensefinder_arg_plot + padloc_arg_plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = 'Figure 5: Defense System - Antibiotic Resistance Gene Correlations',
    subtitle = "Asterisks indicate statistical significance: *** p<0.001, ** p<0.01, * p<0.05 (FDR-corrected)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.margin = margin(t = 15, b = 10)
    )
  )

# Save Figure 5
ggsave(
  filename = file.path(figures_dir, "figure5_arg_correlations_comparison.png"),
  plot     = figure5_arg_correlations,
  width    = 20,
  height   = 12,
  dpi      = 300
)
ggsave(
  filename = file.path(figures_dir, "figure5_arg_correlations_comparison.pdf"),
  plot     = figure5_arg_correlations,
  width    = 20,
  height   = 12
)

# Save summary data
write_tsv(defense_type_counts, 
          file.path(figures_dir, "defensefinder_type_counts.tsv"))
write_tsv(padloc_type_counts, 
          file.path(figures_dir, "padloc_type_counts.tsv"))
write_tsv(defense_systems_per_genome, 
          file.path(figures_dir, "defensefinder_systems_per_genome.tsv"))
write_tsv(padloc_systems_per_genome, 
          file.path(figures_dir, "padloc_systems_per_genome.tsv"))
write_tsv(combined_systems_per_genome, 
          file.path(figures_dir, "combined_systems_per_genome.tsv"))
write_tsv(defense_fisher_df, 
          file.path(figures_dir, "defensefinder_fisher_test_results.tsv"))
write_tsv(padloc_fisher_df, 
          file.path(figures_dir, "padloc_fisher_test_results.tsv"))
write_tsv(defensefinder_arg_results$results, 
          file.path(figures_dir, "defensefinder_arg_correlation_results.tsv"))
write_tsv(padloc_arg_results$results, 
          file.path(figures_dir, "padloc_arg_correlation_results.tsv"))
write_tsv(top_args_counts, 
          file.path(figures_dir, "top_args_prevalence.tsv"))

# Clean up temporary files
if (file.exists(file.path(figures_dir, "temp_defensefinder_circos.png"))) {
  file.remove(file.path(figures_dir, "temp_defensefinder_circos.png"))
}
if (file.exists(file.path(figures_dir, "temp_padloc_circos.png"))) {
  file.remove(file.path(figures_dir, "temp_padloc_circos.png"))
}

cat("Defense distribution analysis completed successfully!\n")
cat("Analyzed both DefenseFinder and PADLOC data:\n")
cat("- DefenseFinder systems:", nrow(defense_systems_data), "\n")
cat("- PADLOC systems:", nrow(padloc_systems_data), "\n") 
cat("Report saved to:", report_file, "\n")
cat("Figures saved to:", figures_dir, "\n")
cat("Generated comparison figures:\n")
cat("- Figure 1 - Side-by-side histograms: figure1_histograms_comparison.png\n")
cat("- Figure 2 - Side-by-side bar charts: figure2_barplots_comparison.png\n")
cat("- Figure 3 - Side-by-side circos networks: figure3_circos_networks_comparison.png\n")
cat("- Figure 4 - Side-by-side Fisher matrices: figure4_fisher_matrices_comparison.png\n")
cat("- Figure 5 - Side-by-side ARG correlations: figure5_arg_correlations_comparison.png\n")
cat("- Individual files also available for detailed analysis\n")
cat("- Top", length(top_args_list), "ARGs analyzed for correlations\n")
cat("- ARG correlation results saved for both DefenseFinder and PADLOC\n")
