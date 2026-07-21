library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)  

# Data root: set MLOMICS_DATA_ROOT or default to repo Main_Dataset
data_root <- Sys.getenv(
  "MLOMICS_DATA_ROOT",
  unset = normalizePath(file.path("..", "..", "Main_Dataset"), mustWork = FALSE)
)

# Construct file paths using sprintf for string interpolation
cancer_type <- 'GBM'
expression_data_path <- file.path(data_root, sprintf('Classification_datasets/GS-%s/Top/%s_mRNA_top.csv', cancer_type, cancer_type))
labels_path <- file.path(data_root, sprintf('Classification_datasets/GS-%s/Top/%s_label_num.csv', cancer_type, cancer_type))

# Load your data
expression_data <- read.csv(expression_data_path, row.names = 1)
# Load the labels
labels <- read.csv(labels_path)

# Add a column for Sample names to the labels data frame
labels$Sample <- colnames(expression_data)

# Transpose expression data for easier manipulation
expression_data_t <- as.data.frame(t(expression_data))
expression_data_t$Sample <- rownames(expression_data_t)

# Merge expression data with labels
expression_data_merged <- merge(expression_data_t, labels, by = "Sample")

# Filter for groups 0, 1, and 2
filtered_data <- expression_data_merged %>%
  filter(Label %in% c(0, 1, 2))

# Define the genes of interest
genes_of_interest <- c("NPR3", "FBXO45", "TRAM1L1", "PI4KB", "KIAA0922", "WBP1")

# Create a list to store plots for each gene
plots <- lapply(genes_of_interest, function(gene) {
  gene_expression_data <- filtered_data %>%
    select(Sample, !!gene, Label)
  
  # Calculate average expression level for each gene and each group
  avg_expression <- gene_expression_data %>%
    group_by(Label) %>%
    summarise(avg_expression = mean(!!sym(gene), na.rm = TRUE))
  
  # Plot using ggplot with jittered points on x-axis
  p <- ggplot(gene_expression_data, aes(x = as.factor(Label), y = !!sym(gene), shape = as.factor(Label), color = as.factor(Label))) +
    geom_jitter(position = position_jitter(width = 0.2), size = 3) +  # Jitter points on x-axis
    scale_shape_manual(values = c(16, 17, 18)) +  # Set different shapes (16 = circle, 17 = triangle, 18 = square)
    labs(x = "Group", y = "Expression Level", color = "Group", shape = "Group") +
    ggtitle(paste("Gene Expression of", gene)) +
    theme_minimal() +
    theme(legend.position = "none") +  # Remove legend for individual plots
    geom_hline(data = avg_expression, aes(yintercept = avg_expression, color = as.factor(Label)), linetype = "dashed")  # Add average expression line per group
  
  return(p)
})

# Arrange plots in a 2x3 grid
multiplot <- cowplot::plot_grid(plotlist = plots, nrow = 2, align = 'hv')

# Display the grid of plots
multiplot
