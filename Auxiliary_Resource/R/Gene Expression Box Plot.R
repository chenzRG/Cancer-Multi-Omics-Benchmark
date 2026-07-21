# Load necessary libraries
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
cancer_type <- 'Pan-cancer'
expression_data_path <- file.path(data_root, sprintf('Classification_datasets/%s/Original/%s_mRNA.csv', cancer_type, cancer_type))
labels_path <- file.path(data_root, sprintf('Classification_datasets/%s/Original/%s_label_num.csv', cancer_type, cancer_type))

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

# Gather data into long format
expression_long <- gather(expression_data_merged, Gene, Expression, -Sample, -Label)

# Calculate log2 fold changes and p-values for each gene
results <- expression_long %>%
  group_by(Gene) %>%
  summarize(
    mean_expr_1 = mean(Expression[Label == 1], na.rm = TRUE),
    mean_expr_0 = mean(Expression[Label == 2], na.rm = TRUE),
    log2FoldChange = ifelse(mean_expr_0 > 0, log2(mean_expr_1 / mean_expr_0),
                            ifelse(mean_expr_0 < 0, log2(mean_expr_1 / abs(mean_expr_0)), 0)),  # Handle negative mean_expr_0
    pValue = t.test(Expression[Label == 1], Expression[Label == 2])$p.value
  )

# Filter out genes with NA log2FoldChange
results <- results %>%
  filter(!is.na(log2FoldChange))

# Create a DEGAll data frame with necessary columns for GO enrichment
DEGAll <- results %>%
  mutate(
    logFC = log2FoldChange,
    PValue = pValue
  )

# Filter DEGAll for genes with pValue < 0.05 and abs(log2FoldChange) > 0.5
selected_genes <- DEGAll %>%
  filter(pValue < 0.05 & abs(log2FoldChange) > 0.5) %>%
  dplyr::select(Gene)

# Filter DEGAll for genes with pValue < 0.05 and abs(log2FoldChange) > 0.5
selected_genes <- DEGAll %>%
  filter(pValue < 0.05 & abs(logFC) > 0.5) %>%
  arrange(desc(abs(logFC))) %>%  # Arrange by descending absolute logFC
  head(15)  # Select top 15 genes

# Filter expression data for the selected genes and groups 0 and 1
filtered_data <- expression_long %>%
  filter(Gene %in% selected_genes$Gene, Label %in% c(0, 1))

# Function to generate boxplot with overlaid dots for each gene
plot_gene <- function(gene) {
  gene_data <- filtered_data %>%
    filter(Gene == gene)
  
  p <- ggplot(gene_data, aes(x = as.factor(Label), y = Expression, fill = as.factor(Label))) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = c("#999999", "#E69F00")) +  # Adjust to match the number of levels in Label
    labs(x = "Group", y = "Expression Level", fill = "Group") +
    ggtitle(paste("Gene:", gene)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
  
  return(p)
}

# Create a list to store plots for each gene
plots <- lapply(selected_genes$Gene, plot_gene)

# Arrange plots in a 3x5 grid using cowplot
multiplot <- cowplot::plot_grid(plotlist = plots, ncol = 5)

# Display the grid of plots
print(multiplot)