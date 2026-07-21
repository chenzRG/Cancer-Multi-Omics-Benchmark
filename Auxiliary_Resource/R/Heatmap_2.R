# Load necessary libraries
library(tidyr)
library(dplyr)
library(pheatmap)

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
labels <- read.csv(labels_path)

# Add a column for Sample names to the labels data frame
labels$Sample <- colnames(expression_data)

# Transpose expression data for easier manipulation
expression_data_t <- as.data.frame(t(expression_data))
expression_data_t$Sample <- rownames(expression_data_t)

# Merge expression data with labels
expression_data_merged <- merge(expression_data_t, labels, by = "Sample")

# Filter for patient groups 0, 1, and 2
filtered_data <- expression_data_merged %>%
  filter(Label %in% c(0, 1, 2))

# Select the top genes based on variance
top_genes <- filtered_data %>%
  select(-Sample, -Label) %>%
  summarise_all(~ var(.)) %>%
  gather(key = "Gene", value = "Variance") %>%
  arrange(desc(Variance)) %>%
  slice(1:100) %>%
  pull(Gene)

# Filter the expression data for these top 100 genes
top_expression_data <- expression_data_merged %>%
  filter(Label %in% c(0, 1, 2)) %>%
  select(Sample, one_of(top_genes))

# Transpose data for pheatmap
top_expression_matrix <- as.matrix(top_expression_data[, -1])  # Remove Sample column for matrix

#Heatmap plot
pheatmap(
  top_expression_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of Top 100 Genes Expression (Patient Groups 0, 1, 2)"
)