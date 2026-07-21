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
    mean_expr_0 = mean(Expression[Label == 0], na.rm = TRUE),
    log2FoldChange = ifelse(mean_expr_0 != 0, log2(mean_expr_1 / mean_expr_0), NA),
    pValue = t.test(Expression[Label == 1], Expression[Label == 0])$p.value
  )

# Filter out genes with NA log2FoldChange
results <- results %>%
  filter(!is.na(log2FoldChange))

# Create a DEGAll data frame with necessary columns
DEGAll <- results %>%
  mutate(
    logFC = log2FoldChange,
    PValue = pValue
  )

# Filter DEGAll for genes with pValue < 0.05 and abs(log2FoldChange) > 0.5
selected_genes <- DEGAll %>%
  filter(pValue < 0.05 & abs(log2FoldChange) > 0.5) %>%
  dplyr::select(Gene)

# Filter the expression data for the selected genes
selected_expression_data <- expression_data_t %>%
  dplyr::select(Sample, one_of(selected_genes$Gene))

# Add labels to the expression data
selected_expression_data <- selected_expression_data %>%
  mutate(Label = expression_data_merged$Label[match(selected_expression_data$Sample, expression_data_merged$Sample)])

# Set row names to Sample and remove Sample column
rownames(selected_expression_data) <- selected_expression_data$Sample
selected_expression_data <- selected_expression_data %>%
  dplyr::select(-Sample)

# Filter to include only valid labels (0 and 1)
selected_expression_data <- selected_expression_data %>%
  filter(Label %in% c(0, 1))

# Order the data by Label
selected_expression_data <- selected_expression_data[order(selected_expression_data$Label),]

# Ensure all values are finite
selected_expression_data_matrix <- as.matrix(selected_expression_data[, -ncol(selected_expression_data)])
selected_expression_data_matrix[!is.finite(selected_expression_data_matrix)] <- NA
selected_expression_data_matrix <- selected_expression_data_matrix[rowSums(is.na(selected_expression_data_matrix)) == 0, ]

# Check annotation_row
annotation_row <- data.frame(Label = selected_expression_data$Label)
rownames(annotation_row) <- rownames(selected_expression_data)

#Heatmap plot
pheatmap(
  selected_expression_data_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of DEGs")


