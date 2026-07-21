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

# Filter the data to include only cancer types labeled from 0 to 13
expression_data_filtered <- expression_data_merged %>% filter(Label >= 0 & Label <= 13)

# Plot boxplot for CD161 expression across all cancer types from 0 to 13
gene_of_interest <- "CD161"

p <- ggplot(expression_data_filtered, aes(x = factor(Label), y = .data[[gene_of_interest]])) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2, notch = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.6) +
  labs(title = sprintf("Expression of %s Across Selected Cancer Types", gene_of_interest),
       x = "Cancer Type",
       y = sprintf("Expression of %s", gene_of_interest)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text_repel(aes(label = ifelse(.data[[gene_of_interest]] > quantile(.data[[gene_of_interest]], 0.75), Sample, "")), 
                  size = 3, box.padding = 0.5)

# Print the plot
print(p)

# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Define the gene of interest as a character vector
gene_of_interest <- c("KLRB1")
# Perform gene ID conversion with bitr() to get Entrez IDs
hg <- bitr(gene_of_interest, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")

# Load necessary libraries
library(ggplot2)
library(ggrepel)

# Perform GO enrichment analysis using enrichGO()
go <- enrichGO(hg$ENTREZID,
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.1, 
               qvalueCutoff = 0.1,
               keyType = 'ENTREZID')

# Plot enriched GO terms using ggplot2 and dotplot
if (!is.null(go) && !is.na(go) && nrow(go) > 0) {
  ggplot(go, split = "ONTOLOGY", showCategory = 5) +
    geom_point(aes(x = GeneRatio,
                   y = Description,
                   color = p.adjust,
                   size = GeneRatio)) +
    scale_color_gradient(low = "#1763a3", high = "#800000") +
    facet_grid(ONTOLOGY ~ ., scale = "free") +
    theme_bw(base_size = 18) +   # Set plot theme
    theme(text = element_text(size = 17))  # Adjust text size
} else {
  print("No enrichment GO terms found!")
}

# Load necessary libraries
library(pheatmap)  # pheatmap for heatmap visualization

# Define cancer types of interest
cancer_types_of_interest <- 0:3

# Filter the data to include only cancer types 0-3
expression_data_filtered <- expression_data_merged[expression_data_merged$Label %in% cancer_types_of_interest, ]

# Identify top genes based on variance (select the top 10 genes with highest variance)
# Calculate variance of each gene
variances <- apply(expression_data_filtered[, !names(expression_data_filtered) %in% c("Sample", "Label")], 2, var)
top_genes <- names(sort(variances, decreasing = TRUE))[1:10]

# Filter expression data to include only top genes
expression_data_top_genes <- expression_data_filtered[, c("Sample", "Label", top_genes)]

# Convert the data to matrix format for heatmap
# Set row names and remove the Label and Sample columns
expression_data_top_genes_matrix <- expression_data_top_genes
rownames(expression_data_top_genes_matrix) <- expression_data_top_genes_matrix$Sample
expression_data_top_genes_matrix <- expression_data_top_genes_matrix[, !names(expression_data_top_genes_matrix) %in% c("Sample", "Label")]
expression_data_top_genes_matrix <- as.matrix(expression_data_top_genes_matrix)

# Create heatmap
pheatmap(
  expression_data_top_genes_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",  # Scaling the rows to better visualize patterns
  main = "Heatmap of Top Genes in Cancer Types 0-3",
  show_rownames = FALSE  # Hide sample names
)
