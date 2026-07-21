# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(cluster)
library(factoextra)
library(survival)
library(survminer)
library(gridExtra)
library(org.Hs.eg.db)

# Data root: set MLOMICS_DATA_ROOT or default to repo Main_Dataset
data_root <- Sys.getenv(
  "MLOMICS_DATA_ROOT",
  unset = normalizePath(file.path("..", "..", "Main_Dataset"), mustWork = FALSE)
)

# Construct file paths using sprintf for string interpolation
cancer_type <- 'LUAD'
expression_data_path <- file.path(data_root, sprintf('Clustering_datasets/%s/Top/%s_mRNA_top.csv', cancer_type, cancer_type))
survival_path <- file.path(data_root, sprintf('Clustering_datasets/%s/Top/survival_%s.csv', cancer_type, cancer_type))

# Load the expression data
expression_data <- read.csv(expression_data_path, row.names = 1)
# Load the survival data
survival <- read.csv(survival_path)

# Add columns for Sample names to the survival data frame
survival$Sample <- colnames(expression_data)

# Transpose expression data for easier manipulation
expression_data_t <- as.data.frame(t(expression_data))
expression_data_t$Sample <- rownames(expression_data_t)

# Merge expression data with survival data
expression_data_merged <- merge(expression_data_t, survival, by = "Sample")

# List of genes of interest
genes_of_interest <- c("UCN2", "RIMS2", "CAVIN2", "GRIA1", "PKHD1L1", "PGM5", "CLIC6")

# Perform gene ID conversion with bitr() to get Entrez IDs
gene_symbols <- genes_of_interest 
hg <- bitr(gene_symbols, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

# Perform GO enrichment analysis using enrichGO()
go <- enrichGO(hg$ENTREZID,
               OrgDb = org.Hs.eg.db, 
               ont = 'ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.1, 
               qvalueCutoff = 0.1,
               keyType = 'ENTREZID')

# Plot enriched GO terms using ggplot2 and dotplot
if (!is.null(go) && !is.na(go) && nrow(go) > 0) {
  ggplot(go, aes(x = GeneRatio, y = Description, color = p.adjust, size = GeneRatio)) +
    geom_point() +
    scale_color_gradient(low = "#1763a3", high = "#800000") +
    facet_grid(ONTOLOGY ~ ., scale = "free") +
    theme_bw(base_size = 18) +   # Set plot theme
    theme(text = element_text(size = 17))  # Adjust text size
} else {
  print("No enrichment GO terms found!")
}