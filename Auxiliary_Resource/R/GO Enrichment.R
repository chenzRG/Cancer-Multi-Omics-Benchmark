# Load necessary libraries
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)

# Data root: set MLOMICS_DATA_ROOT or default to repo Main_Dataset
data_root <- Sys.getenv(
  "MLOMICS_DATA_ROOT",
  unset = normalizePath(file.path("..", "..", "Main_Dataset"), mustWork = FALSE)
)

# Construct file paths using sprintf for string interpolation
cancer_type <- 'BRCA'
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

# Perform gene ID conversion with bitr() to get Entrez IDs
hg<-bitr(selected_genes$Gene,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")

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
