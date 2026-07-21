library(tidyr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(png)
library(grid)

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

#Merge the data for pathway analysis
colnames(hg) <- c("Gene", "ENTREZID")
info_merge <- merge(DEGAll,hg,by='Gene')
PW_input <- info_merge$logFC
names(PW_input) = info_merge$ENTREZID
PW_input = sort(PW_input, decreasing = TRUE)

#Pathway analysis with 'pathview'
pv.out <- pathview(gene.data = PW_input,
                   pathway.id = "04064",
                   species = "hsa",
                   out.suffix = "pathway",
                   kegg.native = TRUE,
                   gene.idtype = "ENTREZID", 
                   same.layer = FALSE,
                   low = list(cpd = "#1763a3"),   
                   high = list(cpd = "#800000"),
                   kegg.dir = getwd()  # Use current working directory for output
)

# Construct the full path to the saved image
output_file <- file.path(getwd(), paste0("hsa04064.pathway.png"))

# Read the PNG file
img <- readPNG(output_file)

# Display the image using grid functions
grid::grid.raster(img)
