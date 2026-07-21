# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

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

# Create a DEGAll data frame with necessary columns for the volcano plot
DEGAll <- results %>%
  mutate(
    logFC = log2FoldChange,
    PValue = pValue,
    color = ifelse(pValue < 0.05 & abs(log2FoldChange) > 0.5,
                   ifelse(log2FoldChange > 0.5, "red", "blue"), "gray")
  )

# Define colors
color <- c(red = "#800000", gray = "#A9A9A9", blue = "#4682B433")

# Define the number of top genes to label
top_n <- 10

# Filter top genes by p-value and log2 fold change
top_genes <- DEGAll %>%
  arrange(pValue) %>%
  slice_head(n = top_n)

# Plot the volcano plot
ggplot(DEGAll, aes(logFC, -log10(PValue), col = color)) +
  geom_point() +
  geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.5,
    segment.color = 'grey50',
    color = 'black'  # Set gene label color to black
  ) +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x = "log2 (fold change)", y = "-log10 (q-value)") +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "darkgreen", lwd = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), lty = 4, col = "darkgreen", lwd = 0.6) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )