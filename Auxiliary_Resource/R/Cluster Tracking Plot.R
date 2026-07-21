# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(cluster)
library(factoextra)

# Data root: set MLOMICS_DATA_ROOT or default to repo Main_Dataset
data_root <- Sys.getenv(
  "MLOMICS_DATA_ROOT",
  unset = normalizePath(file.path("..", "..", "Main_Dataset"), mustWork = FALSE)
)

# Construct file paths using sprintf for string interpolation
cancer_type <- 'KIRP'
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
genes_of_interest <- c("IGF2BP3", "KIAA1429", "HNRNPC")

# Filter genes of interest to those present in the expression data
genes_of_interest <- genes_of_interest[genes_of_interest %in% colnames(expression_data_merged)]

# Filter the expression data for these genes
filtered_data <- expression_data_merged %>%
  dplyr::select(Sample, all_of(genes_of_interest))

# Remove samples with missing gene data
filtered_data <- filtered_data %>%
  filter(complete.cases(.))

# Remove the Sample column for clustering
clustering_data <- filtered_data %>%
  dplyr::select(-Sample)

# Perform k-means clustering for k = 2 to 9 and store results
cluster_results <- lapply(2:9, function(k) {
  kmeans(clustering_data, centers = k, nstart = 25)
})

# Extract cluster assignments for each k into a data frame
cluster_df <- lapply(seq_along(cluster_results), function(i) {
  data.frame(Sample = filtered_data$Sample, k = i, cluster = cluster_results[[i]]$cluster)
})

# Combine all results into a single data frame
cluster_df <- do.call(rbind, cluster_df)

# Use tidyr::gather to reshape data
melted_cluster_df <- cluster_df %>%
  tidyr::gather(variable, value, -Sample, -k)

# Generate a color palette with enough colors
num_clusters <- length(unique(melted_cluster_df$value))
color_palette <- rainbow(num_clusters)

# Determine the order of samples within each cluster k
sample_order <- melted_cluster_df %>%
  arrange(k, value) %>%
  group_by(k) %>%
  mutate(order = row_number()) %>%
  ungroup() %>%
  select(Sample, k, order)

# Merge the order information back to melted_cluster_df
melted_cluster_df <- merge(melted_cluster_df, sample_order, by = c("Sample", "k"))

# Plot the tracking plot with corrected color scale and ordered samples
ggplot(melted_cluster_df, aes(x = interaction(variable, order, lex.order = TRUE), y = k, fill = as.factor(value))) +
  geom_tile() +
  scale_fill_manual(values = color_palette) +
  labs(title = "Tracking Plot for K-means Clustering",
       x = "Samples",
       y = "Number of Clusters (k)",
       fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  
        axis.title.x = element_blank(),  
        legend.position = "bottom")    