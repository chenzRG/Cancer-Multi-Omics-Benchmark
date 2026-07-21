# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(cluster)
library(factoextra)
library(survival)
library(survminer)

# Data root: set MLOMICS_DATA_ROOT or default to repo Main_Dataset
data_root <- Sys.getenv(
  "MLOMICS_DATA_ROOT",
  unset = normalizePath(file.path("..", "..", "Main_Dataset"), mustWork = FALSE)
)

# Construct file paths using sprintf for string interpolation
cancer_type <- 'THYM'
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
genes_of_interest <- c("CELF5", "ODZ1", "CD1C", "DRP2", "PTCRA", "TSHR", "HKDC1", "KCTD19", "RFX8", "ZNF600", "UGT3A2", "PRKCG")

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

# Determine the optimal number of clusters using the Elbow method
fviz_nbclust(clustering_data, kmeans, method = "wss") + 
  geom_vline(xintercept = 2, linetype = 2) +
  labs(subtitle = "Elbow method")

# Set the seed for reproducibility
set.seed(121)

# Perform k-means clustering with k = 2
kmeans_result <- kmeans(clustering_data, centers = 2, nstart = 25)

# Add the cluster assignment to the data
filtered_data$Cluster <- kmeans_result$cluster

# Use the clustering results as labels
labels <- as.data.frame(kmeans_result$cluster)

# Rename the first column to 'Label'
colnames(labels)[1] <- 'Label'

# Add a column for Sample names to the labels data frame
labels$Sample <- colnames(expression_data)

# Merge expression data with labels
expression_data_merged <- merge(expression_data_t, labels, by = "Sample")

# PCA plot for visualization
pca <- prcomp(clustering_data, scale. = TRUE)
pca_data <- as.data.frame(pca$x)
pca_data$Cluster <- as.factor(filtered_data$Cluster)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "PCA of mRNA Expression Data",
       x = "Principal Component 1",
       y = "Principal Component 2")

# Merge the labeled expression with the survival data
final_data <- merge(expression_data_merged, survival, by = "Sample")

# Create a survival object
surv_obj <- with(final_data, Surv(survival_times, event_observed))

# Fit Kaplan-Meier curves for each cluster
km_fit <- survfit(surv_obj ~ Label, data = final_data)

# Plot the Kaplan-Meier curves
ggsurvplot(km_fit, data = final_data, pval = TRUE, conf.int = TRUE, legend.labs = c("Cluster 1", "Cluster 2"),
           title = "Kaplan-Meier Curves by Clusters",
           xlab = "Time (days)", ylab = "Survival Probability",
           palette = c("#E41A1C", "#377EB8", "#4DAF4A"))  # Adjust colors as needed
