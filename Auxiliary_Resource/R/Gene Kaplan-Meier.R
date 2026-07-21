# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(cluster)
library(factoextra)
library(survival)
library(survminer)
library(gridExtra)  

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

# Filter genes of interest to those present in the expression data
genes_of_interest <- genes_of_interest[genes_of_interest %in% colnames(expression_data_merged)]

# Filter the expression data for these genes
filtered_data <- expression_data_merged %>%
  dplyr::select(Sample, all_of(genes_of_interest), survival_times, event_observed) %>%
  filter(complete.cases(.))

# Initialize an empty list to store ggplot objects
plot_list <- list()

# Loop through each gene of interest and plot Kaplan-Meier survival curves
for (gene in genes_of_interest) {
  
  # Normalize and categorize gene expression levels into "low" and "high"
  filtered_data_gene <- filtered_data_gene %>%
    mutate(expression_level = ifelse(.data[[gene]] <= median(.data[[gene]], na.rm = TRUE), "Low", "High"))
  
  # Create survival object
  surv_object <- tryCatch({
    Surv(filtered_data_gene$survival_times, filtered_data_gene$event_observed)
  })
  
  # Fit Kaplan-Meier model
  fit <- tryCatch({
    survfit(surv_object ~ expression_level, data = filtered_data_gene)
  })
  
  # Plot the survival curves using ggsurvplot and store in plot_list
  p <- ggsurvplot(fit, data = filtered_data_gene, pval = TRUE, 
                  title = sprintf("Kaplan-Meier Survival Curve for %s", gene),
                  legend.title = "Expression Level",
                  legend.labs = c("Low", "High"),
                  xlab = "Time (days)",
                  ylab = "Survival Probability")
  
  # Extract the ggplot object from ggsurvplot and store in plot_list
  plot_list[[gene]] <- p$plot
}

# Arrange all ggplot objects in a grid using cowplot::plot_grid
multiplot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)

# Display the grid of plots
print(multiplot)