# Load necessary libraries
library(tidyr)
library(dplyr)
library(survival)
library(ggplot2)

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

# Define survival time and status columns
time_col <- "survival_times" 
status_col <- "event_observed" 

# Function to sanitize gene names
sanitize_gene_name <- function(gene_name) {
  make.names(gene_name, unique = TRUE)
}

# Sanitize column names in the merged data frame
colnames(expression_data_merged) <- sapply(colnames(expression_data_merged), sanitize_gene_name)

# Initialize an empty list to store Cox regression results
cox_results <- list()

# Perform univariate Cox regression for each gene
for (gene in colnames(expression_data_t)[-ncol(expression_data_t)]) {
  sanitized_gene <- sanitize_gene_name(gene)
  if(sanitized_gene %in% colnames(expression_data_merged)) {
    formula <- as.formula(paste("Surv(", time_col, ",", status_col, ") ~", sanitized_gene))
    cox_model <- tryCatch({
      coxph(formula, data = expression_data_merged)
    }, error = function(e) {
      message(paste("Error with gene", gene, ":", e))
      NULL
    })
    if (!is.null(cox_model)) {
      cox_summary <- summary(cox_model)
      cox_results[[gene]] <- cox_summary$coefficients
    }
  } else {
    message(paste("Gene", gene, "not found in the merged data frame"))
  }
}

# Convert the results to a data frame
cox_results_df <- do.call(rbind, cox_results)
cox_results_df <- as.data.frame(cox_results_df)
names(cox_results_df) <- c("coef", "exp_coef", "se_coef", "z", "p")

# Add a column for gene names
cox_results_df$gene <- rownames(cox_results_df)

# Sort the results by p-value and select the top 15
top_genes_df <- cox_results_df %>% arrange(p) %>% head(15)

# Create a forest plot for the top 15 genes
ggplot(top_genes_df, aes(x = exp_coef, y = gene)) +
  geom_point(shape = 15, color = "green", size = 3) +
  geom_errorbarh(aes(xmin = exp_coef - 1.96 * se_coef, xmax = exp_coef + 1.96 * se_coef), height = 0.2, color = "blue") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Hazard ratio", y = "Gene", title = "Forest Plot of Top 15 Cox Regression Results") +
  theme_minimal()

# Display the top 15 results
print(top_genes_df)