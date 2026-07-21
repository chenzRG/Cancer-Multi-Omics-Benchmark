#!/usr/bin/env Rscript
# run_moCluster.R — Wrapper for moCluster (mogsa/mbpca) multi-omics clustering on MLOmics data
#
# mogsa package must be installed:
#   BiocManager::install("mogsa")
#
# Usage:
#   Rscript run_moCluster.R --dataset ACC --version Top --data_path /path/to/Clustering_datasets --k 2
#
# Outputs:
#   results_{dataset}_{version}_moCluster.csv  (in working directory)

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(mogsa))
suppressPackageStartupMessages(library(cluster))

# ---------- Argument parsing ----------
parser <- ArgumentParser(description = "moCluster multi-omics clustering")
parser$add_argument("--dataset",   required = TRUE,  help = "Dataset name (e.g. ACC)")
parser$add_argument("--version",   required = TRUE,  help = "Data version: Top, Original, Aligned")
parser$add_argument("--data_path", required = TRUE,  help = "Root path to Clustering_datasets/")
parser$add_argument("--k",          type = "integer", default = 2,  help = "Number of clusters (default: 2)")
parser$add_argument("--output_dir", default = ".",                  help = "Directory to write results CSV (default: .)")
args <- parser$parse_args()

dataset    <- args$dataset
version    <- args$version
data_path  <- args$data_path
k          <- args$k
output_dir <- args$output_dir

# ---------- File suffix per version ----------
suffix_map <- list(Original = "", Top = "_top", Aligned = "_aligned")
if (!version %in% names(suffix_map)) {
  stop(sprintf("Unknown version '%s'. Must be one of: %s", version, paste(names(suffix_map), collapse = ", ")))
}
suffix <- suffix_map[[version]]

# ---------- Build file paths ----------
base_dir <- file.path(data_path, dataset, version)

make_path <- function(omics) {
  file.path(base_dir, paste0(dataset, "_", omics, suffix, ".csv"))
}

omics_names <- c("mRNA", "miRNA", "CNV", "Methy")
paths <- setNames(lapply(omics_names, make_path), omics_names)

for (nm in omics_names) {
  if (!file.exists(paths[[nm]])) {
    stop(sprintf("File not found: %s", paths[[nm]]))
  }
}

# ---------- Read data ----------
# mogsa mbpca expects a named list of matrices: features x samples
# (mbpca transposes internally: x <- lapply(x, t))
cat(sprintf("Loading %s %s data...\n", dataset, version))
mats <- lapply(paths, function(p) {
  df <- read.csv(p, row.names = 1, check.names = FALSE)
  as.matrix(df)   # features x samples — as expected by mogsa
})

sample_names <- colnames(mats[[1]])
cat(sprintf("  Samples: %d\n", length(sample_names)))
cat(sprintf("  Features per omics: %s\n",
    paste(sapply(mats, nrow), collapse = " / ")))

# ---------- Run moCluster via mbpca ----------
# ncomp: number of components; k controls sparse loading sparsity.
# Use 'globalScore' method which integrates omics via shared sample scores.
n_comp <- min(k + 4, length(sample_names) - 1, 20)
cat(sprintf("Running mbpca with %d components...\n", n_comp))

moa_result <- mbpca(
  x       = mats,
  ncomp   = n_comp,
  method  = "globalScore",
  k       = "all",
  center  = TRUE,
  scale   = FALSE,
  option  = "lambda1",
  moa     = TRUE,
  verbose = FALSE,
  svd.solver = "fast.svd"
)

# Extract sample scores (n_samples x n_comp)
scores <- moaScore(moa_result)

# ---------- Hierarchical clustering on scores ----------
cat(sprintf("Clustering samples with k = %d...\n", k))
d <- dist(scores)
hcl <- hclust(d, method = "ward.D2")
labels <- cutree(hcl, k = k)

# ---------- Save results ----------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
out_file <- file.path(output_dir, sprintf("results_%s_%s_moCluster.csv", dataset, version))
result_df <- data.frame(
  sample  = sample_names,
  cluster = as.integer(labels),
  row.names = NULL
)
write.csv(result_df, file = out_file, row.names = FALSE, quote = FALSE)
cat(sprintf("Results saved to: %s\n", out_file))
cat("Cluster distribution:\n")
print(table(labels))
