#!/usr/bin/env Rscript
# run_CIMLR.R — Wrapper for CIMLR multi-omics clustering on MLOmics data
#
# CIMLR R package must be installed:
#   devtools::install_github("danro9685/CIMLR", ref = "R")
#
# Usage:
#   Rscript run_CIMLR.R --dataset ACC --version Top --data_path /path/to/Clustering_datasets --k 2
#
# Outputs:
#   results_{dataset}_{version}_CIMLR.csv  (in working directory)

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(CIMLR))

# ---------- Argument parsing ----------
parser <- ArgumentParser(description = "CIMLR multi-omics clustering")
parser$add_argument("--dataset",   required = TRUE,  help = "Dataset name (e.g. ACC)")
parser$add_argument("--version",   required = TRUE,  help = "Data version: Top, Original, Aligned")
parser$add_argument("--data_path", required = TRUE,  help = "Root path to Clustering_datasets/")
parser$add_argument("--k",         type = "integer", default = 2, help = "Number of clusters (default: 2)")
args <- parser$parse_args()

dataset   <- args$dataset
version   <- args$version
data_path <- args$data_path
k         <- args$k

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

# ---------- Read and transpose data ----------
# Input CSVs: features x samples  →  CIMLR expects samples x features
cat(sprintf("Loading %s %s data...\n", dataset, version))
mats <- lapply(paths, function(p) {
  df <- read.csv(p, row.names = 1, check.names = FALSE)
  t(as.matrix(df))   # samples x features
})

sample_names <- rownames(mats[[1]])
cat(sprintf("  Samples: %d\n", length(sample_names)))
cat(sprintf("  Features per omics: %s\n",
    paste(sapply(mats, ncol), collapse = " / ")))

# ---------- Concatenate omics ----------
# CIMLR takes a single concatenated data matrix (samples x total_features)
X <- do.call(cbind, mats)
cat(sprintf("  Combined features: %d\n", ncol(X)))

# ---------- Run CIMLR ----------
cat(sprintf("Running CIMLR with k = %d...\n", k))
result <- CIMLR(X, k, nstart = 10)

labels <- result$y$cluster

# ---------- Save results ----------
out_file <- sprintf("results_%s_%s_CIMLR.csv", dataset, version)
result_df <- data.frame(
  sample  = sample_names,
  cluster = as.integer(labels),
  row.names = NULL
)
write.csv(result_df, file = out_file, row.names = FALSE, quote = FALSE)
cat(sprintf("Results saved to: %s\n", out_file))
cat("Cluster distribution:\n")
print(table(labels))
