#!/usr/bin/env Rscript
# run_NEMO.R — Wrapper for NEMO multi-omics clustering on MLOmics data
#
# Usage:
#   Rscript run_NEMO.R --dataset ACC --version Top --data_path /path/to/Clustering_datasets --k 2
#
# Outputs:
#   results_{dataset}_{version}_NEMO.csv  (in working directory)

suppressPackageStartupMessages(library(argparse))

# ---------- Argument parsing ----------
parser <- ArgumentParser(description = "NEMO multi-omics clustering")
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

# ---------- Load NEMO from local source ----------
# NEMO depends on SNFtool for affinityMatrix, dist2, spectralClustering
suppressPackageStartupMessages(library(SNFtool))

# Resolve script directory reliably when called via Rscript (handles paths with spaces)
script_dir <- tryCatch({
  args_all  <- commandArgs(trailingOnly = FALSE)
  file_arg  <- grep("--file=", args_all, value = TRUE)
  if (length(file_arg) > 0) {
    raw_path  <- sub("--file=", "", file_arg[1])
    raw_path  <- gsub("~\\+~", " ", raw_path)   # R encodes spaces as ~+~ in --file
    dirname(normalizePath(raw_path, mustWork = FALSE))
  } else {
    dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE))
  }
}, error = function(e) getwd())

nemo_r <- file.path(script_dir, "NEMO", "R", "NEMO.R")
if (file.exists(nemo_r)) {
  source(nemo_r)
} else {
  stop(sprintf("NEMO.R not found at expected path: %s", nemo_r))
}

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
# NEMO nemo.affinity.graph expects a list of matrices: features x samples
# (it transposes internally via t() before computing distances)
cat(sprintf("Loading %s %s data...\n", dataset, version))
mats <- lapply(paths, function(p) {
  df <- read.csv(p, row.names = 1, check.names = FALSE)
  as.matrix(df)   # features x samples — NEMO uses colnames as patient IDs
})

sample_names <- colnames(mats[[1]])
cat(sprintf("  Samples: %d\n", length(sample_names)))
cat(sprintf("  Features per omics: %s\n",
    paste(sapply(mats, nrow), collapse = " / ")))

# ---------- Run NEMO ----------
cat("Building NEMO affinity graph...\n")
affinity <- nemo.affinity.graph(mats)

cat(sprintf("Spectral clustering with k = %d...\n", k))
labels <- spectralClustering(affinity, k)
names(labels) <- rownames(affinity)

# Align to sample order from data
labels_ordered <- labels[sample_names]

# ---------- Save results ----------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
out_file <- file.path(output_dir, sprintf("results_%s_%s_NEMO.csv", dataset, version))
result_df <- data.frame(
  sample  = sample_names,
  cluster = as.integer(labels_ordered),
  row.names = NULL
)
write.csv(result_df, file = out_file, row.names = FALSE, quote = FALSE)
cat(sprintf("Results saved to: %s\n", out_file))
cat("Cluster distribution:\n")
print(table(labels_ordered))
