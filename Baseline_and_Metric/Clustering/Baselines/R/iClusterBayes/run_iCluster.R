#!/usr/bin/env Rscript
# run_iCluster.R — Wrapper for iClusterBayes multi-omics clustering on MLOmics data
#
# iClusterPlus package must be installed:
#   BiocManager::install("iClusterPlus")
#
# Usage:
#   Rscript run_iCluster.R --dataset ACC --version Top --data_path /path/to/Clustering_datasets --k 2
#
# Outputs:
#   results_{dataset}_{version}_iCluster.csv  (in working directory)

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(iClusterPlus))

# ---------- Argument parsing ----------
parser <- ArgumentParser(description = "iClusterBayes multi-omics clustering")
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

# ---------- Read and transpose data ----------
# iClusterBayes expects matrices: samples x features (n x p)
cat(sprintf("Loading %s %s data...\n", dataset, version))
mats <- lapply(paths, function(p) {
  df <- read.csv(p, row.names = 1, check.names = FALSE)
  t(as.matrix(df))   # samples x features
})

sample_names <- rownames(mats[[1]])
cat(sprintf("  Samples: %d\n", length(sample_names)))
cat(sprintf("  Features per omics: %s\n",
    paste(sapply(mats, ncol), collapse = " / ")))

# iClusterBayes can be slow with thousands of features; use top variance features
select_top_var <- function(X, max_feats = 2000) {
  if (ncol(X) <= max_feats) return(X)
  vars <- apply(X, 2, var, na.rm = TRUE)
  X[, order(vars, decreasing = TRUE)[1:max_feats], drop = FALSE]
}

cat("Selecting top-variance features (max 2000 per omics)...\n")
mats_sel <- lapply(mats, select_top_var, max_feats = 2000)
cat(sprintf("  Features after selection: %s\n",
    paste(sapply(mats_sel, ncol), collapse = " / ")))

# ---------- Run iClusterBayes ----------
# k clusters requires K = k-1 latent factors
K_latent <- k - 1
cat(sprintf("Running iClusterBayes with k = %d (K_latent = %d)...\n", k, K_latent))

# All four omics treated as continuous (Gaussian)
fit <- iClusterBayes(
  dt1   = mats_sel[["mRNA"]],
  dt2   = mats_sel[["miRNA"]],
  dt3   = mats_sel[["CNV"]],
  dt4   = mats_sel[["Methy"]],
  type  = c("gaussian", "gaussian", "gaussian", "gaussian"),
  K     = K_latent,
  n.burnin = 1800,
  n.draw   = 1200,
  prior.gamma = c(0.5, 0.5, 0.5, 0.5),
  sdev  = 0.05,
  beta.var.scale = 1,
  thin  = 3
)

# Extract cluster assignments (1-based)
labels <- fit$clusters

# ---------- Save results ----------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
out_file <- file.path(output_dir, sprintf("results_%s_%s_iCluster.csv", dataset, version))
result_df <- data.frame(
  sample  = sample_names,
  cluster = as.integer(labels),
  row.names = NULL
)
write.csv(result_df, file = out_file, row.names = FALSE, quote = FALSE)
cat(sprintf("Results saved to: %s\n", out_file))
cat("Cluster distribution:\n")
print(table(labels))
