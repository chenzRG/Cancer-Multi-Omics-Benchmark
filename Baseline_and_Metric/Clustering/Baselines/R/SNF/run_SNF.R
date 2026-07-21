#!/usr/bin/env Rscript
# run_SNF.R — Wrapper for SNF-based multi-omics clustering on MLOmics data
#
# Usage:
#   Rscript run_SNF.R --dataset ACC --version Top --data_path /path/to/Clustering_datasets --k 2
#
# Outputs:
#   results_{dataset}_{version}_SNF.csv  (in working directory)

suppressPackageStartupMessages(library(argparse))

# ---------- Argument parsing ----------
parser <- ArgumentParser(description = "SNF multi-omics clustering")
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

# ---------- Load SNF package (installed or from local source) ----------
if (!requireNamespace("SNFtool", quietly = TRUE)) {
  # Fall back to sourcing local R files
  snf_dir <- tryCatch(
    dirname(normalizePath(sys.frames()[[1]]$ofile)),
    error = function(e) getwd()
  )
  for (f in list.files(file.path(snf_dir, "R"), pattern = "\\.R$|\\.r$", full.names = TRUE)) {
    source(f)
  }
} else {
  suppressPackageStartupMessages(library(SNFtool))
}

# ---------- File suffix per version ----------
suffix_map <- list(Original = "", Top = "_top", Aligned = "_aligned")
if (!version %in% names(suffix_map)) {
  stop(sprintf("Unknown version '%s'. Must be one of: %s", version, paste(names(suffix_map), collapse = ", ")))
}
suffix <- suffix_map[[version]]

# ---------- Build file paths ----------
# Data files live directly under {data_path}/{dataset}/
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
# Input CSVs: features x samples  →  transpose to samples x features
cat(sprintf("Loading %s %s data...\n", dataset, version))
mats <- lapply(paths, function(p) {
  df <- read.csv(p, row.names = 1, check.names = FALSE)
  t(as.matrix(df))   # now: samples x features
})

# All must have same row names (samples)
sample_names <- rownames(mats[[1]])
for (nm in omics_names[-1]) {
  if (!identical(rownames(mats[[nm]]), sample_names)) {
    stop(sprintf("Sample mismatch between mRNA and %s", nm))
  }
}
cat(sprintf("  Samples: %d\n", length(sample_names)))
cat(sprintf("  Features per omics: %s\n",
    paste(sapply(mats, ncol), collapse = " / ")))

# ---------- Build per-omics affinity matrices ----------
cat("Building affinity matrices...\n")
K_snf   <- 20    # number of neighbors for SNF kernel
sigma   <- 0.5   # bandwidth

aff_list <- lapply(mats, function(X) {
  D <- dist2(X, X)
  affinityMatrix(D, K = K_snf, sigma = sigma)
})

# ---------- SNF fusion ----------
cat("Running SNF...\n")
W_fused <- SNF(aff_list, K = K_snf, t = 20)

# ---------- Spectral clustering ----------
cat(sprintf("Spectral clustering with k = %d...\n", k))
labels <- spectralClustering(W_fused, K = k)

# ---------- Save results ----------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
out_file <- file.path(output_dir, sprintf("results_%s_%s_SNF.csv", dataset, version))
result_df <- data.frame(
  sample  = sample_names,
  cluster = labels,
  row.names = NULL
)
write.csv(result_df, file = out_file, row.names = FALSE, quote = FALSE)
cat(sprintf("Results saved to: %s\n", out_file))
cat(sprintf("Cluster distribution:\n"))
print(table(labels))
