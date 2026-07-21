#!/usr/bin/env Rscript
# run_DeepCC.R — Wrapper for DeepCC classification on MLOmics data
#
# DeepCC requires gene expression with Entrez gene ID column names.
# This wrapper converts MLOmics gene symbols → Entrez IDs via org.Hs.eg.db,
# then runs getFunctionalSpectra + train_DeepCC_model + get_DeepCC_label.
#
# Usage:
#   Rscript run_DeepCC.R --dataset GS-BRCA --version Top \
#       --data_path /path/to/Classification_datasets \
#       --output_dir /path/to/output [--epochs 5]
#
# Output:
#   {output_dir}/results_{dataset}_{version}_DeepCC.csv

# Point reticulate to the Python that has TensorFlow installed.
# Set the RETICULATE_PYTHON environment variable before calling this script,
# or ensure the Python with TF is first on PATH.
py_env <- Sys.getenv("RETICULATE_PYTHON", unset = "")
if (nchar(py_env) == 0) {
  py_env <- Sys.which("python3")
  if (nchar(py_env) > 0) Sys.setenv(RETICULATE_PYTHON = py_env)
}

suppressPackageStartupMessages({
  library(argparse)
  library(DeepCC)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(foreach)
  library(doParallel)
  # Use keras3 for training (R keras 2.x is incompatible with Python keras 3.x)
  library(keras3)
})

# ---------- Argument parsing ----------
parser <- ArgumentParser(description = "DeepCC cancer subtype classifier")
parser$add_argument("--dataset",    required = TRUE,  help = "Dataset name (e.g. GS-BRCA)")
parser$add_argument("--version",    required = TRUE,  help = "Data version: Top, Original, Aligned")
parser$add_argument("--data_path",  required = TRUE,  help = "Root path to Classification_datasets/")
parser$add_argument("--output_dir", default  = ".",   help = "Directory to write results CSV (default: .)")
parser$add_argument("--epochs",     type = "integer", default = 5, help = "Training epochs (default: 5)")
args <- parser$parse_args()

dataset    <- args$dataset
version    <- args$version
data_path  <- args$data_path
output_dir <- args$output_dir
epochs     <- args$epochs

# ---------- File suffix ----------
suffix_map <- list(Original = "", Top = "_top", Aligned = "_aligned")
if (!version %in% names(suffix_map)) stop(sprintf("Unknown version '%s'", version))
suffix <- suffix_map[[version]]

# Extract cancer type code (e.g. GS-BRCA → BRCA)
cancer <- sub("^GS-", "", dataset)

# ---------- Load mRNA expression ----------
# DeepCC primarily uses mRNA (gene expression) data
mrna_file <- file.path(data_path, dataset, version,
                       paste0(cancer, "_mRNA", suffix, ".csv"))
if (!file.exists(mrna_file)) stop(sprintf("mRNA file not found: %s", mrna_file))

cat(sprintf("Loading mRNA: %s\n", mrna_file))
mrna_df <- read.csv(mrna_file, row.names = 1, check.names = FALSE)
# mrna_df is features × samples; transpose to samples × features
eps <- t(mrna_df)  # samples × genes

cat(sprintf("Expression matrix: %d samples × %d genes\n", nrow(eps), ncol(eps)))

# ---------- Convert gene symbols → Entrez IDs ----------
cat("Converting gene symbols to Entrez IDs...\n")
gene_symbols <- colnames(eps)

entrez_map <- suppressMessages(
  AnnotationDbi::mapIds(org.Hs.eg.db,
                        keys    = gene_symbols,
                        column  = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
)

# Keep only genes with valid Entrez IDs
valid_idx <- !is.na(entrez_map)
cat(sprintf("  Mapped %d / %d gene symbols to Entrez IDs\n",
            sum(valid_idx), length(gene_symbols)))

eps_entrez <- eps[, valid_idx, drop = FALSE]
colnames(eps_entrez) <- entrez_map[valid_idx]

# Convert to data.frame as DeepCC expects
eps_df <- as.data.frame(eps_entrez)

# ---------- Load labels ----------
# Label file is always {cancer}_label_num.csv in the version folder (no omics suffix)
label_file <- file.path(data_path, dataset, version,
                        paste0(cancer, "_label_num.csv"))
if (!file.exists(label_file)) stop(sprintf("Label file not found: %s", label_file))

# Label file has one column "Label" with no sample ID column; order matches mRNA columns
label_df <- read.csv(label_file, check.names = FALSE)
labels <- as.character(label_df[[1]])
cat(sprintf("Labels: %d samples, %d classes\n",
            length(labels), length(unique(labels[!is.na(labels)]))))

# ---------- Compute Functional Spectra ----------
cat("Computing functional spectra (this may take a while)...\n")

# Try MSigDBv7 first; fall back to synthetic gene sets for smoke testing.
# data() issues a warning (not error) when the dataset is missing, so check
# whether the object was actually loaded rather than catching errors.
use_msigdb <- withCallingHandlers(tryCatch({
  e <- new.env()
  data("MSigDBv7", package = "DeepCC", envir = e)
  !is.null(e$MSigDBv7)
}, error = function(err) FALSE),
warning = function(w) invokeRestart("muffleWarning"))

if (use_msigdb) {
  fs <- getFunctionalSpectra(eps_df, geneSets = "MSigDBv7",
                             cores = max(1, parallel::detectCores() - 2))
} else {
  cat("  MSigDBv7 not bundled — computing spectra with synthetic gene sets\n")
  all_entrez <- colnames(eps_df)
  set.seed(42)
  n_sets <- 50
  set_size <- max(5, floor(length(all_entrez) / 10))
  gene_sets_obj <- lapply(seq_len(n_sets), function(i) {
    sample(all_entrez, min(set_size, length(all_entrez)))
  })
  names(gene_sets_obj) <- paste0("GeneSet_", seq_len(n_sets))

  # Replicate getFunctionalSpectra internals with a custom list
  eps_scaled <- scale(as.matrix(eps_df), scale = FALSE)
  n_cores <- max(1, parallel::detectCores() - 2)
  doParallel::registerDoParallel(n_cores)
  fs <- foreach::foreach(idx = 1:nrow(eps_scaled), .combine = rbind) %dopar% {
    geneList <- DeepCC:::preprocessGeneList(eps_scaled[idx, ])
    sapply(gene_sets_obj, function(x) DeepCC:::calcEnrichmentScore(geneList, x))
  }
  rownames(fs) <- rownames(eps_scaled)
}
cat(sprintf("Functional spectra: %d samples × %d gene sets\n", nrow(fs), ncol(fs)))

# ---------- Train DeepCC model ----------
# Use keras3 directly: R keras 2.x (2.16.1) is incompatible with Python keras 3.x.
# This replicates train_DeepCC_model + get_DeepCC_label logic with keras3 API.
cat(sprintf("Training DeepCC model (%d epochs)...\n", epochs))

x_train <- as.matrix(fs)
valid_idx <- !is.na(labels)
x_train_v <- x_train[valid_idx, , drop = FALSE]
y_fac <- factor(labels[valid_idx])
levs  <- levels(y_fac)
n_class <- length(levs)
y_int <- as.integer(y_fac) - 1L
y_oh  <- keras3::to_categorical(y_int, num_classes = n_class)

model <- keras3::keras_model_sequential(input_shape = ncol(x_train_v)) |>
  keras3::layer_dense(units = 256, activation = "relu") |>
  keras3::layer_batch_normalization() |>
  keras3::layer_dense(units = 64,  activation = "relu") |>
  keras3::layer_batch_normalization() |>
  keras3::layer_dense(units = n_class, activation = "softmax")

model |> keras3::compile(
  loss      = "categorical_crossentropy",
  optimizer = keras3::optimizer_adam(learning_rate = 0.001),
  metrics   = c("accuracy")
)
model |> keras3::fit(x_train_v, y_oh,
                     epochs           = epochs,
                     batch_size       = 64,
                     validation_split = 0.1,
                     verbose          = 0)

# ---------- Predict ----------
cat("Predicting labels...\n")
probs  <- model(x_train, training = FALSE)
pred_idx    <- apply(as.matrix(probs), 1, which.max)
pred_labels <- levs[pred_idx]

# ---------- Save results ----------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
out_file <- file.path(output_dir,
                      sprintf("results_%s_%s_DeepCC.csv", dataset, version))
result_df <- data.frame(
  sample    = rownames(fs),
  true      = labels,
  predicted = pred_labels,
  row.names = NULL
)
write.csv(result_df, file = out_file, row.names = FALSE, quote = FALSE)
cat(sprintf("Results saved to: %s\n", out_file))

# Print accuracy
acc <- mean(result_df$true == result_df$predicted, na.rm = TRUE)
cat(sprintf("Training accuracy: %.3f\n", acc))
