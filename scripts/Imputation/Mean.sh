#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATASET=${1:-"Imp-BRCA"}
OMICS=${2:-"mRNA"}
MISS_RATE=${3:-0.3}
DATA_PATH=${4:-"$ROOT_DIR/Main_Dataset/Imputation_datasets"}
OUTPUT_DIR=${5:-"$ROOT_DIR/results/imputation"}
mkdir -p "$OUTPUT_DIR"
python3 "$ROOT_DIR/Baseline_and_Metric/Imputation/run_classical_impute.py" \
  --dataset "$DATASET" --omics "$OMICS" \
  --miss_rate "$MISS_RATE" --method "mean" \
  --data_path "$DATA_PATH" --output_dir "$OUTPUT_DIR"
