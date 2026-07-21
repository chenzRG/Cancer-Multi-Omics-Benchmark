#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATASET=${1:-"GS-BRCA"}
VERSION=${2:-"Top"}
DATA_PATH=${3:-"$ROOT_DIR/Main_Dataset/Classification_datasets"}
OUTPUT_DIR=${4:-"$ROOT_DIR/results/classification"}
EPOCHS=${5:-10}
mkdir -p "$OUTPUT_DIR"
Rscript "$ROOT_DIR/Baseline_and_Metric/Classification/Baselines/Python/DeepCC/run_DeepCC.R" \
  --dataset "$DATASET" --version "$VERSION" \
  --data_path "$DATA_PATH" --output_dir "$OUTPUT_DIR" --epochs "$EPOCHS"
