#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATASET=${1:-"ACC"}
VERSION=${2:-"Top"}
DATA_PATH=${3:-"$ROOT_DIR/Main_Dataset/Clustering_datasets"}
OUTPUT_DIR=${4:-"$ROOT_DIR/results/clustering"}
K=${5:-2}
MCLUST_DIR="$ROOT_DIR/Baseline_and_Metric/Clustering/Baselines/Python/MCluster-VAEs"
mkdir -p "$OUTPUT_DIR"
cd "$MCLUST_DIR"
python3 run_mlomics.py \
  --dataset "$DATASET" --version "$VERSION" \
  --data_path "$DATA_PATH" --k "$K" --output_dir "$OUTPUT_DIR"
