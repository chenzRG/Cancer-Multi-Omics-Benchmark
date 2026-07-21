#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATASET=${1:-"GS-BRCA"}
VERSION=${2:-"Top"}
# MLOmicsLoader expects Main_Dataset root (also accepts Classification_datasets and normalizes)
DATA_PATH=${3:-"$ROOT_DIR/Main_Dataset"}
OUTPUT_DIR=${4:-"$ROOT_DIR/results/classification"}
CUST_DIR="$ROOT_DIR/Baseline_and_Metric/Classification/Baselines/Python/CustOmics"
mkdir -p "$OUTPUT_DIR"
python3 "$CUST_DIR/main.py" \
  --cohorts "$DATASET" --version "$VERSION" \
  --data_directory "$DATA_PATH" --result_directory "$OUTPUT_DIR" \
  --task classification --sources mRNA,miRNA,CNV,Methy
