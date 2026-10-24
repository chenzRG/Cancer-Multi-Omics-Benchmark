#!/usr/bin/env bash
# Run an imputation baseline on an MLOmics imputation dataset.
#
# Usage:
#   ./scripts/run_imputation.sh --dataset Imp-BRCA --omics mRNA --method GAIN --missing 0.3
#   ./scripts/run_imputation.sh --dataset Imp-BRCA --omics mRNA --method Mean --missing 0.3
#   ./scripts/run_imputation.sh --dataset Imp-COAD --omics CNV  --method GRAPE --iterations 100
#
# Supported --method values:
#   Mean  KNN  MICE  SVD  Spectral  GAIN  GRAPE
#
# Options:
#   --dataset     Imp-BRCA|Imp-COAD|Imp-GBM|Imp-LGG|Imp-OV
#   --omics       mRNA | miRNA | CNV | Methy  (default: mRNA)
#   --method      Mean | KNN | MICE | SVD | Spectral | GAIN | GRAPE
#   --missing     Missing rate 0.3|0.5|0.7  (default: 0.3)
#   --iterations  GAIN training iterations (default: 100)
#   --batch_size  GAIN batch size before clamp (default: 128; clamped to n_rows)
#   --data_root   Path to Main_Dataset (default: auto)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

DATASET="Imp-BRCA"
OMICS="mRNA"
METHOD="GAIN"
MISSING="0.3"
ITERATIONS="100"
BATCH_SIZE="128"
DATA_ROOT="${MLOMICS_DATA_ROOT:-$PROJECT_ROOT/Main_Dataset}"
OUTPUT_DIR="$PROJECT_ROOT/results/imputation"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)    DATASET="$2";    shift 2 ;;
    --omics)      OMICS="$2";      shift 2 ;;
    --method)     METHOD="$2";     shift 2 ;;
    --missing)    MISSING="$2";    shift 2 ;;
    --iterations) ITERATIONS="$2"; shift 2 ;;
    --batch_size) BATCH_SIZE="$2"; shift 2 ;;
    --data_root)  DATA_ROOT="$2";  shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

CANCER="${DATASET#Imp-}"
DATA_FILE="$DATA_ROOT/Imputation_datasets/$DATASET/Top/${CANCER}_${OMICS}.csv"
BASELINES_ROOT="$PROJECT_ROOT/Baseline_and_Metric/Imputation"
mkdir -p "$OUTPUT_DIR"

echo "======================================"
echo "MLOmics Imputation Baseline Runner"
echo "  Dataset : $DATASET"
echo "  Omics   : $OMICS"
echo "  Method  : $METHOD"
echo "  Missing : $MISSING"
echo "  File    : $DATA_FILE"
echo "======================================"

if [[ ! -f "$DATA_FILE" ]]; then
  echo "Error: data file not found: $DATA_FILE"; exit 1
fi

# GAIN: effective batch = min(requested, n_rows in CSV body)
gain_batch_size() {
  python3 -c "
import pandas as pd
n = pd.read_csv(r'''$DATA_FILE''', index_col=0).shape[0]
print(max(1, min(int('$BATCH_SIZE'), int(n))))
"
}

case "$METHOD" in

  Mean|KNN|MICE|SVD|Spectral)
    # map display names to run_classical_impute.py choices
    case "$METHOD" in
      Mean)     CLASSICAL="mean" ;;
      KNN)      CLASSICAL="knn" ;;
      MICE)     CLASSICAL="mice" ;;
      SVD)      CLASSICAL="svd" ;;
      Spectral) CLASSICAL="spectral" ;;
    esac
    echo "Running $METHOD via run_classical_impute.py..."
    python3 "$BASELINES_ROOT/run_classical_impute.py" \
      --dataset "$DATASET" \
      --omics "$OMICS" \
      --miss_rate "$MISSING" \
      --method "$CLASSICAL" \
      --data_path "$DATA_ROOT" \
      --output_dir "$OUTPUT_DIR"
    ;;

  GAIN)
    EFF_BATCH="$(gain_batch_size)"
    echo "Running GAIN (batch_size=$EFF_BATCH, iterations=$ITERATIONS)..."
    cd "$BASELINES_ROOT/GAIN"
    python3 main_letter_spam.py \
      --data_path "$DATA_FILE" \
      --miss_rate "$MISSING" \
      --batch_size "$EFF_BATCH" \
      --hint_rate 0.9 \
      --alpha 100 \
      --iterations "$ITERATIONS"
    ;;

  GRAPE)
    echo "Running GRAPE..."
    cd "$BASELINES_ROOT/GRAPE"
    python3 train_mdi.py \
      --epochs 20000 \
      --known "$(python3 -c "print(1 - $MISSING)")" \
      uci \
      --data_path "$DATA_FILE" \
      --data "$DATASET" \
      --train_edge 0.7
    ;;

  *)
    echo "Error: unknown method '$METHOD'"
    echo "Supported: Mean KNN MICE SVD Spectral GAIN GRAPE"
    exit 1
    ;;
esac

echo ""
echo "Imputation run complete: $METHOD on $DATASET/$OMICS (missing=$MISSING)"
