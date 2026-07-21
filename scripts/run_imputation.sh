#!/usr/bin/env bash
# Run an imputation baseline on an MLOmics imputation dataset.
#
# Usage:
#   ./scripts/run_imputation.sh --dataset Imp-BRCA --omics mRNA --method GAIN --missing 0.3
#   ./scripts/run_imputation.sh --dataset Imp-COAD --omics CNV  --method GRAPE
#
# Supported --method values:
#   GAIN  GRAPE
#
# Options:
#   --dataset    Imp-BRCA|Imp-COAD|Imp-GBM|Imp-LGG|Imp-OV
#   --omics      mRNA | miRNA | CNV | Methy  (default: mRNA)
#   --method     GAIN | GRAPE
#   --missing    Missing rate 0.3|0.5|0.7  (default: 0.3)
#   --data_root  Path to Main_Dataset (default: auto)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

DATASET="Imp-BRCA"
OMICS="mRNA"
METHOD="GAIN"
MISSING="0.3"
DATA_ROOT="${MLOMICS_DATA_ROOT:-$PROJECT_ROOT/Main_Dataset}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)   DATASET="$2";   shift 2 ;;
    --omics)     OMICS="$2";     shift 2 ;;
    --method)    METHOD="$2";    shift 2 ;;
    --missing)   MISSING="$2";   shift 2 ;;
    --data_root) DATA_ROOT="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

CANCER="${DATASET#Imp-}"
DATA_FILE="$DATA_ROOT/Imputation_datasets/$DATASET/Top/${CANCER}_${OMICS}.csv"
BASELINES_ROOT="$PROJECT_ROOT/Baseline_and_Metric/Imputation"

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

case "$METHOD" in

  GAIN)
    cd "$BASELINES_ROOT/GAIN"
    python main_letter_spam.py \
      --data_path "$DATA_FILE" \
      --miss_rate "$MISSING" \
      --batch_size 128 \
      --hint_rate 0.9 \
      --alpha 100 \
      --iterations 10000
    ;;

  GRAPE)
    cd "$BASELINES_ROOT/GRAPE"
    python train_mdi.py \
      --epochs 20000 \
      --known "$(python3 -c "print(1 - $MISSING)")" \
      uci \
      --data_path "$DATA_FILE" \
      --data "$DATASET" \
      --train_edge 0.7
    ;;

  *)
    echo "Error: unknown method '$METHOD'"
    echo "Supported: GAIN GRAPE"
    exit 1
    ;;
esac

echo ""
echo "Imputation run complete: $METHOD on $DATASET/$OMICS (missing=$MISSING)"
