#!/usr/bin/env bash
# Run a classification baseline on an MLOmics dataset.
#
# Usage:
#   ./scripts/run_classification.sh --dataset GS-BRCA --version Top --method XGBoost
#   ./scripts/run_classification.sh --dataset Pan-cancer --version Original --method MAUI
#
# Supported --method values:
#   XGBoost  SVM  RF  LR  MAUI  XOmiVAE  DCAP  CustOmics
#
# Options:
#   --dataset   GS-BRCA | GS-COAD | GS-GBM | GS-LGG | GS-OV | Pan-cancer
#   --version   Original | Top | Aligned  (default: Top)
#   --method    Baseline name (default: XGBoost)
#   --data_root Path to Main_Dataset dir (default: env MLOMICS_DATA_ROOT or ./Main_Dataset)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# ---------- defaults ----------
DATASET="GS-BRCA"
VERSION="Top"
METHOD="XGBoost"
DATA_ROOT="${MLOMICS_DATA_ROOT:-$PROJECT_ROOT/Main_Dataset}"
OUTPUT_DIR="$PROJECT_ROOT/results/classification"

# ---------- parse args ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)   DATASET="$2";   shift 2 ;;
    --version)   VERSION="$2";   shift 2 ;;
    --method)    METHOD="$2";    shift 2 ;;
    --data_root) DATA_ROOT="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

echo "======================================"
echo "MLOmics Classification Baseline Runner"
echo "  Dataset : $DATASET"
echo "  Version : $VERSION"
echo "  Method  : $METHOD"
echo "  Data    : $DATA_ROOT"
echo "======================================"

if [[ ! -d "$DATA_ROOT" ]]; then
  echo "Error: data root not found: $DATA_ROOT"
  echo "Run ./download.sh first (or set --data_root / MLOMICS_DATA_ROOT)."
  exit 1
fi

export MLOMICS_DATA_ROOT="$DATA_ROOT"
BASELINES_ROOT="$PROJECT_ROOT/Baseline_and_Metric/Classification/Baselines/Python"
mkdir -p "$OUTPUT_DIR"

case "$METHOD" in

  XGBoost|SVM|RF|LR)
    echo "Running $METHOD via run_classifiers.py..."
    python3 "$BASELINES_ROOT/XGBoost_SVM_RF_LR/run_classifiers.py" \
      --dataset "$DATASET" \
      --version "$VERSION" \
      --data_path "$DATA_ROOT" \
      --output_dir "$OUTPUT_DIR" \
      --method "$METHOD"
    ;;

  MAUI)
    echo "Running MAUI via run_maui.py..."
    if ! python3 -c "import maui" 2>/dev/null; then
      echo "Error: MAUI is not installed."
      echo "Install with: pip install git+https://github.com/BIMSBbioinfo/maui.git"
      echo "(See notes in requirements.txt.)"
      exit 1
    fi
    python3 "$BASELINES_ROOT/MAUI/run_maui.py" \
      --dataset "$DATASET" \
      --version "$VERSION" \
      --data_path "$DATA_ROOT" \
      --output_dir "$OUTPUT_DIR"
    ;;

  XOmiVAE)
    echo "Running XOmiVAE via run_xomivae.py..."
    python3 "$BASELINES_ROOT/XOmiVAE/run_xomivae.py" \
      --dataset "$DATASET" \
      --version "$VERSION" \
      --data_path "$DATA_ROOT" \
      --output_dir "$OUTPUT_DIR"
    ;;

  DCAP)
    echo "Running DCAP via run_dcap.py..."
    python3 "$BASELINES_ROOT/DCAP/run_dcap.py" \
      --dataset "$DATASET" \
      --version "$VERSION" \
      --data_path "$DATA_ROOT" \
      --output_dir "$OUTPUT_DIR"
    ;;

  CustOmics)
    echo "Running CustOmics via main.py..."
    python3 "$BASELINES_ROOT/CustOmics/main.py" \
      --cohorts "$DATASET" \
      --version "$VERSION" \
      --data_directory "$DATA_ROOT" \
      --result_directory "$OUTPUT_DIR" \
      --task classification \
      --sources mRNA,miRNA,CNV,Methy
    ;;

  *)
    echo "Error: unknown method '$METHOD'"
    echo "Supported: XGBoost SVM RF LR MAUI XOmiVAE DCAP CustOmics"
    exit 1
    ;;
esac

echo ""
echo "Classification run complete: $METHOD on $DATASET/$VERSION"
