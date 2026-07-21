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
#   --data_root Path to Main_Dataset dir (default: reads MLOmics Local Data.md or env var)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# ---------- defaults ----------
DATASET="GS-BRCA"
VERSION="Top"
METHOD="XGBoost"
DATA_ROOT="${MLOMICS_DATA_ROOT:-$PROJECT_ROOT/Main_Dataset}"

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

export MLOMICS_DATA_ROOT="$DATA_ROOT"
BASELINES_ROOT="$PROJECT_ROOT/Baseline_and_Metric/Classification/Baselines/Python"

case "$METHOD" in

  XGBoost|SVM|RF|LR)
    NB="$BASELINES_ROOT/XGBoost_ SVM_RF_LR/XGBoost_ SVM_RF_LR.ipynb"
    echo "Running $METHOD via Jupyter notebook..."
    cd "$BASELINES_ROOT/XGBoost_ SVM_RF_LR"
    jupyter nbconvert --to notebook --execute --inplace \
      --ExecutePreprocessor.timeout=600 \
      "XGBoost_ SVM_RF_LR.ipynb" \
      2>&1 | tail -5
    echo "Done. Output notebook in place."
    ;;

  MAUI)
    echo "Running MAUI vignette..."
    cd "$BASELINES_ROOT/MAUI"
    python -c "
import sys; sys.path.insert(0, '$PROJECT_ROOT')
from utils.data_loader import MLOmicsLoader
loader = MLOmicsLoader('$DATA_ROOT')
omics = loader.load_omics('$DATASET', version='$VERSION')
y     = loader.load_labels('$DATASET', version='$VERSION')
print('Loaded:', {k: v.shape for k, v in omics.items()})
from maui import Maui
maui_model = Maui(n_hidden=[900], n_latent=70, epochs=10)
maui_model.fit(omics)
z = maui_model.transform(omics)
print('Latent factors:', z.shape)
print('MAUI test run complete.')
"
    ;;

  XOmiVAE)
    echo "Running XOmiVAE..."
    cd "$BASELINES_ROOT/XOmiVAE"
    python -c "
import sys; sys.path.insert(0, '$PROJECT_ROOT')
from generalHelperFunctions import loadData_mlomics
class_num, labels, samples, omics = loadData_mlomics('$DATASET', version='$VERSION', data_root='$DATA_ROOT')
print(f'Dataset: $DATASET/$VERSION | classes={class_num} | samples={len(samples)}')
for k, arr in omics.items():
    print(f'  {k}: {arr.shape}')
print('XOmiVAE data loading OK.')
"
    ;;

  DCAP)
    echo "Running DCAP (requires R preprocessing first)..."
    cd "$BASELINES_ROOT/DCAP"
    CANCER="${DATASET#GS-}"
    INPUT_CSV="$DATA_ROOT/Classification_datasets/$DATASET/$VERSION/${CANCER}_mRNA_${VERSION,,}.csv"
    OUTPUT_CSV="/tmp/dcap_preprocessed.csv"
    echo "  Step 1: R preprocessing..."
    Rscript Data_preprocessing.R "$DATA_ROOT/Classification_datasets/$DATASET/$VERSION" \
      "${CANCER}_mRNA_${VERSION,,}.csv" \
      "$OUTPUT_CSV"
    echo "  Step 2: Python autoencoder..."
    python autoencoder_DCAP.py --data_path "$OUTPUT_CSV" --output_path /tmp/dcap_output.csv
    echo "DCAP run complete."
    ;;

  CustOmics)
    echo "Running CustOmics data loading test..."
    cd "$BASELINES_ROOT/CustOmics"
    python -c "
import sys; sys.path.insert(0, '$PROJECT_ROOT')
sys.path.insert(0, 'src')
from src.tools.utils import read_data_mlomics
omics, labels = read_data_mlomics('$DATASET', version='$VERSION', data_root='$DATA_ROOT')
print('CustOmics loaded:')
for k, df in omics.items():
    print(f'  {k}: {df.shape}')
print(f'  labels: {labels.shape}')
print('CustOmics data loading OK.')
"
    ;;

  *)
    echo "Error: unknown method '$METHOD'"
    echo "Supported: XGBoost SVM RF LR MAUI XOmiVAE DCAP CustOmics"
    exit 1
    ;;
esac

echo ""
echo "Classification run complete: $METHOD on $DATASET/$VERSION"
