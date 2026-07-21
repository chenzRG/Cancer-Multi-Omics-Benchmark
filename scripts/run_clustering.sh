#!/usr/bin/env bash
# Run a clustering baseline on an MLOmics clustering dataset.
#
# Usage:
#   ./scripts/run_clustering.sh --dataset ACC --version Top --method SNF --k 3
#   ./scripts/run_clustering.sh --dataset LUAD --version Aligned --method SubtypeGAN
#
# Supported --method values:
#   SNF  NEMO  CIMLR  iCluster  moCluster      (R-based)
#   SubtypeGAN  MCluster                       (Python-based)
#
# Options:
#   --dataset   ACC|KIRC|KIRP|LIHC|LUAD|LUSC|PRAD|THCA|THYM
#   --version   Original | Top | Aligned  (default: Top)
#   --method    Baseline name
#   --k         Number of clusters (default: 3)
#   --data_root Path to Main_Dataset/Clustering_datasets (default: auto)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

DATASET="ACC"
VERSION="Top"
METHOD="SNF"
K=3
DATA_ROOT="${MLOMICS_DATA_ROOT:-$PROJECT_ROOT/Main_Dataset}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)   DATASET="$2";   shift 2 ;;
    --version)   VERSION="$2";   shift 2 ;;
    --method)    METHOD="$2";    shift 2 ;;
    --k)         K="$2";         shift 2 ;;
    --data_root) DATA_ROOT="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

CLUST_DATA="$DATA_ROOT/Clustering_datasets"
BASELINES_ROOT="$PROJECT_ROOT/Baseline_and_Metric/Clustering/Baselines"

echo "======================================"
echo "MLOmics Clustering Baseline Runner"
echo "  Dataset : $DATASET"
echo "  Version : $VERSION"
echo "  Method  : $METHOD"
echo "  k       : $K"
echo "  Data    : $CLUST_DATA"
echo "======================================"

case "$METHOD" in

  SNF)
    Rscript "$BASELINES_ROOT/R/SNF/run_SNF.R" \
      --dataset "$DATASET" --version "$VERSION" \
      --data_path "$CLUST_DATA" --k "$K"
    ;;

  NEMO)
    Rscript "$BASELINES_ROOT/R/NEMO/run_NEMO.R" \
      --dataset "$DATASET" --version "$VERSION" \
      --data_path "$CLUST_DATA" --k "$K"
    ;;

  CIMLR)
    Rscript "$BASELINES_ROOT/R/CIMLR/run_CIMLR.R" \
      --dataset "$DATASET" --version "$VERSION" \
      --data_path "$CLUST_DATA" --k "$K"
    ;;

  iCluster)
    Rscript "$BASELINES_ROOT/R/iClusterBayes/run_iCluster.R" \
      --dataset "$DATASET" --version "$VERSION" \
      --data_path "$CLUST_DATA" --k "$K"
    ;;

  moCluster)
    Rscript "$BASELINES_ROOT/R/moCluster/run_moCluster.R" \
      --dataset "$DATASET" --version "$VERSION" \
      --data_path "$CLUST_DATA" --k "$K"
    ;;

  SubtypeGAN)
    cd "$BASELINES_ROOT/Python/Subtype-GAN"
    python SubtypeGAN.py \
      --dataset "$DATASET" --version "$VERSION" \
      --data_path "$CLUST_DATA" -n "$K"
    ;;

  MCluster)
    cd "$BASELINES_ROOT/Python/MCluster-VAEs"
    python run_mlomics.py \
      --dataset "$DATASET" --version "$VERSION" \
      --data_path "$CLUST_DATA" --k "$K"
    ;;

  *)
    echo "Error: unknown method '$METHOD'"
    echo "Supported: SNF NEMO CIMLR iCluster moCluster SubtypeGAN MCluster"
    exit 1
    ;;
esac

echo ""
echo "Clustering run complete: $METHOD on $DATASET/$VERSION (k=$K)"
