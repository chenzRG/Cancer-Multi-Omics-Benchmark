#!/usr/bin/env bash
# Smoke-test all baselines with a small representative run.
# Useful for CI or quick sanity check after environment setup.
#
# Usage:
#   ./scripts/run_all.sh [--data_root /path/to/Main_Dataset]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DATA_ROOT="${MLOMICS_DATA_ROOT:-$PROJECT_ROOT/Main_Dataset}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data_root) DATA_ROOT="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

export MLOMICS_DATA_ROOT="$DATA_ROOT"
PASS=0; FAIL=0

run() {
  local label="$1"; shift
  echo -n "  [$label] ... "
  if "$@" > /tmp/mlomics_run.log 2>&1; then
    echo "PASS"
    ((PASS++)) || true
  else
    echo "FAIL  (see /tmp/mlomics_run.log)"
    ((FAIL++)) || true
  fi
}

echo "======================================"
echo "MLOmics Full Baseline Smoke Test"
echo "Data: $DATA_ROOT"
echo "======================================"

echo ""
echo "[Classification]"
run "XGBoost"   bash "$SCRIPT_DIR/run_classification.sh" --dataset GS-BRCA  --version Top     --method CustOmics
run "MAUI"      bash "$SCRIPT_DIR/run_classification.sh" --dataset GS-COAD  --version Top     --method XOmiVAE
run "XOmiVAE"   bash "$SCRIPT_DIR/run_classification.sh" --dataset GS-BRCA  --version Top     --method XOmiVAE
run "CustOmics" bash "$SCRIPT_DIR/run_classification.sh" --dataset GS-BRCA  --version Top     --method CustOmics

echo ""
echo "[Clustering — Python]"
run "SubtypeGAN" bash "$SCRIPT_DIR/run_clustering.sh" --dataset ACC --version Top --method SubtypeGAN --k 3
run "MCluster"   bash "$SCRIPT_DIR/run_clustering.sh" --dataset ACC --version Top --method MCluster   --k 3

echo ""
echo "[Clustering — R]"
run "SNF"       bash "$SCRIPT_DIR/run_clustering.sh" --dataset ACC --version Top --method SNF       --k 3
run "NEMO"      bash "$SCRIPT_DIR/run_clustering.sh" --dataset ACC --version Top --method NEMO      --k 3
run "CIMLR"     bash "$SCRIPT_DIR/run_clustering.sh" --dataset ACC --version Top --method CIMLR     --k 3
run "iCluster"  bash "$SCRIPT_DIR/run_clustering.sh" --dataset ACC --version Top --method iCluster  --k 3
run "moCluster" bash "$SCRIPT_DIR/run_clustering.sh" --dataset ACC --version Top --method moCluster --k 3

echo ""
echo "[Imputation]"
run "GAIN"  bash "$SCRIPT_DIR/run_imputation.sh" --dataset Imp-BRCA --omics mRNA --method GAIN  --missing 0.3
run "GRAPE" bash "$SCRIPT_DIR/run_imputation.sh" --dataset Imp-BRCA --omics mRNA --method GRAPE --missing 0.3

echo ""
echo "======================================"
echo "Results: $PASS passed, $FAIL failed"
echo "======================================"
[[ $FAIL -eq 0 ]]
