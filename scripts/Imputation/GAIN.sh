#!/bin/bash
# Run GAIN on an MLOmics imputation matrix.
# Usage: ./GAIN.sh [CANCER] [OMICS] [MISS_RATE]
# Example: ./GAIN.sh BRCA mRNA 0.3

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

arg1=${1:-"BRCA"}
arg2=${2:-"mRNA"}
arg3=${3:-0.3}

data_name="${arg1}_${arg2}"
miss_rate="$arg3"
data_file="$ROOT_DIR/Main_Dataset/Imputation_datasets/Imp-${arg1}/Top/${arg1}_${arg2}.csv"

if [[ ! -f "$data_file" ]]; then
  echo "Error: data file not found: $data_file"
  exit 1
fi

# batch = min(128, n_rows); iterations match unified runner default (100)
read -r eff_batch <<EOF
$(python3 -c "import pandas as pd; n=pd.read_csv(r'''$data_file''', index_col=0).shape[0]; print(max(1, min(128, int(n))))")
EOF

cd "$ROOT_DIR/Baseline_and_Metric/Imputation/GAIN"
output_file="$SCRIPT_DIR/GAIN.txt"
python3 main_letter_spam.py \
  --data_name "$data_name" \
  --data_path "$data_file" \
  --miss_rate "$miss_rate" \
  --batch_size "$eff_batch" \
  --hint_rate 0.9 \
  --alpha 100 \
  --iterations 100 \
  | tee -a "$output_file"
