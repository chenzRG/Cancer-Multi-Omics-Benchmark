#!/bin/bash

# Save the current directory
current_dir=$(pwd)

# Navigate to the specific folder
cd ../Baseline_and_Metric/Imputation/GRAPE

# Ensure the logs directory exists
if [ ! -d "./logs" ]; then
    mkdir ./logs
fi

if [ ! -d "./logs/baseline" ]; then
    mkdir ./logs/baseline
fi

# Record the start time
echo "Script started at: $(date)" > logs/baseline/missing.log

# Run the Python script and redirect the output
python baseline_uci_mdi_all.py >> logs/baseline/missing.log 2>&1

# Record the end time
echo "Script ended at: $(date)" >> logs/baseline/missing.log

# Navigate back to the original directory
cd "$current_dir"