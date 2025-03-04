#!/bin/bash

# Save the current directory
current_dir=$(pwd)

# Navigate to the specific folder
cd ../../Baseline_and_Metric/Imputation/GAIN

# 引数を取得
arg1=${1:-"BRCA"}  
arg2=${2:-"CNV"}   
arg3=${3:-0.7}     

# data_name を組み立てる
data_name="${arg1}_${arg2}"
miss_rate="$arg3"

# Define the output file
output_file="$current_dir/GAIN.txt"

# Run all commands and append their output to the output file
python3 main_letter_spam.py --data_name "$data_name" --miss_rate "$miss_rate" --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> "$output_file"

# Navigate back to the original directory
cd "$current_dir"
