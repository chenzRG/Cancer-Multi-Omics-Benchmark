#!/bin/bash

# Save the current directory
current_dir=$(pwd)

# Navigate to the specific folder
cd ../../Baseline_and_Metric/Imputation/GAIN

# Define the output file
output_file="$current_dir/GAIN.txt"

# Run all commands and append their output to the output file
python3 main_letter_spam.py --data_name BCRA --miss_rate 0.7 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name BCRA --miss_rate 0.5 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name BCRA --miss_rate 0.3 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file

python3 main_letter_spam.py --data_name COAD --miss_rate 0.7 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name COAD --miss_rate 0.5 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name COAD --miss_rate 0.3 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file

python3 main_letter_spam.py --data_name GBM --miss_rate 0.7 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name GBM --miss_rate 0.5 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name GBM --miss_rate 0.3 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file

python3 main_letter_spam.py --data_name LGG --miss_rate 0.7 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name LGG --miss_rate 0.5 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name LGG --miss_rate 0.3 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file

python3 main_letter_spam.py --data_name OV --miss_rate 0.7 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name OV --miss_rate 0.5 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file
python3 main_letter_spam.py --data_name OV --miss_rate 0.3 --batch_size 128 --hint_rate 0.9 --alpha 100 --iterations 100 >> $output_file


# Navigate back to the original directory
cd "$current_dir"