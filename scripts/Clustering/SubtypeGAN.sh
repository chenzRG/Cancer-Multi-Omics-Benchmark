#!/bin/bash

task=$1

# Save the current directory
current_dir=$(pwd)

# Navigate to the specific folder
cd ../../../Classification_and_Clustering/Python/Subtype-GAN

# Record the start time
echo "Script started at: $(date)" >> results/$2_$3.cc


python SubtypeGAN.py -m SubtypeGAN -n 4 -t $2_$3 # set dummy number of clusters to 4 
python SubtypeGAN.py -m cc -t $2_$3

# Record the end time
echo "Script ended at: $(date)" >> results/$2_$3.cc

# Navigate back to the original directory
cd "$current_dir"