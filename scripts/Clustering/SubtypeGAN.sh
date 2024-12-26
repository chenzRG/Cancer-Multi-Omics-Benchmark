#!/bin/bash


# Save the current directory
current_dir=$(pwd)

# Navigate to the specific folder
cd ../../../Classification_and_Clustering/Python/Subtype-GAN

# Record the start time
echo "Script started at: $(date)" >> results/$1_$2.cc


python SubtypeGAN.py -m SubtypeGAN -n 4 -t $1_$2 # set dummy number of clusters to 4 
python SubtypeGAN.py -m cc -t $1_$2

# Record the end time
echo "Script ended at: $(date)" >> results/$1_$2.cc

# Navigate back to the original directory
cd "$current_dir"