#!/usr/bin/env bash

task=$1

python SubtypeGAN.py -m SubtypeGAN -n 4 -t $2_$3
python SubtypeGAN.py -m cc -t $2_$3