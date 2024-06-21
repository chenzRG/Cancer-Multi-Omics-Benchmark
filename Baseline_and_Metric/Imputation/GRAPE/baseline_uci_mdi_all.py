import time
import argparse
import sys
import os
import os.path as osp

import numpy as np
import torch
import pandas as pd

from training.baseline import baseline_mdi
from uci.uci_data import load_data
from utils.utils import auto_select_gpu

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--domain', type=str, default='uci')
    parser.add_argument('--data', type=str, default='housing')
    parser.add_argument('--train_edge', type=float, default=0.7)
    parser.add_argument('--split_sample', type=float, default=0.)
    parser.add_argument('--split_by', type=str, default='y') # 'y', 'random'
    parser.add_argument('--split_train', action='store_true', default=False)
    parser.add_argument('--split_test', action='store_true', default=False)
    parser.add_argument('--train_y', type=float, default=0.7)
    parser.add_argument('--node_mode', type=int, default=0)  # 0: feature onehot, sample all 1; 1: all onehot

    parser.add_argument('--method', type=str, default='mean')
    parser.add_argument('--level', type=int, default=0)
    parser.add_argument('--seed', type=int, default=0)

    parser.add_argument('--best_level', action='store_true', default=False)
    parser.add_argument('--comment', type=str, default='v1')
    args = parser.parse_args()
    # device is cpu by default
    if torch.cuda.is_available():
        cuda = auto_select_gpu()
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ['CUDA_VISIBLE_DEVICES'] = str(cuda)
        print('Using GPU {}'.format(os.environ['CUDA_VISIBLE_DEVICES']))
        device = torch.device('cuda:{}'.format(cuda))
    else:
        print('Using CPU')
        device = torch.device('cpu')

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    missing_rates = [0.3, 0.5, 0.7]

    best_levels = {'mean':0, 'knn':2, 'svd':2, 'mice':0, 'spectral':1}

    for args.data in ['BCRA','COAD','GBM','LGG','OV']:
        for missing_rate in missing_rates:
            args.train_edge = missing_rate
            print("missing rate = ", missing_rate)
            if args.best_level:
                args.level = best_levels[args.method]
                print("using best level {} for {}".format(args.level,args.method))
            data = load_data(args)
            print("data_prepared")
            log_path = './uci/mdi_results/results/{}_{}/{}/{}/'.format(args.method, args.comment, args.data, args.seed)
            if not os.path.isdir(log_path):
                os.makedirs(log_path)
            for args.method in ['mean', 'knn', 'svd', 'mice', 'spectral']:
                baseline_mdi(data, args, log_path)
                print("NEXT")

if __name__ == '__main__':
    main()

