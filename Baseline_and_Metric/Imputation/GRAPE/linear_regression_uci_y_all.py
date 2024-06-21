import time
import argparse
import sys
import os
import os.path as osp

import numpy as np
import torch
import pandas as pd

from training.linear_regression import linear_regression
from uci.uci_data import load_data

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

    parser.add_argument('--method', type=str, default=None)
    parser.add_argument('--level', type=int, default=0)
    parser.add_argument('--seed', type=int, default=0)

    parser.add_argument('--best_level', action='store_true', default=False)
    parser.add_argument('--comment', type=str, default='v1')
    args = parser.parse_args()
    # device is cpu by default

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    best_levels = {'mean':0, 'knn':2, 'svd':2, 'mice':0, 'spectral':1}

    ## new
    for args.data in ['concrete', 'energy', 'housing', 'kin8nm',
                    'naval', 'power', 'protein', 'wine', 'yacht']:
        for args.method in ['gnn','gnn_mdi', 'mean', 'knn', 'svd', 'mice', 'spectral']:
            data = load_data(args)
            if args.best_level and (not args.method.startswith('gnn')):
                args.level = best_levels[args.method]
                print("using best level {} for {}".format(args.level,args.method))
            log_path = './uci/y_results/results/{}_{}/{}/{}/'.format(args.method, args.comment, args.data, args.seed)
            if not os.path.isdir(log_path):
                os.makedirs(log_path)
            if args.method == 'gnn':
                load_path = './uci/y_results/results/gnn_{}/{}/{}/'.format(args.comment, args.data, args.seed)
            elif args.method == 'gnn_mdi':
                load_path = './uci/mdi_results/results/gnn_mdi_{}/{}/{}/'.format(args.comment, args.data, args.seed)
            else:
                load_path = ''
            linear_regression(data, args, log_path, load_path)

if __name__ == '__main__':
    main()

