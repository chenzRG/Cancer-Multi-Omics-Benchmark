import os
import os.path as osp
import argparse
import numpy as np
import pandas as pd
import torch
from fancyimpute import KNN, NuclearNormMinimization, SoftImpute, BiScaler, IterativeImputer, SimpleFill, IterativeSVD
from training.baseline import baseline_mdi
from utils.utils import auto_select_gpu
import logging
import inspect

def setup_logging(log_path):
    if not os.path.exists(log_path):
        os.makedirs(log_path)
    logging.basicConfig(filename=os.path.join(log_path, 'experiment.log'),
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def load_data(args):
    uci_path = osp.dirname(osp.abspath(inspect.getfile(inspect.currentframe())))
    file_path = osp.join(uci_path, 'raw_data/{}/data.csv'.format(args.data))
    df = pd.read_csv(file_path, skiprows=1)
    df = df.drop(df.columns[0], axis=1)
    X = df.values
    
    return X

def generate_incomplete_data(X, missing_rate, seed):
    np.random.seed(seed)
    X_incomplete = X.copy()
    n_rows, n_cols = X.shape
    missing_mask = np.zeros(X.shape, dtype=bool)

    for col in range(n_cols):
        missing_indices = np.random.choice(n_rows, int(missing_rate * n_rows), replace=False)
        missing_mask[missing_indices, col] = True

    X_incomplete[missing_mask] = np.nan
    return X_incomplete, missing_mask

def baseline_impute(X_incomplete, method='mean', level=0):
    if method == 'mean':
        X_filled_mean = SimpleFill().fit_transform(X_incomplete)
        return X_filled_mean
    elif method == 'knn':
        k = [3, 10, 50][level]
        X_filled_knn = KNN(k=k, verbose=False).fit_transform(X_incomplete)
        return X_filled_knn
    elif method == 'svd':
        rank = [np.ceil((X_incomplete.shape[1] - 1) / 10), np.ceil((X_incomplete.shape[1] - 1) / 5), X_incomplete.shape[1] - 1][level]
        X_filled_svd = IterativeSVD(rank=int(rank), verbose=False).fit_transform(X_incomplete)
        return X_filled_svd
    elif method == 'mice':
        max_iter = [3, 10, 50][level]
        X_filled_mice = IterativeImputer(max_iter=max_iter).fit_transform(X_incomplete)
        return X_filled_mice
    elif method == 'spectral':
        sparsity = [0.5, None, 3][level]
        X_filled_spectral = SoftImpute(shrinkage_value=sparsity).fit_transform(X_incomplete)
        return X_filled_spectral
    else:
        raise NotImplementedError

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

    # Set device
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

    best_levels = {'mean': 0, 'knn': 2, 'svd': 2, 'mice': 0, 'spectral': 1}

    missing_rates = [0.1, 0.3, 0.5, 0.7]

    for args.data in ['BCRA', 'COAD', 'GBM', 'LGG', 'OV']:
        for args.method in ['mean', 'knn', 'svd', 'mice', 'spectral']:
            for missing_rate in missing_rates:
                log_path = './uci/mdi_results/results/{}_{}_{}_{}_{}/'.format(args.method, args.comment, args.data, args.seed, missing_rate)
                setup_logging(log_path)

                if args.best_level:
                    args.level = best_levels[args.method]
                    logging.info("Using best level {} for {}".format(args.level, args.method))

                # Load and process data
                X = load_data(args)
                X_incomplete, missing_mask = generate_incomplete_data(X, missing_rate, args.seed)
                logging.info("Data prepared with missing rate {}".format(missing_rate))

                # Impute missing values
                X_filled = baseline_impute(X_incomplete, method=args.method, level=args.level)

                # Evaluate imputation
                mse = ((X_filled[missing_mask] - X[missing_mask]) ** 2).mean()
                logging.info("{} Impute MSE for missing rate {}: {}".format(args.method, missing_rate, mse))

                # Ensure log path exists
                if not os.path.isdir(log_path):
                    os.makedirs(log_path)

                # Save imputed data and log results
                np.save(osp.join(log_path, 'X_filled.npy'), X_filled)
                np.save(osp.join(log_path, 'X_incomplete.npy'), X_incomplete)
                np.save(osp.join(log_path, 'missing_mask.npy'), missing_mask)

                baseline_mdi((X, X_filled), args, log_path)
                logging.info("NEXT")

if __name__ == '__main__':
    main()
