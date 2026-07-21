#!/usr/bin/env python3
"""
run_classical_impute.py — Classical imputation methods on MLOmics data.

Methods: mean, knn, mice, svd, spectral

Usage:
    python run_classical_impute.py --dataset Imp-BRCA --omics mRNA \
        --miss_rate 0.3 --method knn \
        --data_path /path/to/Imputation_datasets \
        --output_dir /path/to/output
"""
import argparse, os, sys

METHODS = ['mean', 'knn', 'mice', 'svd', 'spectral']

def main():
    parser = argparse.ArgumentParser(description='Classical imputation')
    parser.add_argument('--dataset',    required=True, help='e.g. Imp-BRCA')
    parser.add_argument('--omics',      default='mRNA',
                        choices=['mRNA', 'miRNA', 'CNV', 'Methy'])
    parser.add_argument('--miss_rate',  type=float, default=0.3)
    parser.add_argument('--method',     required=True, choices=METHODS)
    parser.add_argument('--data_path',  required=True)
    parser.add_argument('--output_dir', default='.')
    args = parser.parse_args()

    import numpy as np
    import pandas as pd
    from sklearn.metrics import mean_absolute_error
    from fancyimpute import SimpleFill, KNN, IterativeImputer, IterativeSVD, SoftImpute

    # Load data
    cancer = args.dataset.replace('Imp-', '')
    fpath  = os.path.join(args.data_path, args.dataset, 'Top',
                          f'{cancer}_{args.omics}.csv')
    if not os.path.exists(fpath):
        raise FileNotFoundError(f'Data file not found: {fpath}')

    df = pd.read_csv(fpath, index_col=0)
    X  = df.T.values.astype('float64')       # samples × features
    print(f'Loaded {args.dataset}/{args.omics}: {X.shape}')

    # Introduce missingness
    rng = np.random.default_rng(42)
    mask = rng.random(X.shape) < args.miss_rate
    X_obs = X.copy()
    X_obs[mask] = np.nan
    print(f'Missing rate: {mask.mean():.3f} ({mask.sum()} values)')

    # Impute
    if args.method == 'mean':
        X_imp = SimpleFill().fit_transform(X_obs)
    elif args.method == 'knn':
        X_imp = KNN(k=5, verbose=False).fit_transform(X_obs)
    elif args.method == 'mice':
        X_imp = IterativeImputer(max_iter=10).fit_transform(X_obs)
    elif args.method == 'svd':
        rank = max(1, X_obs.shape[1] // 10)
        X_imp = IterativeSVD(rank=rank, verbose=False).fit_transform(X_obs)
    elif args.method == 'spectral':
        X_imp = SoftImpute(verbose=False).fit_transform(X_obs)

    # Evaluate on masked positions only
    mae  = mean_absolute_error(X[mask], X_imp[mask])
    rmse = np.sqrt(np.mean((X[mask] - X_imp[mask]) ** 2))
    remaining_nan = int(np.isnan(X_imp).sum())
    print(f'MAE={mae:.6f}  RMSE={rmse:.6f}  remaining_nan={remaining_nan}')

    os.makedirs(args.output_dir, exist_ok=True)
    out_file = os.path.join(args.output_dir,
        f'results_{args.dataset}_{args.omics}_miss{args.miss_rate}_{args.method}.csv')
    result_df = pd.DataFrame({
        'dataset': [args.dataset], 'omics': [args.omics],
        'miss_rate': [args.miss_rate], 'method': [args.method],
        'MAE': [mae], 'RMSE': [rmse]
    })
    result_df.to_csv(out_file, index=False)
    print(f'Results saved to: {out_file}')

if __name__ == '__main__':
    main()
