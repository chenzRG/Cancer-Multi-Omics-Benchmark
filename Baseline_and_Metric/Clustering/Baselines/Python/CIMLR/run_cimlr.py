#!/usr/bin/env python3
"""
run_cimlr.py — CIMLR multi-omics clustering on MLOmics data

Python implementation of:
  CIMLR (Cancer Integration via Multikernel LeaRning)
  Rappoport & Shamir, NAR 2018 / Wang et al. Nature Methods 2017

Usage:
    python run_cimlr.py --dataset ACC --version Top \
        --data_path /path/to/Clustering_datasets --k 3 \
        --output_dir /path/to/output

Output:
    {output_dir}/results_{dataset}_{version}_CIMLR.csv
"""

import argparse, os, sys

def main():
    parser = argparse.ArgumentParser(description='CIMLR multi-omics clustering')
    parser.add_argument('--dataset',    required=True, help='Dataset name (e.g. ACC)')
    parser.add_argument('--version',    required=True, help='Data version: Top, Original, Aligned')
    parser.add_argument('--data_path',  required=True, help='Root path to Clustering_datasets/')
    parser.add_argument('--k',          type=int, default=2, help='Number of clusters (default: 2)')
    parser.add_argument('--n_neighbors',type=int, default=10, help='CIMLR k-neighbours (default: 10)')
    parser.add_argument('--max_iter',   type=int, default=30, help='Optimisation iterations (default: 30)')
    parser.add_argument('--output_dir', default='.', help='Output directory (default: .)')
    args = parser.parse_args()

    import numpy as np
    import pandas as pd

    # Locate this script's directory to import cimlr.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    from cimlr import CIMLR

    # ── File suffix ──────────────────────────────────────────────────────────
    suffix_map = {'Original': '', 'Top': '_top', 'Aligned': '_aligned'}
    if args.version not in suffix_map:
        raise ValueError(f"Unknown version '{args.version}'")
    suffix = suffix_map[args.version]

    # ── Load omics data ──────────────────────────────────────────────────────
    base_dir = os.path.join(args.data_path, args.dataset, args.version)
    omics_types = ['mRNA', 'miRNA', 'CNV', 'Methy']
    omics_list  = []
    sample_ids  = None

    for ot in omics_types:
        fpath = os.path.join(base_dir, f'{args.dataset}_{ot}{suffix}.csv')
        if not os.path.exists(fpath):
            # Try without dataset prefix (clustering datasets use cancer code)
            cancer = args.dataset  # e.g. ACC
            fpath  = os.path.join(base_dir, f'{cancer}_{ot}{suffix}.csv')
        if not os.path.exists(fpath):
            print(f'  Skipping {ot}: file not found at {fpath}')
            continue
        df = pd.read_csv(fpath, index_col=0)
        # features × samples → samples × features
        X = df.T.values.astype('float32')
        X = np.nan_to_num(X, nan=0.0)
        if sample_ids is None:
            sample_ids = list(df.columns)
        print(f'  Loaded {ot}: {X.shape}')
        omics_list.append(X)

    if not omics_list:
        raise RuntimeError('No omics files found — check data_path and dataset name')

    print(f'Running CIMLR: {len(omics_list)} omics, {omics_list[0].shape[0]} samples, '
          f'k={args.k}, n_neighbors={args.n_neighbors}, max_iter={args.max_iter}')

    # ── Run CIMLR ────────────────────────────────────────────────────────────
    model = CIMLR(n_clusters=args.k, n_neighbors=args.n_neighbors,
                  max_iter=args.max_iter)
    labels, S, F = model.fit_predict(omics_list)

    # ── Save results ─────────────────────────────────────────────────────────
    os.makedirs(args.output_dir, exist_ok=True)
    out_file = os.path.join(args.output_dir,
                            f'results_{args.dataset}_{args.version}_CIMLR.csv')
    result_df = pd.DataFrame({'sample': sample_ids, 'cluster': labels})
    result_df.to_csv(out_file, index=False)
    print(f'Results saved to: {out_file}')
    print(f'Cluster distribution: {dict(zip(*np.unique(labels, return_counts=True)))}')


if __name__ == '__main__':
    main()
