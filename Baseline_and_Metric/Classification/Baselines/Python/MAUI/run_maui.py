#!/usr/bin/env python3
"""
run_maui.py — MAUI multi-omics classification on MLOmics data.

Usage:
    python run_maui.py --dataset GS-BRCA --version Top \
        --data_path /path/to/Classification_datasets \
        --output_dir /path/to/output [--n_latent 20] [--epochs 50]
"""
import argparse, os, sys

def main():
    parser = argparse.ArgumentParser(description='MAUI multi-omics VAE')
    parser.add_argument('--dataset',    required=True)
    parser.add_argument('--version',    default='Top')
    parser.add_argument('--data_path',  required=True)
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--n_latent',   type=int, default=20)
    parser.add_argument('--epochs',     type=int, default=50)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir   = os.path.normpath(os.path.join(script_dir, '../../../../..'))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)

    import numpy as np
    import pandas as pd
    from sklearn.cluster import KMeans
    from maui import Maui
    from utils.data_loader import MLOmicsLoader

    data_root = args.data_path
    if os.path.basename(data_root) == 'Classification_datasets':
        data_root = os.path.dirname(data_root)
    loader = MLOmicsLoader(data_root)
    omics  = loader.load_omics(args.dataset, version=args.version)
    labels = loader.load_labels(args.dataset, version=args.version)
    sample_ids = list(list(omics.values())[0].columns)

    print(f'Loaded {len(omics)} omics, {len(sample_ids)} samples')
    model = Maui(n_hidden=[128, 64], n_latent=args.n_latent, epochs=args.epochs)
    model.fit(omics)
    z = model.transform(omics)

    os.makedirs(args.output_dir, exist_ok=True)
    out_file = os.path.join(args.output_dir,
                            f'results_{args.dataset}_{args.version}_MAUI.csv')
    result_df = pd.DataFrame({'sample': sample_ids, 'label': labels,
                               **{f'z{i}': z.values[:, i] for i in range(min(5, z.shape[1]))}})
    result_df.to_csv(out_file, index=False)
    print(f'Results saved to: {out_file}')
    print(f'Latent shape: {z.shape}')

if __name__ == '__main__':
    main()
