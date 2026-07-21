#!/usr/bin/env python3
"""
run_dcap.py — DCAP (Denoising Autoencoder for Cancer Prediction) on MLOmics data.

Reimplements DAE-DCAP.py using Keras (TF2), replacing the TF1 Session API.

Usage:
    python run_dcap.py --dataset GS-BRCA --version Top \
        --data_path /path/to/Classification_datasets \
        --output_dir /path/to/output [--epochs 50]
"""
import argparse, os, sys

def main():
    parser = argparse.ArgumentParser(description='DCAP denoising autoencoder')
    parser.add_argument('--dataset',    required=True)
    parser.add_argument('--version',    default='Top')
    parser.add_argument('--data_path',  required=True)
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--epochs',     type=int, default=50)
    parser.add_argument('--batch_size', type=int, default=64)
    parser.add_argument('--noise_std',  type=float, default=0.1)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir   = os.path.normpath(os.path.join(script_dir, '../../../../..'))
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)

    import numpy as np
    import pandas as pd
    import tensorflow as tf
    from sklearn.cluster import KMeans
    from utils.data_loader import MLOmicsLoader

    data_root = args.data_path
    if os.path.basename(data_root) == 'Classification_datasets':
        data_root = os.path.dirname(data_root)
    loader = MLOmicsLoader(data_root)
    omics  = loader.load_omics(args.dataset, version=args.version)
    labels = loader.load_labels(args.dataset, version=args.version)
    sample_ids = list(list(omics.values())[0].columns)

    X = np.concatenate([df.T.values for df in omics.values()], axis=1).astype('float32')
    print(f'Input shape: {X.shape}')

    n_input   = X.shape[1]
    n_hidden1 = min(200, n_input)
    n_hidden2 = 50
    n_latent  = 10

    # Denoising autoencoder: add Gaussian noise during training
    inp   = tf.keras.Input(shape=(n_input,))
    noisy = tf.keras.layers.GaussianNoise(args.noise_std)(inp)
    enc   = tf.keras.layers.Dense(n_hidden1, activation='tanh')(noisy)
    enc   = tf.keras.layers.Dense(n_hidden2, activation='tanh')(enc)
    latent = tf.keras.layers.Dense(n_latent, activation='tanh')(enc)
    dec   = tf.keras.layers.Dense(n_hidden2, activation='tanh')(latent)
    dec   = tf.keras.layers.Dense(n_hidden1, activation='tanh')(dec)
    out   = tf.keras.layers.Dense(n_input)(dec)

    dae = tf.keras.Model(inp, out)
    dae.compile(optimizer=tf.keras.optimizers.Adam(1e-4), loss='mse')
    dae.fit(X, X, epochs=args.epochs, batch_size=args.batch_size,
            validation_split=0.1, verbose=1)

    encoder_model = tf.keras.Model(inp, latent)
    Z = encoder_model.predict(X)

    n_classes = len(np.unique(labels[~pd.isna(labels)]))
    km = KMeans(n_clusters=n_classes, n_init=20, random_state=0)
    clusters = km.fit_predict(Z)

    os.makedirs(args.output_dir, exist_ok=True)
    out_file = os.path.join(args.output_dir,
                            f'results_{args.dataset}_{args.version}_DCAP.csv')
    result_df = pd.DataFrame({'sample': sample_ids, 'true_label': labels,
                               'cluster': clusters,
                               **{f'z{i}': Z[:, i] for i in range(n_latent)}})
    result_df.to_csv(out_file, index=False)
    print(f'Results saved to: {out_file}')
    print(f'Latent shape: {Z.shape}')

if __name__ == '__main__':
    main()
