#!/usr/bin/env python3
"""
run_xomivae.py — XOmiVAE multi-omics VAE on MLOmics data.

Usage:
    python run_xomivae.py --dataset GS-BRCA --version Top \
        --data_path /path/to/Classification_datasets \
        --output_dir /path/to/output [--epochs 50] [--latent_dim 64]
"""
import argparse, os, sys

def main():
    parser = argparse.ArgumentParser(description='XOmiVAE multi-omics VAE')
    parser.add_argument('--dataset',    required=True)
    parser.add_argument('--version',    default='Top')
    parser.add_argument('--data_path',  required=True)
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--epochs',     type=int, default=50)
    parser.add_argument('--latent_dim', type=int, default=64)
    parser.add_argument('--batch_size', type=int, default=64)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir   = os.path.normpath(os.path.join(script_dir, '../../../../..'))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)

    import numpy as np
    import pandas as pd
    import tensorflow as tf
    import keras
    from generalHelperFunctions import loadData_mlomics

    data_root = args.data_path
    if os.path.basename(data_root) == 'Classification_datasets':
        data_root = os.path.dirname(data_root)
    class_num, labels, sample_ids, omics = loadData_mlomics(
        args.dataset, version=args.version, data_root=data_root)
    print(f'Classes: {class_num}, Samples: {len(sample_ids)}')

    X = np.concatenate([v for v in omics.values()], axis=1).astype('float32')
    n_features = X.shape[1]

    class VAE(keras.Model):
        def __init__(self, n_feat, latent_dim, **kw):
            super().__init__(**kw)
            self.encoder_layers = [
                keras.layers.Dense(512, activation='relu'),
                keras.layers.Dense(256, activation='relu'),
            ]
            self.z_mean_layer    = keras.layers.Dense(latent_dim)
            self.z_log_var_layer = keras.layers.Dense(latent_dim)
            self.decoder_layers  = [
                keras.layers.Dense(256, activation='relu'),
                keras.layers.Dense(512, activation='relu'),
                keras.layers.Dense(n_feat),
            ]
            self.n_feat = n_feat

        def encode(self, x):
            h = x
            for layer in self.encoder_layers:
                h = layer(h)
            return self.z_mean_layer(h), self.z_log_var_layer(h)

        def decode(self, z):
            h = z
            for layer in self.decoder_layers:
                h = layer(h)
            return h

        def call(self, x, training=False):
            z_mean, z_log_var = self.encode(x)
            eps = tf.random.normal(tf.shape(z_mean))
            z   = z_mean + tf.exp(0.5 * z_log_var) * eps
            recon = self.decode(z)
            recon_loss = tf.reduce_mean(tf.square(x - recon)) * self.n_feat
            kl_loss = -0.5 * tf.reduce_mean(
                1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
            self.add_loss(recon_loss + kl_loss)
            return recon

    vae = VAE(n_features, args.latent_dim)
    vae.compile(optimizer=keras.optimizers.Adam(1e-3))
    vae.fit(X, X, epochs=args.epochs, batch_size=args.batch_size, verbose=1)

    # Extract latent means
    z_out = np.zeros((len(X), args.latent_dim), dtype='float32')
    batch = args.batch_size
    for i in range(0, len(X), batch):
        xb = X[i:i+batch]
        z_out[i:i+batch], _ = vae.encode(xb)

    os.makedirs(args.output_dir, exist_ok=True)
    out_file = os.path.join(args.output_dir,
                            f'results_{args.dataset}_{args.version}_XOmiVAE.csv')
    result_df = pd.DataFrame({'sample': sample_ids, 'label': labels,
                               **{f'z{i}': z_out[:, i] for i in range(min(5, z_out.shape[1]))}})
    result_df.to_csv(out_file, index=False)
    print(f'Results saved to: {out_file}')
    print(f'Latent shape: {z_out.shape}')

if __name__ == '__main__':
    main()
