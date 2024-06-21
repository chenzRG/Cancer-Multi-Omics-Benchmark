import argparse
import sys
import numpy as np
import random
import time
import os

import tensorflow as tf
from subprocess import check_output
import h5py
import re
import math
import pandas as pd
from os.path import splitext, basename, isfile
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn import preprocessing
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn import mixture
from keras import backend as K
from keras.layers import Input, Dense, Lambda, Layer, Add, BatchNormalization, Dropout, Activation, merge, Conv2D, \
    MaxPooling2D, Activation, LeakyReLU, concatenate
from keras.models import Model, Sequential
from keras.losses import mse, binary_crossentropy
from keras.optimizers import Adam
from sklearn.ensemble import RandomForestClassifier
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from keras import layers, losses, Model

random.seed(1)
np.random.seed(1)
tf.set_random_seed(1)
session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
K.set_session(sess)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"


class SubtypeGAN():
    def __init__(self, datasets, n_latent_dim, model_path='SubtypeGAN.h5', epochs=100, batch_size=64):
        self.latent_dim = n_latent_dim
        optimizer = Adam()
        self.n = len(datasets)
        self.epochs = 100
        self.batch_size = batch_size
        sample_size = 0
        if self.n > 1:
            sample_size = datasets[0].shape[0]
        print(sample_size)
        self.shape = []
        self.disc_w = 1e-4
        self.model_path = model_path
        input = []
        loss = []
        loss_weights = []
        output = []
        for i in range(self.n):
            self.shape.append(datasets[i].shape[1])
            loss.append('mse')
        loss.append('binary_crossentropy')
        self.decoder, self.disc = self.build_decoder_disc()
        self.disc.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['accuracy'])
        self.encoder = self.build_encoder()
        for i in range(self.n):
            input.append(Input(shape=(self.shape[i],)))
            loss_weights.append(1 / self.n)
        loss_weights.append(self.disc_w)
        z_mean, z_log_var, z = self.encoder(input)
        output = self.decoder(z)
        self.gan = Model(input, output)
        self.gan.compile(loss=loss, loss_weights=loss_weights, optimizer=optimizer)
        print(self.gan.summary())
        return

    def build_encoder(self):
        def sampling(args):
            z_mean, z_log_var = args
            return z_mean + K.exp(0.5 * z_log_var) * K.random_normal(K.shape(z_mean), seed=0)

        encoding_dim = self.latent_dim
        X = []
        dims = []
        denses = []
        for i in range(self.n):
            X.append(Input(shape=(self.shape[i],)))
            dims.append(int(encoding_dim / self.n))
        for i in range(self.n):
            denses.append(Dense(dims[i], kernel_initializer="glorot_normal")(X[i]))
        if self.n > 1:
            merged_dense = concatenate(denses, axis=-1)
        else:
            merged_dense = denses[0]
        model = BatchNormalization()(merged_dense)
        model = Activation('relu')(model)
        model = Dense(encoding_dim, kernel_initializer="glorot_normal")(model)
        z_mean = Dense(encoding_dim, kernel_initializer="glorot_normal")(model)
        z_log_var = Dense(encoding_dim, kernel_initializer="glorot_normal")(model)
        z = Lambda(sampling, output_shape=(encoding_dim,), name='z')([z_mean, z_log_var])
        return Model(X, [z_mean, z_log_var, z])

    def build_decoder_disc(self):
        denses = []
        X = Input(shape=(self.latent_dim,))
        model = Dense(self.latent_dim, kernel_initializer="glorot_normal")(X)
        model = BatchNormalization()(model)
        model = Activation('relu')(model)
        for i in range(self.n):
            denses.append(Dense(self.shape[i], kernel_initializer="glorot_normal")(model))
        dec = Dense(1, activation='sigmoid', kernel_initializer="glorot_normal")(model)
        denses.append(dec)
        m_decoder = Model(X, denses)
        m_disc = Model(X, dec)
        return m_decoder, m_disc

    def build_disc(self):
        X = Input(shape=(self.latent_dim,))
        dec = Dense(1, activation='sigmoid', kernel_initializer="glorot_normal")(X)
        output = Model(X, dec)
        return output

    def train(self, X_train, bTrain=True):
        model_path = self.model_path
        log_file = "./gan.log"
        fp = open(log_file, 'w')
        if bTrain:
            # GAN
            valid = np.ones((self.batch_size, 1))
            fake = np.zeros((self.batch_size, 1))
            for epoch in range(self.epochs * self.batch_size):
                #  Train Discriminator
                iter = epoch / self.batch_size
                data = []
                idx = np.random.randint(0, X_train[0].shape[0], self.batch_size)
                for i in range(self.n):
                    data.append(X_train[i][idx])
                latent_fake = self.encoder.predict(data)[2]
                latent_real = np.random.normal(size=(self.batch_size, self.latent_dim))
                d_loss_real = self.disc.train_on_batch(latent_real, valid)
                d_loss_fake = self.disc.train_on_batch(latent_fake, fake)
                d_loss = 0.5 * np.add(d_loss_real, d_loss_fake)
                outs = data + [valid]
                g_loss = self.gan.train_on_batch(data, outs)
                if epoch % self.batch_size == 0:
                    mse = np.mean(g_loss[1:5])
                    print("%d [D loss: %f, acc: %.2f] [mse: %f]" % (
                        iter + 1, d_loss[0], 100 * d_loss[1], mse))
                    fp.write("%f\n" % (mse))
            fp.close()


class MVAE():
    def __init__(self, datasets, n_latent_dim, epochs=100, batch_size=64):
        self.epochs = epochs
        self.batch_size = 32
        self.n = len(datasets)
        self.batch_size = batch_size
        self.shape = []
        self.epochs = 100
        self.n_components = n_latent_dim
        self.original_dim = 0
        for i in range(self.n):
            self.shape.append(datasets[i].shape[1])
            self.original_dim += datasets[i].shape[1]

    def train(self, X):
        def sampling(args):
            z_mean, z_log_var = args
            batch = K.shape(z_mean)[0]
            dim = K.int_shape(z_mean)[1]
            epsilon = K.random_normal(shape=(batch, dim), seed=0)
            return z_mean + K.exp(0.5 * z_log_var) * epsilon

        input = []
        dims = []
        denses = []
        encoding_dim = self.n_components
        output = []
        for i in range(self.n):
            input.append(Input(shape=(self.shape[i],)))
            dims.append(int(encoding_dim * 1 / self.n))
        for i in range(self.n):
            denses.append(Dense(dims[i])(input[i]))
        if self.n > 1:
            merged_dense = concatenate(denses, axis=-1)
        else:
            merged_dense = denses[0]
        encoded = Dense(encoding_dim)(merged_dense)
        encoded = BatchNormalization()(encoded)
        encoded = Activation('relu')(encoded)
        z_mean = Dense(encoding_dim)(encoded)
        z_log_var = Dense(encoding_dim)(encoded)
        z = Lambda(sampling, output_shape=(encoding_dim,), name='z')([z_mean, z_log_var])
        model = Dense(self.n_components)(z)
        model = BatchNormalization()(model)
        model = Activation('relu')(model)
        for i in range(self.n):
            output.append(Dense(self.shape[i])(model))
        vae = Model(input, output)
        encoder = Model(input, z)
        kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var)
        kl_loss = K.sum(kl_loss, axis=-1)
        kl_loss *= -0.5 / np.sum(dims)
        k_mse_loss = 0
        for i in range(self.n):
            k_mse_loss += mse(input[i], output[i]) / self.n
        vae.add_loss(k_mse_loss)
        vae.add_metric(k_mse_loss, name='mse_loss')
        vae.add_loss(kl_loss)
        vae.add_metric(kl_loss, name='kl_loss')
        vae.compile(optimizer=Adam())
        print(vae.summary())
        h = vae.fit(X, epochs=self.epochs, verbose=2)
        log_file = "./mvae.log"
        fp = open(log_file, 'w')
        for hi in h.history['mse_loss']:
            fp.write("%f\n" % (hi))
        fp.close()
        return


class SubtypeGAN_API(object):
    def __init__(self, model_path='./model/', epochs=200):
        self.model_path = model_path
        self.score_path = './score/'
        self.epochs = epochs
        self.batch_size = 16

    # feature extract
    def feature_gan(self, datasets, n_components=100):
        self.encoder_gan(datasets, n_components)
        return True

    def feature_mvae(self, datasets, n_components=100):
        self.encoder_mvae(datasets, n_components)
        return True

    def impute(self, X):
        X.fillna(X.mean())
        return X

    def encoder_gan(self, ldata, n_components=100):
        egan = SubtypeGAN(ldata, n_components)
        return egan.train(ldata)

    def encoder_mvae(self, ldata, n_components=100):
        emvae = MVAE(ldata, n_components)
        return emvae.train(ldata)


def main(argv=sys.argv):
    parser = argparse.ArgumentParser(description='SubtypeGAN v1.0')
    parser.add_argument("-i", dest='file_input', default="./input/input.list",
                        help="file input")
    parser.add_argument("-e", dest='epochs', type=int, default=200, help="Number of iterations")
    parser.add_argument("-m", dest='run_mode', default="feature", help="run_mode: feature, cluster")
    parser.add_argument("-n", dest='cluster_num', type=int, default=-1, help="cluster number")
    parser.add_argument("-o", dest='output_path', default="./score/", help="file output")
    parser.add_argument("-s", dest='surv_path',
                        default="./data/TCGA/clinical_PANCAN_patient_with_followup.tsv",
                        help="surv input")
    parser.add_argument("-t", dest='type', default="BRCA", help="cancer type: BRCA, GBM")
    args = parser.parse_args()
    model_path = './model/' + args.type + '.h5'
    SubtypeGAN = SubtypeGAN_API(model_path, epochs=args.epochs)
    cancer_dict = {'BRCA': 5, 'BLCA': 5, 'KIRC': 4,
                   'GBM': 3, 'LUAD': 3, 'PAAD': 2,
                   'SKCM': 4, 'STAD': 3, 'UCEC': 4, 'UVM': 4}

    if args.run_mode == 'SubtypeGAN':
        cancer_type = args.type
        if cancer_type not in cancer_dict and args.cluster_num == -1:
            print("Please set the number of clusters!")
        elif args.cluster_num == -1:
            args.cluster_num = cancer_dict[cancer_type]
        fea_tmp_file = './fea/' + cancer_type + '.fea'
        tmp_dir = './fea/' + cancer_type + '/'
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        ldata = []
        l = []
        for line in open(args.file_input, 'rt'):
            base_file = splitext(basename(line.rstrip()))[0]
            fea_save_file = tmp_dir + base_file + '.fea'
            if isfile(fea_save_file):
                df_new = pd.read_csv(fea_save_file, sep=',', header=0, index_col=0)
                l = list(df_new)
            df_new = df_new.T
            ldata.append(df_new.values.astype(float))
        start_time_begin = time.time()
        SubtypeGAN.feature_gan(ldata, n_components=100)
        print(time.time() - start_time_begin)

    elif args.run_mode == 'mvae':
        cancer_type = args.type
        if cancer_type not in cancer_dict and args.cluster_num == -1:
            print("Please set the number of clusters!")
        elif args.cluster_num == -1:
            args.cluster_num = cancer_dict[cancer_type]
        fea_tmp_file = './fea/' + cancer_type + '.fea'
        tmp_dir = './fea/' + cancer_type + '/'
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        ldata = []
        l = []
        for line in open(args.file_input, 'rt'):
            base_file = splitext(basename(line.rstrip()))[0]
            fea_save_file = tmp_dir + base_file + '.fea'
            if isfile(fea_save_file):
                df_new = pd.read_csv(fea_save_file, sep=',', header=0, index_col=0)
                l = list(df_new)
            else:
                print(fea_save_file)
                return
            df_new = df_new.T
            ldata.append(df_new.values.astype(float))
        start_time_begin = time.time()
        SubtypeGAN.feature_mvae(ldata, n_components=100)
        print(time.time() - start_time_begin)


if __name__ == "__main__":
    main()
