import argparse
import sys
import numpy as np
import random
import time
import os

import tensorflow.compat.v1 as tf
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
from keras.layers import Input, Dense, Lambda, Layer, Add, BatchNormalization, Dropout, Activation, Conv2D, \
    MaxPooling2D, Activation, LeakyReLU, concatenate
from keras.models import Model, Sequential
from keras.losses import mse, binary_crossentropy
from keras.optimizers import Adam
from sklearn.ensemble import RandomForestClassifier
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from itertools import combinations
import bisect

random.seed(1)
np.random.seed(1)
tf.set_random_seed(1)
session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
K.set_session(sess)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"


class ConsensusCluster:
    def __init__(self, cluster, L, K, H, resample_proportion=0.8):
        self.cluster_ = cluster
        self.resample_proportion_ = resample_proportion
        self.L_ = L
        self.K_ = K
        self.H_ = H
        self.Mk = None
        self.Ak = None
        self.deltaK = None
        self.bestK = None

    def _internal_resample(self, data, proportion):
        ids = np.random.choice(
            range(data.shape[0]), size=int(data.shape[0] * proportion), replace=False)
        return ids, data[ids, :]

    def fit(self, data):
        Mk = np.zeros((self.K_ - self.L_, data.shape[0], data.shape[0]))
        Is = np.zeros((data.shape[0],) * 2)
        for k in range(self.L_, self.K_):
            i_ = k - self.L_
            for h in range(self.H_):
                ids, dt = self._internal_resample(data, self.resample_proportion_)
                Mh = self.cluster_(n_clusters=k).fit_predict(dt)
                ids_sorted = np.argsort(Mh)
                sorted_ = Mh[ids_sorted]
                for i in range(k):
                    ia = bisect.bisect_left(sorted_, i)
                    ib = bisect.bisect_right(sorted_, i)
                    is_ = ids_sorted[ia:ib]
                    ids_ = np.array(list(combinations(is_, 2))).T
                    if ids_.size != 0:
                        Mk[i_, ids_[0], ids_[1]] += 1
                ids_2 = np.array(list(combinations(ids, 2))).T
                Is[ids_2[0], ids_2[1]] += 1
            Mk[i_] /= Is + 1e-8
            Mk[i_] += Mk[i_].T
            Mk[i_, range(data.shape[0]), range(
                data.shape[0])] = 1
            Is.fill(0)
        self.Mk = Mk
        self.Ak = np.zeros(self.K_ - self.L_)
        for i, m in enumerate(Mk):
            hist, bins = np.histogram(m.ravel(), density=True)
            self.Ak[i] = np.sum(h * (b - a)
                                for b, a, h in zip(bins[1:], bins[:-1], np.cumsum(hist)))
        self.deltaK = np.array([(Ab - Aa) / Aa if i > 2 else Aa
                                for Ab, Aa, i in zip(self.Ak[1:], self.Ak[:-1], range(self.L_, self.K_ - 1))])
        self.bestK = np.argmax(self.deltaK) + \
                     self.L_ if self.deltaK.size > 0 else self.L_

    def predict(self):
        return self.cluster_(n_clusters=self.bestK).fit_predict(
            1 - self.Mk[self.bestK - self.L_])

    def predict_data(self, data):
        return self.cluster_(n_clusters=self.bestK).fit_predict(
            data)


class GeLU(Activation):
    def __init__(self, activation, **kwargs):
        super(GeLU, self).__init__(activation, **kwargs)
        self.__name__ = 'gelu'


def gelu(x):
    return 0.5 * x * (1 + tf.tanh(tf.sqrt(2 / np.pi) * (x + 0.044715 * tf.pow(x, 3))))


get_custom_objects().update({'gelu': GeLU(gelu)})


class AE():
    def __init__(self, X_shape, n_components, epochs=100):
        self.epochs = epochs
        sample_size = X_shape[0]
        self.batch_size = 16
        sample_size = X_shape[0]
        self.epochs = 30
        self.n_components = n_components
        self.shape = X_shape[1]

    def train(self, X):
        encoding_dim = self.n_components
        original_dim = X.shape[1]
        input = Input(shape=(original_dim,))
        encoded = Dense(encoding_dim)(input)
        encoded = BatchNormalization()(encoded)
        encoded = Activation('relu')(encoded)
        z = Dense(encoding_dim, activation='relu')(encoded)
        decoded = Dense(encoding_dim, activation='relu')(z)
        output = Dense(original_dim, activation='sigmoid')(decoded)
        ae = Model(input, output)
        encoder = Model(input, z)
        ae_loss = mse(input, output)
        ae.add_loss(ae_loss)
        ae.compile(optimizer=Adam())
        print(len(ae.layers))
        print(ae.count_params())
        ae.fit(X, epochs=self.epochs, batch_size=self.batch_size, verbose=2)
        return encoder.predict(X)


class VAE():
    def __init__(self, X_shape, n_components, epochs=100):
        self.epochs = epochs
        self.batch_size = 16
        sample_size = X_shape[0]
        self.epochs = 30
        self.n_components = n_components
        self.shape = X_shape[1]

    def train(self, X):
        def sampling(args):
            z_mean, z_log_var = args
            batch = K.shape(z_mean)[0]
            dim = K.int_shape(z_mean)[1]
            epsilon = K.random_normal(shape=(batch, dim), seed=0)
            return z_mean + K.exp(0.5 * z_log_var) * epsilon

        encoding_dim = self.n_components
        original_dim = X.shape[1]
        input = Input(shape=(original_dim,))
        encoded = Dense(encoding_dim)(input)
        encoded = BatchNormalization()(encoded)
        encoded = Activation('relu')(encoded)
        z_mean = Dense(encoding_dim)(encoded)
        z_log_var = Dense(encoding_dim)(encoded)
        z = Lambda(sampling, output_shape=(encoding_dim,), name='z')([z_mean, z_log_var])
        decoded = Dense(encoding_dim, activation='relu')(z)
        output = Dense(original_dim, activation='sigmoid')(decoded)
        vae = Model(input, output)
        encoder = Model(input, z)
        reconstruction_loss = mse(input, output)
        reconstruction_loss *= original_dim
        kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var)
        kl_loss = K.sum(kl_loss, axis=-1)
        kl_loss *= -0.5
        vae_loss = K.mean(reconstruction_loss + kl_loss)
        vae.add_loss(vae_loss)
        vae.compile(optimizer=Adam())
        print(len(vae.layers))
        print(vae.count_params())
        vae.fit(X, epochs=self.epochs, batch_size=self.batch_size, verbose=2)
        return encoder.predict(X)

class SubtypeGAN():
    def __init__(self, datasets, n_latent_dim, weight=0.001, model_path='SubtypeGAN.h5', epochs=100, batch_size=64):
        self.latent_dim = n_latent_dim
        optimizer = Adam()
        self.n = len(datasets)
        self.epochs = epochs
        self.batch_size = batch_size
        sample_size = 0
        if self.n > 1:
            sample_size = datasets[0].shape[0]
        print(sample_size)
        if sample_size > 300:
            self.epochs = 11
        else:
            self.epochs = 10
        self.epochs = 30 * batch_size
        self.shape = []
        self.weight = [0.3, 0.1, 0.1, 0.5]
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
            loss_weights.append((1 - self.disc_w) * self.weight[i])
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
            dims.append(int(encoding_dim * self.weight[i]))
        for i in range(self.n):
            denses.append(Dense(dims[i])(X[i]))
        if self.n > 1:
            merged_dense = concatenate(denses, axis=-1)
        else:
            merged_dense = denses[0]
        model = BatchNormalization()(merged_dense)
        model = Activation('gelu')(model)
        model = Dense(encoding_dim)(model)
        z_mean = Dense(encoding_dim)(model)
        z_log_var = Dense(encoding_dim)(model)
        z = Lambda(sampling, output_shape=(encoding_dim,), name='z')([z_mean, z_log_var])
        return Model(X, [z_mean, z_log_var, z])

    def build_decoder_disc(self):
        denses = []
        X = Input(shape=(self.latent_dim,))
        model = Dense(self.latent_dim)(X)
        model = BatchNormalization()(model)
        model = Activation('gelu')(model)
        for i in range(self.n):
            denses.append(Dense(self.shape[i])(model))
        dec = Dense(1, activation='sigmoid')(model)
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
        log_file = "./run.log"
        fp = open(log_file, 'w')
        if bTrain:
            # GAN
            valid = np.ones((self.batch_size, 1))
            fake = np.zeros((self.batch_size, 1))
            for epoch in range(self.epochs):
                #  Train Discriminator
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
                #  Train Encoder_GAN
                g_loss = self.gan.train_on_batch(data, outs)
            fp.close()
            self.encoder.save(model_path)
        else:
            self.encoder = load_model(model_path)
        mat = self.encoder.predict(X_train)[0]
        return mat


class SubtypeGAN_API(object):
    def __init__(self, model_path='./model/', epochs=200, weight=0.001):
        self.model_path = model_path
        self.score_path = './score/'
        self.epochs = epochs
        self.batch_size = 16
        self.weight = weight

    # feature extract
    def feature_gan(self, datasets, index=None, n_components=100, b_decomposition=True, weight=0.001):
        if b_decomposition:
            X = self.encoder_gan(datasets, n_components)
            fea = pd.DataFrame(data=X, index=index, columns=map(lambda x: 'v' + str(x), range(X.shape[1])))
        else:
            fea = np.concatenate(datasets)
        print("feature extract finished!")
        return fea

    def feature_vae(self, df_ori, n_components=100, b_decomposition=True):
        if b_decomposition:
            X = self.encoder_vae(df_ori, n_components)
            print(X)
            fea = pd.DataFrame(data=X, index=df_ori.index,
                               columns=map(lambda x: 'v' + str(x), range(X.shape[1])))
        else:
            fea = df_ori.copy()
        print("feature extract finished!")
        return fea

    def feature_ae(self, df_ori, n_components=100, b_decomposition=True):
        if b_decomposition:
            X = self.encoder_ae(df_ori, n_components)
            print(X)
            fea = pd.DataFrame(data=X, index=df_ori.index,
                               columns=map(lambda x: 'v' + str(x), range(X.shape[1])))
        else:
            fea = df_ori.copy()
        print("feature extract finished!")
        return fea

    def impute(self, X):
        X.fillna(X.mean())
        return X

    def encoder_gan(self, ldata, n_components=100):
        egan = SubtypeGAN(ldata, n_components, self.weight, self.model_path, self.epochs, self.batch_size)
        return egan.train(ldata)

    def encoder_vae(self, df, n_components=100):
        vae = VAE(df.shape, n_components, self.epochs)
        return vae.train(df)

    def encoder_ae(self, df, n_components=100):
        ae = AE(df.shape, n_components, self.epochs)
        return ae.train(df)

    def tsne(self, X):
        model = TSNE(n_components=2)
        return model.fit_transform(X)

    def pca(self, X):
        fea_model = PCA(n_components=200)
        return fea_model.fit_transform(X)

    def gmm(self, n_clusters=28):
        model = mixture.GaussianMixture(n_components=n_clusters, covariance_type='diag')
        return model

    def kmeans(self, n_clusters=28):
        model = KMeans(n_clusters=n_clusters, random_state=0)
        return model

    def spectral(self, n_clusters=28):
        model = SpectralClustering(n_clusters=n_clusters, random_state=0)
        return model

    def hierarchical(self, n_clusters=28):
        model = AgglomerativeClustering(n_clusters=n_clusters)
        return model


def main(argv=sys.argv):
    parser = argparse.ArgumentParser(description='SubtypeGAN v1.0')
    parser.add_argument("-e", dest='epochs', type=int, default=200, help="Number of iterations")
    parser.add_argument("-m", dest='run_mode', default="feature", help="run_mode: feature, cluster")
    parser.add_argument("-n", dest='cluster_num', type=int, default=-1, help="cluster number")
    parser.add_argument("-w", dest='disc_weight', type=float, default=1e-4, help="weight")
    parser.add_argument("-o", dest='output_path', default="./score/", help="file output")
    parser.add_argument("-t", dest='type', default="ACC_Top", help="data type")
    args = parser.parse_args()
    model_path = './model/' + args.type + '.h5'
    SubtypeGAN = SubtypeGAN_API(model_path, epochs=args.epochs, weight=args.disc_weight)
    omics_types = ['CNV', 'Methy', 'mRNA', 'miRNA']
    
    if args.run_mode == 'SubtypeGAN':
        print("Start SubTypeGAN...")
        data_type = args.type
        if args.cluster_num == -1:
            print("Please set the number of clusters!")
        fea_tmp_file = './fea/' + data_type + '.fea'
        tmp_dir = './fea/' + data_type + '/'
        model_dir ='./model'
        os.makedirs(tmp_dir, exist_ok=True) 
        os.makedirs(model_dir, exist_ok=True) 
        ldata = []
        l = []

        dir_path = '../../../....//Main_Dataset/Clustering_datasets/' + data_type.replace('_', '/')
        files = os.listdir(dir_path)
        for omics in omics_types:
            for file_name in files:
                if omics in file_name and file_name.endswith('.csv') and not file_name.startswith('.'):
                    file_path = os.path.join(dir_path, file_name)
                    if os.path.isfile(file_path):
                        print("Processing ", file_path)
                        df_new = pd.read_csv(file_path, sep=',', header=0, index_col=0)
                        l = list(df_new)
                        df_new = df_new.T
                        ldata.append(df_new.values.astype(float))
                    else:
                        print(f"File not found: {file_path}")

        start_time = time.time()
        vec = SubtypeGAN.feature_gan(ldata, index=l, n_components=100, weight=args.disc_weight)
        df = pd.DataFrame(data=[time.time() - start_time])
        vec.to_csv(fea_tmp_file, header=True, index=True, sep='\t')
        out_file = './results/' + data_type + '.SubtypeGAN.time'
        df.to_csv(out_file, header=True, index=False, sep=',')

        if isfile(fea_tmp_file):
            X = pd.read_csv(fea_tmp_file, header=0, index_col=0, sep='\t')
            X['SubtypeGAN'] = SubtypeGAN.gmm(args.cluster_num).fit_predict(X.values) + 1
            X = X.loc[:, ['SubtypeGAN']]
            out_file = './results/' + data_type + '.SubtypeGAN'
            X.to_csv(out_file, header=True, index=True, sep='\t')
        else:
            print(fea_tmp_file, ' file does not exist!')


    elif args.run_mode == 'cc':
        min_cluster = 2
        max_cluster = 8
        iteration = 10

        data_type = args.type
        fea_tmp_file = './fea/' + data_type + '.fea'
        fs = []
        cc_file = './results/cluster_num.cc'
        fp = open(cc_file, 'a')
        if isfile(fea_tmp_file):
            X = pd.read_csv(fea_tmp_file, header=0, index_col=0, sep='\t')
            cc = ConsensusCluster(SubtypeGAN.gmm, min_cluster, max_cluster, iteration)
            cc.fit(X.values)
            X['cc'] = SubtypeGAN.gmm(cc.bestK).fit_predict(X.values) + 1
            X = X.loc[:, ['cc']]
            out_file = './results/' + data_type + '.cc'
            X.to_csv(out_file, header=True, index=True, sep='\t')
            fp.write("%s, k=%d\n" % (data_type, cc.bestK))
        else:
            print('file does not exist!')
        fp.close()
        print(f"clustering results are saved in {out_file}")


if __name__ == "__main__":
    main()
