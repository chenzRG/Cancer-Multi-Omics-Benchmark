"""
This module contains functions that create different autoencoders
"""

import numpy as np
import tensorflow as tf
from keras.models import Model
from keras import backend as K
import tensorflow as _tf_backend
import keras as _keras_ops
# Keras 2 → Keras 3 backend shims (K.* functions removed from keras.backend in Keras 3)
# Use keras.ops where possible so shims work inside Lambda layers (symbolic KerasTensors).
if not hasattr(K, 'variable'):
    K.variable = lambda val, name=None, **kw: _tf_backend.Variable(val, name=name, dtype='float32')
if not hasattr(K, 'shape'):
    K.shape = _keras_ops.ops.shape
if not hasattr(K, 'exp'):
    K.exp = _keras_ops.ops.exp
if not hasattr(K, 'square'):
    K.square = _keras_ops.ops.square
if not hasattr(K, 'sum'):
    K.sum = _keras_ops.ops.sum
if not hasattr(K, 'mean'):
    K.mean = _keras_ops.ops.mean
if not hasattr(K, 'random_normal'):
    K.random_normal = lambda shape, mean=0.0, stddev=1.0, seed=None, **kw: _keras_ops.random.normal(shape=shape, mean=mean, stddev=stddev)
if not hasattr(K, 'get_value'):
    K.get_value = lambda x: float(x.numpy()) if hasattr(x, 'numpy') else float(x)
if not hasattr(K, 'set_value'):
    K.set_value = lambda x, v: x.assign(v)
from keras.callbacks import Callback
from keras import metrics, optimizers
from keras.layers import BatchNormalization  # keras.layers.normalization removed in Keras 3
from keras.layers import Input, Dense, Lambda, Activation


def make_variational_layer(
    input_layer,
    size,
    batch_normalize_intermediaries,
    relu_intermediaries,
    sampling_function,
    name="",
):
    # variational layer for hidden dim
    z_mean_component = Dense(
        size, kernel_initializer="glorot_uniform", name=f"{name}_mean" if name else None
    )
    z_mean_dense_linear = z_mean_component(input_layer)

    if batch_normalize_intermediaries:
        z_mean_dense_batchnorm = BatchNormalization(
            name=f"batchnorm_{name}_mean" if name else None
        )(z_mean_dense_linear)
    else:
        z_mean_dense_batchnorm = z_mean_dense_linear

    if relu_intermediaries:
        z_mean_encoded = Activation("relu", name=f"relu_{name}_mean" if name else None)(
            z_mean_dense_batchnorm
        )
    else:
        z_mean_encoded = z_mean_dense_batchnorm

    z_log_var_dense_linear = Dense(
        size,
        kernel_initializer="glorot_uniform",
        name=f"{name}_log_var" if name else None,
    )(input_layer)

    if batch_normalize_intermediaries:
        z_log_var_dense_batchnorm = BatchNormalization(
            name=f"batchnorm_{name}_log_var" if name else None
        )(z_log_var_dense_linear)
    else:
        z_log_var_dense_batchnorm = z_log_var_dense_linear

    if relu_intermediaries:
        z_log_var_encoded = Activation(
            "relu", name=f"relu_{name}_logvar" if name else None
        )(z_log_var_dense_batchnorm)
    else:
        z_log_var_encoded = z_log_var_dense_batchnorm

    # return the encoded and randomly sampled z vector (hidden layer 1)
    z = Lambda(
        sampling_function, output_shape=(size,), name=f"sample_{name}" if name else None
    )([z_mean_encoded, z_log_var_encoded])
    return z, z_mean_component


def stacked_vae(
    input_dim,
    hidden_dims=None,
    latent_dim=100,
    initial_beta_val=0,
    learning_rate=0.0005,
    epsilon_std=1.0,
    kappa=1.0,
    epochs=50,
    batch_size=50,
    batch_normalize_inputs=True,
    batch_normalize_intermediaries=True,
    batch_normalize_embedding=True,
    relu_intermediaries=True,
    relu_embedding=True,
    max_beta_val=1,
):
    """
    This is a deep, or stacked, vae.
    `hidden_dims` denotes the size of each successive hidden layer,
    until `latent_dim` which is the middle layer. The default `hidden_dims` is [300].
    """
    if hidden_dims is None:
        hidden_dims = [300]

    # Function for reparameterization trick to make model differentiable
    def sampling(args):
        z_mean, z_log_var = args
        # K.shape uses keras.ops.shape — works on both KerasTensors (symbolic) and eager tensors
        epsilon = K.random_normal(shape=K.shape(z_mean), mean=0.0, stddev=epsilon_std)
        z = z_mean + K.exp(z_log_var / 2) * epsilon
        return z

    # Init beta value
    beta = K.variable(initial_beta_val, name="beta")

    # Input place holder for RNAseq data with specific input size
    original_dim = input_dim

    # Input place holder for RNAseq data with specific input size
    rnaseq_input = Input(shape=(original_dim,), name="input")

    if batch_normalize_inputs:
        batchnorm_input = BatchNormalization(name="batchnorm_input")(rnaseq_input)
    else:
        batchnorm_input = rnaseq_input

    prev = batchnorm_input
    encoder_target = batchnorm_input
    if hidden_dims:
        for i, hidden_dim in enumerate(hidden_dims):
            z, z_mean_component = make_variational_layer(
                prev,
                hidden_dim,
                batch_normalize_intermediaries,
                relu_intermediaries,
                sampling,
                name=f"hidden_dim_{i}",
            )
            prev = z
            # the encoder part to have a path that doesn't do sampling or ReLU'ing
            encoder_target = z_mean_component(encoder_target)
    else:
        z = prev

    # variational layer for latent dim
    l_mean_component = Dense(
        latent_dim, kernel_initializer="glorot_uniform", name="latent_mean"
    )
    l_mean_dense_linear = l_mean_component(z)

    if batch_normalize_embedding:
        l_mean_dense_batchnorm = BatchNormalization(name="batchnorm_latent_mean")(
            l_mean_dense_linear
        )
    else:
        l_mean_dense_batchnorm = l_mean_dense_linear

    if relu_embedding:
        l_mean_encoded = Activation("relu", name="relu_latent_mean")(
            l_mean_dense_batchnorm
        )
    else:
        l_mean_encoded = l_mean_dense_batchnorm

    l_log_var_dense_linear = Dense(
        latent_dim, kernel_initializer="glorot_uniform", name="latent_log_var"
    )(z)

    if batch_normalize_embedding:
        l_log_var_dense_batchnorm = BatchNormalization(name="batchnorm_latent_log_var")(
            l_log_var_dense_linear
        )
    else:
        l_log_var_dense_batchnorm = l_log_var_dense_linear

    if relu_embedding:
        l_log_var_encoded = Activation("relu", name="relu_latent_log_var")(
            l_log_var_dense_batchnorm
        )
    else:
        l_log_var_encoded = l_log_var_dense_batchnorm

    l = Lambda(sampling, output_shape=(latent_dim,), name="sample_latent")(
        [l_mean_encoded, l_log_var_encoded]
    )

    # the encoder part's l to come from the path that only considers mean
    encoder_target = l_mean_component(encoder_target)
    if batch_normalize_embedding:
        encoder_target = BatchNormalization(name="batchnorm_encoder_target")(
            encoder_target
        )
    if relu_embedding:
        encoder_target = Activation("relu", name="relu_encoder_target")(encoder_target)

    # decoder latent->hidden
    prev = l
    if hidden_dims:
        for i, hidden_dim in reversed(list(enumerate(hidden_dims))):
            h = Dense(
                hidden_dim,
                kernel_initializer="glorot_uniform",
                activation="relu",
                name=f"decode_hidden_{i}",
            )(prev)
            prev = h
    else:
        h = l
    reconstruction = Dense(
        original_dim,
        kernel_initializer="glorot_uniform",
        activation="sigmoid",
        name="reconstruction",
    )(h)

    adam = optimizers.Adam(learning_rate=learning_rate)  # 'lr' renamed to 'learning_rate' in Keras 3
    _vae_base = Model(rnaseq_input, reconstruction)
    _z_params  = Model(rnaseq_input, [l_mean_encoded, l_log_var_encoded])

    # non-sampling encoder
    encoder = Model(rnaseq_input, encoder_target)

    # sampling encoder
    sampling_encoder = Model(rnaseq_input, l)

    # decoder model — built before wrapping vae, uses _vae_base.layers
    encoded_input = Input(shape=(latent_dim,))
    prev = encoded_input
    if hidden_dims:
        for i in reversed(range(len(hidden_dims) + 1)):
            prev = _vae_base.layers[-(i + 1)](prev)
    decoder = Model(encoded_input, prev)

    # Keras 3 functional models don't support add_loss() outside call().
    # Wrap with a custom train_step / test_step that recomputes the VAE loss.
    _orig_dim = original_dim
    _beta_ref  = beta

    class _VAETrainable(_keras_ops.Model):
        def __init__(self, **kw):
            super().__init__(**kw)
            # Register sub-models so Keras tracks their trainable_variables
            self._vae_m = _vae_base
            self._z_m   = _z_params

        def call(self, x, training=None):
            return self._vae_m(x, training=training)

        def _compute_loss(self, x):
            import tensorflow as _tf_vae
            recon = self._vae_m(x, training=True)
            _zm, _zlv = self._z_m(x, training=True)
            _recon = float(_orig_dim) * metrics.binary_crossentropy(x, recon)
            _kl = -0.5 * K.sum(
                1 + _zlv - K.square(_zm) - K.exp(_zlv), axis=-1)
            return K.mean(_recon + _beta_ref * _kl)

        def train_step(self, data):
            import tensorflow as _tf_vae
            x = data[0] if isinstance(data, (list, tuple)) else data
            with _tf_vae.GradientTape() as tape:
                loss = self._compute_loss(x)
            grads = tape.gradient(loss, self.trainable_variables)
            self.optimizer.apply_gradients(zip(grads, self.trainable_variables))
            return {'loss': loss}

        def test_step(self, data):
            x = data[0] if isinstance(data, (list, tuple)) else data
            loss = self._compute_loss(x)
            return {'loss': loss}

    vae = _VAETrainable()
    vae.compile(optimizer=adam)

    return vae, encoder, sampling_encoder, decoder, beta


class LadderCallback(Callback):
    """
    This class implements ladder autoecoders:
        https://arxiv.org/abs/1602.02282

    A callback on each epoch end, increments beta by kappa
    """

    def __init__(self, beta, kappa, max_val=1):
        self.beta = beta
        self.kappa = kappa
        self.max_val = max_val

    def on_epoch_end(self, *args, **kwargs):
        if K.get_value(self.beta) <= self.max_val:
            K.set_value(self.beta, K.get_value(self.beta) + self.kappa)


def train_model(vae, x_train, epochs, batch_size, x_val, beta, kappa, max_beta_val, verbose=0):
    hist = vae.fit(
        np.array(x_train),
        shuffle=True,
        epochs=epochs,
        verbose=verbose,
        batch_size=batch_size,
        validation_data=(np.array(x_val), None),
        callbacks=[LadderCallback(beta, kappa, max_beta_val)],
    )
    return hist


def deep_vae(
    input_dim,
    hidden_dims=None,
    latent_dim=100,
    initial_beta_val=0,
    learning_rate=0.0005,
    epsilon_std=1.0,
    kappa=1.0,
    epochs=50,
    batch_size=50,
    batch_normalize_inputs=True,
    batch_normalize_embedding=False,
    relu_embedding=False,
    max_beta_val=1,
):
    """
    This is a deep, not stacked, vae.
    `hidden_dims` denotes the size of each successive hidden layer,
    until `latend_dim` which is the middle layer, which will be mean/variance of a normal distribution.
    The default `hidden_dims` is [300].
    """
    if hidden_dims is None:
        hidden_dims = [300]
    # Function for reparameterization trick to make model differentiable
    def sampling(args):
        z_mean, z_log_var = args
        epsilon = K.random_normal(shape=K.shape(z_mean), mean=0.0, stddev=epsilon_std)
        z = z_mean + K.exp(z_log_var / 2) * epsilon
        return z

    # Init beta value
    beta = K.variable(initial_beta_val)

    # Input place holder for RNAseq data with specific input size
    original_dim = input_dim

    # Input place holder for RNAseq data with specific input size
    rnaseq_input = Input(shape=(original_dim,), name="input")

    if batch_normalize_inputs:
        batchnorm_input = BatchNormalization(name="batchnorm_input")(rnaseq_input)
    else:
        batchnorm_input = rnaseq_input

    prev = batchnorm_input
    if hidden_dims:
        for i, hidden_dim in enumerate(hidden_dims):
            z = Dense(hidden_dim, activation="relu", name=f"hidden_{i}")(prev)
            prev = z
    else:
        z = prev

    # variational layer for latent dim
    l_mean_component = Dense(
        latent_dim, kernel_initializer="glorot_uniform", name="latent_mean"
    )
    l_mean_dense_linear = l_mean_component(z)

    if batch_normalize_embedding:
        l_mean_dense_batchnorm = BatchNormalization(name="batchnorm_latent_mean")(
            l_mean_dense_linear
        )
    else:
        l_mean_dense_batchnorm = l_mean_dense_linear

    if relu_embedding:
        l_mean_encoded = Activation("relu", name="relu_latent_mean")(
            l_mean_dense_batchnorm
        )
    else:
        l_mean_encoded = l_mean_dense_batchnorm

    l_log_var_dense_linear = Dense(
        latent_dim, kernel_initializer="glorot_uniform", name="latent_log_var"
    )(z)

    if batch_normalize_embedding:
        l_log_var_dense_batchnorm = BatchNormalization(name="batchnorm_latent_log_var")(
            l_log_var_dense_linear
        )
    else:
        l_log_var_dense_batchnorm = l_log_var_dense_linear

    if relu_embedding:
        l_log_var_encoded = Activation("relu", name="relu_latent_log_var")(
            l_log_var_dense_batchnorm
        )
    else:
        l_log_var_encoded = l_log_var_dense_batchnorm

    l = Lambda(sampling, output_shape=(latent_dim,), name="sample_latent")(
        [l_mean_encoded, l_log_var_encoded]
    )

    # the encoder part's l to come from the path that only considers mean
    encoder_target = l_mean_component(z)
    if batch_normalize_embedding:
        encoder_target = BatchNormalization(name="batchnorm_encoder_target")(
            encoder_target
        )
    if relu_embedding:
        encoder_target = Activation("relu", name="relu_encoder_target")(encoder_target)

    # decoder latent->hidden
    if hidden_dims:
        prev = l
        for i, hidden_dim in reversed(list(enumerate(hidden_dims))):
            h = Dense(
                hidden_dim,
                kernel_initializer="glorot_uniform",
                activation="relu",
                name=f"decode_hidden_{i}",
            )(prev)
            prev = h
    else:
        h = l

    reconstruction = Dense(
        original_dim,
        kernel_initializer="glorot_uniform",
        activation="sigmoid",
        name="reconstruction",
    )(h)

    adam = optimizers.Adam(learning_rate=learning_rate)  # 'lr' renamed to 'learning_rate' in Keras 3
    _vae_base = Model(rnaseq_input, reconstruction)
    _z_params  = Model(rnaseq_input, [l_mean_encoded, l_log_var_encoded])

    # non-sampling encoder
    encoder = Model(rnaseq_input, encoder_target)

    # sampling encoder
    sampling_encoder = Model(rnaseq_input, l)

    # decoder model — built before wrapping vae, uses _vae_base.layers
    encoded_input = Input(shape=(latent_dim,))
    prev = encoded_input
    if hidden_dims:
        for i in reversed(range(len(hidden_dims) + 1)):
            prev = _vae_base.layers[-(i + 1)](prev)
    decoder = Model(encoded_input, prev)

    # Keras 3 functional models don't support add_loss() outside call().
    _orig_dim = original_dim
    _beta_ref  = beta

    class _VAETrainable(_keras_ops.Model):
        def call(self, x, training=None):
            return _vae_base(x, training=training)

        def _compute_loss(self, x):
            recon = _vae_base(x, training=True)
            _zm, _zlv = _z_params(x, training=True)
            _recon = float(_orig_dim) * metrics.binary_crossentropy(x, recon)
            _kl = -0.5 * K.sum(
                1 + _zlv - K.square(_zm) - K.exp(_zlv), axis=-1)
            return K.mean(_recon + _beta_ref * _kl)

        def train_step(self, data):
            import tensorflow as _tf_vae
            x = data[0] if isinstance(data, (list, tuple)) else data
            with _tf_vae.GradientTape() as tape:
                loss = self._compute_loss(x)
            grads = tape.gradient(loss, self.trainable_variables)
            self.optimizer.apply_gradients(zip(grads, self.trainable_variables))
            return {'loss': loss}

        def test_step(self, data):
            x = data[0] if isinstance(data, (list, tuple)) else data
            loss = self._compute_loss(x)
            return {'loss': loss}

    vae = _VAETrainable()
    vae.compile(optimizer=adam)

    return vae, encoder, sampling_encoder, decoder, beta
