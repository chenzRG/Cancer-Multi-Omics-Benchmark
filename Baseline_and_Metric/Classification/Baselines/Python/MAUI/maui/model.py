import warnings
import maui.utils
import numpy as np
import pandas as pd
from functools import partial
from keras import backend as K
from scipy import spatial, cluster
from sklearn.cluster import KMeans
from sklearn.base import BaseEstimator
from .maui_warnings import MauiWarning
from .autoencoders_architectures import stacked_vae, deep_vae, train_model


class Maui(BaseEstimator):
    """Maui (Multi-omics Autoencoder Integration) model.

    Trains a variational autoencoder to find latent factors in multi-modal data.

    Parameters
    ----------
    n_hidden: array (default [1500])
        The sizes of the hidden layers of the autoencoder architecture.
        Each element of the array specifies the number of nodes in successive
        layers of the autoencoder

    n_latent: int (default 80)
        The size of the latent layer (number of latent features)

    batch_size: int (default 100)
        The size of the mini-batches used for training the network

    epochs: int (default 400)
        The number of epoches to use for training the network

    architecture:
        One of 'stacked' or 'deep'. If 'stacked', will use a stacked VAE model, where
        the intermediate layers are also variational. If 'deep', will train a deep VAE
        where the intermediate layers are regular (ReLU) units, and only the middle
        (latent) layer is variational.

    verbose:
         Integer. 0, 1, or 2. Verbosity mode. 0 = silent, 1 = progress bar, 2 = one line per epoch.
         default: 0
    """

    def __init__(
        self,
        n_hidden=None,
        n_latent=80,
        batch_size=100,
        epochs=400,
        architecture="stacked",
        initial_beta_val=0,
        kappa=1.0,
        max_beta_val=1,
        learning_rate=0.0005,
        epsilon_std=1.0,
        batch_normalize_inputs=True,
        batch_normalize_intermediaries=True,
        batch_normalize_embedding=True,
        relu_intermediaries=True,
        relu_embedding=True,
        input_dim=None,
        verbose=0,
    ):
        if n_hidden is None:
            n_hidden = [1500]
        self.init_args = {k: v for k, v in locals().items() if k != "self"}
        self.n_hidden = n_hidden
        self.n_latent = n_latent
        self.batch_size = batch_size
        self.epochs = epochs

        if architecture == "stacked":
            self.architecture = partial(
                stacked_vae,
                initial_beta_val=initial_beta_val,
                kappa=kappa,
                max_beta_val=max_beta_val,
                learning_rate=learning_rate,
                epsilon_std=epsilon_std,
                batch_normalize_inputs=batch_normalize_inputs,
                batch_normalize_intermediaries=batch_normalize_intermediaries,
                batch_normalize_embedding=batch_normalize_embedding,
                relu_intermediaries=relu_intermediaries,
                relu_embedding=relu_embedding,
            )
        elif architecture == "deep":
            self.architecture = partial(
                deep_vae,
                initial_beta_val=initial_beta_val,
                kappa=kappa,
                max_beta_val=max_beta_val,
                learning_rate=learning_rate,
                epsilon_std=epsilon_std,
                batch_normalize_inputs=batch_normalize_inputs,
                relu_embedding=relu_embedding,
            )
        else:
            raise ValueError("architecture must be one of 'stacked' or 'deep'")

        # typically, a user will not give input_dim
        # then, self.architecture() will be called on first fit()
        # if input_dim given (e.g. when load()ing from disk), self.architecture() will be called here.
        if input_dim is not None:
            vae, encoder, sampling_encoder, decoder, beta = self.architecture(
                input_dim,
                hidden_dims=self.n_hidden,
                latent_dim=self.n_latent,
                batch_size=self.batch_size,
                epochs=self.epochs,
            )
            self.beta = beta
            self.vae = vae
            self.encoder = encoder
            self.sampling_encoder = sampling_encoder
            self.decoder = decoder

        self.training_fn = partial(
            train_model,
            epochs=epochs,
            batch_size=batch_size,
            kappa=kappa,
            max_beta_val=max_beta_val,
            verbose=verbose,
        )

        self.feature_names = None

    def fit(self, X, y=None, X_validation=None, *args, **kwargs):
        """Train autoencoder model

        Parameters
        ----------
        X:  dict with multi-modal dataframes, containing training data, e.g.
            {'mRNA': df1, 'SNP': df2},
            df1, df2, etc. are (n_features, n_samples) pandas.DataFrame's.
            The sample names must match, the feature names need not.
        X_validation: optional, dict with multi-modal dataframes, containing validation data
            will be used to compute validation loss under training
        y:  Not used.

        Returns
        -------
        self : Maui object
        """
        x_ = self._dict2array(X)
        self._validate_indices(x_)
        self.x_ = x_

        self.nwp_ = self.feature_correlations_ = self.w_ = None

        x_test = self._dict2array(X_validation) if X_validation else self.x_

        if self.feature_names != x_test.columns.tolist():
            raise ValueError("Feature mismatch between X and X_val!")

        if not hasattr(self, "vae"):
            vae, encoder, sampling_encoder, decoder, beta = self.architecture(
                self.x_.shape[1],
                hidden_dims=self.n_hidden,
                latent_dim=self.n_latent,
                batch_size=self.batch_size,
                epochs=self.epochs,
            )
            self.init_args["input_dim"] = self.x_.shape[1]
            self.beta = beta
            self.vae = vae
            self.encoder = encoder
            self.sampling_encoder = sampling_encoder
            self.decoder = decoder
        hist = self.training_fn(
            vae=self.vae, x_train=self.x_, x_val=x_test, beta=self.beta, **kwargs
        )
        self.hist = pd.DataFrame(hist.history)
        return self

    def transform(self, X, encoder="mean"):
        """Transform X into the latent space that was previously learned using
        `fit` or `fit_transform`, and return the latent factor representation.

        Parameters
        ----------
        X:          dict with multi-modal dataframes, containing training data, e.g.
                    {'mRNA': df1, 'SNP': df2},
                    df1, df2, etc. are (n_features, n_samples) pandas.DataFrame's.
        encoder:    the mode of the encoder to be used. one of 'mean' or 'sample',
                    where 'mean' indicates the encoder network only uses the mean
                    estimates for each successive layer. 'sample' indicates the
                    encoder should sample from the distribution specified from each
                    successive layer, and results in non-reproducible embeddings.

        Returns
        -------
        z:  DataFrame (n_samples, n_latent_factors)
            Latent factors representation of the data X.
        """
        if encoder == "mean":
            the_encoder = self.encoder
        elif encoder == "sample":
            the_encoder = self.sampling_encoder
        else:
            raise ValueError("`encoder` must be one of 'mean' or 'sample'")

        x_ = self._dict2array(X)
        self._validate_indices(x_)
        self.x_ = x_
        self.z_ = pd.DataFrame(
            the_encoder.predict(self.x_),
            index=self.x_.index,
            columns=[f"LF{i}" for i in range(1, self.n_latent + 1)],
        )

        self.feature_correlations_ = None
        self.w_ = None
        return self.z_

    def fit_transform(self, X, y=None, X_validation=None, encoder="mean"):
        """Train autoencoder model, and return the latent factor representation
        of the data X.

        Parameters
        ----------
        X:  dict with multi-modal dataframes, containing training data, e.g.
            {'mRNA': df1, 'SNP': df2},
            df1, df2, etc. are (n_samples, n_features) pandas.DataFrame's.
            The sample names must match, the feature names need not.
        X_validation: optional, dict with multi-modal dataframes, containing validation data
            will be used to compute validation loss under training
        y:  Not used.

        Returns
        -------
        z:  DataFrame (n_samples, n_latent_factors)
            Latent factors representation of the data X.
        """
        self.fit(X, X_validation=X_validation, y=y)
        return self.transform(X, encoder=encoder)

    def fine_tune(
        self,
        X,
        batch_size=100,
        epochs=100,
        initial_beta_val=1,
        kappa=1.0,
        max_beta_val=1,
        learning_rate=0.0005,
        verbose=0,
        **kwargs,
    ):
        """Fine-tune autoencoder model. This is meant to be used with a pre-trained model.

        Parameters
        ----------
        X:      dict with multi-modal dataframes, containing training data, e.g.
                    {'mRNA': df1, 'SNP': df2},
                    df1, df2, etc. are (n_features, n_samples) pandas.DataFrame's.
                    The sample names must match, the feature names need not.

        batch_size: int (default 100)
                    The size of the mini-batches used for training the network

        epochs:     int (default 100)
                    The number of epoches to use for training the network

        Returns
        -------
        self : Maui object
        """
        if not hasattr(self, "vae"):
            raise Exception(
                "Cannot call fine_tune() on a model that hasn't been trained yet!"
            )
        x_ = self._dict2array(X)
        self._validate_indices(x_)
        self.x_ = x_

        self.nwp_ = self.feature_correlations_ = self.w_ = None

        fine_tunine_fn = partial(
            train_model,
            epochs=epochs,
            batch_size=batch_size,
            kappa=kappa,
            max_beta_val=max_beta_val,
            verbose=verbose,
        )
        K.set_value(self.beta, initial_beta_val)
        hist = fine_tunine_fn(
            vae=self.vae, x_train=self.x_, x_val=self.x_, beta=self.beta, **kwargs
        )
        self.hist = pd.DataFrame(hist.history)
        return self

    def cluster(
        self,
        k=None,
        optimal_k_method="ami",
        optimal_k_range=range(3, 10),
        ami_y=None,
        kmeans_kwargs=None,
    ):
        """Cluster the samples using k-means based on the latent factors.

        Parameters
        ----------
        k:                  optional, the number of clusters to find.
                            if not given, will attempt to find optimal k.
        optimal_k_method:   supported methods are 'ami' and 'silhouette'. Otherwise, callable.
                            if 'ami', will pick K which gives the best AMI
                            (adjusted mutual information) with external labels.
                            if 'silhouette' will pick the K which gives the best
                            mean silhouette coefficient.
                            if callable, should have signature ``scorer(yhat)``
                            and return a scalar score.
        optimal_k_range:    array-like, range of Ks to try to find optimal K among
        ami_y:              array-like (n_samples), the ground-truth labels to use
                            when picking K by "best AMI against ground-truth" method.
        kmeans_kwargs:      optional, kwargs for initialization of sklearn.cluster.KMeans
                            Default: {"n_init": 1000, "n_jobs": 2}

        Returns
        -------
        yhat:   Series (n_samples) cluster labels for each sample
        """
        if kmeans_kwargs is None:
            kmeans_kwargs = {"n_init": 1000, "n_jobs": 2}
        if k is not None:
            return pd.Series(
                KMeans(k, **kmeans_kwargs).fit_predict(self.z_), index=self.z_.index
            )
        else:
            if optimal_k_method == "ami":
                from sklearn.metrics import adjusted_mutual_info_score

                if ami_y is None:
                    raise Exception(
                        "Must provide ``ami_y`` if using 'ami' to select optimal K."
                    )
                z_to_use = self.z_.loc[ami_y.index]
                scorer = lambda yhat: adjusted_mutual_info_score(ami_y, yhat)
            elif optimal_k_method == "silhouette":
                from sklearn.metrics import silhouette_score

                z_to_use = self.z_
                scorer = lambda yhat: silhouette_score(z_to_use, yhat)
            else:
                z_to_use = self.z_
                scorer = optimal_k_method
            yhats = {
                k: pd.Series(
                    KMeans(k, **kmeans_kwargs).fit_predict(z_to_use),
                    index=z_to_use.index,
                )
                for k in optimal_k_range
            }
            score_name = (
                optimal_k_method
                if isinstance(optimal_k_method, str)
                else optimal_k_method.__name__
            )
            self.kmeans_scores = pd.Series(
                [scorer(yhats[k]) for k in optimal_k_range],
                index=optimal_k_range,
                name=score_name,
            )
            self.kmeans_scores.index.name = "K"
            opt_k_index = np.argmax(self.kmeans_scores)
            self.optimal_k_ = self.kmeans_scores.index[opt_k_index]
            self.yhat_ = yhats[opt_k_index]
            return self.yhat_

    def compute_roc(self, y, **kwargs):
        """Compute Receiver Operating Characteristics curve for SVM prediction
        of labels ``y`` from the latent factors. Computes both the ROC curves
        (true positive rate, true negative rate), and the area under the roc (auc).
        ROC and auROC computed for each class (the classes are inferred from ``y``),
        as well as a "mean" ROC, computed by averaging the class ROCs. Only samples
        in the index of ``y`` will be considered.

        Parameters
        ----------
        y:          array-like (n_samples,), the labels of the samples to predict
        **kwargs:   arguments for ``utils.compute_roc``

        Returns
        -------
        roc_curves: dict, one key per class as well as "mean", each value is a dataframe
                    containing the tpr (true positive rate) and fpr (falce positive rate)
                    defining that class (or the mean) ROC.
        """
        self.roc_curves_ = maui.utils.compute_roc(self.z_, y, **kwargs)
        return self.roc_curves_

    def compute_auc(self, y, **kwargs):
        """Compute area under the ROC curve for predicting the labels in y using the
        latent features previously inferred.

        Parameters
        ----------
        y:          labels to predict
        **kwargs:   arguments for ``compute_roc``

        Returns:
        --------
        aucs:   pd.Series, auc per class as well as mean
        """
        self.compute_roc(y, **kwargs)
        self.aucs_ = {
            k: maui.utils.auc(self.roc_curves_[k].FPR, self.roc_curves_[k].TPR)
            for k in self.roc_curves_
        }
        return self.aucs_

    def select_clinical_factors(
        self,
        survival,
        duration_column="duration",
        observed_column="observed",
        alpha=0.05,
        cox_penalizer=0,
    ):
        """Select latent factors which are predictive of survival. This is
        accomplished by fitting a Cox Proportional Hazards (CPH) model to each
        latent factor, while controlling for known covariates, and only keeping
        those latent factors whose coefficient in the CPH is nonzero (adjusted
        p-value < alpha).

        Parameters
        ----------
        survival:           pd.DataFrame of survival information and relevant covariates
                            (such as sex, age at diagnosis, or tumor stage)
        duration_column:    the name of the column in ``survival`` containing the
                            duration (time between diagnosis and death or last followup)
        observed_column:    the name of the column in ``survival`` containing
                            indicating whether time of death is known
        alpha:              threshold for p-value of CPH coefficients to call a latent
                            factor clinically relevant (p < alpha)
        cox_penalizer:      penalty coefficient in Cox PH solver (see ``lifelines.CoxPHFitter``)

        Returns
        -------
        z_clinical: pd.DataFrame, subset of the latent factors which have been
                    determined to have clinical value (are individually predictive
                    of survival, controlling for covariates)
        """
        self.z_clinical_ = maui.utils.select_clinical_factors(
            self.z_,
            survival,
            duration_column=duration_column,
            observed_column=observed_column,
            alpha=alpha,
            cox_penalizer=cox_penalizer,
        )
        return self.z_clinical_

    def c_index(
        self,
        survival,
        clinical_only=True,
        duration_column="duration",
        observed_column="observed",
        cox_penalties=None,
        cv_folds=5,
        sel_clin_alpha=0.05,
        sel_clin_penalty=0,
    ):
        """Compute's Harrell's c-Index for a Cox Proportional Hazards regression modeling
        survival by the latent factors in z.

        Parameters
        ----------
        z:                  pd.DataFrame (n_samples, n_latent factors)
        survival:           pd.DataFrame of survival information and relevant covariates
                            (such as sex, age at diagnosis, or tumor stage)
        clinical_only:      Compute the c-Index for a model containing only
                            individually clinically relevant latent factors
                            (see ``select_clinical_factors``)
        duration_column:    the name of the column in ``survival`` containing the
                            duration (time between diagnosis and death or last followup)
        observed_column:    the name of the column in ``survival`` containing
                            indicating whether time of death is known
        cox_penalties:      penalty coefficient in Cox PH solver (see ``lifelines.CoxPHFitter``)
                            to try. Returns the best c given by the different penalties
                            (by cross-validation). Defaults to [0.1, 1, 10, 100, 1000, 10000].
        cv_folds:           number of cross-validation folds to compute C
        sel_clin_penalty:   CPH penalizer to use when selecting clinical factors
        sel_clin_alpha:     significance level when selecting clinical factors

        Returns
        -------
        cs: array, Harrell's c-Index, an auc-like metric for survival prediction accuracy.
            one value per cv_fold

        """
        if cox_penalties is None:
            cox_penalties = [0.1, 1, 10, 100, 1000, 10000]
        if clinical_only:
            z = self.select_clinical_factors(
                survival,
                duration_column,
                observed_column,
                sel_clin_alpha,
                sel_clin_penalty,
            )
        else:
            z = self.z_
        return maui.utils.compute_harrells_c(
            z, survival, duration_column, observed_column, cox_penalties, cv_folds
        )

    def get_linear_weights(self):
        """Get linear model coefficients obtained from fitting linear models
        predicting feature values from latent factors. One model is fit per latent
        factor, and the coefficients are stored in the matrix.

        Returns
        -------
        W:  (n_features, n_latent_factors) DataFrame
            w_{ij} is the coefficient associated with feature `i` in a linear model
            predicting it from latent factor `j`.
        """
        if not hasattr(self, "w_") or self.w_ is None:
            self.w_ = maui.utils.map_factors_to_feaures_using_linear_models(
                self.z_, self.x_.T
            )
        return self.w_

    def get_feature_correlations(self):
        """Get correlation coefficients between input features and latent factors.

        Returns
        -------
        r:  (n_features, n_latent_factors) DataFrame
            r_{ij} is the correlation coefficient between feature `i`
            and latent factor `j`.
        """
        if (
            not hasattr(self, "feature_correlations_")
            or self.feature_correlations_ is None
        ):
            self.feature_correlations_ = maui.utils.correlate_factors_and_features(
                self.z_, self.x_
            )
        return self.feature_correlations_

    def get_neural_weight_product(self):
        """Get the product of neural weights from feature input to latent factors.

        Returns
        -------
        nwp:    (n_features, n_latent_factors) DataFrame
                nwp_{ij} is the product of the weights in the encoder network
                along the path leading from input feature `i` to latent
                factor `j`.

        """
        if not hasattr(self, "nwp_") or self.nwp_ is None:
            self.nwp_ = maui.utils.neural_path_weight_product(self)
        return self.nwp_

    def drop_unexplanatory_factors(self, threshold=0.02):
        """Drops factors which have a low R^2 score in a univariate linear model
        predicting the features `x` from a column of the latent factors `z`.

        Parameters
        ----------
        threshold:  threshold for R^2, latent factors below this threshold
                    are dropped.

        Returns
        -------
        z_filt:     (n_samples, n_factors) DataFrame of latent factor values,
                    with only those columns from the input `z` which have an R^2
                    above the threshold when using that column as an input
                    to a linear model predicting `x`.
        """
        if not hasattr(self, "z_") or self.z_ is None:
            raise Exception("Must first transform some data in order to drop columns.")
        self.z_ = maui.utils.filter_factors_by_r2(self.z_, self.x_.T, threshold)
        return self.z_

    def merge_similar_latent_factors(
        self,
        distance_in="z",
        distance_metric="correlation",
        linkage_method="complete",
        distance_threshold=0.17,
        merge_fn=np.mean,
        plot_dendrogram=True,
        plot_dendro_ax=None,
    ):
        """Merge latent factorz in z whose distance is below a certain threshold.
        Used to squeeze down latent factor representations if there are many co-linear
        latent factors.

        Parameters
        ----------
        distance_in:        If 'z', latent factors will be merged based on their distance
                            to each other in 'z'. If 'w', favtors will be merged based
                            on their distance in 'w' (see :func:`get_linear_weights`)
        distance_metric:    The distance metric based on which to merge latent factors.
                            One which is supported by :func:`scipy.spatial.distance.pdist`
        linkage_method:     The linkage method used to cluster latent factors. One which
                            is supported by :func:`scipy.cluster.hierarchy.linkage`.
        distance_threshold: Latent factors with distance below this threshold
                            will be merged
        merge_fn:           Function used to determine value of merged latent factor.
                            The default is :func:`numpy.mean`, meaning the merged
                            latent factor will have the mean value of the inputs.
        plot_dendrogram:    Boolean. If true, a dendrogram will be plotted showing
                            which latent factors are merged and the threshold.
        plot_dendro_ax:     A matplotlib axis object to plot the dendrogram on (optional)

        Returns
        -------
        z:                  (n_samples, n_factors) pd.DataFrame of latent factors
                            where some have been merged
        """
        if not hasattr(self, "z_") or self.z_ is None:
            raise Exception(
                "Cannot merge latent factors before fitting/transforming some."
            )

        if distance_in == "z":
            self.z_ = maui.utils.merge_factors(
                self.z_,
                threshold=distance_threshold,
                merge_fn=merge_fn,
                metric=distance_metric,
                linkage=linkage_method,
                plot_dendro=plot_dendrogram,
                plot_dendro_ax=plot_dendro_ax,
            )
        elif distance_in == "w":
            w = self.get_linear_weights()
            d = spatial.distance.pdist(w.T, distance_metric)
            l = cluster.hierarchy.linkage(d, linkage_method)
            self.z_ = maui.utils.merge_factors(
                self.z_,
                l=l,
                threshold=distance_threshold,
                merge_fn=merge_fn,
                plot_dendro=plot_dendrogram,
                plot_dendro_ax=plot_dendro_ax,
            )
        else:
            raise Exception(
                "Only 'z' and 'w' currently supported. See ``maui.utils.merge_factors for more flexibility.``"
            )
        return self.z_

    def save(self, destdir):
        """Save a maui model to disk, so that it may be reloaded later using ``load()``

        Parameters
        ----------

        destdir:    destination directory in which to save model files
        """
        import os
        import json

        with open(os.path.join(destdir, "maui_args.json"), "wt") as outfile:
            json.dump(self.init_args, outfile)
        self.vae.save_weights(os.path.join(destdir, "maui_weights.h5"))
        with open(os.path.join(destdir, "maui_feature_names.txt"), "wt") as outfile:
            for fn in self.feature_names:
                outfile.write(fn)
                outfile.write("\n")

    @classmethod
    def load(cls, directory):
        """Load a maui model from disk, which was previously saved using ``save()``

        Parameters
        ----------

        directory:    The directory from which to load the maui model

        Returns
        -------

        maui_model: a maui model that was previously saved to disk
        """
        import os
        import json

        with open(os.path.join(directory, "maui_args.json"), "rt") as infile:
            init_args = json.load(infile)
        maui_model = cls(**init_args)
        maui_model.vae.load_weights(os.path.join(directory, "maui_weights.h5"))
        if os.path.isfile(os.path.join(directory, "maui_feature_names.txt")):
            with open(
                os.path.join(directory, "maui_feature_names.txt"), "rt"
            ) as infile:
                maui_model.feature_names = [l.strip() for l in infile]
        else:
            warnings.warn(
                MauiWarning(
                    "Feature names not found in save dir. Loaded model will not perform feature name validation on first call!"
                )
            )
        return maui_model

    def _validate_X(self, X):
        if not isinstance(X, dict):
            raise ValueError("data must be a dict")

        df1 = X[list(X.keys())[0]]
        if any(df.columns.tolist() != df1.columns.tolist() for df in X.values()):
            raise ValueError("All dataframes must have same samples (columns)")

        if any(len(df.index) == 0 for df in X.values()):
            raise ValueError("One of the DataFrames was empty.")

        return True

    def _validate_indices(self, x):
        if self.feature_names is None:
            self.feature_names = x.columns.tolist()
        else:
            if self.feature_names != x.columns.tolist():
                raise ValueError("Feature mismatch with previously fit data!")

    def _dict2array(self, X):
        self._validate_X(X)
        new_feature_names = [f"{k}: {c}" for k in sorted(X.keys()) for c in X[k].index]
        sample_names = X[list(X.keys())[0]].columns
        ret = pd.DataFrame(
            np.vstack([X[k] for k in sorted(X.keys())]).T,
            index=sample_names,
            columns=new_feature_names,
        )
        if self.feature_names is not None:
            if set(self.feature_names) != set(ret.columns):
                raise ValueError("Feature mismatch!")
            ret = ret.loc[:, self.feature_names]
        return ret
