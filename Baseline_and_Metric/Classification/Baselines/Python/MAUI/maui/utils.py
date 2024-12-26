"""
The maui.utils model contains utility functions for multi-omics analysis
using maui.
"""
import keras
import numpy as np
import pandas as pd
from scipy import stats
from scipy import interp
from scipy import spatial
from scipy import cluster
from collections import Counter
from sklearn.svm import LinearSVC
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_predict


def merge_factors(
    z,
    l=None,
    threshold=0.17,
    merge_fn=np.mean,
    metric="correlation",
    linkage="single",
    plot_dendro=True,
    plot_dendro_ax=None,
):
    """Merge latent factors in `z` which form clusters, as defined by hierarchical
    clustering where a cluster is formed by cutting at a pre-set threshold, i.e.
    merge factors if their distance to one-another is below `threshold`.

    Parameters
    ----------
    z:              (n_samples, n_factors) DataFrame of latent factor values, output of a maui model
    metric:         Distance metric to merge factors by, one which is supported by
                    :func:`scipy.spatial.distance.pdist`
    linkage:        The kind of linkage to form hierarchical clustering, one which is
                    supported by :func:`scipy.cluster.hierarchy.linkage`
    l:              As an alternative to supplying `metric` and `linkage`, supply a
                    linkage matrix of your own choice, such as one computed by
                    :func:`scipy.cluster.hierarchy.linkage`
    threshold:      The distance threshold. latent factors with similarity below the
                    threshold will be merged to form single latent facator
    merge_fn:       A function which will be used to merge latent factors. The default
                    is :func:`numpy.mean`, i.e. the newly formed (merged) latent factor
                    will be the mean of the merged ones. Supply any function here which
                    has the same interface, i.e. takes a matrix and an axis.
    plot_dendro:    Boolean. If True, the function will plot a dendrogram showing
                    which latent factors are merged and the threshold.
    """
    if l is None:
        d = spatial.distance.pdist(z.T, metric)
        l = cluster.hierarchy.linkage(d, linkage)

    if plot_dendro:
        try:
            import matplotlib.pyplot as plt
        except:
            raise Exception("`plot_dendro` require matplotlib to be installed.")

        if plot_dendro_ax is None:
            fig, plot_dendro_ax = plt.subplots(figsize=(25, 10))
        dendro = cluster.hierarchy.dendrogram(
            l, leaf_font_size=18, color_threshold=threshold, ax=plot_dendro_ax
        )
        plot_dendro_ax.hlines(threshold, *plot_dendro_ax.get_xlim(), "red", "dashed")

    cl_labels = cluster.hierarchy.cut_tree(l, height=threshold).T[0]
    factors_to_delete = set()
    new_factors = list()
    for cl, ct in Counter(cl_labels).items():
        if ct < 2:
            continue
        tomerge = (cl_labels == cl).nonzero()[0]
        factors_to_delete.update(tomerge)
        new_factors.append(
            merge_fn(z.iloc[:, tomerge], axis=1).rename(
                "_".join(f"{i}" for i in tomerge)
            )
        )

    new_z = z.copy()
    new_z = new_z.loc[:, [c for c in new_z.columns if c not in factors_to_delete]]
    new_z = new_z.iloc[
        :, [i for i in range(new_z.shape[1]) if i not in factors_to_delete]
    ]
    new_z = pd.concat([new_z] + new_factors, axis=1)
    return new_z

def r2_score_factors(z, x):
    """Calculate R^2 of a linear model predicting features x from
    latent factors z.

    Parameters
    ----------
    z:  (n_samples, n_factors) DataFrame of latent factor values, output of a maui model
    x:  (n_samples, n_features) DataFrame of concatenated multi-omics data

    Returns
    -------
    scores: (n_factors,) Series of R^2 scores for each latent factor in `z`
            obtained from fitting a linear model predicting `x` from that latent
            factor alone.
    """
    scores = list()
    for i in range(z.shape[1]):
        regressor = LinearRegression()
        regressor.fit(z.values[:, i].reshape(-1, 1), scale(x).T)
        scores.append(regressor.score(z.values[:, i].reshape(-1, 1), scale(x).T))
    return scores


def filter_factors_by_r2(z, x, threshold=0.02):
    """Filter latent factors by the R^2 of a linear model predicting features x
    from latent factors z.

    Parameters
    ----------
    z:  (n_samples, n_factors) DataFrame of latent factor values, output of a maui model
    x:  (n_samples, n_features) DataFrame of concatenated multi-omics data

    Returns
    -------
    z_filtered: (n_samples, n_factors) DataFrame of latent factor values,
                with only those columns from the input `z` which have an R^2
                above the threshold when using that column as an input
                to a linear model predicting `x`.
    """
    scores = r2_score_factors(z, x)
    return z.iloc[:, pd.Series(scores)[pd.Series(scores) > threshold].index]


def neural_path_weight_product(model):
    """The product of weights along the path from from input
    features to latent factors.

    Parameters
    ----------
    model: a maui model

    Returns
    -------
    nwp: (n_features, n_latent_factors) DataFrame of neural weight products
    """
    path_ws = list()
    for layer in model.encoder.layers:
        if isinstance(layer, keras.layers.core.Dense):
            path_ws.append(layer.get_weights()[0])

    return pd.DataFrame(
        np.linalg.multi_dot(path_ws), index=model.x_.columns, columns=model.z_.columns
    )


def map_factors_to_feaures_using_linear_models(z, x):
    """Get feature <-> latent factors mapping from linear models.
    Runs one univariate (multi-output) linear model per latent factor in `z`,
    predicting the values of the features `x`, in order to get weights
    between inputs and outputs.

    Parameters
    ----------
    z:  (n_samples, n_factors) DataFrame of latent factor values, output of a maui model
    x:  (n_samples, n_features) DataFrame of concatenated multi-omics data

    Returns
    -------
    W:  (n_features, n_latent_factors) DataFrame
        w_{ij} is the coefficient associated with feature `i` in a linear model
        predicting it from latent factor `j`.
    """
    ws = list()
    for i in range(z.shape[1]):
        regressor = LinearRegression()
        regressor.fit(z.values[:, i].reshape(-1, 1), scale(x).T)
        ws.append(regressor.coef_)
    W = pd.DataFrame(np.hstack(ws), index=x.index, columns=z.columns)
    return W


def correlate_factors_and_features(z, concatenated_data, pval_threshold=0.001):
    """Compute pearson correlation of latent factors with input features.

    Parameters
    ----------
    z:                  (n_samples, n_factors) DataFrame of latent factor values, output of maui model
    concatenated_data:  (n_samples, n_features) DataFrame of concatenated multi-omics data

    Returns
    -------
    feature_s:  DataFrame (n_features, n_latent_factors)
                Latent factors representation of the data X.
    """
    feature_scores = list()
    for j in range(z.shape[1]):
        corrs = pd.DataFrame(
            [
                stats.pearsonr(concatenated_data.iloc[:, i], z.iloc[:, j])
                for i in range(concatenated_data.shape[1])
            ],
            index=concatenated_data.columns,
            columns=["r", "pval"],
        )
        corrs.loc[corrs.pval > pval_threshold, "r"] = 0
        feature_scores.append(corrs.r)

    feature_s = pd.concat(feature_scores, axis=1)
    feature_s.columns = [i for i in range(z.shape[1])]
    return feature_s


def compute_roc(z, y, classifier=LinearSVC(C=0.001), cv_folds=10):
    """Compute the ROC (false positive rate, true positive rate) using cross-validation.

    Parameters
    ----------
    z:          DataFrame (n_samples, n_latent_factors) of latent factor values
    y:          Series (n_samples,) of ground-truth labels to try to predict
    classifier: Classifier object to use, default ``LinearSVC(C=.001)``

    Returns
    -------
    roc_curves: dict, one key per class as well as "mean", each value is a dataframe
                containing the tpr (true positive rate) and fpr (false positive rate)
                defining that class (or the mean) ROC.
    """
    class_names = sorted(y.unique())
    z_to_use = z.loc[y.index]
    y_true_bin = label_binarize(y, classes=class_names)
    y_proba = cross_val_predict(
        classifier, z_to_use, y, cv=cv_folds, method="decision_function"
    )

    # Compute ROC curve and ROC area for each class
    roc_curves = dict()
    for i, cl_name in enumerate(class_names):
        fpr, tpr, thresholds = roc_curve(y_true_bin[:, i], y_proba[:, i])
        roc_curves[cl_name] = pd.concat(
            [pd.Series(fpr, name="FPR"), pd.Series(tpr, name="TPR")], axis=1
        )

    mean_fpr = np.unique(
        np.concatenate([roc_curves[cl_name].FPR for cl_name in class_names])
    )

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(mean_fpr)
    for cl_name in class_names:
        mean_tpr += interp(mean_fpr, roc_curves[cl_name].FPR, roc_curves[cl_name].TPR)

    # Finally average it
    mean_tpr /= len(class_names)

    roc_curves["mean"] = pd.concat(
        [pd.Series(mean_fpr, name="FPR"), pd.Series(mean_tpr, name="TPR")], axis=1
    )
    return roc_curves


def estimate_kaplan_meier(
    y, survival, duration_column="duration", observed_column="observed"
):
    """Estimate survival curves for groups defined in y based on survival data in ``survival``

    Parameters
    ----------
    y:                  pd.Series, groups (clusters, subtypes). the index is
                        the sample names
    survival:           pd.DataFrame with the same index as y, with columns for
                        the duration (survival time for each patient) and whether
                        or not the death was observed. If the death was not
                        observed (censored), the duration is the time of the last
                        followup.
    duration_column:        the name of the column in  ``survival`` with the duration
    observed_column:    the name of the column in ``survival`` with True/False values
                        for whether death was observed or not

    Returns
    -------
    km_estimates:       pd.DataFrame, index is the timeline, columns are survival
                        functions (estimated by Kaplan-Meier) for each class, as
                        defined in ``y``.
    """
    try:
        import lifelines
    except ImportError:
        raise ImportError(
            "The module ``lifelines`` was not found. It is required for this functionality. You may install it using `pip install lifelines`."
        )
    kmf = lifelines.KaplanMeierFitter()
    sfs = dict()
    for cl in y.unique():
        ixs = list(set(y[y == cl].index) & set(survival.index))
        kmf.fit(
            survival.loc[ixs][duration_column],
            survival.loc[ixs][observed_column],
            label=cl,
        )
        sfs[cl] = kmf.survival_function_
    return pd.concat([sfs[k] for k in sorted(y.unique())], axis=1).interpolate()


def multivariate_logrank_test(
    y, survival, duration_column="duration", observed_column="observed"
):
    """Compute the multivariate log-rank test for differential survival
    among the groups defined by ``y`` in the survival data in ``survival``,
    under the null-hypothesis that all groups have the same survival function
    (i.e. test whether at least one group has different survival rates)

    Parameters
    ----------
    y:                  pd.Series, groups (clusters, subtypes). the index is
                        the sample names
    survival:           pd.DataFrame with the same index as y, with columns for
                        the duration (survival time for each patient) and whether
                        or not the death was observed. If the death was not
                        observed (sensored), the duration is the time of the last
                        followup.
    duration_column:        the name of the column in  ``survival`` with the duration
    observed_column:    the name of the column in ``survival`` with True/False values
                        for whether death was observed or not

    Returns
    -------
    test_statistic:     the test statistic (chi-square)
    p_value:            the associated p_value
    """
    try:
        import lifelines
    except ImportError:
        raise ImportError(
            "The module ``lifelines`` was not found. It is required for this functionality. You may install it using `pip install lifelines`."
        )
    ixs = list(set(y.index) & set(survival.index))
    mlr = lifelines.statistics.multivariate_logrank_test(
        survival.loc[ixs][duration_column],
        y.loc[ixs],
        survival.loc[ixs][observed_column],
    )
    return mlr.test_statistic, mlr.p_value


def select_clinical_factors(
    z,
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
    cox_coefficients = _cph_coefs(
        z, survival, duration_column, observed_column, penalizer=cox_penalizer
    )
    signif_cox_coefs = cox_coefficients.T[cox_coefficients.T.p < alpha]
    return z.loc[:, signif_cox_coefs.index]


def _cph_coefs(z, survival, duration_column, observed_column, penalizer=0):
    """Compute one CPH model for each latent factor (column) in z.
    Return summaries (beta values, p values, confidence intervals)
    """
    try:
        import lifelines
    except ImportError:
        raise ImportError(
            "The module ``lifelines`` was not found. It is required for this functionality. You may install it using `pip install lifelines`."
        )
    return pd.concat(
        [
            lifelines.CoxPHFitter(penalizer=penalizer)
            .fit(
                survival.assign(LF=z.loc[:, i]).dropna(),
                duration_column,
                observed_column,
            )
            .summary.loc["LF"]
            .rename(i)
            for i in z.columns
        ],
        axis=1,
    )


def compute_harrells_c(
    z,
    survival,
    duration_column="duration",
    observed_column="observed",
    cox_penalties=None,
    cv_folds=5,
):
    """Compute's Harrell's c-Index for a Cox Proportional Hazards regression modeling
    survival by the latent factors in z.

    Parameters
    ----------
    z:                  pd.DataFrame (n_samples, n_latent factors)
    survival:           pd.DataFrame of survival information and relevant covariates
                        (such as sex, age at diagnosis, or tumor stage)
    duration_column:    the name of the column in ``survival`` containing the
                        duration (time between diagnosis and death or last followup)
    observed_column:    the name of the column in ``survival`` containing
                        indicating whether time of death is known
    cox_penalties:      penalty coefficient in Cox PH solver (see ``lifelines.CoxPHFitter``)
                        to try. Returns the best c given by the different penalties
                        (by cross-validation). Defualt: [0.1, 1, 10, 100, 1000, 10000]
    cv_folds:           number of cross-validation folds to compute C

    Returns
    -------
    cs: array, Harrell's c-Index, an auc-like metric for survival prediction accuracy.
        one value per cv_fold

    """
    if cox_penalties is None:
        cox_penalties = [0.1, 1, 10, 100, 1000, 10000]
    cvcs = [
        _cv_coxph_c(z, survival, p, duration_column, observed_column, cv_folds)
        for p in cox_penalties
    ]
    return cvcs[np.argmax([np.median(e) for e in cvcs])]


def _cv_coxph_c(
    z,
    survival,
    penalty,
    duration_column="duration",
    observed_column="observed",
    cv_folds=5,
):
    try:
        import lifelines
        import lifelines.utils
    except ImportError:
        raise ImportError(
            "The module ``lifelines`` was not found. It is required for this functionality. You may install it using `pip install lifelines`."
        )

    cph = lifelines.CoxPHFitter(penalizer=penalty)
    survdf = pd.concat([survival, z], axis=1, sort=False).dropna()

    kfold = KFold(cv_folds)
    scores = list()

    for train_index, test_index in kfold.split(survdf):
        x_train, x_test = survdf.iloc[train_index], survdf.iloc[test_index]

        cph.fit(x_train, duration_column, observed_column)
        cindex = lifelines.utils.concordance_index(
            x_test[duration_column],
            -cph.predict_partial_hazard(x_test),
            x_test[observed_column],
        )
        scores.append(cindex)

    return scores


def scale(df):
    """Scale and center data

    Parameters
    ----------
    df:     pd.DataFrame (n_features, n_samples) non-scaled data

    Returns
    -------
    scaled: pd.DataFrame (n_features, n_samples) scaled data
    """
    df_scaled = StandardScaler().fit_transform(df.T)
    return pd.DataFrame(df_scaled, columns=df.index, index=df.columns).T
