import pytest
import numpy as np
import pandas as pd
from scipy import spatial, cluster

from maui import utils


def test_merge_factors():
    z = pd.DataFrame(
        [
            [1, 1, 1, 0, 0, 0, 1, 0, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 1, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 1, 1, 1, 0],
        ],
        index=[f"sample {i}" for i in range(11)],
        columns=[f"LF{i}" for i in range(9)],
        dtype=float,
    )  # expect 0,1,2 to be merged, and 3,7 to be merged

    z_merged = utils.merge_factors(z, metric="euclidean", plot_dendro=False)

    assert z_merged.shape[1] == 6
    assert "0_1_2" in z_merged.columns
    assert "3_7" in z_merged.columns


def test_merge_factors_with_custom_linkage():
    z = pd.DataFrame(
        [
            [1, 1, 1, 0, 0, 0, 1, 0, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 1, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 1, 1, 1, 0],
        ],
        index=[f"sample {i}" for i in range(11)],
        columns=[f"LF{i}" for i in range(9)],
        dtype=float,
    )  # expect 0,1,2 to be merged, and 3,7 to be merged

    l = cluster.hierarchy.linkage(spatial.distance.pdist(z.T, "minkowski"), "average")

    z_merged = utils.merge_factors(z, l=l, plot_dendro=False)

    assert z_merged.shape[1] == 6
    assert "0_1_2" in z_merged.columns
    assert "3_7" in z_merged.columns


def test_filter_factors_by_r2():
    dummy_z = pd.DataFrame(
        [[0, 1, 2], [1, 0, 2]],
        index=["sample 1", "sample 2"],
        columns=["LF1", "LF2", "LF3"],
    )

    dummy_x = pd.DataFrame(
        [[1.0, 1.0, 1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]],
        columns=[f"feature{i}" for i in range(6)],
        index=["sample 1", "sample 2"],
    ).T

    z_filt = utils.filter_factors_by_r2(dummy_z, dummy_x)
    assert z_filt.columns.tolist() == ["LF1", "LF2"]

def test_map_factors_to_feaures_using_linear_models():
    dummy_z = pd.DataFrame(
        [[0, 1], [1, 0]], index=["sample 1", "sample 2"], columns=["LF1", "LF2"]
    )

    dummy_x = pd.DataFrame(
        [[1.0, 1.0, 1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]],
        columns=[f"feature{i}" for i in range(6)],
        index=["sample 1", "sample 2"],
    ).T

    expected_w = np.array(
        [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0], [2.0, -2.0], [2.0, -2.0], [2.0, -2.0]]
    )

    w = utils.map_factors_to_feaures_using_linear_models(dummy_z, dummy_x)

    assert np.allclose(w, expected_w)


def test_correlate_factors_and_features():
    dummy_z = pd.DataFrame(
        [[0, 1], [1, 0]], index=["sample 1", "sample 2"], columns=["LF1", "LF2"]
    )

    dummy_x = pd.DataFrame(
        [[1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1]],
        columns=[f"feature{i}" for i in range(6)],
        index=["sample 1", "sample 2"],
    )

    expected_corrs = np.array(
        [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0], [1.0, -1.0], [1.0, -1.0], [1.0, -1.0]]
    )

    corrs = utils.correlate_factors_and_features(dummy_z, dummy_x)

    assert np.allclose(corrs, expected_corrs)


def test_compute_roc():
    np.random.seed(0)
    dummy_z = pd.DataFrame(
        [
            [0, 1, 1, 1, 0, 1, 1, 0, 0],
            [1, 0, 0, 0, 0, 0, 1, 1, 0],
            [1, 0, 1, 0, 0, 0, 1, 1, 0],
            [1, 0, 0, 1, 0, 0, 1, 1, 0],
            [1, 0, 0, 0, 1, 1, 1, 1, 0],
            [1, 1, 1, 0, 0, 0, 1, 1, 1],
        ],
        index=[f"sample {i}" for i in range(6)],
        columns=[f"LF{i}" for i in range(9)],
    )
    dummy_y = pd.Series(["a", "b", "a", "c", "b", "c"], index=dummy_z.index)

    roc_curves = utils.compute_roc(dummy_z, dummy_y, cv_folds=2)
    assert np.allclose(roc_curves["a"].FPR, [0.0, 0.5, 0.5, 0.75, 1.0])


def test_compute_auc():
    fpr = [0.0, 0.0, 0.5, 0.5, 1.0]
    tpr = [0.0, 0.5, 0.5, 1.0, 1.0]
    roc = utils.auc(fpr, tpr)
    assert roc - 0.75 < 1e-6


def test_estimate_km():
    yhat = pd.Series(
        ["a", "a", "a", "b", "b", "b"], index=[f"Sample {i}" for i in range(6)]
    )
    durations = np.random.poisson(6, 6)
    observed = np.random.randn(6) > 0.1
    survival = pd.DataFrame(
        dict(duration=durations, observed=observed),
        index=[f"Sample {i}" for i in range(6)],
    )
    km = utils.estimate_kaplan_meier(yhat, survival)

    assert "a" in km.columns
    assert "b" in km.columns


def test_multivariate_logrank_test():
    yhat = pd.Series(
        ["a", "a", "a", "b", "b", "b"], index=[f"Sample {i}" for i in range(6)]
    )
    durations = np.random.poisson(6, 6)
    observed = np.random.randn(6) > 0.1
    survival = pd.DataFrame(
        dict(duration=durations, observed=observed),
        index=[f"Sample {i}" for i in range(6)],
    )
    test_stat, p_val = utils.multivariate_logrank_test(yhat, survival)
    assert p_val <= 1.0


def test_select_clinical_factors():
    dummy_z = pd.DataFrame(
        [
            [1, 1, 1, 0, 0, 0, 1, 0, 1],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 1, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 1, 1, 1, 1],
        ],
        index=[f"sample {i}" for i in range(11)],
        columns=[f"LF{i}" for i in range(9)],
    )  # here the first 3 factors separate the groups and the last 6 do not

    durations = [
        1,
        2,
        3,
        4,
        5,
        6,
        1000,
        2000,
        3000,
        4000,
        5000,
    ]  # here the first 3 have short durations, the last 3 longer ones
    observed = [True] * 11  # all events observed
    survival = pd.DataFrame(
        dict(duration=durations, observed=observed),
        index=[f"sample {i}" for i in range(11)],
    )

    z_clinical = utils.select_clinical_factors(dummy_z, survival, cox_penalizer=1, alpha=.1)
    assert "LF0" in z_clinical.columns
    assert "LF1" in z_clinical.columns
    assert "LF2" in z_clinical.columns

    assert "LF3" not in z_clinical.columns
    assert "LF4" not in z_clinical.columns
    assert "LF5" not in z_clinical.columns


def test_compute_harrells_c():
    dummy_z = pd.DataFrame(
        [
            [1, 1, 1, 0, 0, 0, 1, 0, 1],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 0],
            [1, 1, 1, 1, 1, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 0, 1, 1, 1, 1],
        ],
        index=[f"sample {i}" for i in range(11)],
        columns=[f"LF{i}" for i in range(9)],
    )  # here the first 3 factors separate the groups and the last 6 do not

    durations = [
        1,
        2,
        3,
        4,
        5,
        6,
        1000,
        2000,
        3000,
        4000,
        5000,
    ]  # here the first 3 have short durations, the last 3 longer ones
    observed = [True] * 11  # all events observed
    survival = pd.DataFrame(
        dict(duration=durations, observed=observed),
        index=[f"sample {i}" for i in range(11)],
    )
    z_clinical = utils.select_clinical_factors(dummy_z, survival, cox_penalizer=1, alpha=.1)

    np.random.seed(0)
    c = utils.compute_harrells_c(z_clinical, survival, cv_folds=3)
    assert np.allclose(c, [0.5, 0.8, 0.5], atol=.05)
