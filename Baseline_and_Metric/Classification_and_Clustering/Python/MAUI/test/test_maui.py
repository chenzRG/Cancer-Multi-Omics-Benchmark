import os
import pytest
import tempfile
from unittest import mock

import numpy as np
import pandas as pd

from maui import Maui
from maui.maui_warnings import MauiWarning


samples = [f"Sample_{i}" for i in range(10)]

n_features_1 = 20
df1 = pd.DataFrame(
    np.random.randn(n_features_1, len(samples)),
    columns=samples,
    index=[f"data1_feature_{i}" for i in range(n_features_1)],
)
n_features_2 = 6
df2 = pd.DataFrame(
    np.random.randn(n_features_2, len(samples)),
    columns=samples,
    index=[f"data2_feature_{i}" for i in range(n_features_2)],
)

df_empty = pd.DataFrame(
    np.random.randn(0, len(samples)),
    columns=samples,
    index=[f"data0_feature_{i}" for i in range(0)],
)


def test_validate_X_fails_if_not_dict():
    maui_model = Maui()
    with pytest.raises(ValueError):
        maui_model._validate_X([1, 2, 3])


def test_validate_X_fails_if_samples_mismatch():
    maui_model = Maui()
    with pytest.raises(ValueError):
        df2_bad = df2.iloc[:, :2]
        data_with_mismatching_samples = {"a": df1, "b": df2_bad}
        maui_model._validate_X(data_with_mismatching_samples)


def test_validate_X_fails_if_some_data_empty():
    maui_model = Maui()
    with pytest.raises(ValueError):
        maui_model._validate_X({"a": df1, "e": df_empty})


def test_validate_X_returns_true_on_valid_data():
    maui_model = Maui()
    valid_data = {"a": df1, "b": df2}
    assert maui_model._validate_X(valid_data)


def test_dict2array():
    maui_model = Maui()
    arr = maui_model._dict2array({"data1": df1, "data2": df2})
    assert arr.shape[0] == len(df1.columns)
    assert arr.shape[1] == len(df1.index) + len(df2.index)


def test_maui_saves_feature_correlations():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    z = maui_model.fit_transform({"d1": df1, "d2": df2})
    r = maui_model.get_feature_correlations()
    assert r is not None
    assert hasattr(maui_model, "feature_correlations_")


def test_maui_saves_w():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    z = maui_model.fit_transform({"d1": df1, "d2": df2})
    w = maui_model.get_linear_weights()
    assert w is not None
    assert hasattr(maui_model, "w_")


def test_maui_saves_neural_weight_product():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    z = maui_model.fit_transform({"d1": df1, "d2": df2})
    nwp = maui_model.get_neural_weight_product()
    assert nwp is not None
    assert hasattr(maui_model, "nwp_")

    print(maui_model.encoder.summary())

    w1 = maui_model.encoder.layers[2].get_weights()[0]
    w2 = maui_model.encoder.layers[3].get_weights()[0]

    nwp_11 = np.dot(w1[0, :], w2[:, 0])
    assert np.allclose(nwp_11, nwp.iloc[0, 0])


def test_maui_updates_neural_weight_product_when_training():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)

    z_before = maui_model.fit_transform({"d1": df1, "d2": df2})
    nwp_before_fine_tuning = maui_model.get_neural_weight_product()

    maui_model.fine_tune({"d1": df1, "d2": df2})
    z_after = maui_model.transform({"d1": df1, "d2": df2})
    nwp_after_fine_tuning = maui_model.get_neural_weight_product()

    assert not np.allclose(z_before, z_after)
    assert not np.allclose(nwp_before_fine_tuning, nwp_after_fine_tuning)


def test_maui_clusters_with_single_k():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
        np.random.randn(10, 2),
        index=[f"sample {i}" for i in range(10)],
        columns=["LF1", "LF2"],
    )
    maui_model.x_ = pd.DataFrame(
        np.random.randn(20, 10),
        index=[f"feature {i}" for i in range(20)],
        columns=[f"sample {i}" for i in range(10)],
    )

    yhat = maui_model.cluster(5)
    assert yhat.shape == (10,)


def test_maui_clusters_picks_optimal_k_by_ami():
    ami_mock = mock.Mock()
    ami_mock.side_effect = [
        2,
        3,
        1,
    ]  # the optimal AMI will be given at the second trial
    with mock.patch("sklearn.metrics.adjusted_mutual_info_score", ami_mock):
        maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
        maui_model.z_ = pd.DataFrame(
            np.random.randn(10, 2),
            index=[f"sample {i}" for i in range(10)],
            columns=["LF1", "LF2"],
        )
        maui_model.x_ = pd.DataFrame(
            np.random.randn(20, 10),
            index=[f"feature {i}" for i in range(20)],
            columns=[f"sample {i}" for i in range(10)],
        )

        the_y = pd.Series(np.arange(10), index=maui_model.z_.index)

        maui_model.cluster(
            ami_y=the_y, optimal_k_range=[1, 2, 3]
        )  # the second trial is k=2
        print(maui_model.kmeans_scores)
        assert maui_model.optimal_k_ == 2


def test_maui_clusters_picks_optimal_k_by_silhouette():
    silhouette_mock = mock.Mock()
    silhouette_mock.side_effect = [
        2,
        3,
        1,
    ]  # the optimal silhouette will be given at the second trial
    with mock.patch("sklearn.metrics.silhouette_score", silhouette_mock):
        maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
        maui_model.z_ = pd.DataFrame(
            np.random.randn(10, 2),
            index=[f"sample {i}" for i in range(10)],
            columns=["LF1", "LF2"],
        )
        maui_model.x_ = pd.DataFrame(
            np.random.randn(20, 10),
            index=[f"feature {i}" for i in range(20)],
            columns=[f"sample {i}" for i in range(10)],
        )
        maui_model.cluster(
            optimal_k_method="silhouette", optimal_k_range=[1, 2, 3]
        )  # the second trial is k=2

        assert maui_model.optimal_k_ == 2


def test_maui_clusters_picks_optimal_k_with_custom_scoring():
    scorer = mock.Mock()
    scorer.side_effect = [2, 3, 1]  # the optimal AMI will be given at the second trial
    scorer.__name__ = "mock_scorer"

    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
        np.random.randn(10, 2),
        index=[f"sample {i}" for i in range(10)],
        columns=["LF1", "LF2"],
    )
    maui_model.x_ = pd.DataFrame(
        np.random.randn(20, 10),
        index=[f"feature {i}" for i in range(20)],
        columns=[f"sample {i}" for i in range(10)],
    )
    maui_model.cluster(
        optimal_k_method=scorer, optimal_k_range=[1, 2, 3]
    )  # the second trial is k=2

    assert maui_model.optimal_k_ == 2


def test_maui_computes_roc_and_auc():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
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
    y = pd.Series(["a", "b", "a", "c", "b", "c"], index=maui_model.z_.index)
    rocs = maui_model.compute_roc(y, cv_folds=2)
    assert rocs == maui_model.roc_curves_
    assert "a" in rocs
    assert "b" in rocs
    assert "c" in rocs
    assert "mean" in rocs

    aucs = maui_model.compute_auc(y, cv_folds=2)
    assert aucs == maui_model.aucs_


def test_maui_clusters_only_samples_in_y_index_when_optimizing():
    np.random.seed(0)
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
        np.random.randn(10, 2),
        index=[f"sample {i}" for i in range(10)],
        columns=["LF1", "LF2"],
    )
    maui_model.x_ = pd.DataFrame(
        np.random.randn(20, 10),
        index=[f"feature {i}" for i in range(20)],
        columns=[f"sample {i}" for i in range(10)],
    )

    y = pd.Series(
        ["a", "a", "a", "b", "b", "b"], index=[f"sample {i}" for i in range(6)]
    )

    yhat = maui_model.cluster(ami_y=y, optimal_k_range=[1, 2, 3])
    assert set(yhat.index) == set(y.index)


def test_select_clinical_factors():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
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

    z_clin = maui_model.select_clinical_factors(survival, cox_penalizer=1, alpha=0.1)
    assert "LF0" in z_clin.columns
    assert "LF5" not in z_clin.columns


def test_maui_computes_harrells_c():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
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
    cs = maui_model.c_index(
        survival,
        clinical_only=True,
        duration_column="duration",
        observed_column="observed",
        cox_penalties=[0.1, 1, 10, 100, 1000, 10000],
        cv_folds=3,
        sel_clin_alpha=0.1,
        sel_clin_penalty=1,
    )
    print(cs)
    assert np.allclose(cs, [0.5, 0.8, 0.5], atol=0.05)


def test_maui_produces_same_prediction_when_run_twice():
    """This is to show the maui encoder model picks the mean of
    the distribution, not a sample."""
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    z1 = maui_model.transform({"d1": df1, "d2": df2})
    z2 = maui_model.transform({"d1": df1, "d2": df2})
    assert np.allclose(z1, z2)


def test_maui_produces_different_prediction_when_run_twice_with_sampling():
    """This is to show the maui encoder model picks the mean of
    the distribution, not a sample."""
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    z1 = maui_model.transform({"d1": df1, "d2": df2}, encoder="sample")
    z2 = maui_model.transform({"d1": df1, "d2": df2}, encoder="sample")
    assert not np.allclose(z1, z2)


def test_maui_produces_nonnegative_zs_if_relu_embedding_true():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1, relu_embedding=True)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    z1 = maui_model.transform({"d1": df1, "d2": df2})
    assert np.all(z1 >= 0)


def test_maui_produces_pos_and_neg_zs_if_relu_embedding_false():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1, relu_embedding=False)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    z1 = maui_model.transform({"d1": df1, "d2": df2})
    assert not np.all(z1 >= 0)


def test_maui_runs_with_deep_not_stacked_vae():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1, architecture="deep")
    z = maui_model.fit_transform({"d1": df1, "d2": df2})


def test_maui_complains_if_wrong_architecture():
    with pytest.raises(ValueError):
        maui_model = Maui(
            n_hidden=[10], n_latent=2, epochs=1, architecture="wrong value"
        )


def test_maui_supports_single_layer_vae():
    maui_model = Maui(n_hidden=None, n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    z1 = maui_model.transform({"d1": df1, "d2": df2})


def test_maui_supports_not_deep_deep_vae():
    maui_model = Maui(n_hidden=None, n_latent=2, epochs=1, architecture="deep")
    z = maui_model.fit_transform({"d1": df1, "d2": df2})


def test_maui_drops_unexplanatody_factors_by_r2():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
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
    )  # here the first 8 latent factors have R2 above threshold, the last does not
    maui_model.x_ = pd.DataFrame(
        [[1], [1], [1], [1], [1], [1], [0], [0], [0], [0], [0]],
        index=[f"sample {i}" for i in range(11)],
        columns=["Feature 1"],
        dtype=float,
    )

    z_filt = maui_model.drop_unexplanatory_factors()

    assert z_filt.shape[1] == 8


def test_maui_merges_latent_factors():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
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

    z_merged = maui_model.merge_similar_latent_factors(distance_metric="euclidean")
    assert z_merged.shape[1] == 6
    assert "0_1_2" in z_merged.columns
    assert "3_7" in z_merged.columns


def test_maui_merges_latent_factors_by_w():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
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
    )
    maui_model.x_ = pd.DataFrame(
        [[1], [1], [1], [1], [1], [1], [0], [0], [0], [0], [0]],
        index=[f"sample {i}" for i in range(11)],
        columns=["Feature 1"],
        dtype=float,
    )
    # with these z and x, expect 0,1,2 and 4,5 and 3,6,7
    z_merged = maui_model.merge_similar_latent_factors(
        distance_in="w", distance_metric="euclidean"
    )
    assert z_merged.shape[1] == 4
    assert "0_1_2" in z_merged.columns
    assert "3_6_7" in z_merged.columns
    assert "4_5" in z_merged.columns


def test_maui_merge_latent_factors_complains_if_unknown_merge_by():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.z_ = pd.DataFrame(
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

    with pytest.raises(Exception):
        z_merged = maui_model.merge_similar_latent_factors(
            distance_in="xxx", distance_metric="euclidean"
        )


def test_maui_can_save_to_folder():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    with tempfile.TemporaryDirectory() as tmpdirname:
        maui_model.save(tmpdirname)
        assert os.path.isfile(os.path.join(tmpdirname, "maui_weights.h5"))
        assert os.path.isfile(os.path.join(tmpdirname, "maui_args.json"))


def test_maui_can_load_from_folder():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    with tempfile.TemporaryDirectory() as tmpdirname:
        maui_model.save(tmpdirname)
        maui_model_from_disk = Maui.load(tmpdirname)

    assert maui_model_from_disk.n_latent == maui_model.n_latent
    assert np.allclose(
        maui_model.vae.get_weights()[0], maui_model_from_disk.vae.get_weights()[0]
    )
    assert np.allclose(
        maui_model.transform({"d1": df1, "d2": df2}),
        maui_model_from_disk.transform({"d1": df1, "d2": df2}),
    )


def test_maui_can_print_verbose_training(capsys):
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})

    stdout, stderr = capsys.readouterr()
    assert stdout == ""

    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1, verbose=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})

    stdout, stderr = capsys.readouterr()
    assert "Epoch" in stdout


def test_maui_model_makes_2_layer_vae():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1, input_dim=10)
    layers_names = [l.name for l in maui_model.vae.layers]

    assert "hidden_dim_0_mean" in layers_names
    assert "latent_mean" in layers_names
    assert "decode_hidden_0" in layers_names
    assert "reconstruction" in layers_names

    assert "decode_hidden_1" not in layers_names


def test_maui_model_makes_one_layer_vae():
    maui_model = Maui(n_hidden=[], n_latent=2, epochs=1, input_dim=10)
    layers_names = [l.name for l in maui_model.vae.layers]

    print(layers_names)

    assert layers_names[-1] == "reconstruction"

    assert not any(
        "decode_hidden" in name for name in layers_names
    ), "Has a decode hidden..."
    assert not any("hidden_dim" in name for name in layers_names), "Has a hidden dim..."


def test_maui_model_validates_feature_names_on_predict_after_fit():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model.fit({"d1": df1, "d2": df2})

    z = maui_model.transform({"d1": df1, "d2": df2})

    df1_wrong_features = df1.reindex(df1.index[: len(df1.index) - 1])
    with pytest.raises(ValueError):
        z = maui_model.transform({"df1": df1_wrong_features, "df2": df2})


def test_maui_model_saves_feature_names_to_disk():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    with tempfile.TemporaryDirectory() as tmpdirname:
        maui_model.save(tmpdirname)
        maui_model_from_disk = Maui.load(tmpdirname)
    assert maui_model.feature_names == maui_model_from_disk.feature_names


def test_maui_model_loads_model_without_feature_names_from_disk_and_warns():
    maui_model = Maui(n_hidden=[10], n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    with tempfile.TemporaryDirectory() as tmpdirname:
        maui_model.save(tmpdirname)
        os.remove(os.path.join(tmpdirname, "maui_feature_names.txt"))
        with pytest.warns(MauiWarning):
            maui_model_from_disk = Maui.load(tmpdirname)
        assert maui_model_from_disk.feature_names is None


def test_maui_can_fine_tune():
    maui_model = Maui(n_hidden=[], n_latent=2, epochs=1)
    maui_model = maui_model.fit({"d1": df1, "d2": df2})
    maui_model.fine_tune({"d1": df1, "d2": df2}, epochs=1)


def test_maui_complains_if_fine_tune_with_wrong_features():
    maui_model = Maui(n_hidden=[], n_latent=2, epochs=1)
    maui_model.fit({"d1": df1, "d2": df2})

    df1_wrong_features = df1.reindex(df1.index[: len(df1.index) - 1])
    with pytest.raises(ValueError):
        z = maui_model.fine_tune({"df1": df1_wrong_features, "df2": df2})
