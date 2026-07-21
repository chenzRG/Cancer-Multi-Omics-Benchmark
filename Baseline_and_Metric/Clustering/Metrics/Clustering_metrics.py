"""Clustering metrics for MLOmics (SIL / LPS)."""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence, Union

import numpy as np
from lifelines.statistics import multivariate_logrank_test
from sklearn.metrics import silhouette_score

ArrayLike = Union[Sequence[Any], np.ndarray]


def clustering_metrics(
    y_pred: ArrayLike,
    X: Optional[ArrayLike] = None,
    survival_times: Optional[ArrayLike] = None,
    event_observed: Optional[ArrayLike] = None,
) -> Dict[str, Any]:
    """Return silhouette (SIL) and log-rank p-value (LPS) when applicable."""
    metrics: Dict[str, Any] = {}
    y_pred_arr = np.asarray(y_pred)
    n_groups = len(np.unique(y_pred_arr))

    if X is not None and n_groups > 1:
        metrics["sil"] = float(silhouette_score(np.asarray(X), y_pred_arr))

    if (
        survival_times is not None
        and event_observed is not None
        and n_groups > 1
    ):
        results = multivariate_logrank_test(
            survival_times, y_pred_arr, event_observed
        )
        metrics["lps"] = float(results.p_value)
    elif survival_times is not None and event_observed is not None:
        metrics["lps"] = "N/A"

    return metrics


if __name__ == "__main__":
    y_pred = [0, 1, 0, 1, 2, 2, 1, 0, 2, 1]
    X = np.random.default_rng(0).normal(size=(10, 2))
    survival_times = [10, 15, 10, 20, 30, 25, 35, 10, 40, 30]
    event_observed = [1, 0, 1, 1, 0, 1, 0, 1, 0, 1]
    print(clustering_metrics(y_pred, X, survival_times, event_observed))
