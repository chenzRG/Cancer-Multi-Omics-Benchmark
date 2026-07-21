"""Classification metrics for MLOmics (PREC / NMI / ARI)."""
from __future__ import annotations

from typing import Any, Dict, Sequence, Union

import numpy as np
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    precision_score,
)

ArrayLike = Union[Sequence[Any], np.ndarray]


def classification_metrics(y_true: ArrayLike, y_pred: ArrayLike) -> Dict[str, float]:
    """Return weighted precision, NMI, and ARI."""
    return {
        "precision": float(
            precision_score(y_true, y_pred, average="weighted", zero_division=0)
        ),
        "nmi": float(normalized_mutual_info_score(y_true, y_pred)),
        "ari": float(adjusted_rand_score(y_true, y_pred)),
    }


if __name__ == "__main__":
    y_true = [0, 1, 1, 0, 1, 0, 1]
    y_pred = [0, 1, 0, 0, 1, 1, 1]
    metrics = classification_metrics(y_true, y_pred)
    print(f"Precision: {metrics['precision']}")
    print(f"NMI: {metrics['nmi']}")
    print(f"ARI: {metrics['ari']}")
