from __future__ import annotations

import numpy as np

from lineageresolver.metrics import brier_score_multiclass, expected_calibration_error


def test_metrics_ece_brier() -> None:
    y_true = np.array(["A", "B", "A", "B"])
    probs = np.array(
        [
            [0.9, 0.1],
            [0.2, 0.8],
            [0.6, 0.4],
            [0.7, 0.3],
        ]
    )
    labels = ["A", "B"]

    ece, rel = expected_calibration_error(y_true=y_true, probs=probs, class_labels=labels, n_bins=4)
    brier = brier_score_multiclass(y_true=y_true, probs=probs, class_labels=labels)

    assert ece >= 0.0
    assert np.isfinite(ece)
    assert np.isfinite(brier)
    assert brier >= 0.0
    assert {"bin_left", "bin_right", "count", "accuracy", "confidence"}.issubset(rel.columns)
