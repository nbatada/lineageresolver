from __future__ import annotations

import numpy as np

from lineageresolver.metrics import error_coverage_curve


def test_metrics_error_coverage() -> None:
    y_true = np.array(["A", "B", "A", "B", "A"])
    probs = np.array(
        [
            [0.95, 0.05],
            [0.20, 0.80],
            [0.40, 0.60],
            [0.30, 0.70],
            [0.55, 0.45],
        ]
    )
    labels = ["A", "B"]

    curve = error_coverage_curve(y_true=y_true, probs=probs, class_labels=labels)

    assert list(curve.columns) == ["coverage", "error_rate", "confidence_threshold"]
    assert np.all(np.diff(curve["coverage"].to_numpy()) > 0)
    assert curve["coverage"].iloc[-1] == 1.0
    assert np.isfinite(curve["error_rate"]).all()
