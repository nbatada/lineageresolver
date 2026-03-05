"""Model metrics utilities."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def expected_calibration_error(
    y_true: np.ndarray,
    probs: np.ndarray,
    class_labels: list[Any],
    n_bins: int = 10,
) -> tuple[float, pd.DataFrame]:
    """Compute ECE and per-bin reliability summary."""
    y_true = np.asarray(y_true)
    probs = np.asarray(probs, dtype=float)
    preds = np.asarray(class_labels, dtype=object)[np.argmax(probs, axis=1)]
    conf = np.max(probs, axis=1)
    correct = (preds == y_true).astype(float)

    bins = np.linspace(0.0, 1.0, n_bins + 1)
    rows = []
    ece = 0.0
    n = len(y_true)
    for i in range(n_bins):
        left = bins[i]
        right = bins[i + 1]
        if i == n_bins - 1:
            mask = (conf >= left) & (conf <= right)
        else:
            mask = (conf >= left) & (conf < right)
        count = int(mask.sum())
        if count == 0:
            rows.append({"bin_left": left, "bin_right": right, "count": 0, "accuracy": np.nan, "confidence": np.nan})
            continue
        bin_acc = float(correct[mask].mean())
        bin_conf = float(conf[mask].mean())
        ece += (count / n) * abs(bin_acc - bin_conf)
        rows.append(
            {
                "bin_left": left,
                "bin_right": right,
                "count": count,
                "accuracy": bin_acc,
                "confidence": bin_conf,
            }
        )

    return float(ece), pd.DataFrame(rows)


def brier_score_multiclass(
    y_true: np.ndarray,
    probs: np.ndarray,
    class_labels: list[Any],
) -> float:
    """Compute multiclass Brier score."""
    y_true = np.asarray(y_true)
    probs = np.asarray(probs, dtype=float)
    class_to_idx = {label: i for i, label in enumerate(class_labels)}
    one_hot = np.zeros_like(probs)
    for row, label in enumerate(y_true):
        one_hot[row, class_to_idx[label]] = 1.0
    return float(np.mean(np.sum((probs - one_hot) ** 2, axis=1)))


def error_coverage_curve(
    y_true: np.ndarray,
    probs: np.ndarray,
    class_labels: list[Any],
) -> pd.DataFrame:
    """Compute error-vs-coverage points from confidence-ranked predictions."""
    y_true = np.asarray(y_true)
    probs = np.asarray(probs, dtype=float)

    conf = np.max(probs, axis=1)
    pred_idx = np.argmax(probs, axis=1)
    preds = np.asarray(class_labels, dtype=object)[pred_idx]
    is_error = (preds != y_true).astype(float)

    order = np.argsort(-conf)
    sorted_errors = is_error[order]
    n = len(y_true)

    coverage = np.arange(1, n + 1) / n
    cum_error_rate = np.cumsum(sorted_errors) / np.arange(1, n + 1)
    return pd.DataFrame(
        {
            "coverage": coverage,
            "error_rate": cum_error_rate,
            "confidence_threshold": conf[order],
        }
    )
