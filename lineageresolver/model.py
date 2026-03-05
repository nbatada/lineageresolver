"""Model scoring and posterior inference utilities."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def stable_softmax(logits: np.ndarray) -> np.ndarray:
    """Compute row-wise softmax with log-sum-exp stabilization."""
    logits = np.asarray(logits, dtype=float)
    shifted = logits - np.max(logits, axis=1, keepdims=True)
    exp_values = np.exp(shifted)
    sums = np.sum(exp_values, axis=1, keepdims=True)
    return exp_values / sums


def compute_class_scores(
    obs: pd.DataFrame,
    classes: list[str],
    beta_by_class: dict[str, dict[str, float]],
) -> np.ndarray:
    """Compute linear class scores from module/evidence features."""
    n_obs = len(obs)
    logits = np.zeros((n_obs, len(classes)), dtype=float)
    k_j = _get_obs_column(obs, "lineageresolver_tcr_junction_umis")
    amb_ll = _get_obs_column(obs, "lineageresolver_amb_ll_tcr")
    qc = _get_obs_column(obs, "percent_mito")

    for idx, class_name in enumerate(classes):
        beta = beta_by_class.get(class_name, {})
        b0 = float(beta.get("b0", 0.0))
        b_module = float(beta.get("b_module", 1.0))
        b_kj = float(beta.get("b_kJ", 0.0))
        b_ambll = float(beta.get("b_ambll", 0.0))
        b_qc = float(beta.get("b_qc", 0.0))

        module = _get_obs_column(obs, f"lineageresolver_module_{class_name}")
        logits[:, idx] = b0 + b_module * module + b_kj * k_j - b_ambll * amb_ll + b_qc * qc

    return logits


def infer_posteriors(
    obs: pd.DataFrame,
    classes: list[str],
    beta_by_class: dict[str, dict[str, float]],
    tau: float,
) -> dict[str, Any]:
    """Compute class probabilities and abstention-aware calls."""
    logits = compute_class_scores(obs=obs, classes=classes, beta_by_class=beta_by_class)
    probs = stable_softmax(logits)
    max_idx = np.argmax(probs, axis=1)
    max_p = probs[np.arange(probs.shape[0]), max_idx]
    label_map = np.asarray(classes, dtype=object)[max_idx]
    label_call = np.where(max_p >= tau, label_map, "uncertain")
    entropy = -np.sum(probs * np.log(np.clip(probs, 1e-12, 1.0)), axis=1)

    posterior_df = pd.DataFrame(
        {f"lineageresolver_p_{cls}": probs[:, i] for i, cls in enumerate(classes)},
        index=obs.index,
    )
    return {
        "logits": logits,
        "posterior_df": posterior_df,
        "label_map": label_map,
        "label_call": label_call,
        "max_p": max_p,
        "entropy": entropy,
    }


def _get_obs_column(obs: pd.DataFrame, column: str) -> np.ndarray:
    if column not in obs.columns:
        return np.zeros(len(obs), dtype=float)
    return pd.to_numeric(obs[column], errors="coerce").fillna(0.0).to_numpy(dtype=float)
