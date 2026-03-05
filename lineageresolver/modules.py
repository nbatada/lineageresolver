"""Module score feature utilities."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import sparse


def compute_weighted_module_scores(
    adata: Any,
    modules_cfg: dict[str, dict[str, Any]],
    classes: list[str] | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Compute sparse-safe weighted module scores for each class."""
    if classes is None:
        classes = list(modules_cfg.keys())

    var_index = {str(gene): idx for idx, gene in enumerate(adata.var_names)}
    scores: dict[str, np.ndarray] = {}
    diagnostics: dict[str, Any] = {"missing_genes": {}, "classes_with_zero_features": []}

    for class_name in classes:
        module = modules_cfg.get(class_name, {})
        genes = [str(g) for g in module.get("genes", [])]
        weights_raw = module.get("weights")
        if weights_raw is None:
            weights_raw = [1.0] * len(genes)
        if len(weights_raw) != len(genes):
            raise ValueError(
                f"Module weights length mismatch for class '{class_name}': "
                f"{len(weights_raw)} weights for {len(genes)} genes."
            )

        present_indices: list[int] = []
        present_weights: list[float] = []
        missing_genes: list[str] = []

        for gene, weight in zip(genes, weights_raw):
            idx = var_index.get(gene)
            if idx is None:
                missing_genes.append(gene)
                continue
            present_indices.append(idx)
            present_weights.append(float(weight))

        diagnostics["missing_genes"][class_name] = missing_genes

        if not present_indices:
            diagnostics["classes_with_zero_features"].append(class_name)
            scores[class_name] = np.zeros(len(adata.obs_names), dtype=float)
            continue

        x_sub = adata.X[:, present_indices]
        weight_vec = np.asarray(present_weights, dtype=float)
        denominator = float(np.abs(weight_vec).sum())
        if denominator == 0.0:
            denominator = float(len(weight_vec))

        if sparse.issparse(x_sub):
            class_score = np.asarray(x_sub @ weight_vec).ravel()
        else:
            class_score = np.asarray(x_sub).dot(weight_vec)

        scores[class_name] = class_score / denominator

    return scores, diagnostics
