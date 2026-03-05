"""Public API entrypoints."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from lineageresolver import ambient as ambient_mod
from lineageresolver.io import resolve_candidate_mask, validate_adata_for_infer


def infer(
    adata: Any,
    task_config: Any,
    raw_10x_path: str | None = None,
    ambient_profile: Any | None = None,
    evidence_table: str | None = None,
    candidate_set: Any | None = None,
    tau: float = 0.9,
    return_posteriors: bool = True,
    inplace: bool = True,
):
    """Run v1 infer() skeleton with candidate routing and placeholder outputs."""
    del raw_10x_path, ambient_profile, evidence_table, return_posteriors

    validate_adata_for_infer(adata)
    if not (0.0 < float(tau) <= 1.0):
        raise ValueError("tau must be in the interval (0, 1].")

    adata_out = adata if inplace else adata.copy()
    candidate_mask = resolve_candidate_mask(adata_out, candidate_set)

    label_map = np.full(len(adata_out.obs_names), "not_evaluated", dtype=object)
    label_call = np.full(len(adata_out.obs_names), "not_evaluated", dtype=object)
    max_p = np.full(len(adata_out.obs_names), np.nan, dtype=float)

    label_map[candidate_mask] = "unresolved"
    label_call[candidate_mask] = "uncertain"
    max_p[candidate_mask] = 0.0

    adata_out.obs["lineageresolver_label_map"] = label_map
    adata_out.obs["lineageresolver_label_call"] = label_call
    adata_out.obs["lineageresolver_max_p"] = max_p

    adata_out.uns["lineageresolver"] = {
        "mode": "placeholder",
        "candidate_count": int(candidate_mask.sum()),
        "n_obs": int(len(adata_out.obs_names)),
        "task_config_type": type(task_config).__name__,
        "task_config_source": str(task_config) if isinstance(task_config, (str, Path)) else "in-memory",
    }
    adata_out.uns["lineageresolver_candidate_mask"] = candidate_mask.tolist()

    return None if inplace else adata_out


def estimate_ambient_profile(*args: Any, **kwargs: Any):
    """Ambient estimation API placeholder."""
    return ambient_mod.estimate_ambient_profile(*args, **kwargs)
