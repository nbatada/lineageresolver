"""I/O helpers for LineageResolver."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import numpy as np
import pandas as pd

DEFAULT_CLUSTER_COLUMNS = ("cluster", "clusters", "leiden", "seurat_clusters")


def validate_adata_for_infer(adata: Any) -> None:
    """Validate that an AnnData-like object has fields required by infer()."""
    required_attrs = ("X", "obs", "obs_names", "var_names", "uns")
    for attr in required_attrs:
        if not hasattr(adata, attr):
            raise ValueError(f"adata is missing required attribute '{attr}'.")

    if len(adata.obs_names) == 0:
        raise ValueError("adata.obs_names must not be empty.")
    if len(adata.var_names) == 0:
        raise ValueError("adata.var_names must not be empty.")
    if getattr(adata, "obs", None) is None:
        raise ValueError("adata.obs must be a pandas DataFrame.")
    if not isinstance(adata.obs, pd.DataFrame):
        raise ValueError("adata.obs must be a pandas DataFrame.")
    if len(adata.obs) != len(adata.obs_names):
        raise ValueError("adata.obs length must match adata.obs_names length.")


def resolve_candidate_mask(adata: Any, candidate_set: Any | None) -> np.ndarray:
    """Resolve user-provided candidate_set input to a boolean mask."""
    n_obs = len(adata.obs_names)
    obs_names = pd.Index(adata.obs_names)

    if candidate_set is None:
        return np.ones(n_obs, dtype=bool)

    if isinstance(candidate_set, str):
        return _resolve_from_string(adata, candidate_set)

    if isinstance(candidate_set, np.ndarray):
        return _coerce_bool_vector(candidate_set, n_obs)

    if isinstance(candidate_set, pd.Series):
        series = candidate_set.reindex(obs_names) if candidate_set.index.is_unique else candidate_set
        return _coerce_bool_vector(series.to_numpy(), n_obs)

    if isinstance(candidate_set, Iterable) and not isinstance(candidate_set, (bytes, bytearray)):
        values = list(candidate_set)
        if not values:
            return np.zeros(n_obs, dtype=bool)

        as_str = [str(v) for v in values]
        if all(v in obs_names for v in as_str):
            return obs_names.isin(as_str)

        cluster_col = _find_cluster_column(adata.obs)
        if cluster_col is None:
            raise ValueError(
                "candidate_set list was not recognized as barcodes and no cluster column "
                f"was found among {DEFAULT_CLUSTER_COLUMNS}."
            )
        return adata.obs[cluster_col].astype(str).isin(as_str).to_numpy(dtype=bool)

    raise ValueError(
        "Unsupported candidate_set type. Expected None, str, bool vector, or iterable of identifiers."
    )


def _resolve_from_string(adata: Any, value: str) -> np.ndarray:
    obs_names = pd.Index(adata.obs_names)
    if value in adata.obs.columns:
        column = adata.obs[value]
        if pd.api.types.is_bool_dtype(column):
            return column.to_numpy(dtype=bool)
        unique = set(column.dropna().unique().tolist())
        if unique.issubset({0, 1}):
            return column.astype(bool).to_numpy()
        raise ValueError(
            f"candidate_set column '{value}' exists but is not boolean. "
            "Use a boolean column, a barcode list, or a cluster-id list."
        )

    if value in obs_names:
        return obs_names.isin([value])

    raise ValueError(
        f"candidate_set string '{value}' did not match a boolean obs column or an obs barcode."
    )


def _coerce_bool_vector(values: np.ndarray, n_obs: int) -> np.ndarray:
    arr = np.asarray(values)
    if arr.shape[0] != n_obs:
        raise ValueError(
            f"candidate_set boolean vector length {arr.shape[0]} does not match n_obs ({n_obs})."
        )
    if arr.dtype == bool:
        return arr

    unique = set(np.unique(arr).tolist())
    if unique.issubset({0, 1}):
        return arr.astype(bool)
    raise ValueError("candidate_set vector must be boolean (or integer 0/1).")


def _find_cluster_column(obs: pd.DataFrame) -> str | None:
    for column in DEFAULT_CLUSTER_COLUMNS:
        if column in obs.columns:
            return column
    return None
