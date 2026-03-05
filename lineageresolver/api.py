"""Public API entrypoints."""

from __future__ import annotations

from pathlib import Path
from typing import Any
import warnings

import numpy as np
import pandas as pd

from lineageresolver import ambient as ambient_mod
from lineageresolver.config import load_task_config
from lineageresolver.evidence import derive_tcr_features
from lineageresolver.io import resolve_candidate_mask, validate_adata_for_infer
from lineageresolver.modules import compute_weighted_module_scores


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
    del return_posteriors

    validate_adata_for_infer(adata)
    if not (0.0 < float(tau) <= 1.0):
        raise ValueError("tau must be in the interval (0, 1].")

    loaded_config = load_task_config(task_config)
    classes = [str(c) for c in loaded_config.get("classes", [])]
    modules_cfg = loaded_config.get("modules", {})

    adata_out = adata if inplace else adata.copy()
    candidate_mask = resolve_candidate_mask(adata_out, candidate_set)

    ambient_result = None
    ambient_mode: str
    if ambient_profile is not None:
        ambient_mode = "provided"
    elif raw_10x_path is not None:
        ambient_result = ambient_mod.estimate_ambient_profile(raw_10x_path)
        ambient_mode = "raw10x"
    else:
        warnings.warn(
            "raw_10x_path was not provided; using fallback ambient estimation from filtered cells.",
            UserWarning,
            stacklevel=2,
        )
        ambient_result = ambient_mod.estimate_fallback_ambient_profile(adata_out)
        ambient_mode = "fallback"

    label_map = np.full(len(adata_out.obs_names), "not_evaluated", dtype=object)
    label_call = np.full(len(adata_out.obs_names), "not_evaluated", dtype=object)
    max_p = np.full(len(adata_out.obs_names), np.nan, dtype=float)

    label_map[candidate_mask] = "unresolved"
    label_call[candidate_mask] = "uncertain"
    max_p[candidate_mask] = 0.0

    adata_out.obs["lineageresolver_label_map"] = label_map
    adata_out.obs["lineageresolver_label_call"] = label_call
    adata_out.obs["lineageresolver_max_p"] = max_p

    ambient_profile_df = _resolve_ambient_profile(ambient_profile, ambient_result)
    tcr_features, evidence_diagnostics = derive_tcr_features(
        adata=adata_out,
        evidence_table=evidence_table,
        ambient_profile=ambient_profile_df,
        tcr_marker_genes=loaded_config.get("tcr_marker_set_genes_for_alpha_tcr"),
    )
    for column in tcr_features.columns:
        adata_out.obs[column] = tcr_features[column].to_numpy()

    module_scores, module_diagnostics = compute_weighted_module_scores(
        adata_out,
        modules_cfg=modules_cfg,
        classes=classes if classes else list(modules_cfg.keys()),
    )
    for class_name, values in module_scores.items():
        adata_out.obs[f"lineageresolver_module_{class_name}"] = values

    adata_out.uns["lineageresolver_ambient_estimation_mode"] = ambient_mode
    if ambient_result is not None:
        adata_out.uns["lineageresolver_ambient_report"] = ambient_result.report

    adata_out.uns["lineageresolver"] = {
        "mode": "placeholder",
        "candidate_count": int(candidate_mask.sum()),
        "n_obs": int(len(adata_out.obs_names)),
        "task_config_type": type(task_config).__name__,
        "task_config_source": str(task_config) if isinstance(task_config, (str, Path)) else "in-memory",
        "ambient_mode": ambient_mode,
        "module_diagnostics": module_diagnostics,
        "evidence_diagnostics": evidence_diagnostics,
    }
    adata_out.uns["lineageresolver_candidate_mask"] = candidate_mask.tolist()

    return None if inplace else adata_out


def estimate_ambient_profile(*args: Any, **kwargs: Any):
    """Ambient estimation API placeholder."""
    return ambient_mod.estimate_ambient_profile(*args, **kwargs)


def _resolve_ambient_profile(
    ambient_profile: Any | None,
    ambient_result: Any | None,
) -> pd.DataFrame | None:
    if ambient_profile is not None:
        if isinstance(ambient_profile, pd.DataFrame):
            return ambient_profile
        if isinstance(ambient_profile, (str, Path)):
            return pd.read_csv(ambient_profile, sep="\t")
        raise ValueError("ambient_profile must be a DataFrame or path to TSV when provided.")

    if ambient_result is not None and hasattr(ambient_result, "profile"):
        return ambient_result.profile
    return None
