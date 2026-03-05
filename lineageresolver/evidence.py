"""Evidence ingestion and feature derivation utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.special import gammaln

from lineageresolver.ambient import TR_MARKER_GENES

EVIDENCE_REQUIRED_KEY = "barcode"


def read_evidence_table(evidence_table: str | Path | pd.DataFrame | None) -> pd.DataFrame | None:
    """Read evidence table from TSV or DataFrame."""
    if evidence_table is None:
        return None
    if isinstance(evidence_table, pd.DataFrame):
        return evidence_table.copy()
    path = Path(evidence_table)
    if not path.exists():
        raise FileNotFoundError(f"Evidence table not found: {path}")
    return pd.read_csv(path, sep="\t")


def poisson_logpmf(k: np.ndarray, lam: np.ndarray) -> np.ndarray:
    """Numerically stable Poisson log PMF."""
    k_arr = np.asarray(k, dtype=float)
    lam_arr = np.asarray(lam, dtype=float)

    safe_lam = np.where(lam_arr <= 0.0, 1e-12, lam_arr)
    logpmf = k_arr * np.log(safe_lam) - safe_lam - gammaln(k_arr + 1.0)
    zero_zero = (lam_arr <= 0.0) & (k_arr == 0.0)
    logpmf = np.where(zero_zero, 0.0, logpmf)
    return logpmf


def poisson_pmf(k: np.ndarray, lam: np.ndarray) -> np.ndarray:
    """Poisson PMF from stable log-PMF."""
    return np.exp(poisson_logpmf(k, lam))


def amb_ll_tcr(k_j: np.ndarray, lambda_amb_tcr: np.ndarray) -> np.ndarray:
    """Ambient negative log-likelihood for TCR evidence."""
    return -poisson_logpmf(k_j, lambda_amb_tcr)


def derive_tcr_features(
    adata: Any,
    evidence_table: str | Path | pd.DataFrame | None,
    ambient_profile: pd.DataFrame | None,
    tcr_marker_genes: list[str] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Align evidence table and derive model-ready TCR features."""
    evidence_df = read_evidence_table(evidence_table)
    obs_names = pd.Index(adata.obs_names.astype(str))
    aligned = _align_evidence_to_obs(evidence_df, obs_names)

    trg = _get_column_or_zeros(aligned, "tcr_trg_junction_umis", len(obs_names))
    trd = _get_column_or_zeros(aligned, "tcr_trd_junction_umis", len(obs_names))
    total = _get_column_or_zeros(aligned, "tcr_total_junction_umis", len(obs_names))
    k_j = np.where(total > 0.0, total, trg + trd)

    n_umi = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
    a_hat = _estimate_a_hat(adata, n_umi)

    alpha_tcr = _derive_alpha_tcr(ambient_profile, tcr_marker_genes)
    lambda_amb = alpha_tcr * n_umi * a_hat
    amb_ll = amb_ll_tcr(k_j, lambda_amb)

    features = pd.DataFrame(
        {
            "lineageresolver_tcr_junction_umis": k_j,
            "lineageresolver_lambda_amb_tcr": lambda_amb,
            "lineageresolver_amb_ll_tcr": amb_ll,
            "lineageresolver_a_hat": a_hat,
        },
        index=obs_names,
    )

    diagnostics = {
        "alpha_tcr": float(alpha_tcr),
        "evidence_rows_input": int(0 if evidence_df is None else len(evidence_df)),
        "evidence_rows_aligned": int(aligned.notna().any(axis=1).sum()) if aligned is not None else 0,
        "evidence_missing_rows": int(0 if aligned is None else aligned.isna().all(axis=1).sum()),
        "evidence_table_missing": evidence_df is None,
    }
    return features, diagnostics


def _align_evidence_to_obs(
    evidence_df: pd.DataFrame | None,
    obs_names: pd.Index,
) -> pd.DataFrame | None:
    if evidence_df is None:
        return None
    if EVIDENCE_REQUIRED_KEY in evidence_df.columns:
        aligned = evidence_df.copy().set_index(EVIDENCE_REQUIRED_KEY)
    else:
        aligned = evidence_df.copy()
        aligned.index = aligned.index.astype(str)

    aligned.index = aligned.index.astype(str)
    return aligned.reindex(obs_names)


def _get_column_or_zeros(df: pd.DataFrame | None, column: str, n_obs: int) -> np.ndarray:
    if df is None or column not in df.columns:
        return np.zeros(n_obs, dtype=float)
    values = pd.to_numeric(df[column], errors="coerce").fillna(0.0).to_numpy(dtype=float)
    return np.clip(values, a_min=0.0, a_max=None)


def _estimate_a_hat(adata: Any, n_umi: np.ndarray) -> np.ndarray:
    if "lineageresolver_a_hat" in adata.obs.columns:
        values = pd.to_numeric(adata.obs["lineageresolver_a_hat"], errors="coerce").fillna(0.0).to_numpy()
        return np.clip(values.astype(float), 0.0, 1.0)
    if "synthetic_ambient_fraction" in adata.obs.columns:
        values = pd.to_numeric(adata.obs["synthetic_ambient_fraction"], errors="coerce").fillna(0.0).to_numpy()
        return np.clip(values.astype(float), 0.0, 1.0)

    # Conservative proxy: lower-depth cells are expected to carry higher ambient proportion.
    proxy = 50.0 / (n_umi + 50.0)
    return np.clip(proxy, 0.02, 0.2)


def _derive_alpha_tcr(
    ambient_profile: pd.DataFrame | None,
    tcr_marker_genes: list[str] | None,
) -> float:
    if ambient_profile is None:
        return 0.0
    marker_set = set(tcr_marker_genes or TR_MARKER_GENES)
    if "gene" not in ambient_profile.columns or "ambient_fraction" not in ambient_profile.columns:
        return 0.0
    mask = ambient_profile["gene"].astype(str).isin(marker_set)
    alpha = float(ambient_profile.loc[mask, "ambient_fraction"].sum())
    return max(alpha, 0.0)
