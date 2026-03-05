"""Ambient profile estimation utilities."""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from lineageresolver.io import write_ambient_artifacts

TR_MARKER_GENES = {"TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2"}


@dataclass
class AmbientEstimationResult:
    """Container for ambient estimation outputs."""

    profile: pd.DataFrame
    report: dict[str, Any]
    selected_barcodes: list[str]


def estimate_ambient_profile(
    raw_10x_path: str | Path,
    min_umi: int = 0,
    max_umi: int = 50,
    bottom_percent: float | None = None,
    output_dir: str | Path | None = None,
) -> AmbientEstimationResult:
    """Estimate ambient profile from unfiltered 10x-style data."""
    matrix, genes, barcodes = read_raw_10x_mtx(raw_10x_path)
    barcode_umis = np.asarray(matrix.sum(axis=0)).ravel()
    selected_mask = select_empty_droplets(
        barcode_umis=barcode_umis,
        min_umi=min_umi,
        max_umi=max_umi,
        bottom_percent=bottom_percent,
    )

    if not selected_mask.any():
        raise ValueError(
            "No droplets matched the ambient selection criteria. "
            "Adjust min_umi/max_umi or bottom_percent."
        )

    selected_barcodes = [barcodes[i] for i in np.where(selected_mask)[0]]
    ambient_counts = np.asarray(matrix[:, selected_mask].sum(axis=1)).ravel().astype(float)
    total_ambient = float(ambient_counts.sum())
    if total_ambient <= 0.0:
        raise ValueError("Selected droplets had zero total counts; cannot estimate ambient profile.")

    ambient_fraction = ambient_counts / total_ambient
    profile = pd.DataFrame(
        {
            "gene": genes,
            "ambient_count": ambient_counts,
            "ambient_fraction": ambient_fraction,
        }
    )

    report = build_ambient_report(
        profile=profile,
        barcode_umis=barcode_umis,
        selected_mask=selected_mask,
        min_umi=min_umi,
        max_umi=max_umi,
        bottom_percent=bottom_percent,
    )

    if output_dir is not None:
        write_ambient_artifacts(profile=profile, report=report, output_dir=output_dir)

    return AmbientEstimationResult(
        profile=profile,
        report=report,
        selected_barcodes=selected_barcodes,
    )


def estimate_fallback_ambient_profile(
    adata: Any,
    bottom_percent: float = 10.0,
    output_dir: str | Path | None = None,
) -> AmbientEstimationResult:
    """Estimate ambient profile from low-UMI cells in filtered AnnData."""
    if not (0.0 < bottom_percent <= 100.0):
        raise ValueError("bottom_percent must be in (0, 100].")

    n_umi = np.asarray(adata.X.sum(axis=1)).ravel()
    cutoff = np.percentile(n_umi, bottom_percent)
    selected_mask = n_umi <= cutoff
    if not selected_mask.any():
        raise ValueError("No cells matched fallback ambient selection criteria.")

    ambient_counts = np.asarray(adata.X[selected_mask, :].sum(axis=0)).ravel().astype(float)
    total_ambient = float(ambient_counts.sum())
    if total_ambient <= 0.0:
        raise ValueError("Fallback selected cells had zero counts.")

    ambient_fraction = ambient_counts / total_ambient
    genes = [str(g) for g in adata.var_names]
    selected_barcodes = [str(bc) for bc in np.asarray(adata.obs_names)[selected_mask]]

    profile = pd.DataFrame(
        {
            "gene": genes,
            "ambient_count": ambient_counts,
            "ambient_fraction": ambient_fraction,
        }
    )
    report = {
        "mode": "fallback",
        "n_barcodes_used": int(selected_mask.sum()),
        "total_umis_used": total_ambient,
        "umi_thresholds": {
            "bottom_percent": float(bottom_percent),
            "cutoff_n_umi": float(cutoff),
        },
        "tr_marker_ambient_fraction": float(
            profile[profile["gene"].isin(TR_MARKER_GENES)]["ambient_fraction"].sum()
        ),
    }

    if output_dir is not None:
        write_ambient_artifacts(profile=profile, report=report, output_dir=output_dir)

    return AmbientEstimationResult(
        profile=profile,
        report=report,
        selected_barcodes=selected_barcodes,
    )


def read_raw_10x_mtx(raw_10x_path: str | Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """Read 10x raw matrix from MTX directory layout."""
    matrix_dir = _resolve_matrix_dir(Path(raw_10x_path))
    matrix_path = matrix_dir / "matrix.mtx.gz"
    features_path = matrix_dir / "features.tsv.gz"
    barcodes_path = matrix_dir / "barcodes.tsv.gz"

    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr()

    genes = _read_features(features_path)
    barcodes = _read_lines(barcodes_path)

    if matrix.shape == (len(genes), len(barcodes)):
        return matrix, genes, barcodes
    if matrix.shape == (len(barcodes), len(genes)):
        return matrix.T.tocsr(), genes, barcodes

    raise ValueError(
        "Raw 10x matrix shape does not match features/barcodes lengths: "
        f"matrix={matrix.shape}, genes={len(genes)}, barcodes={len(barcodes)}"
    )


def select_empty_droplets(
    barcode_umis: np.ndarray,
    min_umi: int = 0,
    max_umi: int = 50,
    bottom_percent: float | None = None,
) -> np.ndarray:
    """Select candidate empty droplets by UMI thresholds or bottom percentile."""
    if bottom_percent is not None:
        if not (0.0 < bottom_percent <= 100.0):
            raise ValueError("bottom_percent must be in (0, 100].")
        cutoff = np.percentile(barcode_umis, bottom_percent)
        return barcode_umis <= cutoff

    if min_umi < 0 or max_umi < min_umi:
        raise ValueError("Require 0 <= min_umi <= max_umi.")
    return (barcode_umis >= min_umi) & (barcode_umis <= max_umi)


def build_ambient_report(
    profile: pd.DataFrame,
    barcode_umis: np.ndarray,
    selected_mask: np.ndarray,
    min_umi: int,
    max_umi: int,
    bottom_percent: float | None,
) -> dict[str, Any]:
    """Build diagnostics report for ambient estimation."""
    selected_umis = barcode_umis[selected_mask]
    top_genes = (
        profile.nlargest(10, "ambient_fraction")[["gene", "ambient_fraction"]].to_dict(orient="records")
    )
    tr_fraction = float(profile[profile["gene"].isin(TR_MARKER_GENES)]["ambient_fraction"].sum())

    return {
        "n_barcodes_used": int(selected_mask.sum()),
        "total_umis_used": float(selected_umis.sum()),
        "umi_thresholds": {
            "min_umi": int(min_umi),
            "max_umi": int(max_umi),
            "bottom_percent": bottom_percent,
        },
        "selected_umi_summary": {
            "min": float(selected_umis.min()),
            "median": float(np.median(selected_umis)),
            "max": float(selected_umis.max()),
        },
        "top_genes": top_genes,
        "tr_marker_ambient_fraction": tr_fraction,
    }


def _resolve_matrix_dir(raw_10x_path: Path) -> Path:
    if raw_10x_path.is_dir():
        if (raw_10x_path / "matrix.mtx.gz").exists():
            return raw_10x_path
        nested = raw_10x_path / "raw_feature_bc_matrix"
        if nested.is_dir():
            return nested
    if raw_10x_path.suffix == ".h5":
        raise NotImplementedError("10x H5 input is not implemented in v1.")
    raise FileNotFoundError(
        "Could not find raw 10x MTX files. Expected either:\n"
        "- <raw_10x_path>/matrix.mtx.gz\n"
        "- <raw_10x_path>/raw_feature_bc_matrix/matrix.mtx.gz"
    )


def _read_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def _read_features(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        genes: list[str] = []
        for line in handle:
            parts = line.strip().split("\t")
            if not parts:
                continue
            if len(parts) >= 2:
                genes.append(parts[1])
            else:
                genes.append(parts[0])
    return genes


def load_manifest(path: str | Path) -> dict[str, Any]:
    """Load JSON manifest helper for tests."""
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)
