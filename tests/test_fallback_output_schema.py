from __future__ import annotations

from pathlib import Path

import anndata as ad

from lineageresolver.ambient import estimate_fallback_ambient_profile


H5AD_PATH = Path("tests/data/filtered_small.h5ad")


def test_fallback_output_schema_and_normalization() -> None:
    adata = ad.read_h5ad(H5AD_PATH)
    result = estimate_fallback_ambient_profile(adata)

    assert set(["gene", "ambient_count", "ambient_fraction"]).issubset(result.profile.columns)
    assert len(result.profile) == adata.n_vars
    assert abs(float(result.profile["ambient_fraction"].sum()) - 1.0) < 1e-8
    assert (result.profile["ambient_count"] >= 0).all()
    assert result.report["mode"] == "fallback"
    assert result.report["n_barcodes_used"] > 0
