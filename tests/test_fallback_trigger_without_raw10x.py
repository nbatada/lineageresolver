from __future__ import annotations

from pathlib import Path

import anndata as ad
import pytest

from lineageresolver.api import infer


H5AD_PATH = Path("tests/data/filtered_small.h5ad")


def test_fallback_trigger_without_raw10x() -> None:
    adata = ad.read_h5ad(H5AD_PATH)

    with pytest.warns(UserWarning, match="fallback ambient"):
        out = infer(
            adata=adata,
            task_config={"classes": ["NK", "gdT", "abT"]},
            raw_10x_path=None,
            inplace=False,
        )

    assert out.uns["lineageresolver_ambient_estimation_mode"] == "fallback"
    assert "lineageresolver_ambient_report" in out.uns
