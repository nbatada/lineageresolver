from __future__ import annotations

from pathlib import Path

import anndata as ad
import pytest

from lineageresolver.api import infer


H5AD_PATH = Path("tests/data/filtered_small.h5ad")


def test_fallback_mode_flag_written() -> None:
    adata = ad.read_h5ad(H5AD_PATH)
    with pytest.warns(UserWarning):
        infer(adata=adata, task_config={"classes": ["A", "B"]}, inplace=True)

    assert adata.uns["lineageresolver_ambient_estimation_mode"] == "fallback"
    assert adata.uns["lineageresolver"]["ambient_mode"] == "fallback"
