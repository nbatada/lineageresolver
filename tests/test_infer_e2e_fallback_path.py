from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from lineageresolver.api import infer


H5AD_PATH = Path("tests/data/filtered_small.h5ad")
EVIDENCE_PATH = Path("tests/data/evidence_small.tsv")
CONFIG_PATH = Path("tests/data/configs/cytotoxic_adjudication_v1.yaml")


def test_infer_e2e_fallback_path() -> None:
    adata = ad.read_h5ad(H5AD_PATH)

    with pytest.warns(UserWarning, match="fallback ambient"):
        out = infer(
            adata=adata,
            task_config=CONFIG_PATH,
            raw_10x_path=None,
            evidence_table=EVIDENCE_PATH,
            inplace=False,
        )

    probs = out.obs[["lineageresolver_p_NK", "lineageresolver_p_gdT", "lineageresolver_p_abT"]].to_numpy()
    assert np.isfinite(probs).all()
    assert np.allclose(probs.sum(axis=1), 1.0, atol=1e-6)
    assert out.uns["lineageresolver_ambient_estimation_mode"] == "fallback"
