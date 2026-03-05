from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from lineageresolver.api import infer


H5AD_PATH = Path("tests/data/filtered_small.h5ad")


def test_missing_evidence_defaults_do_not_crash_infer() -> None:
    adata = ad.read_h5ad(H5AD_PATH)

    with pytest.warns(UserWarning):
        out = infer(
            adata=adata,
            task_config={
                "classes": ["NK", "gdT", "abT"],
                "modules": {
                    "NK": {"genes": ["NKG7"]},
                    "gdT": {"genes": ["TRDC"]},
                    "abT": {"genes": ["TRAC"]},
                },
            },
            evidence_table=None,
            inplace=False,
        )

    assert "lineageresolver_tcr_junction_umis" in out.obs.columns
    assert "lineageresolver_lambda_amb_tcr" in out.obs.columns
    assert "lineageresolver_amb_ll_tcr" in out.obs.columns
    assert "lineageresolver_a_hat" in out.obs.columns

    assert np.allclose(out.obs["lineageresolver_tcr_junction_umis"].to_numpy(), 0.0)
    assert np.isfinite(out.obs["lineageresolver_lambda_amb_tcr"].to_numpy()).all()
    assert np.isfinite(out.obs["lineageresolver_amb_ll_tcr"].to_numpy()).all()
    assert out.uns["lineageresolver"]["evidence_diagnostics"]["evidence_table_missing"] is True
