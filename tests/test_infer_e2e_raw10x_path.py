from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np

from lineageresolver.api import infer


H5AD_PATH = Path("tests/data/filtered_small.h5ad")
EVIDENCE_PATH = Path("tests/data/evidence_small.tsv")
RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")
CONFIG_PATH = Path("tests/data/configs/cytotoxic_adjudication_v1.yaml")


def test_infer_e2e_raw10x_path() -> None:
    adata = ad.read_h5ad(H5AD_PATH)
    out = infer(
        adata=adata,
        task_config=CONFIG_PATH,
        raw_10x_path=RAW10X_PATH,
        evidence_table=EVIDENCE_PATH,
        inplace=False,
    )

    required = [
        "lineageresolver_label_map",
        "lineageresolver_label_call",
        "lineageresolver_max_p",
        "lineageresolver_uncertainty_entropy",
        "lineageresolver_tcr_junction_umis",
        "lineageresolver_lambda_amb_tcr",
        "lineageresolver_amb_ll_tcr",
        "lineageresolver_a_hat",
        "lineageresolver_p_NK",
        "lineageresolver_p_gdT",
        "lineageresolver_p_abT",
    ]
    for col in required:
        assert col in out.obs.columns

    probs = out.obs[["lineageresolver_p_NK", "lineageresolver_p_gdT", "lineageresolver_p_abT"]].to_numpy()
    assert np.isfinite(probs).all()
    assert np.allclose(probs.sum(axis=1), 1.0, atol=1e-6)

    assert out.uns["lineageresolver_ambient_estimation_mode"] == "raw10x"
    assert "lineageresolver_ambient_report" in out.uns
    assert "lineageresolver_task_config" in out.uns
