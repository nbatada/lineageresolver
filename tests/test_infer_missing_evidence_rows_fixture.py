from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np

from lineageresolver.api import infer


H5AD_PATH = Path("tests/data/filtered_small.h5ad")
EVIDENCE_PATH = Path("tests/data/evidence_small.tsv")
RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")
CONFIG_PATH = Path("tests/data/configs/cytotoxic_adjudication_v1.yaml")
MANIFEST_PATH = Path("tests/data/raw10x_small/manifest.json")


def test_infer_missing_evidence_rows_fixture() -> None:
    with open(MANIFEST_PATH, "r", encoding="utf-8") as handle:
        missing_rows = json.load(handle)["missing_evidence_barcodes"]

    adata = ad.read_h5ad(H5AD_PATH)
    out = infer(
        adata=adata,
        task_config=CONFIG_PATH,
        raw_10x_path=RAW10X_PATH,
        evidence_table=EVIDENCE_PATH,
        inplace=False,
    )

    missing_kj = out.obs.loc[missing_rows, "lineageresolver_tcr_junction_umis"].to_numpy()
    assert np.allclose(missing_kj, 0.0)
    assert np.isfinite(out.obs.loc[missing_rows, "lineageresolver_lambda_amb_tcr"].to_numpy()).all()
    assert np.isfinite(out.obs.loc[missing_rows, "lineageresolver_amb_ll_tcr"].to_numpy()).all()
