from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np

from lineageresolver.api import infer


H5AD_PATH = Path("tests/data/filtered_small.h5ad")
EVIDENCE_PATH = Path("tests/data/evidence_small.tsv")
RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")
CONFIG_PATH = Path("tests/data/configs/cytotoxic_adjudication_v1.yaml")


def test_infer_candidate_subset_alignment() -> None:
    adata = ad.read_h5ad(H5AD_PATH)
    adata.obs["is_candidate"] = [i % 2 == 0 for i in range(adata.n_obs)]

    out = infer(
        adata=adata,
        task_config=CONFIG_PATH,
        raw_10x_path=RAW10X_PATH,
        evidence_table=EVIDENCE_PATH,
        candidate_set="is_candidate",
        inplace=False,
    )

    candidate = out.obs["is_candidate"].to_numpy(dtype=bool)
    non_candidate = ~candidate

    assert np.isfinite(out.obs.loc[candidate, "lineageresolver_max_p"].to_numpy()).all()
    assert out.obs.loc[non_candidate, "lineageresolver_label_map"].eq("not_evaluated").all()
    assert out.obs.loc[non_candidate, "lineageresolver_label_call"].eq("not_evaluated").all()

    prob_cols = ["lineageresolver_p_NK", "lineageresolver_p_gdT", "lineageresolver_p_abT"]
    assert np.isnan(out.obs.loc[non_candidate, prob_cols].to_numpy()).all()
    assert np.isfinite(out.obs.loc[candidate, prob_cols].to_numpy()).all()
