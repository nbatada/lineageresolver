from __future__ import annotations

import numpy as np
import pandas as pd

from lineageresolver.evidence import amb_ll_tcr, derive_tcr_features, poisson_pmf
from tests.helpers import make_mini_adata


def test_tcr_feature_derivation_matches_expected_calculation() -> None:
    adata = make_mini_adata(n_cells=3, n_genes=4)
    adata.var_names = pd.Index(["TRAC", "GENE1", "GENE2", "TRBC1"])
    adata.obs["synthetic_ambient_fraction"] = [0.1, 0.2, 0.05]

    evidence = pd.DataFrame(
        {
            "barcode": ["CELL_0001-1", "CELL_0002-1", "CELL_0003-1"],
            "tcr_total_junction_umis": [1, 4, 0],
        }
    )
    ambient = pd.DataFrame(
        {
            "gene": ["TRAC", "TRBC1", "GENE1"],
            "ambient_fraction": [0.03, 0.02, 0.95],
            "ambient_count": [30, 20, 950],
        }
    )

    features, diagnostics = derive_tcr_features(
        adata=adata,
        evidence_table=evidence,
        ambient_profile=ambient,
        tcr_marker_genes=["TRAC", "TRBC1"],
    )

    alpha_tcr = 0.05
    n_umi = np.asarray(adata.X.sum(axis=1)).ravel()
    expected_lambda = alpha_tcr * n_umi * np.array([0.1, 0.2, 0.05])

    np.testing.assert_allclose(features["lineageresolver_tcr_junction_umis"].to_numpy(), [1, 4, 0])
    np.testing.assert_allclose(features["lineageresolver_lambda_amb_tcr"].to_numpy(), expected_lambda)
    assert diagnostics["alpha_tcr"] == alpha_tcr


def test_poisson_and_amb_ll_behavior_stable_and_monotonic() -> None:
    # For k values above lambda, -log PMF should increase with k.
    k = np.array([2, 5, 8], dtype=float)
    lam = np.array([1.5, 1.5, 1.5], dtype=float)
    amb = amb_ll_tcr(k, lam)

    assert amb[0] < amb[1] < amb[2]

    pmf_small = poisson_pmf(np.array([0, 1]), np.array([1e-12, 1e-12]))
    pmf_large = poisson_pmf(np.array([10000, 10050]), np.array([10000.0, 10000.0]))
    assert np.isfinite(pmf_small).all()
    assert np.isfinite(pmf_large).all()
