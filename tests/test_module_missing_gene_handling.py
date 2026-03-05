from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from lineageresolver.modules import compute_weighted_module_scores
from tests.helpers import MiniAnnData


def test_module_missing_gene_handling() -> None:
    x = sparse.csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]]))
    adata = MiniAnnData(
        X=x,
        obs_names=pd.Index(["c1", "c2"]),
        var_names=pd.Index(["A", "B"]),
        obs=pd.DataFrame(index=["c1", "c2"]),
    )
    modules_cfg = {
        "has_missing": {"genes": ["A", "MISSING"], "weights": [1.0, 1.0]},
        "all_missing": {"genes": ["X", "Y"], "weights": [1.0, 1.0]},
    }

    scores, diagnostics = compute_weighted_module_scores(
        adata,
        modules_cfg,
        classes=["has_missing", "all_missing"],
    )

    np.testing.assert_allclose(scores["has_missing"], np.array([1.0, 3.0]))
    np.testing.assert_allclose(scores["all_missing"], np.array([0.0, 0.0]))
    assert diagnostics["missing_genes"]["has_missing"] == ["MISSING"]
    assert diagnostics["missing_genes"]["all_missing"] == ["X", "Y"]
    assert "all_missing" in diagnostics["classes_with_zero_features"]
