from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from lineageresolver.modules import compute_weighted_module_scores
from tests.helpers import MiniAnnData


def test_module_scoring_weighted_toy() -> None:
    x = sparse.csr_matrix(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
    adata = MiniAnnData(
        X=x,
        obs_names=pd.Index(["c1", "c2"]),
        var_names=pd.Index(["A", "B", "C"]),
        obs=pd.DataFrame(index=["c1", "c2"]),
    )
    modules_cfg = {
        "Class1": {
            "genes": ["A", "C"],
            "weights": [1.0, 2.0],
        }
    }

    scores, _ = compute_weighted_module_scores(adata, modules_cfg, classes=["Class1"])
    expected = np.array([(1.0 + 2.0 * 3.0) / 3.0, (4.0 + 2.0 * 6.0) / 3.0])
    np.testing.assert_allclose(scores["Class1"], expected)
