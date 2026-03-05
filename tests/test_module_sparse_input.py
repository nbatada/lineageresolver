from __future__ import annotations

import numpy as np

from lineageresolver.modules import compute_weighted_module_scores
from tests.helpers import make_mini_adata


def test_module_sparse_input_returns_finite_scores() -> None:
    adata = make_mini_adata(n_cells=10, n_genes=6)
    adata.var_names = adata.var_names.map(lambda g: g.replace("GENE", "G"))
    modules_cfg = {
        "k1": {"genes": ["G_001", "G_003"], "weights": [1.0, 0.5]},
        "k2": {"genes": ["G_002", "G_004"]},
    }

    scores, _ = compute_weighted_module_scores(adata, modules_cfg)

    assert set(scores.keys()) == {"k1", "k2"}
    assert scores["k1"].shape == (adata.n_obs,)
    assert scores["k2"].shape == (adata.n_obs,)
    assert np.isfinite(scores["k1"]).all()
    assert np.isfinite(scores["k2"]).all()
