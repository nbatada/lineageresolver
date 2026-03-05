from __future__ import annotations

import numpy as np
import pandas as pd

from lineageresolver.model import infer_posteriors


def test_softmax_probability_normalization() -> None:
    obs = pd.DataFrame(
        {
            "lineageresolver_module_A": [0.2, 1.5, -0.1],
            "lineageresolver_module_B": [1.0, -0.2, 0.3],
            "lineageresolver_tcr_junction_umis": [0.0, 2.0, 1.0],
            "lineageresolver_amb_ll_tcr": [0.1, 1.2, 0.5],
        }
    )
    beta = {
        "A": {"b0": 0.0, "b_module": 1.0, "b_kJ": 0.3, "b_ambll": 0.2},
        "B": {"b0": 0.1, "b_module": 0.8, "b_kJ": -0.2, "b_ambll": 0.1},
    }

    result = infer_posteriors(obs=obs, classes=["A", "B"], beta_by_class=beta, tau=0.8)
    probs = result["posterior_df"].to_numpy()

    assert np.allclose(probs.sum(axis=1), 1.0, atol=1e-8)
    assert np.isfinite(probs).all()
