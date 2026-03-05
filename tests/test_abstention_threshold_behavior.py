from __future__ import annotations

import numpy as np
import pandas as pd

from lineageresolver.model import infer_posteriors


def test_abstention_threshold_behavior() -> None:
    obs = pd.DataFrame(
        {
            "lineageresolver_module_A": [1.0, 0.5],
            "lineageresolver_module_B": [1.0, 0.5],
            "lineageresolver_tcr_junction_umis": [0.0, 0.0],
            "lineageresolver_amb_ll_tcr": [0.0, 0.0],
        }
    )
    beta = {
        "A": {"b0": 0.0, "b_module": 1.0, "b_kJ": 0.0, "b_ambll": 0.0},
        "B": {"b0": 0.0, "b_module": 1.0, "b_kJ": 0.0, "b_ambll": 0.0},
    }

    conservative = infer_posteriors(obs=obs, classes=["A", "B"], beta_by_class=beta, tau=0.9)
    permissive = infer_posteriors(obs=obs, classes=["A", "B"], beta_by_class=beta, tau=0.49)

    assert (conservative["label_call"] == "uncertain").all()
    assert (permissive["label_call"] != "uncertain").all()
    assert np.allclose(conservative["max_p"], 0.5)
