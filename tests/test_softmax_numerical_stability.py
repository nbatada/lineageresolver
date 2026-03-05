from __future__ import annotations

import numpy as np

from lineageresolver.model import stable_softmax


def test_softmax_numerical_stability_large_logits() -> None:
    logits = np.array(
        [
            [10000.0, 10001.0, 9999.0],
            [-10000.0, -10001.0, -9999.0],
            [1e6, 1e6 + 1.0, 1e6 - 2.0],
        ]
    )
    probs = stable_softmax(logits)

    assert np.isfinite(probs).all()
    assert not np.isnan(probs).any()
    assert np.allclose(probs.sum(axis=1), 1.0, atol=1e-8)
