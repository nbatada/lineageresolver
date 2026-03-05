from __future__ import annotations

import pytest

from lineageresolver.api import infer
from tests.helpers import make_mini_adata


def test_infer_signature_runs_and_returns_copy_when_not_inplace() -> None:
    adata = make_mini_adata()
    out = infer(adata=adata, task_config={"classes": ["A", "B"]}, inplace=False)

    assert out is not adata
    assert "lineageresolver_label_map" in out.obs.columns
    assert "lineageresolver_label_call" in out.obs.columns
    assert "lineageresolver_max_p" in out.obs.columns
    assert "lineageresolver_label_map" not in adata.obs.columns


def test_infer_raises_clear_error_for_missing_required_adata_fields() -> None:
    class InvalidAdata:
        pass

    with pytest.raises(ValueError, match="missing required attribute"):
        infer(adata=InvalidAdata(), task_config={"classes": ["A", "B"]})


def test_infer_rejects_invalid_tau() -> None:
    adata = make_mini_adata()
    with pytest.raises(ValueError, match="tau"):
        infer(adata=adata, task_config={"classes": ["A", "B"]}, tau=0.0)
