from __future__ import annotations

import math

from lineageresolver.api import infer
from tests.helpers import make_mini_adata


def test_placeholder_output_columns_and_diagnostics_written() -> None:
    adata = make_mini_adata(n_cells=4)
    adata.obs["is_candidate"] = [True, False, True, False]

    infer(adata=adata, task_config={"classes": ["NK", "gdT", "abT"]}, candidate_set="is_candidate")

    assert "lineageresolver_label_map" in adata.obs.columns
    assert "lineageresolver_label_call" in adata.obs.columns
    assert "lineageresolver_max_p" in adata.obs.columns

    assert adata.obs["lineageresolver_label_map"].tolist() == [
        "unresolved",
        "not_evaluated",
        "unresolved",
        "not_evaluated",
    ]
    assert adata.obs["lineageresolver_label_call"].tolist() == [
        "uncertain",
        "not_evaluated",
        "uncertain",
        "not_evaluated",
    ]
    max_p = adata.obs["lineageresolver_max_p"].tolist()
    assert max_p[0] == 0.0
    assert math.isnan(max_p[1])
    assert max_p[2] == 0.0
    assert math.isnan(max_p[3])

    assert adata.uns["lineageresolver"]["mode"] == "placeholder"
    assert adata.uns["lineageresolver"]["candidate_count"] == 2
