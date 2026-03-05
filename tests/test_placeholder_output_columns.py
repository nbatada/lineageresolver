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

    assert adata.obs["lineageresolver_label_map"].tolist()[1] == "not_evaluated"
    assert adata.obs["lineageresolver_label_map"].tolist()[3] == "not_evaluated"
    assert adata.obs["lineageresolver_label_call"].tolist()[1] == "not_evaluated"
    assert adata.obs["lineageresolver_label_call"].tolist()[3] == "not_evaluated"

    assert adata.obs["lineageresolver_label_map"].tolist()[0] in {"NK", "gdT", "abT"}
    assert adata.obs["lineageresolver_label_map"].tolist()[2] in {"NK", "gdT", "abT"}
    assert adata.obs["lineageresolver_label_call"].tolist()[0] in {"uncertain", "NK", "gdT", "abT"}
    assert adata.obs["lineageresolver_label_call"].tolist()[2] in {"uncertain", "NK", "gdT", "abT"}
    max_p = adata.obs["lineageresolver_max_p"].tolist()
    assert 0.0 <= max_p[0] <= 1.0
    assert math.isnan(max_p[1])
    assert 0.0 <= max_p[2] <= 1.0
    assert math.isnan(max_p[3])

    assert adata.uns["lineageresolver"]["mode"] == "softmax_v1"
    assert adata.uns["lineageresolver"]["candidate_count"] == 2
