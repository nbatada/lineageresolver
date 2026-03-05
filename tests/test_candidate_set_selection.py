from __future__ import annotations

import numpy as np
import pytest

from lineageresolver.io import resolve_candidate_mask
from tests.helpers import make_mini_adata


def test_candidate_set_accepts_boolean_obs_column() -> None:
    adata = make_mini_adata(n_cells=5)
    adata.obs["is_candidate"] = [True, False, True, False, True]

    mask = resolve_candidate_mask(adata, "is_candidate")
    assert mask.tolist() == [True, False, True, False, True]


def test_candidate_set_accepts_barcode_list() -> None:
    adata = make_mini_adata(n_cells=5)
    mask = resolve_candidate_mask(adata, ["CELL_0001-1", "CELL_0004-1"])
    assert mask.tolist() == [True, False, False, True, False]


def test_candidate_set_accepts_cluster_id_list() -> None:
    adata = make_mini_adata(n_cells=5)
    adata.obs["leiden"] = ["0", "1", "2", "0", "2"]

    mask = resolve_candidate_mask(adata, ["2"])
    assert mask.tolist() == [False, False, True, False, True]


def test_candidate_set_rejects_mismatched_boolean_vector_length() -> None:
    adata = make_mini_adata(n_cells=5)
    with pytest.raises(ValueError, match="length"):
        resolve_candidate_mask(adata, np.array([True, False]))
