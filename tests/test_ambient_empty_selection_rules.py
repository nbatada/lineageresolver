from __future__ import annotations

from pathlib import Path

from lineageresolver.ambient import estimate_ambient_profile, load_manifest


RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")
MANIFEST_PATH = Path("tests/data/raw10x_small/manifest.json")


def test_ambient_empty_selection_rules_overlap_manifest_empties() -> None:
    result = estimate_ambient_profile(RAW10X_PATH)
    manifest = load_manifest(MANIFEST_PATH)

    selected = set(result.selected_barcodes)
    empties = set(manifest["empty_barcodes"])

    overlap = selected & empties
    precision = len(overlap) / max(len(selected), 1)
    recall = len(overlap) / max(len(empties), 1)

    assert precision >= 0.95
    assert recall >= 0.95
