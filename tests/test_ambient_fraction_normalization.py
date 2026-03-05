from __future__ import annotations

from pathlib import Path

from lineageresolver.ambient import estimate_ambient_profile


RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")


def test_ambient_fraction_normalization_and_nonnegativity() -> None:
    result = estimate_ambient_profile(RAW10X_PATH)

    assert (result.profile["ambient_count"] >= 0).all()
    assert (result.profile["ambient_fraction"] >= 0).all()
    assert abs(float(result.profile["ambient_fraction"].sum()) - 1.0) < 1e-8
    assert abs(float(result.profile["ambient_count"].sum()) - result.report["total_umis_used"]) < 1e-8
