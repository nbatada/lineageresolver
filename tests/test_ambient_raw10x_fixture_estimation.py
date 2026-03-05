from __future__ import annotations

from pathlib import Path

from lineageresolver.ambient import TR_MARKER_GENES, estimate_ambient_profile


RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")


def test_ambient_raw10x_fixture_estimation() -> None:
    result = estimate_ambient_profile(RAW10X_PATH)

    assert set(["gene", "ambient_count", "ambient_fraction"]).issubset(result.profile.columns)
    assert result.report["n_barcodes_used"] > 0
    assert result.report["total_umis_used"] > 0.0
    assert abs(float(result.profile["ambient_fraction"].sum()) - 1.0) < 1e-8

    tr_rows = result.profile[result.profile["gene"].isin(TR_MARKER_GENES)]
    assert not tr_rows.empty
    assert tr_rows["ambient_fraction"].sum() > 0.0


def test_ambient_artifacts_are_written(tmp_path: Path) -> None:
    estimate_ambient_profile(RAW10X_PATH, output_dir=tmp_path)

    assert (tmp_path / "ambient_profile.tsv").exists()
    assert (tmp_path / "ambient_report.json").exists()
