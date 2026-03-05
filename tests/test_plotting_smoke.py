from __future__ import annotations

from pathlib import Path

import pandas as pd

from lineageresolver.plotting import plot_error_coverage_curve, plot_reliability_curve


def test_plotting_smoke(tmp_path: Path) -> None:
    reliability = pd.DataFrame(
        {
            "bin_left": [0.0, 0.5],
            "bin_right": [0.5, 1.0],
            "count": [10, 10],
            "accuracy": [0.4, 0.8],
            "confidence": [0.3, 0.85],
        }
    )
    coverage = pd.DataFrame(
        {
            "coverage": [0.2, 0.4, 0.6, 0.8, 1.0],
            "error_rate": [0.0, 0.1, 0.15, 0.2, 0.25],
            "confidence_threshold": [0.95, 0.9, 0.8, 0.7, 0.5],
        }
    )

    rel_path = plot_reliability_curve(reliability, tmp_path / "reliability.png")
    cov_path = plot_error_coverage_curve(coverage, tmp_path / "coverage.png")

    assert rel_path.exists()
    assert cov_path.exists()
    assert rel_path.stat().st_size > 0
    assert cov_path.stat().st_size > 0
