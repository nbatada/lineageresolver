from __future__ import annotations

from pathlib import Path
import subprocess
import sys


H5AD_PATH = Path("tests/data/filtered_small.h5ad")
EVIDENCE_PATH = Path("tests/data/evidence_small.tsv")
RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")
CONFIG_PATH = Path("tests/data/configs/cytotoxic_adjudication_v1.yaml")


def test_cli_report_smoke(tmp_path: Path) -> None:
    infer_out = tmp_path / "infer_out.h5ad"
    infer_cmd = [
        sys.executable,
        "-m",
        "lineageresolver.cli",
        "infer",
        "--adata",
        str(H5AD_PATH),
        "--task-config",
        str(CONFIG_PATH),
        "--raw-10x-path",
        str(RAW10X_PATH),
        "--evidence-table",
        str(EVIDENCE_PATH),
        "--output",
        str(infer_out),
    ]
    infer_completed = subprocess.run(infer_cmd, check=False, capture_output=True, text=True)
    assert infer_completed.returncode == 0, infer_completed.stderr

    report_dir = tmp_path / "report"
    report_cmd = [
        sys.executable,
        "-m",
        "lineageresolver.cli",
        "report",
        "--adata",
        str(infer_out),
        "--label-column",
        "truth_label",
        "--output-dir",
        str(report_dir),
    ]
    report_completed = subprocess.run(report_cmd, check=False, capture_output=True, text=True)

    assert report_completed.returncode == 0, report_completed.stderr
    assert (report_dir / "metrics.json").exists()
    assert (report_dir / "reliability_bins.csv").exists()
    assert (report_dir / "error_coverage.csv").exists()
    assert (report_dir / "reliability.png").exists()
    assert (report_dir / "coverage.png").exists()
