from __future__ import annotations

from pathlib import Path
import subprocess
import sys


RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")


def test_cli_estimate_ambient_smoke(tmp_path: Path) -> None:
    cmd = [
        sys.executable,
        "-m",
        "lineageresolver.cli",
        "estimate-ambient",
        "--raw-10x-path",
        str(RAW10X_PATH),
        "--output-dir",
        str(tmp_path),
    ]
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True)

    assert completed.returncode == 0, completed.stderr
    assert (tmp_path / "ambient_profile.tsv").exists()
    assert (tmp_path / "ambient_report.json").exists()
