from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import anndata as ad


H5AD_PATH = Path("tests/data/filtered_small.h5ad")
EVIDENCE_PATH = Path("tests/data/evidence_small.tsv")
RAW10X_PATH = Path("tests/data/raw10x_small/raw_feature_bc_matrix")
CONFIG_PATH = Path("tests/data/configs/cytotoxic_adjudication_v1.yaml")


def test_cli_infer_smoke(tmp_path: Path) -> None:
    out_path = tmp_path / "infer_out.h5ad"
    cmd = [
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
        str(out_path),
    ]
    completed = subprocess.run(cmd, check=False, capture_output=True, text=True)

    assert completed.returncode == 0, completed.stderr
    assert out_path.exists()

    out = ad.read_h5ad(out_path)
    assert "lineageresolver_label_call" in out.obs.columns
    assert "lineageresolver_p_NK" in out.obs.columns
