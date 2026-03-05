from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from lineageresolver.evidence import derive_tcr_features


H5AD_PATH = Path("tests/data/filtered_small.h5ad")
EVIDENCE_PATH = Path("tests/data/evidence_small.tsv")


def test_evidence_barcode_alignment() -> None:
    adata = ad.read_h5ad(H5AD_PATH)
    ambient_profile = pd.DataFrame(
        {
            "gene": adata.var_names.astype(str),
            "ambient_count": np.ones(adata.n_vars),
            "ambient_fraction": np.ones(adata.n_vars) / adata.n_vars,
        }
    )

    features, diagnostics = derive_tcr_features(
        adata=adata,
        evidence_table=EVIDENCE_PATH,
        ambient_profile=ambient_profile,
    )

    # Missing evidence rows in fixture should default to zero kJ.
    assert features.loc["CELL_0026-1", "lineageresolver_tcr_junction_umis"] == 0.0
    assert features.loc["CELL_0032-1", "lineageresolver_tcr_junction_umis"] == 0.0

    # A populated row should align and retain known value.
    assert features.loc["CELL_0001-1", "lineageresolver_tcr_junction_umis"] == 5.0
    assert diagnostics["evidence_missing_rows"] > 0
