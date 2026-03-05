from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass
class MiniAnnData:
    X: Any
    obs_names: pd.Index
    var_names: pd.Index
    obs: pd.DataFrame
    uns: dict[str, Any] = field(default_factory=dict)

    @property
    def n_obs(self) -> int:
        return len(self.obs_names)

    def copy(self) -> "MiniAnnData":
        if sparse.issparse(self.X):
            x_copy = self.X.copy()
        else:
            x_copy = np.array(self.X, copy=True)
        return MiniAnnData(
            X=x_copy,
            obs_names=self.obs_names.copy(),
            var_names=self.var_names.copy(),
            obs=self.obs.copy(deep=True),
            uns={k: v for k, v in self.uns.items()},
        )


def make_mini_adata(n_cells: int = 6, n_genes: int = 4) -> MiniAnnData:
    x = sparse.csr_matrix(np.arange(n_cells * n_genes, dtype=float).reshape(n_cells, n_genes))
    obs_names = pd.Index([f"CELL_{i:04d}-1" for i in range(1, n_cells + 1)], dtype="object")
    var_names = pd.Index([f"GENE_{i:03d}" for i in range(1, n_genes + 1)], dtype="object")
    obs = pd.DataFrame(index=obs_names)
    return MiniAnnData(X=x, obs_names=obs_names, var_names=var_names, obs=obs, uns={})
