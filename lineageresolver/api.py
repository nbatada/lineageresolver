"""Public API entrypoints."""

from __future__ import annotations

from typing import Any


def infer(
    adata: Any,
    task_config: Any,
    raw_10x_path: str | None = None,
    ambient_profile: Any | None = None,
    evidence_table: str | None = None,
    candidate_set: Any | None = None,
    tau: float = 0.9,
    return_posteriors: bool = True,
    inplace: bool = True,
):
    """Inference API placeholder; implemented in later tickets."""
    raise NotImplementedError("infer() is not implemented yet.")


def estimate_ambient_profile(*args: Any, **kwargs: Any):
    """Ambient estimation API placeholder; implemented in later tickets."""
    raise NotImplementedError("estimate_ambient_profile() is not implemented yet.")
