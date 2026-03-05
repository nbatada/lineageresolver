"""Plotting utilities."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def plot_reliability_curve(reliability_df: pd.DataFrame, output_path: str | Path) -> Path:
    """Plot reliability diagram and save to file."""
    out = Path(output_path)
    valid = reliability_df.dropna(subset=["accuracy", "confidence"])

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot([0, 1], [0, 1], linestyle="--", color="gray", linewidth=1)
    ax.plot(valid["confidence"], valid["accuracy"], marker="o", linewidth=1.5)
    ax.set_xlabel("Mean predicted confidence")
    ax.set_ylabel("Empirical accuracy")
    ax.set_title("Reliability Diagram")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_error_coverage_curve(curve_df: pd.DataFrame, output_path: str | Path) -> Path:
    """Plot error-rate vs coverage and save to file."""
    out = Path(output_path)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(curve_df["coverage"], curve_df["error_rate"], linewidth=1.5)
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Error rate")
    ax.set_title("Error-Coverage Curve")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out
