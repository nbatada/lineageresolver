"""Command line interface for LineageResolver."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from lineageresolver.api import estimate_ambient_profile, infer
from lineageresolver.metrics import brier_score_multiclass, error_coverage_curve, expected_calibration_error
from lineageresolver.plotting import plot_error_coverage_curve, plot_reliability_curve


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lineageresolver")
    subparsers = parser.add_subparsers(dest="command")

    ambient_parser = subparsers.add_parser(
        "estimate-ambient",
        help="Estimate ambient profile from raw 10x data",
    )
    ambient_parser.add_argument("--raw-10x-path", required=True, help="Path to raw 10x MTX directory")
    ambient_parser.add_argument("--output-dir", required=True, help="Directory for ambient artifacts")
    ambient_parser.add_argument("--min-umi", type=int, default=0)
    ambient_parser.add_argument("--max-umi", type=int, default=50)
    ambient_parser.add_argument("--bottom-percent", type=float, default=None)

    infer_parser = subparsers.add_parser("infer", help="Run one-call inference on AnnData")
    infer_parser.add_argument("--adata", required=True, help="Input AnnData .h5ad path")
    infer_parser.add_argument("--task-config", required=True, help="Task config YAML/JSON path")
    infer_parser.add_argument("--output", required=True, help="Output .h5ad path")
    infer_parser.add_argument("--raw-10x-path", default=None, help="Raw 10x path for ambient estimation")
    infer_parser.add_argument("--ambient-profile", default=None, help="Ambient profile TSV override")
    infer_parser.add_argument("--evidence-table", default=None, help="Evidence table TSV path")
    infer_parser.add_argument("--tau", type=float, default=0.9)
    infer_parser.add_argument("--no-posteriors", action="store_true", help="Skip posterior columns")

    report_parser = subparsers.add_parser("report", help="Generate diagnostics report artifacts")
    report_parser.add_argument("--adata", required=True, help="Input .h5ad with infer outputs")
    report_parser.add_argument("--label-column", default="truth_label", help="Ground-truth label column in obs")
    report_parser.add_argument("--output-dir", required=True, help="Directory for report outputs")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "estimate-ambient":
        return _cmd_estimate_ambient(args)
    if args.command == "infer":
        return _cmd_infer(args)
    if args.command == "report":
        return _cmd_report(args)
    parser.error(f"Unknown command: {args.command}")
    return 2


def _cmd_estimate_ambient(args: argparse.Namespace) -> int:
    estimate_ambient_profile(
        raw_10x_path=args.raw_10x_path,
        min_umi=args.min_umi,
        max_umi=args.max_umi,
        bottom_percent=args.bottom_percent,
        output_dir=args.output_dir,
    )
    return 0


def _cmd_infer(args: argparse.Namespace) -> int:
    try:
        import anndata as ad
    except ModuleNotFoundError as exc:
        raise RuntimeError("anndata is required for CLI infer command.") from exc

    adata = ad.read_h5ad(args.adata)
    out = infer(
        adata=adata,
        task_config=args.task_config,
        raw_10x_path=args.raw_10x_path,
        ambient_profile=args.ambient_profile,
        evidence_table=args.evidence_table,
        tau=args.tau,
        return_posteriors=not args.no_posteriors,
        inplace=False,
    )
    _sanitize_uns_for_h5ad(out)
    out.write_h5ad(args.output)
    return 0


def _cmd_report(args: argparse.Namespace) -> int:
    try:
        import anndata as ad
    except ModuleNotFoundError as exc:
        raise RuntimeError("anndata is required for CLI report command.") from exc

    adata = ad.read_h5ad(args.adata)
    if args.label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{args.label_column}' not found in AnnData obs.")

    posterior_cols = [col for col in adata.obs.columns if col.startswith("lineageresolver_p_")]
    if not posterior_cols:
        raise ValueError("No posterior columns found (expected lineageresolver_p_<class>).")

    classes = [col.replace("lineageresolver_p_", "", 1) for col in posterior_cols]
    probs = adata.obs[posterior_cols].to_numpy(dtype=float)
    y_true = adata.obs[args.label_column].astype(str).to_numpy()
    valid = np.isin(y_true, classes) & np.isfinite(probs).all(axis=1)
    if valid.sum() == 0:
        raise ValueError("No rows with valid ground-truth labels matching posterior classes.")

    y_true = y_true[valid]
    probs = probs[valid]

    ece, rel = expected_calibration_error(y_true=y_true, probs=probs, class_labels=classes)
    brier = brier_score_multiclass(y_true=y_true, probs=probs, class_labels=classes)
    curve = error_coverage_curve(y_true=y_true, probs=probs, class_labels=classes)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics_path = out_dir / "metrics.json"
    rel_csv = out_dir / "reliability_bins.csv"
    cov_csv = out_dir / "error_coverage.csv"
    rel_png = out_dir / "reliability.png"
    cov_png = out_dir / "coverage.png"

    with open(metrics_path, "w", encoding="utf-8") as handle:
        json.dump({"ece": float(ece), "brier": float(brier), "n": int(len(y_true))}, handle, indent=2)
    rel.to_csv(rel_csv, index=False)
    curve.to_csv(cov_csv, index=False)
    plot_reliability_curve(rel, rel_png)
    plot_error_coverage_curve(curve, cov_png)
    return 0


def _sanitize_uns_for_h5ad(adata: Any) -> None:
    """Convert nested uns payloads that may not round-trip through h5ad."""
    for key in list(adata.uns.keys()):
        if not str(key).startswith("lineageresolver"):
            continue
        value = adata.uns[key]
        if isinstance(value, dict):
            adata.uns[key] = json.dumps(value, default=_json_default)


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


if __name__ == "__main__":
    raise SystemExit(main())
