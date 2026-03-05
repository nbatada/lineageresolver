"""Command line interface for LineageResolver."""

from __future__ import annotations

import argparse
from typing import Sequence


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lineageresolver")
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("estimate-ambient", help="Estimate ambient profile from raw 10x data")
    subparsers.add_parser("infer", help="Run one-call inference on AnnData")
    subparsers.add_parser("report", help="Generate diagnostics report artifacts")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    parser.parse_args(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
