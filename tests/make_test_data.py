#!/usr/bin/env python3
"""Generate deterministic synthetic fixtures for LineageResolver CI tests.

Design summary:
- Genes: 200 (required immune markers + filler genes).
- Raw 10x barcodes: 500 total (400 empty droplets + 100 real cells).
- Empty droplets: UMI totals sampled uniformly from [0, 30], dominated by ambient profile.
- Cells: UMI totals sampled uniformly from [300, 1000], generated from class-specific
  expression plus additive ambient contamination where ambient fraction is sampled from
  clipped Beta(2.5, 20.0) in [0.02, 0.20].
- Truth labels: 34 NK, 33 gdT, 33 abT.
- Evidence table: barcode keyed, intentionally omits 8 cell barcodes to test missing evidence.

Outputs:
- raw10x fixture (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`)
- filtered AnnData fixture (`filtered_small.h5ad`)
- evidence fixture (`evidence_small.tsv`)
- task config fixture (`configs/cytotoxic_adjudication_v1.yaml`)
- manifest (`raw10x_small/manifest.json`) with partitions, labels, and contamination metadata

The generation uses a fixed RNG seed and writes small files suitable for CI.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse


DEFAULT_SEED = 20260305
N_GENES = 200
N_EMPTY = 400
N_CELLS = 100


def _build_gene_list() -> List[str]:
    required = [
        "NKG7",
        "GNLY",
        "GZMB",
        "PRF1",
        "TRAC",
        "TRBC1",
        "TRBC2",
        "TRDC",
        "TRGC1",
        "TRGC2",
    ]
    extras = [
        "MALAT1",
        "LTB",
        "ACTB",
        "RPLP0",
        "HBB",
        "HBA1",
        "LST1",
        "TYROBP",
        "MT-ND1",
        "MT-CO1",
        "MT-ATP6",
    ]

    genes: List[str] = []
    for g in required + extras:
        if g not in genes:
            genes.append(g)

    filler_idx = 1
    while len(genes) < N_GENES:
        genes.append(f"GENE_{filler_idx:03d}")
        filler_idx += 1

    return genes


def _normalize(weights: np.ndarray) -> np.ndarray:
    total = float(weights.sum())
    if total <= 0:
        raise ValueError("Weight vector must have positive sum")
    return weights / total


def _ambient_profile_weights(gene_to_idx: Dict[str, int], n_genes: int) -> np.ndarray:
    w = np.ones(n_genes, dtype=float)
    boost = {
        "MALAT1": 120.0,
        "HBB": 80.0,
        "HBA1": 60.0,
        "LST1": 40.0,
        "TYROBP": 25.0,
        "NKG7": 15.0,
        "GNLY": 10.0,
        "TRAC": 1.5,
        "TRBC1": 1.0,
        "TRBC2": 1.0,
        "TRDC": 0.8,
        "TRGC1": 0.8,
        "TRGC2": 0.8,
    }
    for gene, value in boost.items():
        w[gene_to_idx[gene]] = value
    return _normalize(w)


def _cell_class_weights(
    klass: str,
    gene_to_idx: Dict[str, int],
    n_genes: int,
) -> np.ndarray:
    w = np.full(n_genes, 0.5, dtype=float)

    common = {
        "MALAT1": 2.0,
        "ACTB": 2.0,
        "RPLP0": 1.5,
        "LTB": 1.0,
    }
    for gene, value in common.items():
        w[gene_to_idx[gene]] = value

    if klass == "NK":
        boosts = {
            "NKG7": 25.0,
            "GNLY": 22.0,
            "GZMB": 14.0,
            "PRF1": 12.0,
            "TRAC": 0.3,
            "TRBC1": 0.3,
            "TRBC2": 0.3,
            "TRDC": 0.2,
            "TRGC1": 0.2,
            "TRGC2": 0.2,
        }
    elif klass == "gdT":
        boosts = {
            "TRDC": 20.0,
            "TRGC1": 18.0,
            "TRGC2": 14.0,
            "NKG7": 6.0,
            "GNLY": 5.0,
            "GZMB": 4.0,
            "PRF1": 4.0,
            "TRAC": 1.0,
            "TRBC1": 1.0,
            "TRBC2": 1.0,
        }
    elif klass == "abT":
        boosts = {
            "TRAC": 20.0,
            "TRBC1": 16.0,
            "TRBC2": 16.0,
            "NKG7": 5.0,
            "GNLY": 3.0,
            "GZMB": 3.0,
            "PRF1": 2.0,
            "TRDC": 1.0,
            "TRGC1": 1.0,
            "TRGC2": 1.0,
        }
    else:
        raise ValueError(f"Unknown class: {klass}")

    for gene, value in boosts.items():
        w[gene_to_idx[gene]] = value

    return _normalize(w)


def _write_gzip_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename=str(path), mode="wb", mtime=0) as gz:
        gz.write(payload)


def _write_gzip_text(path: Path, text: str) -> None:
    _write_gzip_bytes(path, text.encode("utf-8"))


def _write_mtx_gz(path: Path, matrix: sparse.spmatrix) -> None:
    buf = io.BytesIO()
    spio.mmwrite(buf, matrix.tocoo())
    _write_gzip_bytes(path, buf.getvalue())


def _make_raw10x_and_cell_counts(seed: int) -> Tuple[
    List[str],
    List[str],
    List[str],
    Dict[str, str],
    Dict[str, float],
    np.ndarray,
    np.ndarray,
]:
    rng = np.random.default_rng(seed)

    genes = _build_gene_list()
    n_genes = len(genes)
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    empty_barcodes = [f"EMPTY_{i:04d}-1" for i in range(1, N_EMPTY + 1)]
    cell_barcodes = [f"CELL_{i:04d}-1" for i in range(1, N_CELLS + 1)]
    all_barcodes = empty_barcodes + cell_barcodes

    class_labels = np.array(["NK"] * 34 + ["gdT"] * 33 + ["abT"] * 33, dtype=object)
    rng.shuffle(class_labels)
    truth_label_map = {bc: str(lbl) for bc, lbl in zip(cell_barcodes, class_labels)}

    ambient_probs = _ambient_profile_weights(gene_to_idx, n_genes)
    class_probs = {
        "NK": _cell_class_weights("NK", gene_to_idx, n_genes),
        "gdT": _cell_class_weights("gdT", gene_to_idx, n_genes),
        "abT": _cell_class_weights("abT", gene_to_idx, n_genes),
    }

    counts_all = np.zeros((len(all_barcodes), n_genes), dtype=np.int32)

    for i in range(N_EMPTY):
        total = int(rng.integers(0, 31))
        if total > 0:
            counts_all[i, :] = rng.multinomial(total, ambient_probs)

    contamination_fraction_map: Dict[str, float] = {}
    for j, bc in enumerate(cell_barcodes, start=N_EMPTY):
        label = truth_label_map[bc]
        total = int(rng.integers(300, 1001))
        ambient_frac = float(np.clip(rng.beta(2.5, 20.0), 0.02, 0.20))
        n_amb = int(round(total * ambient_frac))
        n_sig = total - n_amb

        sig_counts = rng.multinomial(n_sig, class_probs[label])
        amb_counts = rng.multinomial(n_amb, ambient_probs)
        counts_all[j, :] = sig_counts + amb_counts

        contamination_fraction_map[bc] = ambient_frac

    cell_counts = counts_all[N_EMPTY:, :]

    return (
        genes,
        empty_barcodes,
        cell_barcodes,
        truth_label_map,
        contamination_fraction_map,
        counts_all,
        cell_counts,
    )


def _make_evidence_table(
    seed: int,
    cell_barcodes: List[str],
    truth_label_map: Dict[str, str],
) -> Tuple[pd.DataFrame, List[str]]:
    rng = np.random.default_rng(seed + 17)

    missing = sorted(rng.choice(cell_barcodes, size=8, replace=False).tolist())

    records = []
    for bc in cell_barcodes:
        if bc in missing:
            continue

        label = truth_label_map[bc]
        if label == "NK":
            total = int(rng.choice([0, 1, 2], p=[0.82, 0.16, 0.02]))
            split = rng.multinomial(total, [0.35, 0.35, 0.30])
        elif label == "gdT":
            total = int(rng.integers(1, 6))
            split = rng.multinomial(total, [0.45, 0.45, 0.10])
        elif label == "abT":
            total = int(rng.integers(1, 6))
            split = rng.multinomial(total, [0.10, 0.10, 0.80])
        else:
            raise ValueError(f"Unknown class label: {label}")

        trg = int(split[0])
        trd = int(split[1])
        records.append(
            {
                "barcode": bc,
                "tcr_trg_junction_umis": trg,
                "tcr_trd_junction_umis": trd,
                "tcr_total_junction_umis": total,
            }
        )

    df = pd.DataFrame.from_records(records).sort_values("barcode").reset_index(drop=True)
    return df, missing


def _write_raw10x_fixture(
    output_root: Path,
    genes: List[str],
    all_barcodes: List[str],
    counts_all: np.ndarray,
) -> None:
    raw_dir = output_root / "raw10x_small" / "raw_feature_bc_matrix"
    raw_dir.mkdir(parents=True, exist_ok=True)

    feature_ids = [f"SYNTH_GENE_{i:04d}" for i in range(1, len(genes) + 1)]
    features_lines = [f"{gid}\t{gname}\tGene Expression" for gid, gname in zip(feature_ids, genes)]
    barcodes_lines = all_barcodes

    matrix = sparse.coo_matrix(counts_all.T)

    _write_mtx_gz(raw_dir / "matrix.mtx.gz", matrix)
    _write_gzip_text(raw_dir / "features.tsv.gz", "\n".join(features_lines) + "\n")
    _write_gzip_text(raw_dir / "barcodes.tsv.gz", "\n".join(barcodes_lines) + "\n")


def _write_filtered_h5ad(
    output_root: Path,
    genes: List[str],
    cell_barcodes: List[str],
    truth_label_map: Dict[str, str],
    contamination_fraction_map: Dict[str, float],
    cell_counts: np.ndarray,
) -> None:
    mito_gene_idx = [i for i, g in enumerate(genes) if g.startswith("MT-")]

    n_counts = cell_counts.sum(axis=1)
    if mito_gene_idx:
        mito_counts = cell_counts[:, mito_gene_idx].sum(axis=1)
        percent_mito = (mito_counts / np.maximum(n_counts, 1)) * 100.0
    else:
        percent_mito = np.zeros_like(n_counts, dtype=float)

    obs = pd.DataFrame(
        {
            "truth_label": [truth_label_map[bc] for bc in cell_barcodes],
            "n_counts": n_counts.astype(np.int64),
            "percent_mito": percent_mito.astype(float),
            "synthetic_ambient_fraction": [contamination_fraction_map[bc] for bc in cell_barcodes],
        },
        index=pd.Index(cell_barcodes, name="barcode"),
    )

    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    var["gene_id"] = [f"SYNTH_GENE_{i:04d}" for i in range(1, len(genes) + 1)]

    adata = ad.AnnData(X=sparse.csr_matrix(cell_counts), obs=obs, var=var)
    adata.write_h5ad(output_root / "filtered_small.h5ad", compression="gzip")


def _write_task_config(output_root: Path) -> None:
    config_dir = output_root / "configs"
    config_dir.mkdir(parents=True, exist_ok=True)

    config_text = """task_name: cytotoxic_adjudication_v1
classes: [NK, gdT, abT]
modules:
  NK:
    genes: [NKG7, GNLY, GZMB, PRF1]
    weights: [1.0, 1.0, 0.8, 0.8]
  gdT:
    genes: [TRDC, TRGC1, TRGC2, NKG7]
    weights: [1.2, 1.0, 1.0, 0.5]
  abT:
    genes: [TRAC, TRBC1, TRBC2, LTB]
    weights: [1.2, 1.0, 1.0, 0.5]
tcr_marker_set_genes_for_alpha_tcr: [TRAC, TRBC1, TRBC2, TRDC, TRGC1, TRGC2]
parameters:
  beta:
    NK:  {b0: 0.0, b_module: 1.0, b_kJ: -0.3, b_ambll: 0.3, b_qc: 0.0}
    gdT: {b0: 0.0, b_module: 1.0, b_kJ: 0.8,  b_ambll: 0.6, b_qc: 0.0}
    abT: {b0: 0.0, b_module: 1.0, b_kJ: 0.5,  b_ambll: 0.4, b_qc: 0.0}
defaults:
  tau: 0.9
  candidate_set_required: false
  warn_if_all_cells: true
"""
    (config_dir / "cytotoxic_adjudication_v1.yaml").write_text(config_text, encoding="utf-8")


def _write_manifest(
    output_root: Path,
    seed: int,
    genes: List[str],
    empty_barcodes: List[str],
    cell_barcodes: List[str],
    truth_label_map: Dict[str, str],
    contamination_fraction_map: Dict[str, float],
    missing_evidence_barcodes: List[str],
) -> None:
    class_counts = {
        "NK": sum(1 for b in cell_barcodes if truth_label_map[b] == "NK"),
        "gdT": sum(1 for b in cell_barcodes if truth_label_map[b] == "gdT"),
        "abT": sum(1 for b in cell_barcodes if truth_label_map[b] == "abT"),
    }

    manifest = {
        "seed": seed,
        "summary": {
            "n_genes": len(genes),
            "n_barcodes_total": len(empty_barcodes) + len(cell_barcodes),
            "n_empty_barcodes": len(empty_barcodes),
            "n_cell_barcodes": len(cell_barcodes),
            "class_counts": class_counts,
        },
        "design": {
            "empty_droplets": "UMI totals sampled uniformly in [0, 30], counts dominated by ambient profile.",
            "cells": "UMI totals sampled uniformly in [300, 1000], with class-specific signal plus ambient contamination.",
            "contamination_regime": {
                "ambient_fraction_distribution": "clipped Beta(2.5, 20.0)",
                "clip_range": [0.02, 0.20],
                "ambient_TR_genes_nonzero": True,
            },
            "evidence_missing_rows": len(missing_evidence_barcodes),
        },
        "genes": genes,
        "empty_barcodes": empty_barcodes,
        "cell_barcodes": cell_barcodes,
        "truth_label_by_barcode": truth_label_map,
        "synthetic_ambient_fraction_by_barcode": {
            bc: round(contamination_fraction_map[bc], 6) for bc in cell_barcodes
        },
        "missing_evidence_barcodes": missing_evidence_barcodes,
    }

    raw_root = output_root / "raw10x_small"
    raw_root.mkdir(parents=True, exist_ok=True)
    (raw_root / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")


def generate_fixtures(output_root: Path, seed: int) -> None:
    (
        genes,
        empty_barcodes,
        cell_barcodes,
        truth_label_map,
        contamination_fraction_map,
        counts_all,
        cell_counts,
    ) = _make_raw10x_and_cell_counts(seed)

    all_barcodes = empty_barcodes + cell_barcodes

    _write_raw10x_fixture(output_root, genes, all_barcodes, counts_all)
    _write_filtered_h5ad(
        output_root,
        genes,
        cell_barcodes,
        truth_label_map,
        contamination_fraction_map,
        cell_counts,
    )

    evidence_df, missing_evidence_barcodes = _make_evidence_table(seed, cell_barcodes, truth_label_map)
    evidence_df.to_csv(output_root / "evidence_small.tsv", sep="\t", index=False)

    _write_task_config(output_root)
    _write_manifest(
        output_root,
        seed,
        genes,
        empty_barcodes,
        cell_barcodes,
        truth_label_map,
        contamination_fraction_map,
        missing_evidence_barcodes,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate deterministic synthetic test fixtures.")
    parser.add_argument(
        "--output-root",
        default="tests/data",
        help="Directory where fixtures are generated (default: tests/data)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"RNG seed (default: {DEFAULT_SEED})",
    )
    args = parser.parse_args()

    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    generate_fixtures(output_root=output_root, seed=args.seed)


if __name__ == "__main__":
    main()
