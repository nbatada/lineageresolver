# LineageResolver

LineageResolver is a Python/Scanpy-first toolkit for adjudicating ambiguous single-cell immune lineage identity using joint probabilistic inference over RNA, orthogonal evidence, and ambient RNA explanations.

## Source of truth
The canonical product and technical specification is in [docs/SOFTWARE_SPEC.md](docs/SOFTWARE_SPEC.md).

## Minimal usage (`infer` one-call workflow)

```python
import scanpy as sc
import lineageresolver as lr

adata = sc.read_h5ad("filtered_cells.h5ad")

lr.infer(
    adata=adata,
    task_config="examples/configs/cytotoxic_adjudication_v1.yaml",
    raw_10x_path="raw_feature_bc_matrix/",  # enables automatic ambient estimation
    evidence_table="evidence.tsv",
    tau=0.9,
    return_posteriors=True,
    inplace=True,
)

# outputs are written to adata.obs / adata.uns
print(adata.obs[["lineageresolver_label_map", "lineageresolver_label_call", "lineageresolver_max_p"]].head())
```

## CLI commands
- `lineageresolver estimate-ambient --raw-10x-path tests/data/raw10x_small/raw_feature_bc_matrix --output-dir out/ambient`
- `lineageresolver infer --adata tests/data/filtered_small.h5ad --task-config tests/data/configs/cytotoxic_adjudication_v1.yaml --raw-10x-path tests/data/raw10x_small/raw_feature_bc_matrix --evidence-table tests/data/evidence_small.tsv --output out/infer.h5ad`
- `lineageresolver report --adata out/infer.h5ad --label-column truth_label --output-dir out/report`

## Output columns
`infer()` writes per-cell outputs to `adata.obs`:
- `lineageresolver_label_map`
- `lineageresolver_label_call`
- `lineageresolver_max_p`
- `lineageresolver_uncertainty_entropy`
- `lineageresolver_p_<class>` (if `return_posteriors=True`)
- `lineageresolver_tcr_junction_umis`
- `lineageresolver_lambda_amb_tcr`
- `lineageresolver_amb_ll_tcr`
- `lineageresolver_a_hat`

Diagnostics are written to `adata.uns`, including:
- `lineageresolver_ambient_estimation_mode`
- `lineageresolver_ambient_report`
- `lineageresolver_task_config`

## Notebook walkthrough
A runnable Scanpy walkthrough is provided at:
- `examples/scanpy_notebook.ipynb`
- `examples/configs/cytotoxic_adjudication_v1.yaml`

## Project status
Early development. Core docs, architecture, and ticketed implementation plan are defined; code implementation is in progress.
