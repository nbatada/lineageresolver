# Project Plan: LineageResolver v1

## Purpose
Deliver a Python/Scanpy-first package that resolves ambiguous immune lineage identity per cell using joint probabilistic inference over RNA evidence, orthogonal evidence (TCR fragments in v1), and ambient RNA explanation.

## Goals
- Provide a one-call user workflow through `lineageresolver.infer(...)`.
- Automatically estimate ambient profile from raw 10x input by default.
- Produce calibrated posteriors and an abstention-aware call (`label` or `uncertain`).
- Ship a usable CLI (`estimate-ambient`, `infer`, `report`).
- Keep implementation sparse-matrix safe (no densification for core paths).

## Non-Goals (v1)
- No R package or Seurat-native implementation.
- No C/C++ acceleration layer.
- No general-purpose atlas annotation framework.
- No replacement for full ambient correction pipelines.

## Milestones

### M0: Scaffolding and Interfaces
- Repository skeleton, packaging, test harness, lint/type config.
- API and CLI entrypoint stubs.
- Ticket(s): `TICKET_00`, `TICKET_01`.

### M1: Ambient Estimation
- Primary raw 10x ambient estimator.
- Fallback ambient estimator from low-UMI cells.
- Diagnostics and mode flagging.
- Ticket(s): `TICKET_02`, `TICKET_03`.

### M2: Feature Engineering
- Weighted RNA module scoring.
- Evidence table ingestion and TCR-derived features (`kJ`, `lambda_amb_tcr`, `amb_ll_tcr`).
- Ticket(s): `TICKET_04`, `TICKET_05`.

### M3: Inference Core
- Class score computation and numerically stable softmax.
- Posterior outputs, entropy, MAP label, abstention thresholding.
- Ticket(s): `TICKET_06`.

### M4: End-to-End Workflow and Reporting
- Integrated one-call inference pipeline with automatic ambient routing.
- Metrics and reporting plots.
- Ticket(s): `TICKET_07`, `TICKET_08`.

### M5: CLI and User Documentation
- CLI wiring for ambient estimation and inference.
- Example notebook and end-user docs.
- Ticket(s): `TICKET_09`, `TICKET_10`.

### Optional M6: Performance Hardening
- Profiling at larger candidate sizes and targeted optimization.
- Ticket(s): `TICKET_11`.

## Ordered Ticket Plan With Dependencies
1. `TICKET_00_repo_scaffold` (no deps)
2. `TICKET_01_io_api_skeleton` (depends on 00)
3. `TICKET_02_ambient_raw10x_estimation` (depends on 00)
4. `TICKET_03_ambient_fallback_estimation` (depends on 01, 02)
5. `TICKET_04_module_scoring_engine` (depends on 01)
6. `TICKET_05_evidence_table_and_tcr_features` (depends on 02, 03, 04)
7. `TICKET_06_softmax_inference_and_abstention` (depends on 04, 05)
8. `TICKET_07_end_to_end_infer_auto_ambient` (depends on 03, 04, 05, 06)
9. `TICKET_08_metrics_and_plots` (depends on 06, 07)
10. `TICKET_09_cli_wiring` (depends on 02, 07, 08)
11. `TICKET_10_docs_and_example_notebook` (depends on 07, 08, 09)
12. `TICKET_11_profiling_and_optimization` (optional; depends on 07, 09)

## Testing and Fixtures

Required committed fixtures (CI-safe):
- `tests/data/raw10x_small/raw_feature_bc_matrix/`
  - `matrix.mtx.gz`
  - `features.tsv.gz`
  - `barcodes.tsv.gz`
- `tests/data/raw10x_small/manifest.json`
- `tests/data/filtered_small.h5ad`
- `tests/data/evidence_small.tsv`
- `tests/data/configs/cytotoxic_adjudication_v1.yaml`

Fixture generation:
- `tests/make_test_data.py` is the canonical generator.
- Fixtures are deterministic via fixed seed and must be reproducible in CI.
- Total fixture footprint must remain small (target under 10 MB committed).

Required tests depending on fixtures:
- Ambient estimator tests using `raw10x_small` for empty-droplet detection and `ambient_fraction` normalization.
- Fallback ambient tests using `filtered_small.h5ad` without `raw_10x_path`.
- E2E infer tests using `filtered_small.h5ad` + `evidence_small.tsv` + task config for:
  - posterior normalization and softmax stability
  - abstention behavior under configured `tau`
  - missing evidence-row handling for barcodes absent from evidence TSV
  - ambient mode routing (`raw10x` vs `fallback`)

## v1 Success Metrics

### Functional Correctness
- `infer()` writes expected `adata.obs` columns for map call, abstention-aware call, posterior max, entropy, and per-class probabilities.
- Ambient estimation mode is explicit in `adata.uns` (`raw10x` or `fallback`).
- Missing evidence table path is handled gracefully and deterministically.

### Statistical Correctness
- Posterior probabilities sum to 1 (within tolerance).
- Abstention threshold behavior exactly follows configured `tau`.
- Poisson likelihood feature (`amb_ll_tcr`) is monotonic in expected directions.
- Optional calibration path works when truth labels are supplied.

### Performance
- Ambient estimation runs with sparse operations and no densification.
- Inference complexity scales roughly `O(N_candidate * K)`.
- Candidate-set run on 100k cells completes in minutes on commodity CPU.

## Key Risks and Mitigations
- Risk: ambient profile quality degrades without raw 10x.
  - Mitigation: explicit fallback diagnostics, warnings, and report outputs.
- Risk: unstable posteriors for low-signal cells.
  - Mitigation: conservative defaults and explicit abstention policy.
- Risk: memory blowups with sparse-to-dense conversions.
  - Mitigation: tests that assert sparse-safe implementations in core paths.

## Definition of Done for v1
- All required tickets (00-10) accepted.
- Docs in `docs/` are consistent with API, config schema, and outputs.
- CLI and Python API both execute full one-call workflow.
- Unit tests and a small e2e test pass in CI.
