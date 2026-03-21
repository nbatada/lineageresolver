# TICKET 07: End-to-End infer() With Automatic Ambient

## Title
Integrate one-call inference pipeline with automatic ambient estimation routing.

## Goal
Deliver functional end-to-end `infer()` behavior that chooses ambient source automatically and writes complete output schema.

## Scope
- Orchestrate ambient estimation, module scoring, evidence features, and model inference.
- Route ambient path: provided profile vs raw10x estimate vs fallback estimate.
- Persist diagnostics and config snapshots.
- Exercise deterministic e2e path with committed fixtures in `tests/data/`.

## Non-Goals
- No CLI wiring in this ticket.
- No benchmark-level optimization.

## Detailed Steps
1. Implement orchestration in `api.infer` with explicit path decision logic.
2. Ensure candidate mask is honored consistently through all stages.
3. Write required output columns to `adata.obs`.
4. Write ambient/mode/config diagnostics to `adata.uns`.
5. Validate behavior for missing optional inputs.
6. Add fixture-based e2e coverage using:
   - `tests/data/raw10x_small/raw_feature_bc_matrix/`
   - `tests/data/filtered_small.h5ad`
   - `tests/data/evidence_small.tsv`
   - `tests/data/configs/cytotoxic_adjudication_v1.yaml`

## Acceptance Criteria
- One-call `infer()` runs with raw10x path and without raw10x path.
- Required obs/uns fields are present and internally consistent.
- Candidate subset behavior does not corrupt full-object alignment.
- Fixture-based e2e test passes for both ambient modes and handles barcodes missing from `evidence_small.tsv` without failure.
- Posterior columns are finite and row-wise probabilities sum to 1 within tolerance on `filtered_small.h5ad`.

## Tests To Add
- `test_infer_e2e_raw10x_path.py` (fixture-based)
- `test_infer_e2e_fallback_path.py` (fixture-based)
- `test_infer_candidate_subset_alignment.py`
- `test_infer_missing_evidence_rows_fixture.py`

## Files/Modules Expected To Change
- `lineageresolver/api.py`
- `lineageresolver/ambient.py`
- `lineageresolver/modules.py`
- `lineageresolver/evidence.py`
- `lineageresolver/model.py`
- `tests/test_infer_e2e_raw10x_path.py`
- `tests/test_infer_e2e_fallback_path.py`
- `tests/test_infer_missing_evidence_rows_fixture.py`
