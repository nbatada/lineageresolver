# TICKET 01: IO and API Skeleton

## Title
Implement AnnData validation and one-call API skeleton.

## Goal
Provide an operational `infer(...)` entrypoint with validated inputs, candidate selection, and placeholder outputs written to `adata.obs`.

## Scope
- Implement the public `infer(...)` signature from spec.
- Validate required AnnData fields.
- Support candidate set routing (mask/obs column/barcode list/cluster ids).
- Write required output columns with placeholder values.

## Non-Goals
- No real ambient estimation or probabilistic model yet.
- No calibration or metrics outputs.

## Detailed Steps
1. Define `infer(...)` signature and runtime argument checks.
2. Implement candidate selection helper in `io.py` or `api.py`.
3. Create consistent output columns for required names.
4. Add diagnostics entry under `adata.uns` for placeholder mode.
5. Ensure `inplace=True/False` behavior is predictable.

## Acceptance Criteria
- `infer()` runs on toy AnnData without errors.
- Required columns exist in `adata.obs` after call.
- Candidate set modes are accepted and validated.
- Missing required AnnData attributes produce clear errors.

## Tests To Add
- `test_api_infer_signature.py`
- `test_candidate_set_selection.py`
- `test_placeholder_output_columns.py`

## Files/Modules Expected To Change
- `lineageresolver/api.py`
- `lineageresolver/io.py`
- `tests/test_api_infer_signature.py`
- `tests/test_candidate_set_selection.py`
