# TICKET 03: Fallback Ambient Estimation

## Title
Implement fallback ambient estimator using low-UMI cells in AnnData.

## Goal
Provide robust fallback behavior when `raw_10x_path` is unavailable while keeping one-call workflow functional.

## Scope
- Estimate ambient profile from low-UMI cells in filtered AnnData.
- Set explicit fallback mode diagnostics and warnings.
- Return profile in same format as raw10x estimator.
- Validate fallback behavior with `tests/data/filtered_small.h5ad`.

## Non-Goals
- No replacement for raw10x estimator quality.
- No external contamination modeling beyond basic proxy strategy.

## Detailed Steps
1. Add fallback selection logic (bottom percent or explicit threshold).
2. Aggregate counts and normalize to ambient fractions.
3. Emit warning and set `lineageresolver_ambient_estimation_mode=fallback`.
4. Record fallback diagnostics in `adata.uns`.
5. Ensure shape and column parity with primary estimator output.
6. Assert fallback path can run deterministically on `tests/data/filtered_small.h5ad` when `raw_10x_path=None`.

## Acceptance Criteria
- Fallback runs when raw10x path is absent.
- Mode flag and warnings are emitted consistently.
- Output schema matches primary ambient estimator.
- Running fallback on `tests/data/filtered_small.h5ad` returns a profile with fraction sum near 1 and expected gene coverage.

## Tests To Add
- `test_fallback_trigger_without_raw10x.py` (uses `tests/data/filtered_small.h5ad`)
- `test_fallback_output_schema.py`
- `test_fallback_mode_flag_written.py`

## Files/Modules Expected To Change
- `lineageresolver/ambient.py`
- `lineageresolver/api.py`
- `tests/test_fallback_trigger_without_raw10x.py`
