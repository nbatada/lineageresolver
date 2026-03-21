# TICKET 08: Metrics and Plots

## Title
Add calibration and uncertainty reporting utilities.

## Goal
Provide quantitative and visual diagnostics for posterior quality and abstention behavior.

## Scope
- Implement ECE, Brier score, and error-coverage calculations.
- Implement plotting helpers for reliability and coverage curves.
- Support report artifact generation for CLI and notebook use.

## Non-Goals
- No dashboard service or UI layer.
- No large-scale benchmarking framework.

## Detailed Steps
1. Add metrics functions with clear input contracts.
2. Add plotting utilities with deterministic outputs.
3. Add serialization helpers for report JSON/CSV.
4. Document expected usage in module docstrings.

## Acceptance Criteria
- Metrics functions produce expected values on toy data.
- Plot functions generate files without runtime errors.
- Error messages are clear when truth labels are absent.

## Tests To Add
- `test_metrics_ece_brier.py`
- `test_metrics_error_coverage.py`
- `test_plotting_smoke.py`

## Files/Modules Expected To Change
- `lineageresolver/metrics.py`
- `lineageresolver/plotting.py`
- `tests/test_metrics_ece_brier.py`
