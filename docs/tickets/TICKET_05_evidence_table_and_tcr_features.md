# TICKET 05: Evidence Table and TCR Features

## Title
Implement evidence table ingestion and TCR-derived features.

## Goal
Integrate orthogonal evidence keyed by barcode and derive model-ready TCR features including ambient-adjusted likelihood terms.

## Scope
- Parse TSV evidence table and align with `adata.obs_names`.
- Build `kJ`, `lambda_amb_tcr`, and `amb_ll_tcr` features.
- Support absent evidence table with neutral defaults.

## Non-Goals
- No plugin system beyond v1 TCR feature path.
- No model coefficient fitting.

## Detailed Steps
1. Implement evidence table loader and barcode join logic.
2. Add handling for missing rows/columns with diagnostics.
3. Compute `kJ` from available TCR junction fields.
4. Compute `alpha_tcr` from ambient profile and configured marker genes.
5. Compute `lambda_amb_tcr = alpha_tcr * nUMI * a_hat` and Poisson-derived `amb_ll_tcr`.
6. Write evidence feature columns to `adata.obs`.

## Acceptance Criteria
- Evidence rows align correctly to AnnData barcodes.
- Missing evidence table does not crash inference.
- Derived features match toy expected calculations.

## Tests To Add
- `test_evidence_barcode_alignment.py`
- `test_tcr_feature_derivation.py`
- `test_missing_evidence_defaults.py`

## Files/Modules Expected To Change
- `lineageresolver/evidence.py`
- `lineageresolver/ambient.py`
- `lineageresolver/api.py`
- `tests/test_evidence_barcode_alignment.py`
