# TICKET 02: Ambient Estimation From Raw 10x

## Title
Implement primary ambient estimator from unfiltered 10x matrix.

## Goal
Estimate ambient profile (`ambient_fraction` and `ambient_count`) from empty/low-UMI droplets in raw 10x data.

## Scope
- Read raw 10x matrix from directory MTX path.
- Optional support for 10x H5 input.
- Select low-UMI droplets using configured thresholds.
- Compute and return normalized ambient profile.
- Emit ambient diagnostics report.
- Use committed synthetic fixture `tests/data/raw10x_small/raw_feature_bc_matrix/` for deterministic validation.

## Non-Goals
- No fallback estimation from filtered AnnData (separate ticket).
- No model integration yet.

## Detailed Steps
1. Add matrix reader for raw10x sparse formats.
2. Compute total UMI per barcode.
3. Implement empty-droplet selection strategy (`min_umi`, `max_umi`, optional bottom percent).
4. Aggregate selected droplet counts per gene.
5. Normalize to ambient fractions and build report structure.
6. Add artifact writing helpers (`ambient_profile.tsv`, `ambient_report.json`).
7. Validate behavior against `tests/data/raw10x_small/manifest.json` barcode partitions.

## Acceptance Criteria
- Function returns non-negative counts/fractions and normalized sum close to 1.
- Ambient report includes barcode count, UMI totals, and threshold metadata.
- Sparse path avoids dense matrix conversion.
- Running on `tests/data/raw10x_small/raw_feature_bc_matrix/` produces nonzero ambient fraction for at least one TR gene in `{TRAC, TRBC1, TRBC2, TRDC, TRGC1, TRGC2}`.
- Empty-droplet selection overlaps strongly with `manifest.json` `empty_barcodes` set.

## Tests To Add
- `test_ambient_raw10x_fixture_estimation.py` (uses `tests/data/raw10x_small/raw_feature_bc_matrix/`)
- `test_ambient_fraction_normalization.py`
- `test_ambient_empty_selection_rules.py` (cross-checks `manifest.json`)

## Files/Modules Expected To Change
- `lineageresolver/ambient.py`
- `lineageresolver/io.py`
- `tests/test_ambient_raw10x_fixture_estimation.py`
