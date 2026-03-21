# Task Configuration Spec (YAML/JSON)

## Purpose
`task_config` defines class set, module features, evidence settings, and model parameters for `lineageresolver.infer(...)`.

## Accepted Formats
- YAML (`.yaml` / `.yml`)
- JSON (`.json`)

Both formats map to the same schema.

## Schema

### Top-Level Fields
- `task_name` (string, required)
- `classes` (array of string, required, length >= 2)
- `modules` (object, required)
- `tcr_marker_set_genes_for_alpha_tcr` (array of string, optional)
- `parameters` (object, required)
- `defaults` (object, optional)

### `modules` object
For each class name in `classes`, `modules.<class>` is required with:
- `genes` (array of string, required, non-empty)
- `weights` (array of number, optional)

Validation rules:
- Every class in `classes` must exist in `modules`.
- `genes` must be unique per class after normalization.
- If `weights` is provided, it must have same length as `genes`.
- If `weights` is missing, default to `1.0` per gene.

### `parameters.beta`
`parameters.beta` must contain one object per class with numeric coefficients:
- `b0`
- `b_module`
- `b_kJ`
- `b_ambll`
- `b_qc` (optional, default `0.0`)

Validation rules:
- Every class in `classes` must have beta parameters.
- All beta values must be finite numeric values.

### `defaults`
Optional defaults for runtime behavior:
- `tau` (float, default `0.9`, valid range `(0, 1]`)
- `candidate_set_required` (bool, default `false`)
- `warn_if_all_cells` (bool, default `true`)

## Runtime Validation Rules
- Unknown top-level keys are warnings by default (strict mode may reject).
- Duplicate class names are invalid.
- Configured class names must be identifier-safe for output column generation.
- Missing genes in `adata.var_names` are allowed but must be reported in diagnostics.
- If all module genes for a class are missing, inference should continue with warning (class module score defaults to 0).

## Cytotoxic Adjudication Example

```yaml
task_name: cytotoxic_adjudication_v1
classes: [NK, gdT, abT]

modules:
  NK:
    genes: [NKG7, GNLY, CTSW, PRF1]
    weights: [1.0, 1.0, 0.8, 0.8]
  gdT:
    genes: [TRDC, TRGC1, TRGC2, CCL5]
    weights: [1.2, 1.0, 1.0, 0.7]
  abT:
    genes: [TRAC, TRBC1, TRBC2, LTB]
    weights: [1.2, 1.0, 1.0, 0.6]

tcr_marker_set_genes_for_alpha_tcr:
  - TRAC
  - TRBC1
  - TRBC2
  - TRDC
  - TRGC1
  - TRGC2

parameters:
  beta:
    NK:  {b0: 0.0, b_module: 1.1, b_kJ: -0.2, b_ambll: 0.4, b_qc: 0.0}
    gdT: {b0: 0.0, b_module: 1.0, b_kJ: 0.9,  b_ambll: 0.6, b_qc: 0.0}
    abT: {b0: 0.0, b_module: 1.0, b_kJ: 0.5,  b_ambll: 0.3, b_qc: 0.0}

defaults:
  tau: 0.9
  candidate_set_required: false
  warn_if_all_cells: true
```

## Defaulting Behavior Summary
- `tau` defaults to `0.9` if not set.
- Missing module weights default to uniform `1.0`.
- `b_qc` defaults to `0.0` when omitted.
- Evidence-free inference is valid; TCR-derived terms are set to neutral defaults with diagnostics.
