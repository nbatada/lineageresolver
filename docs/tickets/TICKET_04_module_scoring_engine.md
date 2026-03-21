# TICKET 04: Module Scoring Engine

## Title
Implement weighted RNA module scoring for each class.

## Goal
Generate sparse-safe per-cell module features (`lineageresolver_module_<class>`) based on config-defined gene sets and weights.

## Scope
- Parse module definitions per class.
- Map genes to `adata.var_names`.
- Compute weighted module scores per cell.
- Handle missing genes with diagnostics.

## Non-Goals
- No classifier logic in this ticket.
- No score calibration in this ticket.

## Detailed Steps
1. Implement gene index resolution with robust missing-gene handling.
2. Compute weighted sum/average module scores on sparse matrix input.
3. Normalize score outputs consistently across classes.
4. Write per-class module score columns to `adata.obs`.
5. Emit diagnostics for missing genes and zero-feature classes.

## Acceptance Criteria
- Expected module scores match toy matrix calculations.
- Implementation avoids densifying full matrix.
- Missing gene behavior is deterministic and documented.

## Tests To Add
- `test_module_scoring_weighted_toy.py`
- `test_module_missing_gene_handling.py`
- `test_module_sparse_input.py`

## Files/Modules Expected To Change
- `lineageresolver/modules.py`
- `lineageresolver/api.py`
- `tests/test_module_scoring_weighted_toy.py`
