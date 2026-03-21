# TICKET 11: Profiling and Optimization (Optional)

## Title
Profile runtime/memory and optimize bottlenecks.

## Goal
Improve runtime and memory behavior for large candidate sets and raw10x ambient estimation while preserving correctness.

## Scope
- Profile key paths (`ambient`, `modules`, `model`, end-to-end infer).
- Identify and implement high-impact optimizations.
- Document performance baselines before and after changes.

## Non-Goals
- No algorithm redesign that changes output semantics.
- No low-level compiled acceleration requirement for v1.

## Detailed Steps
1. Add reproducible profiling scripts for representative inputs.
2. Capture baseline CPU time and memory for critical flows.
3. Eliminate avoidable allocations and dense conversions.
4. Tune sparse operations and I/O chunking where needed.
5. Record benchmark deltas and any tradeoffs.

## Acceptance Criteria
- Performance report produced with reproducible commands.
- At least one meaningful bottleneck improvement demonstrated or justified.
- No regression in unit/e2e correctness tests.

## Tests To Add
- `test_perf_regression_guard_small.py` (optional threshold guard)
- Existing correctness tests rerun post-optimization

## Files/Modules Expected To Change
- `lineageresolver/ambient.py`
- `lineageresolver/modules.py`
- `lineageresolver/model.py`
- `scripts/profile_*.py`
- `docs/PROJECT_PLAN.md` (performance notes)
