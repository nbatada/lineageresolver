# TICKET 09: CLI Wiring

## Title
Implement user-facing CLI commands for ambient estimation and inference.

## Goal
Expose primary workflows through `lineageresolver` command without requiring custom Python scripts.

## Scope
- `lineageresolver estimate-ambient`
- `lineageresolver infer`
- `lineageresolver report`
- Input/output argument parsing and artifact writing.

## Non-Goals
- No distributed execution support.
- No interactive TUI.

## Detailed Steps
1. Wire command group and shared options.
2. Add `estimate-ambient` command and output paths.
3. Add `infer` command for full one-call workflow.
4. Add `report` command backed by metrics/plot modules.
5. Ensure consistent exit codes and error handling.

## Acceptance Criteria
- Each command appears in `lineageresolver --help`.
- Commands run successfully on toy data.
- CLI outputs match Python API output contracts.

## Tests To Add
- `test_cli_estimate_ambient_smoke.py`
- `test_cli_infer_smoke.py`
- `test_cli_report_smoke.py`

## Files/Modules Expected To Change
- `lineageresolver/cli.py`
- `lineageresolver/api.py`
- `lineageresolver/io.py`
- `tests/test_cli_estimate_ambient_smoke.py`
