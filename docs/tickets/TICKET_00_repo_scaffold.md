# TICKET 00: Repository Scaffold

## Title
Repository scaffold for Python/Scanpy-first LineageResolver implementation.

## Goal
Create the baseline package, tooling, and test harness needed to implement v1 features incrementally.

## Scope
- Create canonical package/test/example directory layout.
- Add package metadata and editable install support.
- Add linting, formatting, and test configuration.
- Add CLI entrypoint stub and API stubs.

## Non-Goals
- No feature-complete ambient estimation or model logic.
- No finalized output schema implementation.

## Detailed Steps
1. Add `pyproject.toml` with dependencies and console script entrypoint (`lineageresolver`).
2. Create package files under `lineageresolver/` with import-safe stubs.
3. Add test package and minimal smoke tests.
4. Add formatter/linter/type-check config.
5. Ensure `pip install -e .` and command discovery both work.

## Acceptance Criteria
- Editable install succeeds.
- `lineageresolver --help` returns usage text.
- Test suite runs (`pytest`) with scaffold-level checks.
- Directory tree aligns with `docs/FILE_STRUCTURE.md`.

## Tests To Add
- `test_imports.py`: package import smoke test.
- `test_cli_smoke.py`: CLI executable responds to `--help`.

## Files/Modules Expected To Change
- `pyproject.toml`
- `lineageresolver/__init__.py`
- `lineageresolver/api.py`
- `lineageresolver/cli.py`
- `tests/test_imports.py`
- `tests/test_cli_smoke.py`
