# Canonical File Structure

## Repository Tree

```text
lineageresolver/
  pyproject.toml
  README.md
  LICENSE
  docs/
    SOFTWARE_SPEC.md
    PROJECT_PLAN.md
    ARCHITECTURE.md
    FILE_STRUCTURE.md
    CONFIG_SPEC.md
    DEVELOPMENT_WORKFLOW.md
    tickets/
      PM_0_KICKOFF.md
      TICKET_00_repo_scaffold.md
      TICKET_01_io_api_skeleton.md
      TICKET_02_ambient_raw10x_estimation.md
      TICKET_03_ambient_fallback_estimation.md
      TICKET_04_module_scoring_engine.md
      TICKET_05_evidence_table_and_tcr_features.md
      TICKET_06_softmax_inference_and_abstention.md
      TICKET_07_end_to_end_infer_auto_ambient.md
      TICKET_08_metrics_and_plots.md
      TICKET_09_cli_wiring.md
      TICKET_10_docs_and_example_notebook.md
      TICKET_11_profiling_and_optimization.md
  lineageresolver/
    __init__.py
    api.py
    ambient.py
    modules.py
    evidence.py
    model.py
    metrics.py
    plotting.py
    io.py
    cli.py
    config.py
  tests/
    test_ambient.py
    test_modules.py
    test_evidence.py
    test_model.py
    test_infer_e2e.py
  examples/
    scanpy_notebook.ipynb
    configs/
      cytotoxic_adjudication_v1.yaml
```

## Module Responsibilities
- `lineageresolver/api.py`: public Python entrypoints and orchestration.
- `lineageresolver/ambient.py`: raw 10x and fallback ambient estimators.
- `lineageresolver/modules.py`: weighted module scoring from expression matrix.
- `lineageresolver/evidence.py`: evidence ingestion, barcode alignment, and derived TCR features.
- `lineageresolver/model.py`: score computation, softmax posterior inference, abstention logic.
- `lineageresolver/metrics.py`: ECE, Brier, and error-coverage calculations.
- `lineageresolver/plotting.py`: reliability and coverage plotting.
- `lineageresolver/io.py`: input/output helpers for AnnData and artifact files.
- `lineageresolver/cli.py`: command-line command wiring.
- `lineageresolver/config.py`: schema validation and defaults application.

## Documentation Roles
- `docs/SOFTWARE_SPEC.md`: source-of-truth product and technical spec.
- `docs/PROJECT_PLAN.md`: roadmap, dependencies, milestones, success criteria.
- `docs/ARCHITECTURE.md`: module architecture, interfaces, data flow.
- `docs/CONFIG_SPEC.md`: task config schema and examples.
- `docs/DEVELOPMENT_WORKFLOW.md`: implementation, testing, and release process.
- `docs/tickets/*.md`: implementation tickets and acceptance contracts.

## Test Coverage Mapping
- `test_ambient.py`: raw/fallback ambient correctness and diagnostics.
- `test_modules.py`: module score correctness on toy sparse matrices.
- `test_evidence.py`: evidence alignment and derived feature correctness.
- `test_model.py`: softmax stability, posterior normalization, abstention.
- `test_infer_e2e.py`: one-call integration behavior and output schema.
