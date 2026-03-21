TICKET PM-0 — Project kickoff: create specs, docs, and implementation tickets

Objective
Initialize the LineageResolver repository so that:
1) the project specification is captured in a canonical location,
2) supporting documentation is generated from the spec,
3) a complete set of implementation tickets exists in docs/tickets/,
4) the repo has an agreed file structure and development workflow.

Repository context
Local path:
  ~/git/lineageresolver

Source of truth
The main specification must live at:
  docs/SOFTWARE_SPEC.md

If docs/SOFTWARE_SPEC.md does not exist yet, create it and populate it with the latest software spec
(provided in this conversation).

Deliverables (must create these files in-repo)

A) Spec and docs (all markdown)
1) docs/SOFTWARE_SPEC.md
   - full updated spec (Python/Scanpy-first, automatic ambient estimation from raw 10x, one-call infer API)
   - include: scope, inputs/outputs, ambient estimation, statistical model v1, config schema, performance, tests, repo structure

2) docs/PROJECT_PLAN.md
   - implementation roadmap and milestones
   - ordered ticket list with dependencies
   - success metrics for v1 (functional + statistical correctness + speed)

3) docs/ARCHITECTURE.md
   - module-level architecture diagram (text-based ok)
   - data flow: raw10x -> ambient_profile -> features -> posteriors -> adata.obs
   - interfaces between modules: ambient/modules/evidence/model/api/cli

4) docs/FILE_STRUCTURE.md
   - canonical repository tree
   - role of each module file

5) docs/CONFIG_SPEC.md
   - YAML/JSON schema for task_config
   - required fields, defaults, validation rules
   - example config for cytotoxic adjudication (NK vs gdT vs abT)

6) docs/DEVELOPMENT_WORKFLOW.md
   - implement-test-correct loop guidance
   - testing strategy (unit + small e2e)
   - code style (formatting, linting, type hints)
   - release strategy (versioning, changelog)

B) Tickets
Create directory:
  docs/tickets/

Create one markdown file per ticket with:
  - title
  - goal
  - scope
  - non-goals
  - detailed steps
  - acceptance criteria
  - tests to add
  - files/modules expected to change

Tickets to create (required)
- TICKET_00_repo_scaffold.md
- TICKET_01_io_api_skeleton.md
- TICKET_02_ambient_raw10x_estimation.md
- TICKET_03_ambient_fallback_estimation.md
- TICKET_04_module_scoring_engine.md
- TICKET_05_evidence_table_and_tcr_features.md
- TICKET_06_softmax_inference_and_abstention.md
- TICKET_07_end_to_end_infer_auto_ambient.md
- TICKET_08_metrics_and_plots.md
- TICKET_09_cli_wiring.md
- TICKET_10_docs_and_example_notebook.md
Optional:
- TICKET_11_profiling_and_optimization.md

C) README update
Update README.md to include:
- a short description of LineageResolver
- a link/reference to docs/SOFTWARE_SPEC.md
- a minimal usage snippet for infer() (even if placeholder)
- project status note (early development)

Constraints and priorities
- Primary environment: Python + Scanpy/AnnData
- No R package in v1; outputs must be stable for later R integration
- No C core in v1
- Ambient profile estimation must be automatic by default
- Must support raw/unfiltered 10x input without requiring loading it into R
- Must preserve a one-call user workflow: infer() does everything when raw_10x_path is provided
- Must use sparse matrices; avoid densification
- Must include unit tests for mathematical correctness (Poisson likelihood, softmax stability, abstention behavior)

Acceptance criteria
- All docs exist at the specified paths and are internally consistent
- Tickets are ordered, implementable sequentially, and each has clear acceptance tests
- README references docs/SOFTWARE_SPEC.md and reflects the one-call workflow