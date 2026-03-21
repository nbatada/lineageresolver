# TICKET 10: Documentation and Example Notebook

## Title
Create user documentation and runnable notebook for v1 workflow.

## Goal
Provide clear onboarding docs and an example demonstrating one-call inference with ambient estimation.

## Scope
- Update docs for API/CLI usage and outputs.
- Add notebook walkthrough using toy/small dataset.
- Add example task config file.

## Non-Goals
- No large benchmark publication.
- No production deployment guide.

## Detailed Steps
1. Add or update quickstart documentation and command examples.
2. Create notebook with steps:
   - load AnnData
   - run infer with raw10x path
   - inspect posterior/abstention outputs
3. Add expected outputs and interpretation guidance.
4. Validate notebook execution end-to-end.

## Acceptance Criteria
- Notebook executes without manual edits on reference environment.
- Docs match actual API and CLI flags.
- README links to spec and notebook resources.

## Tests To Add
- `test_docs_snippets_sync.py` (optional smoke check)
- `test_notebook_smoke.py` (optional, if notebook CI enabled)

## Files/Modules Expected To Change
- `README.md`
- `docs/*.md`
- `examples/scanpy_notebook.ipynb`
- `examples/configs/cytotoxic_adjudication_v1.yaml`
