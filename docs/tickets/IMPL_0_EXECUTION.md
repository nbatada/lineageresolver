TICKET IMPL-0 — Execute implementation tickets with implement-test-correct loop

Objective
Implement LineageResolver by completing the existing tickets in docs/tickets/ in strict numeric order.
For each ticket:
  implement -> run tests -> fix until tests pass -> commit -> proceed.

Source of truth
- docs/SOFTWARE_SPEC.md
- ticket markdown files in docs/tickets/

Rules
1) Follow ticket order:
   TICKET_00_* then TICKET_01_* … through TICKET_10_* (and optional TICKET_11_*).
2) Do not start a new ticket until all tests pass for the current ticket.
3) Add or update tests required by each ticket before considering it complete.
4) Do not change public interfaces (API/CLI/config schema) without updating:
   - docs/SOFTWARE_SPEC.md
   - docs/CONFIG_SPEC.md
   - impacted ticket acceptance criteria
5) Every ticket must end with:
   - tests passing locally
   - a commit message: "Complete TICKET_XX: <short title>"

Testing requirements (critical)
- Unit tests for mathematical correctness:
  - PoissonPMF and amb_ll_tcr behavior (monotonicity with kJ; stability for small/large lambda)
  - softmax numerical stability (prob sums to 1; no NaNs for large logits)
  - abstention logic (tau threshold)
- Functional tests:
  - infer() writes required columns to adata.obs
  - estimate_ambient_profile() produces correct s_g on toy raw10x matrices
  - end-to-end infer() runs with:
      (a) raw_10x_path provided
      (b) fallback ambient mode (no raw10x_path)
  - CLI entry points run on example data

Test data policy
- Do NOT commit large real raw 10x datasets to git.
- Create small synthetic raw10x-style matrices in tests (MatrixMarket format) stored under tests/data/.
- Provide scripts to download public example datasets for manual benchmarking under examples/ (optional).
- All automated CI tests must run using only small synthetic fixtures.

Deliverables
- All tickets completed with passing tests.
- A working one-call API:
    lineageresolver.infer(adata, task_config, raw_10x_path=..., evidence_table=...)
- A working CLI:
    lineageresolver estimate-ambient ...
    lineageresolver infer ...