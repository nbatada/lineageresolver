UPDATED SOFTWARE SPEC — LineageResolver
Ambiguous identity adjudication for scRNA-seq with joint probabilistic inference
(Primary implementation in Python/Scanpy; R support deferred)

--------------------------------------------------------------------
0) PRODUCT GOAL (ONE SENTENCE)

Provide a fast, easy-to-use Scanpy/AnnData-first tool that resolves ambiguous immune cell identities at single-cell
resolution by jointly integrating RNA evidence, partial orthogonal evidence (e.g., TCR fragments), and an estimated ambient
RNA profile, producing calibrated posteriors and an abstention-aware label per cell.

--------------------------------------------------------------------
1) SCOPE AND POSITIONING

What the tool is
- A post-processing "adjudicator" for ambiguous identity problems where gene expression alone is insufficient.
- Intended to run on an ambiguous candidate set (recommended), but can run on all cells.

What the tool is not (v1)
- Not an aligner or UMI counter; assumes counts matrix exists.
- Not a general-purpose atlas labeler.
- Not a full SoupX replacement; it estimates an ambient profile sufficient for discounting spurious evidence in adjudication.

Primary v1 flagship task
- Cytotoxic lymphocyte adjudication:
  NK vs γδ T vs cytotoxic αβ T (optional MAIT extension).

Secondary v1 tasks (optional)
- CD4 vs CD8 boundary with ambient-driven double-positive artifacts.
- Provide a general plugin structure but implement only TCR-fragment plugin in v1.

--------------------------------------------------------------------
2) USER EXPERIENCE REQUIREMENTS

2.1 One-call API (primary)
Python function:
  lineageresolver.infer(
      adata,
      task_config,
      raw_10x_path=None,
      ambient_profile=None,
      evidence_table=None,
      candidate_set=None,
      tau=0.9,
      return_posteriors=True,
      inplace=True
  )

Behavior:
- If ambient_profile is None and raw_10x_path is provided:
    estimate ambient profile automatically from unfiltered 10x data.
- If ambient_profile is None and raw_10x_path is not provided:
    use fallback ambient estimation from adata counts and emit a warning + diagnostic flag.

Outputs (written to adata.obs):
- lineageresolver_label_map                 (MAP label: argmax posterior)
- lineageresolver_label_call                (abstention-aware: label or "uncertain")
- lineageresolver_max_p                     (max posterior probability)
- lineageresolver_uncertainty_entropy       (entropy over posterior)
- lineageresolver_p_<class>                 (one column per class if return_posteriors=True)

Plus evidence summary columns (if evidence present):
- lineageresolver_tcr_junction_umis
- lineageresolver_lambda_amb_tcr
- lineageresolver_amb_ll_tcr
- lineageresolver_module_<class>
- lineageresolver_a_hat (estimated ambient fraction proxy; if computed)

2.2 CLI (secondary but required for broad usability)
Executable:
  lineageresolver

Commands:
- estimate-ambient
- infer
- report

--------------------------------------------------------------------
3) INPUTS

3.1 Required for infer
- AnnData object with:
  adata.X (UMI counts; sparse matrix)
  adata.var_names (gene symbols)
  adata.obs_names (cell barcodes)
- task_config (YAML/JSON): class definitions, gene modules, evidence model parameters.

3.2 Optional inputs
A) raw_10x_path (preferred for robust ambient estimation)
- Path to 10x unfiltered output:
  raw_feature_bc_matrix/  OR  raw_feature_bc_matrix.h5

B) ambient_profile (advanced override)
- TSV: gene, ambient_fraction, ambient_count
- If supplied, skip estimation.

C) evidence_table (recommended for lineage tasks)
- TSV keyed by barcode, may include:
  tcr_trg_junction_umis
  tcr_trd_junction_umis
  tcr_total_junction_umis
  tcr_mapq_mean (optional)
  doublet_score (optional)
  other plugin-specific fields (future)

D) candidate_set
One of:
- obs column name containing boolean mask
- list of barcodes
- list of cluster ids (if obs has clustering column)
Default:
- all cells (warn if large; recommend providing candidate_set)

--------------------------------------------------------------------
4) AMBIENT PROFILE ESTIMATION (MANDATORY DEFAULT)

4.1 estimate_ambient_profile(raw_10x_path, ...)
Goal:
Estimate s_g = ambient gene frequency profile from empty/low-UMI droplets.

Algorithm (v1)
1) Read unfiltered 10x matrix as sparse.
2) Compute per-barcode total UMI counts.
3) Select empty droplets:
   - default: barcodes with total UMIs in [min_umi, max_umi] where max_umi is low (e.g., <= 50)
   - or select bottom p% by UMI (configurable).
4) Sum counts across selected droplets to get ambient counts per gene.
5) Normalize to fractions s_g.

Output:
- ambient_profile.tsv:
    gene, ambient_count, ambient_fraction
- ambient_report.json:
    n_barcodes_used
    total_umis_used
    umi_thresholds
    top_genes

Fallback ambient estimation (when raw_10x_path missing)
- Use low-UMI cells within adata as proxy empties (bottom p% by nUMI).
- Mark adata.uns["lineageresolver_ambient_estimation_mode"] = "fallback"
- Provide warning and include diagnostics.

--------------------------------------------------------------------
5) STATISTICAL MODEL (v1: FEATURE-LEVEL JOINT INFERENCE)

Design requirement:
The model must be truly joint in the sense that ambient is a competing explanation for both RNA markers and orthogonal evidence,
not a pre-correction step.

5.1 Candidate identities
Let z_c be identity class for cell c, z_c ∈ {1..K}.

We infer posterior:
  P(z_c | features_c)

5.2 Features per cell (computed on candidate set)
A) RNA module evidence
For each class k:
  module_k(c) = weighted module score from RNA counts (normalized)
Store as:
  lineageresolver_module_<class>

B) Orthogonal evidence (TCR fragments v1)
Let:
  kJ(c) = tcr_total_junction_umis for cell c (UMI-collapsed; junction-supporting)
If missing evidence_table:
  set kJ(c)=0 and mark evidence_missing flag.

C) Ambient likelihood for orthogonal evidence
Compute expected ambient TCR junction UMIs:
  lambda_amb_tcr(c) = alpha_tcr * nUMI(c) * a_hat(c)

Where:
  nUMI(c) is total UMIs in cell c (from adata.X)
  a_hat(c) is a per-cell ambient proxy (v1 simple estimator)
  alpha_tcr is derived from ambient profile:
    alpha_tcr = sum_{g in TCR_marker_set} s_g
or from a user-provided constant if TCR not present in gene space.

Then compute ambient log-likelihood for observed kJ:
  amb_ll_tcr(c) = -log PoissonPMF(kJ(c); lambda_amb_tcr(c))
(or NB alternative later).

D) QC covariates (optional)
- doublet_score (if provided)
- percent_mito (if available)

5.3 Classifier (softmax)
For each class k:
  score_k(c) =
      βk0
    + βk1 * module_k(c)
    + βk2 * kJ(c)
    - βk3 * amb_ll_tcr(c)
    + βk4 * qc_terms(c)

Posterior:
  p_k(c) = softmax_k(score_k(c))

Abstention:
  label_map(c)  = argmax_k p_k(c)
  label_call(c) = label_map(c) if max_k p_k(c) >= tau else "uncertain"

Calibration (v1)
- If user provides truth labels for a subset:
    optional Platt scaling on logits.
- Otherwise:
    ship with pretrained parameters calibrated on internal benchmark datasets (to be created by user later).

NOTE: v1 can ship with default parameters that are intentionally conservative (favor abstention) if no training is provided.

--------------------------------------------------------------------
6) TASK CONFIG FORMAT (YAML/JSON)

task_name: cytotoxic_adjudication_v1
classes: [NK, gdT, abT]
modules:
  NK:
    genes: [NKG7, GNLY, ...]
    weights: [...]
  gdT:
    genes: [TRDC, TRGC1, TRGC2, ...]
    weights: [...]
  abT:
    genes: [TRAC, TRBC1, TRBC2, ...]
    weights: [...]
tcr_marker_set_genes_for_alpha_tcr: [TRAC, TRBC1, TRBC2, TRDC, TRGC1, TRGC2]
parameters:
  beta:
    NK:  {b0:..., b_module:..., b_kJ:..., b_ambll:..., b_qc:...}
    gdT: {...}
    abT: {...}
defaults:
  tau: 0.9
  candidate_set_required: false
  warn_if_all_cells: true

--------------------------------------------------------------------
7) PERFORMANCE TARGETS

- Ambient estimation from raw 10x matrix:
  should handle common raw matrices (hundreds of thousands to millions of barcodes) without loading into R.
  Python implementation must use sparse I/O and avoid densification.

- Inference:
  O(N_candidate * K) where N_candidate is candidate set size.
  Should run on 100k candidate cells in minutes on typical CPU.

--------------------------------------------------------------------
8) QUALITY AND CORRECTNESS REQUIREMENTS

Unit tests must cover:
- ambient estimation correctness on toy matrices
- module scoring correctness on toy matrices
- Poisson likelihood and amb_ll correctness
- softmax numerical stability and posterior sums
- abstention behavior
- missing evidence behavior (evidence_table absent)
- deterministic output with fixed seed

Diagnostics:
- store ambient_report.json in adata.uns
- store model config and parameters in adata.uns
- include “ambient_estimation_mode” flag ("raw10x" vs "fallback")

--------------------------------------------------------------------
9) REPO STRUCTURE (PYTHON-FIRST)

repo/
  pyproject.toml
  README.md
  LICENSE
  lineageresolver/
    __init__.py
    api.py                  # infer(), estimate_ambient_profile()
    ambient.py              # raw10x reading + ambient estimation
    modules.py              # module scoring
    evidence.py             # evidence table parsing + derived features
    model.py                # softmax model + calibration (optional)
    metrics.py              # ECE, Brier, error-coverage
    plotting.py             # reliability, error-coverage
    io.py                   # AnnData I/O helpers
    cli.py                  # click/typer CLI entrypoint
    config.py               # config schema validation
  tests/
    test_ambient.py
    test_modules.py
    test_model.py
    test_infer_e2e.py
  examples/
    scanpy_notebook.ipynb
    configs/
      cytotoxic_adjudication_v1.yaml

R support deferred
- do not build R package in v1
- keep outputs in TSV + AnnData so future R wrapper can load results and attach to Seurat meta.data

--------------------------------------------------------------------
10) UPDATED TICKETS (ORDERED IMPLEMENT-TEST-CORRECT LOOP)

TICKET 0 — Kickoff scaffolding
- Create repo structure above
- Add pyproject, CLI skeleton, CI, formatting, unit test harness
Acceptance:
- `pip install -e .` works
- `lineageresolver --help` works
- tests run (empty tests allowed)

TICKET 1 — AnnData IO + one-call API skeleton
- Implement infer() signature in api.py
- Validate required AnnData fields
- Candidate set selection logic
- Write placeholder columns into adata.obs
Acceptance:
- infer() runs on toy AnnData and writes columns:
  lineageresolver_label_map, lineageresolver_label_call, lineageresolver_max_p

TICKET 2 — Ambient estimation from raw 10x (primary path)
- Implement estimate_ambient_profile(raw_10x_path)
- Support MTX directory and 10x h5 if feasible (optional)
- Write ambient_profile.tsv and ambient_report.json
Acceptance:
- unit test on small synthetic raw10x-style matrix
- returns correct gene fractions

TICKET 3 — Fallback ambient estimation (no raw10x)
- Implement fallback estimator using low-UMI cells from adata
- Set flags in adata.uns and warnings
Acceptance:
- test ensures fallback mode triggers without raw10x_path

TICKET 4 — Module scoring engine
- Implement weighted module scores from config
- Works on sparse matrices without densifying
Acceptance:
- unit test on toy matrix with expected scores

TICKET 5 — Evidence table ingestion + derived TCR features
- Parse evidence_table TSV keyed by barcode
- Compute kJ(c), lambda_amb_tcr(c), amb_ll_tcr(c)
- Handle missing barcodes and missing columns gracefully
Acceptance:
- tests for alignment of evidence to adata.obs_names
- amb_ll increases with higher kJ given fixed lambda

TICKET 6 — Softmax model inference + abstention
- Implement score_k(c), softmax probabilities, label_map/call
- Store p_class columns, entropy, max_p
Acceptance:
- posterior sums to 1 per cell
- abstention threshold works

TICKET 7 — End-to-end infer() with automatic ambient estimation
- infer() calls estimate_ambient_profile if raw_10x_path provided
- otherwise calls fallback
- compute modules + evidence features + model posteriors
Acceptance:
- one-call infer works end-to-end on synthetic dataset
- writes ambient report + config snapshot into adata.uns

TICKET 8 — Metrics + plots (calibration, error-coverage)
- Implement ECE, Brier, reliability data, error vs coverage
- Plotting functions (matplotlib)
Acceptance:
- functions run and produce files on toy data

TICKET 9 — CLI wiring
- `lineageresolver estimate-ambient ...`
- `lineageresolver infer ...`
- `lineageresolver report ...` (optional)
Acceptance:
- CLI commands run on example data

TICKET 10 — Documentation + example notebook
- Provide scanpy tutorial with:
  - ambient estimation from raw 10x
  - inference on filtered AnnData
  - interpret results (max_p, uncertain)
Acceptance:
- notebook runs end-to-end

TICKET 11 (optional) — Performance profiling + optimization
- profile ambient estimation on large raw10x
- optimize via chunking or sparse ops if needed
Acceptance:
- document runtime/memory on a representative dataset

--------------------------------------------------------------------
11) “PM KICKOFF” TICKET (FOR PROJECT MANAGER INSTANCE)

TICKET PM-0 — Project kickoff: create detailed ticket specs and repo plan

Objective
Initialize the project management workstream with:
- the repo structure
- the ordered tickets
- acceptance criteria per ticket
- interfaces between modules
- coding standards and test requirements

Deliverables
1) A PROJECT_PLAN.md that includes:
   - goals and non-goals
   - architecture overview
   - data flow
   - API and CLI specs
   - config schema example
   - output schema (obs columns and uns fields)

2) A TICKETS/ directory containing one markdown file per ticket:
   - TICKET_00_kickoff_scaffolding.md
   - TICKET_01_io_api_skeleton.md
   - ...
   Each ticket file must include:
     - scope
     - out-of-scope
     - detailed steps
     - acceptance criteria
     - tests to write
     - files/modules expected to change

3) A FILE_STRUCTURE.md that matches the repo layout and explains each module.

4) A CONFIG_SPEC.md describing the YAML schema and required fields.

Acceptance criteria
- All ticket docs are internally consistent and reference the same interfaces.
- Each ticket can be implemented independently in sequence.
- The one-call infer() path is preserved as the user-facing default.

--------------------------------------------------------------------
12) RECOMMENDATION ON R SUPPORT

Defer R to v2.

Rationale
- You prefer Scanpy/Python and will benchmark there first.
- Once the algorithm and outputs are stable, an R wrapper is straightforward:
  - export per-cell TSV
  - join into Seurat meta.data by barcode
  - optionally call python via reticulate later

Therefore:
- build Python-first, produce stable TSV/AnnData outputs
- add R wrapper only after v1 benchmarking demonstrates value

--------------------------------------------------------------------
END SPEC