# Architecture: LineageResolver v1

## Architectural Principles
- Python-first with Scanpy/AnnData-native interfaces.
- One-call inference path is the primary product surface.
- Ambient profile acts as a competing explanation in the model (not only preprocessing).
- Sparse-matrix operations only in the critical path.

## Module-Level Architecture

```text
+-------------------+         +--------------------+
|  User / Pipeline  |         |   Task Config      |
| infer(...) / CLI  |         | (YAML or JSON)     |
+---------+---------+         +----------+---------+
          |                               |
          v                               v
    +-----+-------------------------------+------+
    |                 API Layer                  |
    | lineageresolver.api.infer / estimate_*     |
    +-----+-------------------------------+------+
          |                               |
          |                               |
          v                               v
+---------+---------+           +---------+---------+
| Ambient Estimator |           | Config Validator  |
| raw10x/fallback   |           | schema + defaults |
+---------+---------+           +---------+---------+
          |                               |
          +---------------+---------------+
                          |
                          v
                +---------+---------+
                | Feature Builders  |
                | modules + evidence|
                +---------+---------+
                          |
                          v
                +---------+---------+
                |  Model Core       |
                | softmax + tau     |
                +---------+---------+
                          |
                          v
                +---------+---------+
                | Outputs            |
                | adata.obs / uns    |
                +---------+---------+
                          |
                          v
                +---------+---------+
                | Metrics / Plots    |
                | report CLI command |
                +--------------------+
```

## Data Flow
1. Inputs enter via `infer(adata, task_config, raw_10x_path, ambient_profile, evidence_table, candidate_set, tau, ...)`.
2. `config.py` validates task schema, applies defaults, and prepares model parameters.
3. Ambient path selection:
- If `ambient_profile` is provided, use it directly.
- Else if `raw_10x_path` is provided, run raw 10x ambient estimation.
- Else run fallback ambient estimation from low-UMI cells in `adata` and set fallback diagnostics.
4. `modules.py` computes class module features from sparse count matrix.
5. `evidence.py` aligns evidence rows by barcode and derives `kJ`, `lambda_amb_tcr`, and `amb_ll_tcr`.
6. `model.py` computes class scores, posterior probabilities, MAP label, entropy, max probability, and abstention-aware call.
7. `api.py` writes outputs to `adata.obs` and diagnostics/config snapshots to `adata.uns`.
8. `metrics.py` and `plotting.py` optionally compute calibration/coverage artifacts for reports.

## Module Interfaces

### `api.py`
- `infer(...)` orchestrates config validation, ambient selection, feature generation, model inference, and writeback.
- `estimate_ambient_profile(...)` wrapper to generate ambient profile and report from raw 10x.

### `ambient.py`
- `estimate_ambient_from_raw10x(raw_10x_path, min_umi, max_umi, bottom_pct, ...)`.
- `estimate_ambient_from_adata(adata, bottom_pct, ...)` for fallback path.
- Output: `ambient_profile_df`, `ambient_report_dict`.

### `modules.py`
- `compute_module_scores(adata, modules_config, candidate_mask=None)`.
- Output: per-class module score vector aligned to `adata.obs_names`.

### `evidence.py`
- `load_evidence_table(path_or_df, obs_names)`.
- `build_tcr_features(evidence_df, adata, ambient_profile, config)`.
- Output: aligned feature frame and evidence diagnostics.

### `model.py`
- `infer_posteriors(feature_df, beta_config, tau, qc_df=None)`.
- Output: posterior matrix, labels, uncertainty metrics.

### `cli.py`
- `estimate-ambient`: generate ambient profile artifacts.
- `infer`: run full workflow on an AnnData input.
- `report`: produce metrics/plots from inference output.

## Output Contract

### `adata.obs` required outputs
- `lineageresolver_label_map`
- `lineageresolver_label_call`
- `lineageresolver_max_p`
- `lineageresolver_uncertainty_entropy`
- `lineageresolver_p_<class>` for each class when requested

### Optional evidence/diagnostic outputs in `adata.obs`
- `lineageresolver_tcr_junction_umis`
- `lineageresolver_lambda_amb_tcr`
- `lineageresolver_amb_ll_tcr`
- `lineageresolver_module_<class>`
- `lineageresolver_a_hat`

### `adata.uns` outputs
- `lineageresolver_ambient_report`
- `lineageresolver_ambient_estimation_mode`
- `lineageresolver_model_config`
- `lineageresolver_runtime_diagnostics`
