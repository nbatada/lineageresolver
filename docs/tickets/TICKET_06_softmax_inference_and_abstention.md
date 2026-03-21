# TICKET 06: Softmax Inference and Abstention

## Title
Implement class scoring, posterior inference, and abstention policy.

## Goal
Produce calibrated-ready posterior probabilities and abstention-aware calls from class scores.

## Scope
- Compute class scores from module and evidence features.
- Apply numerically stable softmax per cell.
- Derive `label_map`, `max_p`, entropy, and abstention-aware `label_call`.

## Non-Goals
- No external calibration training pipeline (optional follow-on).
- No advanced hierarchical models.

## Detailed Steps
1. Implement score computation using configured beta terms.
2. Implement stable softmax (log-sum-exp pattern).
3. Compute posterior summary columns and uncertainty entropy.
4. Apply abstention threshold `tau`.
5. Add optional hook for calibration transform if labels are provided.

## Acceptance Criteria
- Posterior probabilities sum to 1 for each evaluated cell.
- Abstention behavior matches threshold semantics.
- Numerical stability maintained for large-magnitude logits.

## Tests To Add
- `test_softmax_probability_normalization.py`
- `test_softmax_numerical_stability.py`
- `test_abstention_threshold_behavior.py`

## Files/Modules Expected To Change
- `lineageresolver/model.py`
- `lineageresolver/api.py`
- `tests/test_softmax_probability_normalization.py`
