# 2. Locked PPMI-to-PDBP validation and the cross-sex transfer gap

* Status: Accepted
* Date: 2026-11-23

The gene-level predictive analysis trains and tunes every model in PPMI only, freezes it, and scores PDBP once per model named in a committed spec (`docs/specs/step2-cross-cohort-ml-spec.md`); anything added after PDBP is scored is labelled exploratory. Its primary estimand is the cross-sex transfer gap (within-sex minus cross-sex AUC in PDBP, training matched on N), judged against a pre-specified ±0.05 AUC equivalence margin with a three-way verdict (shared, sex-specific, inconclusive). We chose this over pooled training with random splits, which cannot separate overfitting from cohort shift, and over gene × sex interaction terms or stratified-vs-pooled comparisons, which are underpowered at gene level and cannot support a claim of *no* difference. PDBP was already used to replicate the MOFA factors but not to evaluate any classifier, so a locked plan keeps "externally validated" defensible.

## Considered options

- **Pooled PPMI + PDBP with leave-one-cohort-out**: rejected as primary; uses both cohorts for tuning, so PDBP is no longer independent.
- **Penalized gene × sex interaction terms**: rejected; few interactions would be selected and none could be tested for equivalence.
- **Stratified vs pooled model AUC**: kept as a secondary comparison only.

## Consequences

- PDBP scoring code reads only frozen model files; a test enforces that no PDBP sample influences feature selection, scaling parameters or hyperparameters.
- A null result is reported as "shared" only when the 90% interval lies inside the margin; otherwise it is "inconclusive", not "no difference".
