# 07: Variance-explained fixes, proteomics zero handling, retrain and stability (record)

**What was built:** corrections found while evaluating PDBP replication, a retrained discovery model, and a real stability check. Details and results are in `docs/mofa-step1-validation-report.md` (Oct 7 version).

**Blocked by:** None

**Status:** Done (2026-10-05 to 2026-10-07, verified on Terra)

- [x] Fix variance-explained extraction and comparison. The old code lost view names, summed R² across views and sexes (giving Factor 1 ≈ 120%), and used different R² formulas for PPMI and PDBP (`extract_factors.py`, `evaluate_factors.py`).
- [x] Make the PDBP projection match training: centering within each sex, view scaling, only observed features (no zero-filling), and the saved Tau and AlphaZ (`project_pdbp.py`).
- [x] Diagnose proteomics zeros as non-detects that dominated plasma variance (`plasma_diagnostics.py`). Treat zeros as missing, and add zero-based QC dropping samples with > 50% zeros and proteins with > 30% (`proteomics_qc.py`, `preprocess_pdbp.py`, notebook cell 25).
- [x] Train to convergence: `n_iter` 1,000 (was 20), with Tau and AlphaZ saved (`train_mofa.py`). All seeds converged at 220 iterations.
- [x] Seed checks: report convergence and seed agreement (`select_best_model.py`). Found that mofapy2 seeds are identical (PCA initialization).
- [x] Add subsampling stability: 10 × 80% subsets drawn within each sex; stable means reproduced with |r| ≥ 0.8 in ≥ 90% of subsets (`subsample_stability.py`).
- [x] Drop plasma from the main model. Even after the zero fix it explains 1.2–1.5% in PPMI and below 0 in PDBP. The main model is RNA-seq + CSF (`ppmi_mofa_outputs_noplasma/`), with the three-view model kept as a sensitivity analysis (`prep_utils.py` with plasma optional, plus new notebook cells).
- [x] Tests: the suite grew from 2 to 20 tests (`tests/`).

**Outcome:** CSF Factor 2 replicates closely (24–27% of CSF in both cohorts). No factor is specific to one sex. Factor 7 (PD axis) and Factor 10 (aging axis) replicate and are stable, and they are the inputs to issue 04's interaction test.
