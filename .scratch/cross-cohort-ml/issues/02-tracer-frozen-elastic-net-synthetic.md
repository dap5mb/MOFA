# 02: Tracer bullet: frozen elastic-net model, PPMI → PDBP, on synthetic cohorts

**What to build:** The thinnest complete path of the analysis, run on two synthetic cohorts (one PPMI-like, one PDBP-like) with a planted PD signal: elastic-net logistic regression tuned by nested cross-validation in "PPMI" (outer 5, inner 10, class weights), refit on all "PPMI" and saved as a **frozen model**; a scorer that only reads frozen models scores "PDBP" and returns per-participant scores and AUC. A minimal Snakemake rule runs it end to end. This establishes both test seams from the spec.

**Blocked by:** None (can start immediately)

**Status:** done (2026-11-23)

- [x] The frozen model records the gene list, scaling rule, covariate definition, hyperparameters, coefficients and the internal nested-CV AUC.
- [x] The scorer cannot fit: it takes a frozen model plus a matrix and returns scores only.
- [x] Seam 2: changing PDBP values or labels leaves every frozen parameter unchanged; passing PDBP participants to the trainer or feature selector raises an error.
- [x] Seam 1 (first case): a signal planted in both synthetic cohorts gives PDBP AUC clearly above 0.5; a null synthetic dataset gives AUC near 0.5.
- [x] One Snakemake rule runs synthetic data → frozen model → PDBP scores → AUC table.
