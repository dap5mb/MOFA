# 02: Multi-group MOFA training on PPMI discovery cohort

**What to build:** Implement and execute multi-group multi-omics factor analysis (`MOFA2` / `mofapy2`) on the PPMI complete-case baseline dataset with `groups = sex`. Train factor representations to isolate variation shared between sexes from sex-specific variation, verifying convergence and stability across random seeds.

**Blocked by:** 01: Baseline cohort and multi-omic matrix harmonization

**Status:** Done via `mofa-tasks` 02–04 and the Oct 6–7 retrain (closed 2026-10-07)

- [x] Instantiate MOFA model with RNAseq and proteomics views using sex as the grouping variable.
- [x] Train the multi-group MOFA model on the PPMI baseline matrix. All participants were used rather than complete cases only (see decisions doc).
- [x] Evaluate model convergence across multiple random seeds and compute variance explained per view and per sex. Seeds are identical in mofapy2, so stability comes from subsampling instead.
- [x] Characterize factors into shared versus sex-specific categories based on group variance decomposition. Result: all factors are shared, with similar R² in women and men in both cohorts. The sex question moves to the sex × diagnosis interaction test (`pdbp-projection` 04).
- [x] Export the trained MOFA model object and factor weight matrices.
