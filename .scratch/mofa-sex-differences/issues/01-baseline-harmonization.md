# 01: Baseline cohort and multi-omic matrix harmonization

**What to build:** Process raw AMP-PD clinical metadata and biospecimen assay files to extract baseline participants (`visit_month = 0.0`) for PPMI and PDBP cohorts. Normalize whole blood RNAseq counts with variance-stabilizing transformation, normalize plasma and CSF proteomics abundances, and generate aligned complete-case matrices along with single-modality feature tables.

**Blocked by:** None (can start immediately)

**Status:** Done via `mofa-tasks` 01 and `pdbp-projection` 01–02 (closed 2026-10-07)

- [x] Filter clinical metadata for PPMI and PDBP participants at baseline visit month 0.0.
- [x] Align participant IDs across RNAseq and proteomics assay tables.
- [x] Apply variance-stabilizing transformation to RNAseq count matrices and select top variable genes. Implemented as size-factor normalization + log1p, then the top 5,000 variable genes, with batch correction in `adjust_batch_effects.R`.
- [x] Normalize plasma and CSF proteomics abundance values. Zeros are treated as missing, with zero-based QC (2026-10-06). Plasma was later dropped from the main model.
- [x] ~~Export harmonized complete-case training matrices and single-modality projection tables.~~ Superseded: training uses all participants with MOFA's native missing-data handling (`step1_mofa_tasks_and_decisions.md`), so no separate projection tables are needed.
- [x] Verify harmonization outputs against sample count expectations. Actual PPMI counts: 1,329 participants; 1,308 RNA-seq, 354 CSF, 165 plasma; 158 complete cases for all three views.
