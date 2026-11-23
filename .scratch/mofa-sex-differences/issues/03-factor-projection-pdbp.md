# 03: Factor projection on single-modality and validation cohorts

**What to build:** Project the single-modality PPMI participants ($N \approx 1,200$) and the full baseline PDBP cohort ($N \approx 1,468$) into the pretrained PPMI factor space. Assess factor alignment, correlation, and reproducibility between the discovery and validation cohorts.

**Blocked by:** 02: Multi-group MOFA training on PPMI discovery cohort

**Status:** Done via `pdbp-projection` 01–03 (closed 2026-10-07)

- [x] ~~Project single-modality PPMI participants into the pretrained factor space.~~ Not needed: every PPMI participant is in training, so all have factor scores.
- [x] Project PDBP baseline multi-omics matrices into the PPMI factor space (1,248 participants).
- [x] Measure factor score consistency and variance explained across the PDBP validation cohort (per view and sex; validation report, Section 1).
- [x] Export consolidated participant-by-factor tables for discovery and validation cohorts (`latent_factors_Z_with_meta.csv` for PPMI, `pdbp_latent_factors_Z.csv` for PDBP).
