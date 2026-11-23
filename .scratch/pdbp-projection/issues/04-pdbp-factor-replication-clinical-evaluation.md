# 04: PDBP Factor Replication & Clinical Association Evaluation

**What to build:** An evaluation script (`evaluate_pdbp_factors.py`) that computes view-wise variance explained ($R^2$), evaluates structural factor replication, and runs linear/logistic regression between projected PDBP latent factors and preprocessed clinical metadata (diagnosis, MDS-UPDRS, and audited secondary scores).

**Blocked by:** 03: PDBP Latent Factor Projection

**Status:** Done (2026-10-07, verified on Terra)

- [x] ~~Implement `evaluate_pdbp_factors.py`~~. `evaluate_factors.py` evaluates both cohorts (associations plus variance-explained replication), so no separate script is needed.
- [x] Compute view-wise variance explained ($R^2$) and compare replication against PPMI discovery metrics (per view x sex group; implemented in `project_pdbp.py` / `evaluate_factors.py` on 2026-10-05; rerun on Terra for the main RNA-seq + CSF model on 2026-10-07).
- [x] Deduplicate the PPMI metadata merge before any regression (1,329 participants currently become 1,432 rows). `evaluate_factors.dedupe_metadata` keeps the first row per participant and reports which columns differ between a participant's rows. It's used by `evaluate_factors.py` and `clinical_models.py` (2026-10-07). On Terra: 103 PPMI participants had 2 rows, and only `mds_updrs_part_iii_summary_score` differed (90 participants); PDBP had no duplicates.
- [x] Confirm the `diagnosis` categories in PPMI and PDBP. PPMI: PD 753, Ctrl 555, Other 21; PDBP: PD 767, Ctrl 481. Models use PD vs Ctrl, with Other excluded (2026-10-07).
- [x] (Run 2026-10-07; no replicated interaction, see report Section 3) Run sex x diagnosis interaction tests per factor (`factor ~ diagnosis * sex + age`, `factor ~ mds_updrs_part_iii * sex + age`), one row per participant, BH-corrected, in PPMI and projected PDBP. See `docs/mofa-step1-validation-report.md` Section 3.
- [x] (Run 2026-10-07; Factor 7 is the replicated PD axis) Run linear and logistic regression models testing projected factors against verified clinical measures (diagnosis, MDS-UPDRS, and audited secondary scores).
- [x] Generate comparative replication plots and summary audit reports for validated factors (`figures/*_noplasma.png`, `docs/mofa-step1-validation-report.md`). Extend with the interaction-test results.

**Note (2026-10-07):** start with Factor 7 (PD axis) and Factor 10 (aging axis), the stable, replicated clinical factors of the main model, and include age as a covariate (Factor 7 also tracks age). Factor 6 isn't stable, and Factors 3 and 9 don't replicate.

**Implementation (2026-10-07):** `scripts/clinical_models.py`. Factors are z-scored per cohort, and clinical scores are z-scored within PD cases. Models per factor:
- `factor ~ pd * male + age`, plus the PD effect within each sex;
- `factor ~ score * male + age` and the within-sex score effect, for MDS-UPDRS III and MoCA, in PD cases only;
- `logit(pd) ~ factor + age + male`;
- `score ~ factor + age + male`, in PD cases.

p-values are BH-corrected across factors per cohort, analysis and term. A result counts as replicated when PPMI q < 0.05, PDBP p < 0.05 and the signs agree. Outputs: `clinical_models_long.csv`, `clinical_models_replication.csv`, `within_sex_diagnosis_effects.png`. Tests: `tests/test_clinical_models.py` (planted female-only and shared PD effects are recovered and replicated). Needs `statsmodels` on Terra.

- [x] (Rerun 2026-10-07: 162 PPMI participants had repeated exams, all on medication; kept off-medication for 137, state not recorded for 19, on medication for 6) Fix the repeated baseline MDS-UPDRS III exams: 90 PPMI participants have two different scores. `query_comprehensive_metadata.py` now keeps one exam per participant (off medication, then untreated, then state not recorded, then on medication). On Terra: rerun the query for PPMI, rebuild `metadata_m0_PPMI.csv` from notebook cell 4 then cell 5 (cell 5 alone won't overwrite the existing UPDRS column), then rerun `clinical_models.py`.
- [x] (Rerun 2026-10-07 for PPMI; PDBP had no duplicates; report Section 2 updated) Rerun `evaluate_factors.py` for PPMI (and PDBP, for consistency) on the main model after the exam fix. The PPMI associations in report Section 2 were computed on the 1,432-row merge, before deduplication, so their significance is inflated. Then update Section 2.

**Outcome (2026-10-07):** no factor's association with diagnosis, MDS-UPDRS III or MoCA differs by sex in a way that replicates (smallest PPMI interaction q = 0.21). Factor 7 is the replicated PD axis (lower in PD; PPMI −0.19, PDBP −0.46 log-odds per SD), with near-identical effects in women and men. Within PD, Factor 7 tracks MDS-UPDRS III and MoCA in PPMI only. Details: `docs/mofa-step1-validation-report.md` Sections 2–3.
