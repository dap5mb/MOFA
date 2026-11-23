# 06: Single-group MOFA sensitivity check

**What to build:** A second PPMI MOFA model trained with a single group (no `groups = sex`), same views and preprocessing as the sex-grouped model, used as a sensitivity check for sex-differential PD factors. Rationale: the sex-grouped model centers features within each sex, which removes average sex differences; a single-group model keeps them, so sex x diagnosis effects in mean level can appear directly. See `docs/mofa-step1-validation-report.md` Section 3.

**Blocked by:** None (04's interaction test is implemented; use `clinical_models.py` so both models are evaluated identically)

**Status:** closed (2026-11-19)

- [x] Decide and document sex-chromosome feature handling (XIST, RPS4Y1, other chrY genes) before training, since they will form a pure sex factor in a single-group model. **Decision (2026-11-19): keep them**, consistent with not regressing out biological drivers before MOFA (decisions doc) and with no gene-to-chromosome annotation needed on Terra. The sex factor is identified afterwards as the factor with no counterpart in the sex-grouped model (`compare_models.py`) and a strong sex association.
- [x] (Trained 2026-11-19: converged at 415 iterations; subsampling stability still running) Train single-group PPMI MOFA on the main model's inputs (RNA-seq + CSF, `ppmi_mofa_inputs_noplasma/`), with the long table's `group` column set to one value, written to `ppmi_mofa_inputs_singlegroup/`. One seed (seeds are identical); `subsample_stability.py` for stability (subsets drawn from the whole cohort, since there's one group).
- [x] (Run 2026-11-19: no replicated interaction; the PD axis replicates in both sexes) Run the same `factor ~ diagnosis * sex + age` and `factor ~ mds_updrs_part_iii * sex + age` interaction tests (`clinical_models.py`, same settings as 04). PDBP is projected from `pdbp_mofa_inputs_noplasma/`; `project_pdbp.py` now ignores sex for grouping when the model has one group.
- [x] (Run 2026-11-19: single-group Factor 3 is the sex factor, with no counterpart; main Factor 7 (PD) ↔ single Factor 8, |r| = 0.94; main Factor 10 (aging) ↔ single Factor 10, |r| = 0.93; the main CSF factor isn't recovered, best |r| = 0.27; see report Section 3) Match factors to the sex-grouped model and report which sex-differential PD axes agree across both models. `compare_models.py` matches on factor-score |r| (both models share the same participants) and reports weight |r| and top-5% overlap per view, in both directions.

**Implementation (2026-11-19):** `scripts/compare_models.py` (new); `project_pdbp.py` (single-group models put every sample in their one group). Tests: `tests/test_compare_models.py`. Toy check with mofapy2: shared factors matched across the two models (|r| ≥ 0.98), and a planted sex-gene block formed a single-group-only factor (|r| = 0.06 with any grouped factor).

**Note (2026-10-07):** triage to `ready-for-agent` once 04's interaction test exists, so both models are evaluated with the same code.

- [x] Subsampling stability for the single-group model (2026-11-19): all 10 factors reproduced in 10 of 10 subsets (minimum |r| 0.88–1.00). In the report, Section 3.
- [x] (2026-11-19) PDBP rerun after the library-size normalization fix: no replicated interaction; single-group Factor 5 (↔ main Factor 4) now replicates as a PD axis, like main Factor 4. Results in `results/single_group/`.
