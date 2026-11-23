# 03: Analysis set and PDBP field-coverage report

**What to build:** The baseline analysis set for both cohorts from AMP-PD 2023_v4, matching the MOFA participant definition and restricted to PD vs control, with the flags the later analyses need: genetic carrier status, disease duration (age at baseline minus age at diagnosis), PD medication use at baseline, and MDS-UPDRS III selected by the MOFA exam rule (off-medication first), plus MoCA. A coverage report says how complete duration and medication are in PDBP, and the duration bands and **PPMI-like subset** thresholds are fixed from it before any PDBP classifier score exists.

**Blocked by:** 01

**Status:** done (2026-11-23)

- [x] Participant counts by cohort, sex and diagnosis match the MOFA write-up, less PPMI's "other" diagnoses. *(PDBP matches exactly: 1,248, i.e. 767 PD and 481 controls. PPMI: 1,308 RNA-seq participants, 753 PD and 555 controls, the same as the write-up; none of the 21 "other" had RNA-seq.)*
- [x] Carrier flags come from WGS-confirmed PD mutation and enrollment genetic status; how they were combined is documented. *(`analysis_set.py` docstring. The flag rates look implausible; see the notes below before the carrier sensitivity analysis.)*
- [x] Coverage report: share of PDBP PD cases with duration and with medication status, by sex. *(Duration 98.6% of women, 97.5% of men; medication 96.1% and 98.8%.)*
- [x] Duration bands and the PPMI-like-subset definition are written into the step-2 spec (dated addendum) before any PDBP scoring; the subset's size by sex is reported. *(Addendum 2026-11-23: bands 0–1, 2–4, 5–9 and ≥ 10 years; PPMI-like = PD with duration ≤ 2 years and no medication criterion, which gives 99 women and 144 men, plus all controls. Coded as `analysis_set.FIXED_DEFINITIONS`, `duration_band()` and `ppmi_like()`.)*
- [x] Runs on Terra; only aggregate counts leave Terra.

## Terra run (from the repository root; participant-level files stay under `results/`)

```
python3 ml/scripts/query_ml_fields.py --cohort PPMI --output results/ml/analysis_set/ppmi_ml_fields.csv
python3 ml/scripts/query_ml_fields.py --cohort PDBP --output results/ml/analysis_set/pdbp_ml_fields.csv
python3 ml/scripts/analysis_set.py --cohort PPMI --mofa_metadata <mofa>/data/metadata_m0_PPMI.csv \
    --ml_fields results/ml/analysis_set/ppmi_ml_fields.csv --rnaseq_matrix <mofa>/data/ppmi_rnaseq_adjusted.csv \
    --output results/ml/analysis_set/ppmi_participants.csv
python3 ml/scripts/analysis_set.py --cohort PDBP --mofa_metadata <mofa>/data/metadata_m0_PDBP.csv \
    --ml_fields results/ml/analysis_set/pdbp_ml_fields.csv --rnaseq_matrix <mofa>/data/pdbp_rnaseq_adjusted.csv \
    --output results/ml/analysis_set/pdbp_participants.csv
python3 ml/scripts/coverage_report.py \
    --participants results/ml/analysis_set/ppmi_participants.csv results/ml/analysis_set/pdbp_participants.csv \
    --output_dir results/ml/analysis_set/coverage_report
```

Only `coverage_report/` (counts, percentages, percentiles) leaves Terra.

## Notes from the build (2026-11-23)

- Field locations come from the AMP-PD schema (`data/amp_pd_schema.json`, 2020_v2; to confirm in 2023_v4): `age_at_diagnosis` and medication fields in tier-1 `PD_Medical_History`; `genetic_status_enrollment` in tier-2 `Clinically_Reported_Genetic_Status`. `has_known_PD_Mutation_in_WGS`, MDS-UPDRS III (off-medication rule) and MoCA come from MOFA's own query, which is already merged into `metadata_m0_<cohort>.csv`.
- Rules decided in the build, open to change before the addendum: ages are whole years, so a duration of −1 is read as 0 and anything below −1 is set to missing. Medication comes from the medical-history fields first, then the MDS-UPDRS 3a answer of the selected exam. A carrier is a participant positive in either source.
- Any genetic-status or yes/no value the code does not recognize stops the run with the list of values, so the first Terra run may need its vocabulary extended.
- PPMI reference: the write-up's 753 PD / 555 control / 21 other counts cover all 1,329 participants with any view; the RNA-seq analysis set (1,308) is compared against them as is.

## Findings from the Terra report (2026-11-23)

- **Hardly any PDBP PD cases are unmedicated:** 13 women and 12 men. Requiring "no medication" leaves 13 women and 9 men in the PPMI-like subset at ≤ 2 years, too few for a per-sex AUC. Reporting AUC by medication use within PDBP (story 33) has the same problem.
- **PDBP duration bands:** the candidate bands (0–1, 2–4, 5–9 and ≥ 10 whole years) each hold at least 57 PD cases per sex.
- **Carrier flags look wrong:**
  - `has_known_PD_Mutation_in_WGS` is "yes" for 31% of PDBP participants, with almost the same rate in controls (28%) as in PD (33%).
  - It is also "yes" for 67% of PPMI participants, including 243 of 290 female controls.
  - PDBP has no `genetic_status_enrollment` at all.
  - PPMI's genetic cohorts may explain part of the PPMI figure, since unaffected carriers count as controls, but not the PDBP figure.
  - Before issue 09 (carrier-excluded sensitivity), check on Terra: the raw value counts of `has_known_PD_Mutation_in_WGS` and `genetic_status_enrollment`; the field's definition in the 2023_v4 schema; and how it agrees with `genetic_status_wgs` and the per-gene `has_known_{GBA,LRRK2,SNCA}_mutation_in_WGS` flags. If PPMI controls are mostly unaffected carriers, the write-up should say so.
- **PPMI MoCA is mostly missing at baseline:** it is recorded for 58 of 316 PD women and 64 of 437 PD men, against 99% in PDBP. PPMI probably records MoCA at screening rather than at `visit_month = 0`. This affects issue 11 (MoCA regression trained in PPMI), and MOFA's query uses the same rule.
- **PPMI baseline medication:** 44% of PD women and 27% of PD men are medicated, mostly according to the MDS-UPDRS 3a answer.

## Follow-up 2026-11-23: screening-visit MoCA

`query_ml_fields.py` now returns MoCA from baseline, else from the latest screening visit before baseline (`visit_month < 0`), with the visit recorded. `analysis_set.py` uses this score in place of MOFA's baseline-only one, and the coverage report counts how many scores came from screening (`moca_from_screening`). MOFA's own query is unchanged. To finish: rerun the query, analysis-set and report steps on Terra and check PPMI MoCA coverage. If coverage stays low, screening visits may be stored without a `visit_month`; check `visit_name` values in `MOCA`.
