# 03: Analysis set and PDBP field-coverage report

**What to build:** The baseline analysis set for both cohorts from AMP-PD 2023_v4, matching the MOFA participant definition and restricted to PD vs control, with the flags the later analyses need: genetic carrier status, disease duration (age at baseline minus age at diagnosis), PD medication use at baseline, and MDS-UPDRS III selected by the MOFA exam rule (off-medication first), plus MoCA. A coverage report says how complete duration and medication are in PDBP, and the duration bands and **PPMI-like subset** thresholds are fixed from it before any PDBP classifier score exists.

**Blocked by:** 01

**Status:** in progress (code built and tested locally on synthetic data, 2026-11-23; waiting for the Terra run)

- [ ] Participant counts by cohort, sex and diagnosis match the MOFA write-up, less PPMI's "other" diagnoses. *(Report compares against the write-up; check on Terra.)*
- [x] Carrier flags come from WGS-confirmed PD mutation and enrollment genetic status; how they were combined is documented. *(`analysis_set.py` docstring; to confirm on Terra: the enrollment vocabulary.)*
- [ ] Coverage report: share of PDBP PD cases with duration and with medication status, by sex. *(`coverage_report.py`; needs the Terra run.)*
- [ ] Duration bands and the PPMI-like-subset definition are written into the step-2 spec (dated addendum) before any PDBP scoring; the subset's size by sex is reported. *(Report lists candidate subset sizes; `analysis_set.FIXED_DEFINITIONS` stays `None`, and `require_fixed_definitions()` refuses to run until the addendum exists.)*
- [ ] Runs on Terra; only aggregate counts leave Terra.

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
