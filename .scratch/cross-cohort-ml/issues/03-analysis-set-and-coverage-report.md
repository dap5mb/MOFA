# 03: Analysis set and PDBP field-coverage report

**What to build:** The baseline analysis set for both cohorts from AMP-PD 2023_v4, matching the MOFA participant definition and restricted to PD vs control, with the flags the later analyses need: genetic carrier status, disease duration (age at baseline minus age at diagnosis), PD medication use at baseline, and MDS-UPDRS III selected by the MOFA exam rule (off-medication first), plus MoCA. A coverage report says how complete duration and medication are in PDBP, and the duration bands and **PPMI-like subset** thresholds are fixed from it before any PDBP classifier score exists.

**Blocked by:** 01

**Status:** ready-for-agent

- [ ] Participant counts by cohort, sex and diagnosis match the MOFA write-up, less PPMI's "other" diagnoses.
- [ ] Carrier flags come from WGS-confirmed PD mutation and enrollment genetic status; how they were combined is documented.
- [ ] Coverage report: share of PDBP PD cases with duration and with medication status, by sex.
- [ ] Duration bands and the PPMI-like-subset definition are written into the step-2 spec (dated addendum) before any PDBP scoring; the subset's size by sex is reported.
- [ ] Runs on Terra; only aggregate counts leave Terra.
