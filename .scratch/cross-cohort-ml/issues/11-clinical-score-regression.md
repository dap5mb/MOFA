# 11: Clinical-score regression

**What to build:** Elastic-net regression of MDS-UPDRS III and of MoCA in PD cases only, trained and frozen in PPMI and scored in PDBP with the same preprocessing and leakage guards, to see whether the signal tracks severity or only diagnosis.

**Blocked by:** 04

**Status:** ready-for-agent

- [ ] Uses the MDS-UPDRS III exam rule from ticket 03.
- [ ] Reports Spearman ρ with bootstrap intervals, internal (PPMI) and external (PDBP), overall and by sex.
