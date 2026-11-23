# 04: Label-free preprocessing and the first real-data run (PPMI internal AUC only)

**What to build:** Real RNA-seq for both cohorts, prepared with the leakage guards: within-cohort size-factor normalization and log transform, plate correction without diagnosis, the 5,000 most variable genes chosen in PPMI only and applied unchanged to PDBP, and standardization within cohort and sex. The tracer model from 02 is run on real PPMI and reports its internal nested-CV AUC. PDBP is prepared but **not scored** in this ticket.

**Blocked by:** 01, 02, 03

**Status:** ready-for-agent

- [ ] The gene list is derived from PPMI alone; a check confirms no PDBP sample entered its selection.
- [ ] Diagnosis does not appear in any preprocessing design.
- [ ] The overlap with MOFA's 5,000-gene list is reported (expected to be high).
- [ ] Internal PPMI nested-CV AUC is reported, overall and per sex; no PDBP classifier scores exist yet.
