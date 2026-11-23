# 01: Prefactor: shared preprocessing callable from the ML module

**What to build:** The ML module can reuse the MOFA metadata query, within-cohort RNA-seq normalization and PDBP preparation without copying them. Plate correction gains an option to leave diagnosis out of the design (needed so the classifier label never enters the features). MOFA's default behaviour and outputs are unchanged. Spec: step 2 (cross-cohort ML); ADR-0002.

**Blocked by:** None (can start immediately)

**Status:** done (2026-11-23)

- [x] The ML module imports the shared metadata, normalization and PDBP-prep functions from one place; no copies of that code exist under the ML module.
- [x] Plate correction can run with or without diagnosis in the design; the default (with diagnosis) reproduces current MOFA outputs exactly.
- [x] A test shows that with the label-free option, permuting diagnosis labels leaves the corrected matrix unchanged.
- [x] All MOFA tests pass, run one file at a time (packages missing locally, such as gseapy and matplotlib, excepted as before).
