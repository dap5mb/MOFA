# 01: PDBP Metadata Harmonization & Clinical Preprocessing

**What to build:** A unified master sample manifest and clinical preprocessing pipeline for PDBP baseline participants (`visit_month = 0.0`), isolating biological sex, demographic variables, and available clinical scores (diagnosis, MDS-UPDRS, and auditing secondary clinical scales).

**Blocked by:** None (can start immediately)

**Status:** Completed (Verified on Terra & Local Tests) 

- [x] Extract and filter PDBP baseline records for `visit_month = 0.0`.
- [x] Harmonize demographic metadata (biological sex, age at baseline) and establish universal `patient_id` mapping.
- [x] Preprocess and harmonize clinical metadata (diagnosis `Ctrl` vs. `PD`, MDS-UPDRS scores, and auditing availability of secondary clinical measures like MoCA).
- [x] Output a clean master sample and clinical manifest ready for multi-omic alignment.
