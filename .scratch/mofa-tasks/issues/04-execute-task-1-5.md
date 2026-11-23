# 04: Execute Task 1.5 (Validate and Extract Factors)

**What to build:**
Evaluate convergence stability, select the optimal model, and extract latent factors ($Z$) and feature weights ($W$).

**Blocked by:** 
03: Execute Task 1.4 (Train Discovery MOFA Model)

**Status:** Done (Validated on Terra)

- [x] Audit ELBO plots to confirm model convergence.
- [x] Select the optimal model seed based on highest ELBO and stability.
- [x] Extract feature weights and latent factor coordinates for downstream analysis.
- [x] Perform metadata association audit to identify and isolate technical plate artifacts (documented in `docs/mofa-step1-validation-report.md`).

**End of Day Note:** Step 1 discovery phase (Tasks 1.2 – 1.5) is fully complete, tested on Terra, and documented. Next session will begin by defining tickets for Tasks 1.6 & 1.7 (PDBP Validation & Projection).

**Superseded (2026-10-07):** the Sep 25 model was trained for only 20 iterations, and its plate finding no longer applies after batch correction. The current main model is `ppmi_mofa_outputs_noplasma/model_seed_0.hdf5` (RNA-seq + CSF); see `docs/mofa-step1-validation-report.md`.
