# 01: Finalize Task 1.2 (Preprocess Discovery Matrices)

**What to build:**
Finalize the preprocessing of PPMI omic matrices by implementing the pre-flight integrity check and sample intersection logic to ensure metadata/omic ID alignment.

**Blocked by:** 
None (can start immediately)

**Status:** Done (Verified on Terra)

- [x] Implement `run_preflight_check` to filter missing 'sex' samples and align IDs across views.
- [x] Update `mofa_wf.ipynb` to use the refactored `prep_utils.py` workflow.
- [x] Validate alignment output (logs confirm intersection counts).

**Closing Note:** Verified on Terra. Intersection size: 1329 participants (160 complete cases). Artifacts successfully exported.
