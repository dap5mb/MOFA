# 03: PDBP Latent Factor Projection

**What to build:** A dedicated Python projection script (`project_pdbp.py`) that loads the trained PPMI HDF5 model and preprocessed PDBP omic matrices, leveraging `mofapy2`'s native missing data handling to compute and save validation latent factor coordinates ($Z_{\text{PDBP}}$) across sex groups.

**Blocked by:** 02: PDBP Omic Preprocessing & Feature Alignment

**Status:** Completed (Verified on Terra)

- [x] Implement `project_pdbp.py` to load the trained PPMI MOFA HDF5 model and PDBP preprocessed matrices.
- [x] Utilize `mofapy2` native missing data handling to project single-view and multi-view PDBP subjects into the PPMI factor space.
- [x] Extract, validate, and save PDBP sample latent factor coordinates ($Z_{\text{PDBP}}$) for male and female groups.
- [x] Verify output factor shapes and participant ID alignment against the master manifest.
