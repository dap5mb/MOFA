# 03: Execute Task 1.4 (Train Discovery MOFA Model)

**What to build:**
Execute the discovery MOFA model training loop across multiple random initializations to ensure stability.

**Blocked by:** 
02: Execute Task 1.3 (Configure Multi-Group MOFA)

**Status:** Done, with the scope changed (2026-10-06)

- [x] ~~Run n=20 random seed loop training.~~ Ran 5 seeds. More would add nothing: mofapy2 0.7.5 starts from a PCA of the data, so seeds give identical solutions (|r| = 1.0, ELBOs within 0.07). Stability is assessed by subsampling instead (`subsample_stability.py`).
- [x] Verify HDF5 artifact generation for each seed.
- [x] Generate ELBO convergence diagnostic plots for every seed. All seeds converged at 220 iterations; `select_best_model.py` now reports convergence and seed agreement.
