# 02: Execute Task 1.3 (Configure Multi-Group MOFA)

**What to build:**
Configure the `mofapy2` architecture for the PPMI discovery cohort, defining the multi-group (sex) variable and view specifications.

**Blocked by:** 
01: Finalize Task 1.2 (Preprocess Discovery Matrices)

**Status:** Done (2026-10-06)

- [x] Configure `train_mofa.py` to correctly map sex groups and scale views.
- [x] Define model options (n_factors, seeds, iter) in `train_mofa.py`.

**Note (2026-10-07):** `train_mofa.build_and_train()` holds the settings: `groups = sex`, `scale_views=True`, 10 factors, ARD on factors and weights, `n_iter=1000` (was 20), and it saves Tau and AlphaZ. `config.yaml` still says 15 factors; see the validation report, Section 2.
