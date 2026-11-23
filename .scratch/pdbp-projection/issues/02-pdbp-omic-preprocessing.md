# 02: PDBP Omic Preprocessing & Feature Alignment

**What to build:** Preprocessed PDBP whole blood RNA-seq (VST and filtering) and CSF/plasma proteomics (log2 transformation and plate median centering), strictly subset and aligned to the exact feature set (rows) of the trained PPMI MOFA model weights ($W$).

**Blocked by:** 01: PDBP Metadata Harmonization & Clinical Preprocessing

**Status:** Completed (uses preprocess_pdbp.py) 

- [x] Apply variance-stabilizing transformation (VST) and lowly-expressed gene filtering to PDBP whole blood RNA-seq data.
- [x] Apply log2 transformation and plate batch/normalization to PDBP CSF and plasma proteomics data.
- [x] Subset and align PDBP omic matrices to precisely match the feature set (genes and proteins) present in the trained PPMI MOFA model weights ($W$).
- [x] Export validated, MOFA-ready PDBP omic matrices and long-format table aligned with the master participant manifest.
