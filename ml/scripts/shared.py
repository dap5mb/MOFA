"""
ml/scripts/shared.py

The one place the cross-cohort ML module gets MOFA's preprocessing from. Nothing here is
reimplemented: each name is loaded from mofa/scripts on first use, so a missing optional
dependency (google-cloud-bigquery, pydeseq2, h5py) only matters for the function that needs it.

    from shared import correct_plate, preprocess_rnaseq
"""

import importlib
import os
import sys

MOFA_SCRIPTS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "mofa", "scripts"))
if MOFA_SCRIPTS not in sys.path:
    sys.path.insert(0, MOFA_SCRIPTS)

# exported name -> module in mofa/scripts that defines it
_EXPORTS = {
    # baseline metadata (BigQuery on Terra; MDS-UPDRS III off-medication exam rule)
    "query_master_clinical_metadata": "query_comprehensive_metadata",
    "create_master_metadata": "create_master_metadata",
    "process_pdbp_metadata": "prepare_pdbp_metadata",
    # within-cohort RNA-seq normalization (size factors + log1p; n_top_genes=None keeps all genes)
    "preprocess_rnaseq": "preprocess_rnaseq",
    # plate correction (include_diagnosis=False for the label-free design)
    "correct_plate": "plate_correction",
    # PDBP preparation
    "preprocess_pdbp_view": "preprocess_pdbp",
    "align_features_to_model": "preprocess_pdbp",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module 'shared' has no attribute {name!r}")
    return getattr(importlib.import_module(_EXPORTS[name]), name)
