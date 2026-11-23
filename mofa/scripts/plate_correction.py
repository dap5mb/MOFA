"""
mofa/scripts/plate_correction.py

Python entry point to adjust_batch_effects.R (limma removeBatchEffect on plate), so code
outside the notebook can call the same correction without shelling out by hand.

include_diagnosis=True is MOFA's design (~ sex + diagnosis). include_diagnosis=False is the
label-free design (~ sex) used by the cross-cohort ML module: the diagnosis column is
dropped before the metadata is written, so the R script never sees the label.
"""

import os
import subprocess
import tempfile

import pandas as pd

R_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "adjust_batch_effects.R")


def correct_plate(expr: pd.DataFrame, metadata: pd.DataFrame,
                  include_diagnosis: bool = True, rscript: str = "Rscript") -> pd.DataFrame:
    """
    expr: genes x participants (log scale). metadata: indexed by participant_id, with
    'plate' and 'sex' (and 'diagnosis' when include_diagnosis). Returns the corrected
    matrix for the participants found in both, in metadata order.
    """
    required = ["plate", "sex"] + (["diagnosis"] if include_diagnosis else [])
    missing = [c for c in required if c not in metadata.columns]
    if missing:
        raise ValueError(f"Metadata is missing columns needed for plate correction: {missing}")
    meta = metadata[required].copy()
    meta.index.name = "participant_id"

    with tempfile.TemporaryDirectory() as tmp:
        expr_path = os.path.join(tmp, "expr.csv")
        meta_path = os.path.join(tmp, "meta.csv")
        out_path = os.path.join(tmp, "adjusted.csv")
        expr.to_csv(expr_path)
        meta.to_csv(meta_path)
        cmd = [rscript, R_SCRIPT, "-e", expr_path, "-m", meta_path, "-o", out_path]
        if not include_diagnosis:
            cmd.append("--no_diagnosis")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return pd.read_csv(out_path, index_col=0)
