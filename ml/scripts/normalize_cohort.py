#!/usr/bin/env python3
"""
ml/scripts/normalize_cohort.py

Label-free normalization of one cohort's whole-blood RNA-seq (spec stories 8-9; issue 04).
Each cohort is normalized on its own, with MOFA's code (imported through shared.py):

  1. Size factors (DESeq2 median-of-ratios, estimated within the cohort) and log1p, by
     preprocess_rnaseq. PPMI keeps every gene that passes MOFA's low-count filter (counts >= 10
     in >= 20% of PPMI participants), which is the "all expressed genes" set; no variance
     selection happens here. PDBP keeps all genes, because its genes come from the PPMI list.
  2. Plate correction by limma removeBatchEffect with the label-free design (~ sex), by
     correct_plate(include_diagnosis=False). --no_plate_correction skips it (sensitivity,
     story 10).

Diagnosis never reaches either step: the participants file is read without its diagnosis
columns, and the metadata handed to preprocess_rnaseq and correct_plate holds only
participant_id, study, sex and plate. The participants are the cohort's analysis set (PD or
control, from analysis_set.py); restricting to them is the cohort definition, not a design.

Output: genes x participants (log scale, centered per gene), stays on Terra.

    python3 ml/scripts/normalize_cohort.py --cohort PPMI --counts <mofa>/data/ppmi_rnaseq.csv \\
        --participants results/ml/analysis_set/ppmi_participants.csv \\
        --output results/ml/preprocessed/ppmi_expr.csv
"""

import argparse
import os
import tempfile

import pandas as pd

import shared

COHORTS = ("PPMI", "PDBP")
# The only participant columns preprocessing reads; diagnosis is not among them
PARTICIPANT_COLUMNS = ["participant_id", "cohort", "sex", "plate", "in_analysis_set"]
DESIGN_COLUMNS = ["sex", "plate"]


def read_participants(path: str) -> pd.DataFrame:
    """The analysis-set participants of a participants file, without any diagnosis column."""
    df = pd.read_csv(path, usecols=lambda c: c in PARTICIPANT_COLUMNS)
    missing = [c for c in PARTICIPANT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing columns {missing}")
    df = df[df["in_analysis_set"].astype(bool)].set_index("participant_id")
    return df.drop(columns="in_analysis_set")


def _count_genes(counts_path: str) -> int:
    return len(pd.read_csv(counts_path, usecols=[0]))


def normalize_cohort(counts_path: str, participants: pd.DataFrame, cohort: str,
                     plate_correction: bool = True) -> pd.DataFrame:
    """
    counts_path: raw counts, genes x participants (MOFA's aligned matrix). participants:
    indexed by participant_id with cohort, sex and plate (read_participants()). Returns
    genes x participants for the participants found in both.
    """
    if cohort not in COHORTS:
        raise ValueError(f"cohort must be one of {COHORTS}")
    if "diagnosis" in participants.columns:
        raise ValueError("Preprocessing must not see diagnosis; use read_participants()")
    other = participants.loc[participants["cohort"] != cohort, "cohort"]
    if len(other):
        raise ValueError(f"{cohort} normalization got {len(other)} participants from {sorted(other.unique())}")
    design = participants[DESIGN_COLUMNS]
    if design.isna().any().any():
        raise ValueError(f"{int(design.isna().any(axis=1).sum())} participants have no sex or plate")

    with tempfile.TemporaryDirectory() as tmp:
        meta_path = os.path.join(tmp, "metadata.csv")
        out_path = os.path.join(tmp, "log1p.csv")
        participants[["sex"]].assign(study=cohort).rename_axis("participant_id").to_csv(meta_path)
        # PPMI: n_top_genes = number of genes keeps every gene passing the low-count filter,
        # ordered by variance (MOFA's own path, without cutting at 5,000). PDBP: all genes.
        n_top = _count_genes(counts_path) if cohort == "PPMI" else None
        shared.preprocess_rnaseq(counts_file=counts_path, metadata_file=meta_path, output_file=out_path,
                                 visit_filter=0.0, cohort_name=cohort, n_top_genes=n_top)
        expr = pd.read_csv(out_path, index_col=0)

    expr = expr.loc[:, [p for p in participants.index if p in expr.columns]]
    if plate_correction:
        expr = shared.correct_plate(expr, design.loc[expr.columns], include_diagnosis=False)
    expr.index.name = "gene"
    return expr


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cohort", required=True, choices=COHORTS)
    p.add_argument("--counts", required=True, help="Raw counts, genes x participants")
    p.add_argument("--participants", required=True, help="analysis_set.py output for this cohort")
    p.add_argument("--output", required=True, help="Normalized matrix, genes x participants (stays on Terra)")
    p.add_argument("--no_plate_correction", action="store_true", help="Sensitivity analysis (story 10)")
    a = p.parse_args()

    participants = read_participants(a.participants)
    expr = normalize_cohort(a.counts, participants, a.cohort, plate_correction=not a.no_plate_correction)
    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    expr.to_csv(a.output)
    print(f"{a.cohort}: {expr.shape[0]} genes x {expr.shape[1]} participants "
          f"({len(participants)} in the analysis set); plate correction: {not a.no_plate_correction}")
    print(f"Saved: {a.output}")


if __name__ == "__main__":
    main()
