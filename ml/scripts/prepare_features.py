#!/usr/bin/env python3
"""
ml/scripts/prepare_features.py

Applies the PPMI gene list (select_genes.py) unchanged to one cohort's normalized matrix
(normalize_cohort.py) and standardizes each gene within cohort and sex (spec story 13; the
frozen model's "within_cohort_sex" scaling rule). Uses no labels. Run once per cohort; a
cohort missing any listed gene stops the run, because the list must not be changed to fit it.

Output: participants x genes, in gene-list order (stays on Terra).

    python3 ml/scripts/prepare_features.py --expr results/ml/preprocessed/pdbp_expr.csv \\
        --participants results/ml/analysis_set/pdbp_participants.csv \\
        --genes results/ml/preprocessed/genes_top5000.txt \\
        --output results/ml/preprocessed/pdbp_features_top5000.csv
"""

import argparse
import os

import pandas as pd

from features import standardize
from select_genes import read_expr, read_gene_list


def prepare_features(expr: pd.DataFrame, meta: pd.DataFrame, genes: list,
                     scaling_rule: str = "within_cohort_sex") -> pd.DataFrame:
    """expr: participants x genes; meta: cohort and sex per participant."""
    missing = [g for g in genes if g not in expr.columns]
    if missing:
        raise ValueError(f"{len(missing)} of {len(genes)} listed genes are missing from the matrix, "
                         f"e.g. {missing[:5]}")
    cohorts = meta.loc[expr.index, "cohort"].unique()
    if len(cohorts) != 1:
        raise ValueError(f"Prepare one cohort at a time; got {sorted(cohorts)}")
    return standardize(expr[genes], meta, scaling_rule)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--expr", required=True, help="normalize_cohort.py output (genes x participants)")
    p.add_argument("--participants", required=True, help="analysis_set.py output for the same cohort")
    p.add_argument("--genes", required=True, help="select_genes.py gene list (from PPMI)")
    p.add_argument("--output", required=True, help="Features, participants x genes (stays on Terra)")
    a = p.parse_args()

    expr = read_expr(a.expr)
    meta = pd.read_csv(a.participants, usecols=["participant_id", "cohort", "sex"], index_col="participant_id")
    out = prepare_features(expr, meta, read_gene_list(a.genes))
    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    out.to_csv(a.output)
    print(f"{meta.loc[out.index, 'cohort'].iloc[0]}: {out.shape[0]} participants x {out.shape[1]} genes")
    print(f"Saved: {a.output}")


if __name__ == "__main__":
    main()
