#!/usr/bin/env python3
"""
ml/scripts/train_ppmi.py

Trains the tracer model (issue 02: elastic-net logistic regression, nested CV with outer 5 and
inner 10 folds, class weights, no covariates) on real PPMI features (prepare_features.py) and
writes the frozen model and its internal nested-CV AUC, overall and per sex (issue 04).
It reads PPMI only; the trainer refuses any other participant. Nothing here scores PDBP.

Outputs
  --model   frozen model (JSON): gene list, scaling rule, hyperparameters, coefficients, internal AUC
  --auc     internal AUC table (TSV, aggregate): all, Female, Male. auc_pooled_oof is the AUC of
            the out-of-fold scores pooled over the 5 outer folds; auc_fold_mean / _sd (overall
            row only) summarize the per-fold AUCs.
  --oof     out-of-fold PPMI scores (CSV, participant level, stays on Terra)

    python3 ml/scripts/train_ppmi.py --features results/ml/preprocessed/ppmi_features_top5000.csv \\
        --participants results/ml/analysis_set/ppmi_participants.csv \\
        --model results/ml/ppmi_elastic_net/frozen_model.json \\
        --auc results/ml/ppmi_elastic_net/internal_auc.tsv \\
        --oof results/ml/ppmi_elastic_net/ppmi_oof_scores.csv --n_jobs -1
"""

import argparse
import os

import pandas as pd

from train import train_elastic_net


def internal_auc_table(model) -> pd.DataFrame:
    cv = model.internal_cv_auc
    rows = [{"cohort": "PPMI", "evaluation": "internal_nested_cv", "sex": "all",
             "n": model.training["n"], "n_pd": model.training["n_pd"],
             "auc_pooled_oof": cv["pooled_oof"], "auc_fold_mean": cv["mean"], "auc_fold_sd": cv["sd"]}]
    for sex, r in cv["by_sex"].items():
        rows.append({"cohort": "PPMI", "evaluation": "internal_nested_cv", "sex": sex,
                     "n": r["n"], "n_pd": r["n_pd"], "auc_pooled_oof": r["auc_pooled_oof"],
                     "auc_fold_mean": float("nan"), "auc_fold_sd": float("nan")})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--features", required=True, help="PPMI prepare_features.py output")
    p.add_argument("--participants", required=True, help="PPMI analysis_set.py output")
    p.add_argument("--model", required=True)
    p.add_argument("--auc", required=True)
    p.add_argument("--oof", required=True)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n_jobs", type=int, default=None, help="Parallel inner grid search (-1 = all cores)")
    a = p.parse_args()

    expr = pd.read_csv(a.features, index_col=0)
    meta = pd.read_csv(a.participants, usecols=["participant_id", "cohort", "sex", "diagnosis"],
                       index_col="participant_id")
    model, oof = train_elastic_net(expr, meta, list(expr.columns), seed=a.seed, n_jobs=a.n_jobs,
                                   return_oof=True)

    for path in (a.model, a.auc, a.oof):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    model.save(a.model)
    table = internal_auc_table(model)
    table.to_csv(a.auc, sep="\t", index=False)
    pd.concat([meta.loc[oof.index, ["sex"]], oof], axis=1).to_csv(a.oof)
    print(table.to_string(index=False))
    print(f"Hyperparameters: {dict(model.hyperparameters)}; "
          f"non-zero coefficients: {sum(c != 0 for c in model.coefficients)} of {len(model.genes)}")
    print(f"Saved: {a.model}, {a.auc}, {a.oof}")


if __name__ == "__main__":
    main()
