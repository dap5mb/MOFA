#!/usr/bin/env python3
"""
ml/scripts/run_tracer.py

Thinnest complete path: PPMI -> feature selection -> nested-CV elastic net -> frozen model
-> PDBP scores -> AUC table. run() is the entry point; the CLI runs it on synthetic cohorts
and writes the frozen model (JSON), PDBP scores (CSV) and the AUC table (TSV).
"""

import argparse
import json
import os

import pandas as pd
from sklearn.metrics import roc_auc_score

from features import select_top_variance
from frozen_model import FrozenModel
from score import score
from synthetic import SCENARIOS, make_synthetic_cohorts
from train import labels, train_elastic_net


def run(ppmi_expr, ppmi_meta, pdbp_expr, pdbp_meta, n_top: int = 5000, grid: dict = None,
        seed: int = 0):
    """Returns (frozen model, PDBP scores, AUC table). PDBP reaches only the scorer."""
    genes = select_top_variance(ppmi_expr, ppmi_meta, n_top)
    model = train_elastic_net(ppmi_expr, ppmi_meta, genes, grid=grid, seed=seed)
    # Freeze for real: what is scored is what was written
    model = FrozenModel(**json.loads(json.dumps(model.to_dict())))

    scores = score(model, pdbp_expr, pdbp_meta[["cohort", "sex"]])

    y = labels(pdbp_meta, scores.index)
    auc = pd.DataFrame([
        {"cohort": "PPMI", "evaluation": "internal_nested_cv", "n": model.training["n"],
         "n_pd": model.training["n_pd"], "auc": model.internal_cv_auc["mean"]},
        {"cohort": "PDBP", "evaluation": "external_frozen_model", "n": len(y),
         "n_pd": int(y.sum()), "auc": float(roc_auc_score(y, scores))},
    ])
    return model, scores, auc


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario", choices=SCENARIOS, default="signal")
    p.add_argument("--model", required=True, help="Output frozen model (JSON)")
    p.add_argument("--scores", required=True, help="Output PDBP scores (CSV)")
    p.add_argument("--auc", required=True, help="Output AUC table (TSV)")
    p.add_argument("--n_ppmi", type=int, default=300)
    p.add_argument("--n_pdbp", type=int, default=300)
    p.add_argument("--n_genes", type=int, default=300)
    p.add_argument("--n_top", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()

    ppmi_expr, ppmi_meta, pdbp_expr, pdbp_meta = make_synthetic_cohorts(
        a.scenario, n_ppmi=a.n_ppmi, n_pdbp=a.n_pdbp, n_genes=a.n_genes, seed=a.seed)
    model, scores, auc = run(ppmi_expr, ppmi_meta, pdbp_expr, pdbp_meta, n_top=a.n_top, seed=a.seed)

    for path in (a.model, a.scores, a.auc):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    model.save(a.model)
    pd.concat([pdbp_meta.loc[scores.index, ["sex"]], scores], axis=1).to_csv(a.scores)
    auc.assign(scenario=a.scenario).to_csv(a.auc, sep="\t", index=False)
    print(auc.to_string(index=False))


if __name__ == "__main__":
    main()
