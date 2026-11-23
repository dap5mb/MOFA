"""
ml/scripts/train.py

Model trainer: elastic-net logistic regression tuned by nested cross-validation in PPMI
(outer 5 folds for the internal AUC, inner 10 folds for hyperparameters, class weights),
then refit on all PPMI and returned as a FrozenModel. Refuses any non-PPMI participant.

Hyperparameters use scikit-learn's names: l1_ratio is glmnet's alpha (L1 share), and
C is the inverse penalty strength (glmnet lambda = 1 / (C * n)).
"""

import numpy as np
import pandas as pd
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_score

from features import require_training_cohort, standardize
from frozen_model import FrozenModel

DEFAULT_GRID = {
    "C": [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0],
    "l1_ratio": [0.1, 0.5, 0.9],
}


def _elastic_net(seed: int) -> LogisticRegression:
    kw = dict(solver="saga", class_weight="balanced", max_iter=5000, random_state=seed)
    # scikit-learn < 1.8 needs penalty="elasticnet" for l1_ratio to apply; 1.8+ deprecates it
    if tuple(int(x) for x in sklearn.__version__.split(".")[:2]) < (1, 8):
        kw["penalty"] = "elasticnet"
    return LogisticRegression(**kw)


def labels(meta: pd.DataFrame, index) -> np.ndarray:
    """1 = PD, 0 = control."""
    dx = meta.loc[index, "diagnosis"]
    unknown = set(dx.unique()) - {"PD", "Control"}
    if unknown:
        raise ValueError(f"diagnosis must be 'PD' or 'Control'; found {sorted(map(str, unknown))}")
    return (dx == "PD").to_numpy().astype(int)


def train_elastic_net(expr: pd.DataFrame, meta: pd.DataFrame, genes: list,
                      scaling_rule: str = "within_cohort_sex", grid: dict = None,
                      outer_folds: int = 5, inner_folds: int = 10, seed: int = 0) -> FrozenModel:
    require_training_cohort(expr, meta, "Trainer")
    grid = {k: list(v) for k, v in (grid or DEFAULT_GRID).items()}

    X = standardize(expr[genes], meta, scaling_rule).to_numpy()
    y = labels(meta, expr.index)

    inner = StratifiedKFold(inner_folds, shuffle=True, random_state=seed)
    outer = StratifiedKFold(outer_folds, shuffle=True, random_state=seed + 1)
    search = GridSearchCV(_elastic_net(seed), grid, scoring="roc_auc", cv=inner)

    fold_auc = cross_val_score(search, X, y, cv=outer, scoring="roc_auc")
    search.fit(X, y)
    est = search.best_estimator_

    return FrozenModel(
        model_type="elastic_net_logistic",
        genes=list(genes),
        scaling_rule=scaling_rule,
        covariates=[],
        hyperparameters={k: float(v) for k, v in search.best_params_.items()},
        tuning={"grid": grid, "outer_folds": outer_folds, "inner_folds": inner_folds,
                "class_weight": "balanced", "scoring": "roc_auc", "seed": seed,
                "solver": "saga", "sklearn_version": sklearn.__version__},
        coefficients=[float(c) for c in est.coef_[0]],
        intercept=float(est.intercept_[0]),
        internal_cv_auc={"mean": float(np.mean(fold_auc)), "sd": float(np.std(fold_auc, ddof=1)),
                         "folds": [float(a) for a in fold_auc]},
        training={"cohort": "PPMI", "n": int(len(y)), "n_pd": int(y.sum()),
                  "n_control": int(len(y) - y.sum())},
    )
