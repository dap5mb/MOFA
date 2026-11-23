"""
ml/scripts/score.py

Scorer: the only entry point that applies a classifier to PDBP. It takes a frozen model and
a matrix and returns per-participant scores. It reads no labels and imports no estimator,
so it cannot fit anything.
"""

import numpy as np
import pandas as pd

from features import standardize
from frozen_model import FrozenModel


def score(model: FrozenModel, expr: pd.DataFrame, meta: pd.DataFrame) -> pd.Series:
    """
    Predicted probability of PD for each participant (rows of expr). meta needs only
    'cohort' and 'sex' (for the scaling rule); any other column is ignored.
    """
    missing = [g for g in model.genes if g not in expr.columns]
    if missing:
        raise ValueError(f"{len(missing)} model genes are missing from the matrix, e.g. {missing[:5]}")
    groups = meta.loc[expr.index, ["cohort", "sex"]]
    X = standardize(expr[list(model.genes)], groups, model.scaling_rule).to_numpy()
    logit = X @ np.asarray(model.coefficients) + model.intercept
    return pd.Series(1.0 / (1.0 + np.exp(-logit)), index=expr.index, name="score")
