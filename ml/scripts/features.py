"""
ml/scripts/features.py

Label-free feature selection (PPMI only) and the scaling rule shared by trainer and scorer.

Matrices are participants x genes, indexed by participant_id. Metadata is indexed by
participant_id with 'cohort' ("PPMI" or "PDBP") and 'sex'.
"""

import numpy as np
import pandas as pd

TRAINING_COHORT = "PPMI"
SCALING_RULES = ("within_cohort_sex",)


def require_training_cohort(expr: pd.DataFrame, meta: pd.DataFrame, what: str) -> None:
    """Raise unless every participant in expr is a PPMI participant (the leakage guard)."""
    missing = expr.index.difference(meta.index)
    if len(missing):
        raise ValueError(f"{what}: {len(missing)} participants have no metadata")
    cohorts = meta.loc[expr.index, "cohort"]
    other = cohorts[cohorts != TRAINING_COHORT]
    if len(other):
        raise ValueError(f"{what} accepts {TRAINING_COHORT} participants only; got "
                         f"{len(other)} from {sorted(other.unique())}")


def select_top_variance(expr: pd.DataFrame, meta: pd.DataFrame, n_top: int) -> list:
    """The n_top most variable genes in PPMI. Uses no labels; refuses any non-PPMI participant."""
    require_training_cohort(expr, meta, "Feature selector")
    var = expr.var(axis=0)
    return list(var.nlargest(min(n_top, expr.shape[1])).index)


def standardize(expr: pd.DataFrame, meta: pd.DataFrame, rule: str = "within_cohort_sex") -> pd.DataFrame:
    """
    within_cohort_sex: z-score each gene within each cohort x sex group, so cohort offsets and
    sex mean differences do not reach the model (as in the MOFA projection). Genes constant
    within a group are set to 0.
    """
    if rule not in SCALING_RULES:
        raise ValueError(f"Unknown scaling rule {rule!r}; known: {SCALING_RULES}")
    if expr.isna().any().any():
        raise ValueError("Expression matrix has missing values")
    groups = meta.loc[expr.index, ["cohort", "sex"]]
    if groups.isna().any().any():
        raise ValueError("Every participant needs cohort and sex to be standardized")
    out = expr.groupby([groups["cohort"], groups["sex"]], group_keys=False).transform(
        lambda g: (g - g.mean()) / g.std(ddof=1))
    return out.replace([np.inf, -np.inf], np.nan).fillna(0.0).loc[expr.index]
