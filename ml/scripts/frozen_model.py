"""
ml/scripts/frozen_model.py

A frozen model: everything needed to score a cohort, fixed after training in PPMI.
It is plain data (no estimator object), saved as JSON, so whatever reads it can only
apply it, never refit it.
"""

import json
from dataclasses import asdict, dataclass, field
from types import MappingProxyType


@dataclass(frozen=True)
class FrozenModel:
    model_type: str                 # "elastic_net_logistic"
    genes: tuple                    # feature order; coefficients follow it
    scaling_rule: str               # "within_cohort_sex": each cohort standardized within sex
    covariates: tuple               # covariate definition (none in the tracer)
    hyperparameters: dict           # chosen by inner CV on all PPMI (C, l1_ratio)
    tuning: dict                    # grid, folds, class weights, seed
    coefficients: tuple             # one per gene, same order as genes
    intercept: float
    internal_cv_auc: dict           # nested-CV estimate in PPMI: mean, sd, per outer fold
    training: dict                  # cohort and counts
    notes: dict = field(default_factory=dict)

    def __post_init__(self):
        if len(self.genes) != len(self.coefficients):
            raise ValueError("genes and coefficients differ in length")
        # Freeze the nested containers too
        for name in ("genes", "covariates", "coefficients"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        for name in ("hyperparameters", "tuning", "internal_cv_auc", "training", "notes"):
            object.__setattr__(self, name, MappingProxyType(dict(getattr(self, name))))

    def to_dict(self) -> dict:
        d = {k: (dict(v) if isinstance(v, MappingProxyType) else v)
             for k, v in ((f, getattr(self, f)) for f in self.__dataclass_fields__)}
        for k in ("genes", "covariates", "coefficients"):
            d[k] = list(d[k])
        return d

    def save(self, path: str) -> None:
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=2, sort_keys=True)

    @classmethod
    def load(cls, path: str) -> "FrozenModel":
        with open(path) as fh:
            return cls(**json.load(fh))
