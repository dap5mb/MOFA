"""Seam 2: the train/score boundary (leakage guard)."""

import dataclasses
import inspect
import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts")))
import score as score_module
from features import select_top_variance
from frozen_model import FrozenModel
from run_tracer import run
from score import score
from synthetic import make_synthetic_cohorts
from train import train_elastic_net

SMALL_GRID = {"C": [0.01, 0.1, 1.0], "l1_ratio": [0.5]}


class TestTrainScoreBoundary(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ppmi_expr, cls.ppmi_meta, cls.pdbp_expr, cls.pdbp_meta = make_synthetic_cohorts("signal", seed=0)
        cls.model, cls.scores, _ = run(cls.ppmi_expr, cls.ppmi_meta, cls.pdbp_expr, cls.pdbp_meta,
                                       n_top=200, grid=SMALL_GRID)

    def test_changing_pdbp_leaves_frozen_model_unchanged(self):
        rng = np.random.default_rng(99)
        other_expr = pd.DataFrame(rng.normal(size=self.pdbp_expr.shape) * 5 + 3,
                                  index=self.pdbp_expr.index, columns=self.pdbp_expr.columns)
        other_meta = self.pdbp_meta.copy()
        other_meta["diagnosis"] = rng.permutation(other_meta["diagnosis"].to_numpy())
        model, scores, _ = run(self.ppmi_expr, self.ppmi_meta, other_expr, other_meta,
                               n_top=200, grid=SMALL_GRID)
        self.assertEqual(model.to_dict(), self.model.to_dict())
        self.assertFalse(np.allclose(scores, self.scores))  # the new PDBP was in fact scored

    def test_feature_selector_refuses_pdbp(self):
        mixed = pd.concat([self.ppmi_expr, self.pdbp_expr.iloc[:1]])
        meta = pd.concat([self.ppmi_meta, self.pdbp_meta])
        with self.assertRaisesRegex(ValueError, "PPMI participants only"):
            select_top_variance(mixed, meta, 100)
        with self.assertRaisesRegex(ValueError, "PPMI participants only"):
            select_top_variance(self.pdbp_expr, self.pdbp_meta, 100)

    def test_trainer_refuses_pdbp(self):
        mixed = pd.concat([self.ppmi_expr, self.pdbp_expr.iloc[:1]])
        meta = pd.concat([self.ppmi_meta, self.pdbp_meta])
        genes = list(self.model.genes)
        with self.assertRaisesRegex(ValueError, "PPMI participants only"):
            train_elastic_net(mixed, meta, genes, grid=SMALL_GRID)
        with self.assertRaisesRegex(ValueError, "PPMI participants only"):
            train_elastic_net(self.pdbp_expr, self.pdbp_meta, genes, grid=SMALL_GRID)

    def test_scorer_takes_no_labels_and_cannot_fit(self):
        self.assertEqual(list(inspect.signature(score).parameters), ["model", "expr", "meta"])
        self.assertFalse(hasattr(FrozenModel, "fit"))
        # No estimator is reachable from the scorer module
        self.assertFalse(any(getattr(v, "__module__", "").startswith("sklearn")
                             for v in vars(score_module).values()))
        # Labels are never read: dropping or scrambling diagnosis gives the same scores
        no_dx = self.pdbp_meta.drop(columns="diagnosis")
        pd.testing.assert_series_equal(score(self.model, self.pdbp_expr, no_dx), self.scores)

    def test_scoring_does_not_change_the_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "model.json")
            self.model.save(path)
            before = open(path, "rb").read()
            score(FrozenModel.load(path), self.pdbp_expr, self.pdbp_meta)
            self.assertEqual(open(path, "rb").read(), before)
        with self.assertRaises(dataclasses.FrozenInstanceError):
            self.model.intercept = 0.0
        with self.assertRaises(TypeError):
            self.model.hyperparameters["C"] = 1e6

    def test_scorer_refuses_matrix_missing_model_genes(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            score(self.model, self.pdbp_expr.drop(columns=self.model.genes[0]), self.pdbp_meta)


if __name__ == "__main__":
    unittest.main()
