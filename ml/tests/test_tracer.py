"""Seam 1 (first case): the whole tracer on synthetic cohorts, through run()."""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts")))
from frozen_model import FrozenModel
from run_tracer import run
from synthetic import make_synthetic_cohorts

SMALL_GRID = {"C": [0.01, 0.1, 1.0], "l1_ratio": [0.5]}


def auc_of(table, cohort):
    return float(table.loc[table["cohort"] == cohort, "auc"].iloc[0])


class TestTracerSignal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = make_synthetic_cohorts("signal", seed=0)
        cls.model, cls.scores, cls.auc = run(*cls.data, n_top=200, grid=SMALL_GRID)

    def test_planted_signal_transfers_to_pdbp(self):
        self.assertGreater(auc_of(self.auc, "PDBP"), 0.7)
        self.assertGreater(auc_of(self.auc, "PPMI"), 0.7)

    def test_scores_every_pdbp_participant(self):
        pdbp_meta = self.data[3]
        self.assertEqual(list(self.scores.index), list(pdbp_meta.index))
        self.assertTrue(((self.scores > 0) & (self.scores < 1)).all())

    def test_frozen_model_records_what_scoring_needs(self):
        m = self.model
        self.assertEqual(len(m.genes), 200)
        self.assertEqual(len(m.coefficients), len(m.genes))
        self.assertEqual(m.scaling_rule, "within_cohort_sex")
        self.assertEqual(tuple(m.covariates), ())
        self.assertEqual(set(m.hyperparameters), {"C", "l1_ratio"})
        self.assertEqual(m.tuning["outer_folds"], 5)
        self.assertEqual(m.tuning["inner_folds"], 10)
        self.assertEqual(m.tuning["class_weight"], "balanced")
        self.assertEqual(len(m.internal_cv_auc["folds"]), 5)
        self.assertEqual(m.training["cohort"], "PPMI")

    def test_frozen_model_round_trips_through_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "model.json")
            self.model.save(path)
            self.assertEqual(FrozenModel.load(path).to_dict(), self.model.to_dict())


class TestTracerNull(unittest.TestCase):
    def test_null_data_gives_chance_auc(self):
        data = make_synthetic_cohorts("null", n_pdbp=1000, seed=1)
        _, _, auc = run(*data, n_top=200, grid=SMALL_GRID)
        self.assertLess(abs(auc_of(auc, "PDBP") - 0.5), 0.08)


if __name__ == "__main__":
    unittest.main()
