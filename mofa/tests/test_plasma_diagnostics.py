import unittest
import pandas as pd
import numpy as np
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))

from plasma_diagnostics import process, zero_share_of_variance, pca, variance_captured


class TestPlasmaDiagnostics(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.raw = pd.DataFrame(2 ** rng.normal(17, 1, size=(20, 50)),
                                index=[f"P{i}" for i in range(20)], columns=[f"S{j}" for j in range(50)])
        self.raw.iloc[0, :5] = 0.0  # non-detects

    def test_zeros_become_outliers_only_in_current_treatment(self):
        current = process(self.raw, "current")
        masked = process(self.raw, "zero_as_nan")
        self.assertLess(current.iloc[0, :5].max(), -10)
        self.assertTrue(masked.iloc[0, :5].isna().all())
        self.assertGreater(zero_share_of_variance(self.raw, current), 10)

    def test_own_pcs_capture_their_variance(self):
        x = process(self.raw, "zero_as_nan")
        var_pct, loadings, _ = pca(x, 3)
        self.assertAlmostEqual(variance_captured(x, loadings), float(np.sum(var_pct)), places=6)


if __name__ == "__main__":
    unittest.main()
