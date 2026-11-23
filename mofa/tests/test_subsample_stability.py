import unittest
import numpy as np
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))

from subsample_stability import draw_subset, top_jaccard


class TestSubsampleStability(unittest.TestCase):
    def test_draw_subset_keeps_group_balance_and_is_reproducible(self):
        samples = [f"S{i}" for i in range(100)]
        groups = ["Female"] * 30 + ["Male"] * 70
        data = pd.DataFrame({"sample": samples * 2, "group": groups * 2, "view": ["A"] * 100 + ["B"] * 100})
        a = draw_subset(data, 0.8, np.random.default_rng(0))
        b = draw_subset(data, 0.8, np.random.default_rng(0))
        self.assertEqual(a, b)
        self.assertEqual(len(a), 80)
        self.assertEqual(len(set(a)), 80)
        sex = dict(zip(samples, groups))
        self.assertEqual(sum(sex[s] == "Female" for s in a), 24)

    def test_top_jaccard_is_sign_invariant(self):
        w = pd.Series(np.arange(1.0, 21.0), index=[f"f{i}" for i in range(20)])
        self.assertEqual(top_jaccard(w, -w, 0.1), 1.0)


if __name__ == "__main__":
    unittest.main()
