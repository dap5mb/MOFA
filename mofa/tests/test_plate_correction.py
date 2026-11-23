import os
import shutil
import subprocess
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts")))
from plate_correction import correct_plate


def _limma_available():
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(["Rscript", "-e", "suppressPackageStartupMessages({library(limma); library(dplyr); library(ggplot2)})"],
                       capture_output=True)
    return r.returncode == 0


@unittest.skipUnless(_limma_available(), "Rscript with limma, dplyr and ggplot2 not available")
class TestPlateCorrection(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        n, g = 48, 25
        ids = [f"P{i}" for i in range(n)]
        self.meta = pd.DataFrame({
            "plate": rng.choice(["A", "B", "C"], n),
            "sex": rng.choice(["Male", "Female"], n),
            "diagnosis": rng.choice(["Case", "Control"], n),
        }, index=pd.Index(ids, name="participant_id"))
        expr = rng.normal(5, 1, (g, n))
        expr += (self.meta["plate"] == "B").to_numpy() * 1.5          # plate shift
        expr[:5] += (self.meta["diagnosis"] == "Case").to_numpy() * 1.0  # diagnosis signal
        self.expr = pd.DataFrame(expr, index=[f"G{i}" for i in range(g)], columns=ids)

        self.permuted = self.meta.copy()
        self.permuted["diagnosis"] = rng.permutation(self.meta["diagnosis"].to_numpy())

    def test_label_free_ignores_diagnosis(self):
        a = correct_plate(self.expr, self.meta, include_diagnosis=False)
        b = correct_plate(self.expr, self.permuted, include_diagnosis=False)
        pd.testing.assert_frame_equal(a, b, check_exact=True)

    def test_default_design_uses_diagnosis(self):
        a = correct_plate(self.expr, self.meta)
        b = correct_plate(self.expr, self.permuted)
        self.assertFalse(np.allclose(a.to_numpy(), b.to_numpy()))

    def test_label_free_differs_from_default(self):
        a = correct_plate(self.expr, self.meta)
        b = correct_plate(self.expr, self.meta, include_diagnosis=False)
        self.assertFalse(np.allclose(a.to_numpy(), b.to_numpy()))

    def test_label_free_runs_without_diagnosis_column(self):
        out = correct_plate(self.expr, self.meta.drop(columns="diagnosis"), include_diagnosis=False)
        self.assertEqual(out.shape, self.expr.shape)


if __name__ == "__main__":
    unittest.main()
