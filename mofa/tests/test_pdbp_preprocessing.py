import unittest
import pandas as pd
import numpy as np
import io
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))

from preprocess_pdbp import align_features_to_model

class TestPDBPPreprocessing(unittest.TestCase):
    def test_align_features_to_model(self):
        # Mock PDBP matrix with extra features and missing features
        pdbp_df = pd.DataFrame({
            "PDBP_subj1": [10.0, 20.0, 30.0, 40.0],
            "PDBP_subj2": [12.0, 22.0, 32.0, 42.0]
        }, index=["GeneA", "GeneB", "GeneC", "GeneX"])
        
        # Target model features (from PPMI training)
        model_features = ["GeneA", "GeneC", "GeneY"]
        
        aligned = align_features_to_model(pdbp_df, model_features, view_type="RNA_seq")
        
        # Should contain GeneA, GeneC, and GeneY (with NaN for GeneY since it's missing in PDBP)
        self.assertIn("GeneA", aligned.index)
        self.assertIn("GeneC", aligned.index)
        self.assertIn("GeneY", aligned.index)
        self.assertNotIn("GeneX", aligned.index)
        self.assertEqual(len(aligned), 3)

    def test_proteomics_zeros_become_missing_and_bad_samples_dropped(self):
        from preprocess_pdbp import preprocess_pdbp_view
        rng = np.random.default_rng(0)
        raw = pd.DataFrame(2 ** rng.normal(17, 1, size=(4, 6)),
                           index=["P1", "P2", "P3", "P4"], columns=[f"S{i}" for i in range(6)])
        raw.loc["P1", "S0"] = 0.0           # isolated non-detect
        raw.loc[["P1", "P2", "P3"], "S5"] = 0.0  # 75% zeros -> fails sample QC

        out = preprocess_pdbp_view(raw, "Plasma_Proteomics", ["P1", "P2", "P3", "P4"])

        self.assertNotIn("S5", out.columns)
        self.assertTrue(np.isnan(out.loc["P1", "S0"]))
        self.assertGreater(out.min().min(), -5)  # no log2(0 + 1) outliers

        kept = preprocess_pdbp_view(raw, "Plasma_Proteomics", ["P1", "P2", "P3", "P4"],
                                    zeros_as_missing=False, max_sample_zero_frac=None)
        self.assertIn("S5", kept.columns)
        self.assertLess(kept.loc["P1", "S0"], -10)


class TestProteomicsQC(unittest.TestCase):
    def test_filter_drops_samples_then_proteins(self):
        from proteomics_qc import filter_proteomics
        raw = pd.DataFrame(1000.0, index=["A", "B", "C"], columns=["s1", "s2", "s3", "s4"])
        raw.loc[["A", "B"], "s4"] = 0.0  # s4: 67% zeros -> dropped
        raw.loc["C", ["s1", "s2"]] = 0.0  # C: 50% zeros among s1-s3 -> > 30%, dropped
        out = filter_proteomics(raw, max_sample_zero_frac=0.50, max_protein_zero_frac=0.30)
        self.assertEqual(list(out.columns), ["s1", "s2", "s3"])
        self.assertEqual(list(out.index), ["A", "B"])

if __name__ == "__main__":
    unittest.main()
