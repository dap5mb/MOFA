import unittest
import pandas as pd
import numpy as np
import tempfile
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))

from evaluate_factors import load_and_merge_factors_metadata, evaluate_factor_associations

class TestEvaluateFactors(unittest.TestCase):
    def test_evaluate_factors(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create dummy Z csv with 20 samples to meet the minimum threshold
            z_df = pd.DataFrame({
                "participant_id": [f"P{i}" for i in range(20)],
                "Factor1": np.random.randn(20),
                "Factor2": np.random.randn(20)
            })
            z_path = os.path.join(tmpdir, "z.csv")
            z_df.to_csv(z_path, index=False)
            
            # Create dummy metadata csv with 20 samples
            meta_df = pd.DataFrame({
                "participant_id": [f"P{i}" for i in range(20)],
                "diagnosis": ["PD"]*10 + ["Ctrl"]*10,
                "age_at_baseline": np.random.randint(50, 80, 20)
            })
            meta_path = os.path.join(tmpdir, "meta.csv")
            meta_df.to_csv(meta_path, index=False)
            
            merged = load_and_merge_factors_metadata(z_path, meta_path)
            self.assertEqual(len(merged), 20)
            self.assertIn("diagnosis", merged.columns)
            self.assertIn("Factor1", merged.columns)
            
            out_dir = os.path.join(tmpdir, "output")
            evaluate_factor_associations(merged, out_dir, ["diagnosis", "age_at_baseline"])
            
            assoc_csv = os.path.join(out_dir, "factor_metadata_associations.csv")
            self.assertTrue(os.path.exists(assoc_csv))
            
            res_df = pd.read_csv(assoc_csv)
            self.assertGreater(len(res_df), 0)
            self.assertIn("Factor", res_df.columns)
            self.assertIn("Metadata", res_df.columns)
            self.assertIn("P-value", res_df.columns)

    def test_variance_replication_compares_per_view_and_group(self):
        from evaluate_factors import compare_variance_replication
        with tempfile.TemporaryDirectory() as tmpdir:
            ppmi = pd.DataFrame({
                "Factor1": [60.0, 58.0, 1.0, 2.0], "Factor2": [1.0, 0.5, 12.0, 11.0],
                "group": ["Female", "Male", "Female", "Male"],
                "view": ["RNA_seq", "RNA_seq", "CSF_Proteomics", "CSF_Proteomics"],
            })
            pdbp = ppmi.copy()
            pdbp[["Factor1", "Factor2"]] *= 0.9
            pdbp["view"] = pdbp["view"].str.lower()
            ppmi_path, pdbp_path = os.path.join(tmpdir, "ppmi.csv"), os.path.join(tmpdir, "pdbp.csv")
            ppmi.to_csv(ppmi_path, index=False)
            pdbp.to_csv(pdbp_path, index=False)

            compare_variance_replication(ppmi_path, pdbp_path, tmpdir)

            comp = pd.read_csv(os.path.join(tmpdir, "cohort_variance_replication_comparison.csv"))
            self.assertEqual(len(comp), 8)  # 2 views x 2 groups x 2 factors, nothing summed
            self.assertLessEqual(comp["PPMI_R2"].max(), 100)
            row = comp[(comp.view == "RNA_seq") & (comp.group == "Female") & (comp.Factor == "Factor1")].iloc[0]
            self.assertAlmostEqual(row["PDBP_R2"], 54.0)

    def test_variance_replication_rejects_mismatched_views(self):
        from evaluate_factors import compare_variance_replication
        with tempfile.TemporaryDirectory() as tmpdir:
            ppmi = pd.DataFrame({"Factor1": [1.0], "group": ["Female"], "view": ["View_0"]})
            pdbp = pd.DataFrame({"Factor1": [1.0], "group": ["Female"], "view": ["RNA_seq"]})
            ppmi.to_csv(os.path.join(tmpdir, "a.csv"), index=False)
            pdbp.to_csv(os.path.join(tmpdir, "b.csv"), index=False)
            with self.assertRaises(ValueError):
                compare_variance_replication(os.path.join(tmpdir, "a.csv"), os.path.join(tmpdir, "b.csv"), tmpdir)

if __name__ == "__main__":
    unittest.main()
