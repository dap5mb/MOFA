import unittest
import pandas as pd
import numpy as np
import sys
import os

# Add mofa/scripts to path to import
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts")))
from prep_utils import run_preflight_check

class TestMOFAPrep(unittest.TestCase):
    def setUp(self):
        self.meta = pd.DataFrame({
            "participant_id": ["P1", "P2", "P3"],
            "sex": ["M", "F", None]
        }).set_index("participant_id")
        
        self.rna = pd.DataFrame(np.random.rand(5, 3), columns=["P1", "P2", "P3"], index=[f"G{i}" for i in range(5)])
        self.csf = pd.DataFrame(np.random.rand(5, 2), columns=["P1", "P2"], index=[f"P{i}" for i in range(5)])
        self.pla = pd.DataFrame(np.random.rand(5, 2), columns=["P1", "P2"], index=[f"P{i}" for i in range(5)])

    def test_preflight_filters_sex(self):
        cleaned_meta = run_preflight_check(self.rna, self.csf, self.pla, self.meta)
        self.assertNotIn("P3", cleaned_meta.index)
        self.assertEqual(len(cleaned_meta), 2)

    def test_validate_without_plasma_view(self):
        import tempfile
        from prep_utils import validate_and_intersect_views
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, long_df = validate_and_intersect_views(self.rna, self.csf, None, self.meta, output_dir=tmpdir)
            self.assertEqual(sorted(long_df["view"].unique()), ["CSF_Proteomics", "RNA_seq"])
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "ppmi_plasma_proteomics_mofa_ready.csv")))
            self.assertTrue(os.path.exists(os.path.join(tmpdir, "ppmi_mofa_input_long.csv")))

if __name__ == "__main__":
    unittest.main()
