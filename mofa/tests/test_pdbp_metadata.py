import unittest
import pandas as pd
import io
import os
import sys

# Add mofa/scripts to path if needed
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))

from prepare_pdbp_metadata import process_pdbp_metadata

class TestPDBPMetadata(unittest.TestCase):
    def test_process_pdbp_metadata(self):
        csv_data = """study,visit_month,participant_id,sex,age_at_baseline,diagnosis,mds_updrs_3,moca
PDBP,0.0,1001,Female,65.5,PD,22.0,26.0
PDBP,6.0,1001,Female,66.0,PD,25.0,25.0
PDBP,0.0,1002,Male,70.1,Ctrl,2.0,29.0
PPMI,0.0,2001,Male,62.0,PD,18.0,28.0
"""
        df = pd.read_csv(io.StringIO(csv_data))
        processed = process_pdbp_metadata(df)
        
        # Verify baseline filtering for PDBP only
        self.assertEqual(len(processed), 2)
        self.assertTrue((processed["study"] == "PDBP").all())
        self.assertTrue((processed["visit_month"] == 0.0).all())
        
        # Verify columns
        self.assertIn("participant_id", processed.columns)
        self.assertIn("sex", processed.columns)
        self.assertIn("diagnosis", processed.columns)
        self.assertIn("age_at_baseline", processed.columns)
        
        # Verify missing sex drop if any
        self.assertFalse(processed["sex"].isna().any())

if __name__ == "__main__":
    unittest.main()
