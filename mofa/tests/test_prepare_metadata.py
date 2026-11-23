import unittest
import pandas as pd
import io

def filter_baseline_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters clinical metadata for PPMI and PDBP participants at visit_month == 0.0.
    """
    valid_studies = {"PPMI", "PDBP"}
    filtered = df[(df["study"].isin(valid_studies)) & (df["visit_month"] == 0.0)].copy()
    return filtered

class TestPrepareMetadata(unittest.TestCase):
    def test_filter_baseline_metadata(self):
        csv_data = """study,visit_month,sex,participant_id
BioFIND,0.5,Female,92
PPMI,0.0,Male,101
PPMI,6.0,Male,101
PDBP,0.0,Female,201
HBS,0.0,Male,301
"""
        df = pd.read_csv(io.StringIO(csv_data))
        res = filter_baseline_metadata(df)
        self.assertEqual(len(res), 2)
        self.assertSetEqual(set(res["study"]), {"PPMI", "PDBP"})
        self.assertTrue((res["visit_month"] == 0.0).all())

if __name__ == "__main__":
    unittest.main()
