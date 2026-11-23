import unittest
import numpy as np
import pandas as pd
import tempfile
import subprocess
import os
import sys

SCRIPTS = os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts'))
sys.path.insert(0, SCRIPTS)

from evaluate_factors import dedupe_metadata


def make_cohort(prefix, n, rng, female_only_effect=1.0):
    """Factor1: PD effect in women only. Factor2: PD effect in both sexes. Factor3: noise."""
    ids = [f"{prefix}{i}" for i in range(n)]
    sex = rng.choice(["Female", "Male"], n)
    dx = rng.choice(["Parkinson's Disease", "Control", "Prodromal"], n, p=[0.5, 0.4, 0.1])
    age = rng.normal(65, 8, n)
    pdx = (dx == "Parkinson's Disease").astype(float)
    female = (sex == "Female").astype(float)
    updrs = np.where(pdx == 1, rng.normal(20, 8, n), rng.normal(2, 2, n))
    z = pd.DataFrame({
        "participant_id": ids,
        "Factor1": female_only_effect * pdx * female + rng.normal(size=n),
        "Factor2": 0.8 * pdx + rng.normal(size=n),
        "Factor3": rng.normal(size=n),
        "group": sex,
    })
    meta = pd.DataFrame({"participant_id": ids, "sex": sex, "diagnosis": dx, "age_at_baseline": age,
                         "mds_updrs_part_iii_summary_score": updrs, "moca_total_score": rng.normal(27, 2, n)})
    # Duplicate rows for a few participants, with a conflicting UPDRS value (e.g. two exam records)
    extra = meta.iloc[:5].copy()
    extra["mds_updrs_part_iii_summary_score"] += 10
    return z, pd.concat([meta, extra], ignore_index=True)


class TestDedupe(unittest.TestCase):
    def test_keeps_first_row_per_participant(self):
        meta = pd.DataFrame({"participant_id": ["a", "a", "b"], "x": [1, 2, 3]})
        out = dedupe_metadata(meta)
        self.assertEqual(list(out["participant_id"]), ["a", "b"])
        self.assertEqual(int(out.loc[out.participant_id == "a", "x"].iloc[0]), 1)


class TestClinicalModels(unittest.TestCase):
    def test_detects_female_only_effect_and_replication(self):
        rng = np.random.default_rng(0)
        with tempfile.TemporaryDirectory() as tmp:
            paths = {}
            for cohort, n in (("ppmi", 800), ("pdbp", 600)):
                z, meta = make_cohort(cohort.upper(), n, rng)
                paths[cohort] = (os.path.join(tmp, f"{cohort}_z.csv"), os.path.join(tmp, f"{cohort}_meta.csv"))
                z.to_csv(paths[cohort][0], index=False)
                meta.to_csv(paths[cohort][1], index=False)
            out = os.path.join(tmp, "out")
            cmd = [sys.executable, os.path.join(SCRIPTS, "clinical_models.py"),
                   "--ppmi_z", paths["ppmi"][0], "--ppmi_metadata", paths["ppmi"][1],
                   "--pdbp_z", paths["pdbp"][0], "--pdbp_metadata", paths["pdbp"][1],
                   "--case_levels", "Parkinson's Disease", "--control_levels", "Control", "--output", out]
            run = subprocess.run(cmd, capture_output=True, text=True, cwd=SCRIPTS)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("more than one row", run.stdout)

            res = pd.read_csv(os.path.join(out, "clinical_models_long.csv"))
            n_ppmi = res[(res.cohort == "PPMI") & (res.analysis == "interaction_dx") & (res.factor == "Factor1")]["n"].iloc[0]
            self.assertLess(n_ppmi, 800)  # prodromal excluded, duplicates removed

            rep = pd.read_csv(os.path.join(out, "clinical_models_replication.csv")).set_index(["analysis", "factor", "sex"])
            # Factor1: PD effect only in women -> interaction replicates (negative: smaller in men)
            f1 = rep.loc[("interaction_dx", "Factor1", "all")]
            self.assertTrue(f1["replicated"])
            self.assertLess(f1["PPMI_estimate"], 0)
            self.assertTrue(rep.loc[("within_sex_dx", "Factor1", "Female"), "replicated"])
            self.assertFalse(rep.loc[("within_sex_dx", "Factor1", "Male"), "replicated"])
            # Factor2: same PD effect in both sexes -> main effect replicates, interaction doesn't
            self.assertTrue(rep.loc[("main_dx", "Factor2", "all"), "replicated"])
            self.assertFalse(rep.loc[("interaction_dx", "Factor2", "all"), "replicated"])
            # Factor3: noise
            self.assertFalse(rep.loc[("main_dx", "Factor3", "all"), "replicated"])

    def test_lists_labels_without_levels(self):
        rng = np.random.default_rng(1)
        with tempfile.TemporaryDirectory() as tmp:
            z, meta = make_cohort("P", 100, rng)
            z.to_csv(os.path.join(tmp, "z.csv"), index=False)
            meta.to_csv(os.path.join(tmp, "m.csv"), index=False)
            run = subprocess.run([sys.executable, os.path.join(SCRIPTS, "clinical_models.py"),
                                  "--ppmi_z", os.path.join(tmp, "z.csv"), "--ppmi_metadata", os.path.join(tmp, "m.csv"),
                                  "--output", os.path.join(tmp, "out")], capture_output=True, text=True, cwd=SCRIPTS)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("Prodromal", run.stdout)
            self.assertIn("Pass --case_levels", run.stdout)
            self.assertFalse(os.path.exists(os.path.join(tmp, "out")))


if __name__ == "__main__":
    unittest.main()
