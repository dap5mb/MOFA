"""Issue 03: analysis set and coverage report, through build_participants() and coverage_report()."""

import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts")))
import analysis_set
from analysis_set import analysis_set as pd_vs_control
from analysis_set import build_participants
from coverage_report import coverage_report, load, write_report


def mofa_meta(rows):
    """MOFA metadata rows: (participant_id, sex, diagnosis label, extra dict)."""
    base = {"study": "PPMI", "visit_month": 0.0, "age_at_baseline": 60, "rnaseq_sample_id": "s",
            "has_known_PD_Mutation_in_WGS": "No", "mds_updrs_part_iii_summary_score": 20,
            "moca_total_score": 27, "updrs3_medication_for_pd": np.nan, "plate": "P1",
            "diagnosis_at_baseline": "Parkinson's Disease"}
    out = []
    for pid, sex, dx, extra in rows:
        r = dict(base, participant_id=pid, sex=sex, diagnosis=dx)
        r.update(extra)
        out.append(r)
    return pd.DataFrame(out)


def ml_fields(rows):
    """ML-field rows: (participant_id, extra dict)."""
    base = {"age_at_diagnosis": 59, "use_of_pd_medication": "No", "on_levodopa": np.nan,
            "on_dopamine_agonist": np.nan, "on_other_pd_medications": np.nan,
            "medication_visit_name": "M0", "genetic_status_enrollment": np.nan, "genetic_status_wgs": np.nan}
    return pd.DataFrame([dict(base, participant_id=pid, **extra) for pid, extra in rows])


class TestParticipantDefinition(unittest.TestCase):
    def setUp(self):
        self.meta = mofa_meta([
            ("A", "Female", "PD", {}),
            ("A", "Female", "PD", {"plate": "P2"}),              # duplicate: first kept
            ("B", "Male", "Control", {}),
            ("C", "Male", "Prodromal", {}),                      # Other: kept but out of the set
            ("D", "Female", "PD", {"visit_month": 12.0}),        # not baseline
            ("E", np.nan, "PD", {}),                             # no recorded sex
            ("F", "Male", "PD", {"rnaseq_sample_id": np.nan}),   # no RNA-seq
            ("G", "Female", "PD", {"study": "PDBP"}),            # other cohort
        ])
        self.fields = ml_fields([(p, {}) for p in "ABCDEFG"])

    def test_mofa_rnaseq_participants_one_row_each(self):
        df = build_participants(self.meta, self.fields, "PPMI")
        self.assertEqual(sorted(df.index), ["A", "B", "C"])
        self.assertEqual(df.loc["A", "plate"], "P1")
        self.assertEqual(df.loc["C", "diagnosis"], "Other")
        self.assertEqual(sorted(pd_vs_control(df).index), ["A", "B"])
        self.assertEqual(set(pd_vs_control(df)["diagnosis"]), {"PD", "Control"})

    def test_rnaseq_matrix_ids_define_membership(self):
        df = build_participants(self.meta, self.fields, "PPMI", rnaseq_ids=["B", "A"])
        self.assertEqual(list(df.index), ["B", "A"])
        with self.assertRaises(ValueError):
            build_participants(self.meta, self.fields, "PPMI", rnaseq_ids=["A", "Z"])

    def test_pdbp_diagnosis_comes_from_diagnosis_at_baseline(self):
        meta = mofa_meta([
            ("P1", "Female", "Ctrl", {"study": "PDBP"}),
            ("P2", "Male", "PD", {"study": "PDBP", "diagnosis_at_baseline": "No PD Nor Other Neurological Disorder"}),
            ("P3", "Male", "PD", {"study": "PDBP", "diagnosis_at_baseline": "Essential Tremor"}),
        ])
        df = build_participants(meta, ml_fields([("P1", {}), ("P2", {}), ("P3", {})]), "PDBP")
        self.assertEqual(df["diagnosis"].to_dict(), {"P1": "PD", "P2": "Control", "P3": "Other"})


class TestDerivedFlags(unittest.TestCase):
    def build(self, meta_extra, field_extra, dx="PD"):
        meta = mofa_meta([(f"X{i}", "Female", dx, m) for i, m in enumerate(meta_extra)])
        fields = ml_fields([(f"X{i}", f) for i, f in enumerate(field_extra)])
        return build_participants(meta, fields, "PPMI")

    def test_carrier_combines_wgs_and_enrollment(self):
        df = self.build(
            [{"has_known_PD_Mutation_in_WGS": "Yes"}, {"has_known_PD_Mutation_in_WGS": "No"},
             {"has_known_PD_Mutation_in_WGS": "No"}, {"has_known_PD_Mutation_in_WGS": np.nan},
             {"has_known_PD_Mutation_in_WGS": np.nan}],
            [{"genetic_status_enrollment": "No Mutation"}, {"genetic_status_enrollment": "LRRK2+"},
             {}, {}, {"genetic_status_enrollment": "Not tested"}])
        self.assertEqual([None if pd.isna(v) else bool(v) for v in df["carrier"]],
                         [True, True, False, None, None])

    def test_unrecognized_vocabulary_stops_the_run(self):
        with self.assertRaisesRegex(ValueError, "genetic_status_enrollment"):
            self.build([{}], [{"genetic_status_enrollment": "pending review"}])
        with self.assertRaisesRegex(ValueError, "use_of_pd_medication"):
            self.build([{}], [{"use_of_pd_medication": "sometimes"}])

    def test_duration_pd_only_with_rounding_rule(self):
        df = self.build(
            [{"age_at_baseline": 65}, {"age_at_baseline": 65}, {"age_at_baseline": 65}, {"age_at_baseline": 65}],
            [{"age_at_diagnosis": 60}, {"age_at_diagnosis": 66}, {"age_at_diagnosis": 68},
             {"age_at_diagnosis": np.nan}])
        d = df["disease_duration_years"]
        self.assertEqual(d.iloc[0], 5)
        self.assertEqual(d.iloc[1], 0)
        self.assertTrue(np.isnan(d.iloc[2]) and np.isnan(d.iloc[3]))
        self.assertEqual(df["duration_status"].tolist(), ["ok", "rounded_to_zero", "negative_set_missing", "missing"])
        controls = self.build([{}], [{"age_at_diagnosis": 50}], dx="Control")
        self.assertTrue(np.isnan(controls["disease_duration_years"].iloc[0]))

    def test_medication_history_first_then_updrs3_3a(self):
        df = self.build(
            [{}, {}, {"updrs3_medication_for_pd": "Yes"}, {"updrs3_medication_for_pd": "Yes"}, {}],
            [{"use_of_pd_medication": "No", "on_levodopa": "Yes"},          # any yes wins
             {"use_of_pd_medication": "No"},
             {"use_of_pd_medication": np.nan},                              # falls back to 3a
             {"use_of_pd_medication": "No"},                                # history beats 3a
             {"use_of_pd_medication": np.nan}])
        self.assertEqual([None if pd.isna(v) else bool(v) for v in df["pd_medication"]],
                         [True, False, True, False, None])
        self.assertEqual([None if pd.isna(v) else v for v in df["pd_medication_source"]],
                         ["medical_history", "medical_history", "updrs3_3a", "medical_history", None])

    def test_bands_and_subset_refuse_until_fixed(self):
        self.assertIsNone(analysis_set.FIXED_DEFINITIONS)
        with self.assertRaises(RuntimeError):
            analysis_set.require_fixed_definitions()


class TestCoverageReport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(0)
        rows, fields = [], []
        for i in range(40):
            pid = f"PD{i:02d}"
            sex = "Female" if i % 2 else "Male"
            dx = "Control" if i < 10 else "Parkinson's Disease"
            rows.append((pid, sex, "x", {"study": "PDBP", "diagnosis_at_baseline": dx}))
            dur = int(rng.integers(0, 12))
            fields.append((pid, {"age_at_diagnosis": 60 - dur if i % 5 else np.nan,
                                 "use_of_pd_medication": ("Yes" if dur > 2 else "No") if i % 3 else np.nan}))
        meta = mofa_meta(rows)
        cls.participants = build_participants(meta, ml_fields(fields), "PDBP")
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "pdbp.csv")
            cls.participants.to_csv(path)
            cls.tables = coverage_report(load([path]))
            with open(write_report(cls.tables, os.path.join(tmp, "out"))) as f:
                cls.report = f.read()

    def test_coverage_matches_participant_level_flags(self):
        cov = self.tables["coverage_pd"].set_index("sex")
        for sex in ("Female", "Male"):
            g = self.participants[(self.participants["sex"] == sex) & (self.participants["diagnosis"] == "PD")]
            self.assertEqual(cov.loc[sex, "n_pd"], len(g))
            self.assertEqual(cov.loc[sex, "duration_n"], g["disease_duration_years"].notna().sum())
            self.assertEqual(cov.loc[sex, "medication_n"], g["pd_medication"].notna().sum())

    def test_ppmi_like_candidates_count_short_unmedicated_pd(self):
        cand = self.tables["ppmi_like_candidates"]
        row = cand[(cand["max_duration_years"] == 2) & cand["unmedicated_required"]].iloc[0]
        p = self.participants
        expected = ((p["sex"] == "Female") & (p["diagnosis"] == "PD") & (p["disease_duration_years"] <= 2)
                    & p["pd_medication"].eq(False).fillna(False)).sum()
        self.assertEqual(row["Female_pd"], expected)
        self.assertEqual(row["Female_controls"], ((p["sex"] == "Female") & (p["diagnosis"] == "Control")).sum())

    def test_reference_check_reports_differences(self):
        ref = self.tables["reference_check"]
        total = ref[(ref["cohort"] == "PDBP") & (ref["count"] == "RNA-seq participants")].iloc[0]
        self.assertEqual(total["observed"], 40)
        self.assertEqual(total["difference"], 40 - 1248)

    def test_no_participant_ids_leave_in_the_report(self):
        for pid in self.participants.index:
            self.assertNotIn(pid, self.report)
        for name, t in self.tables.items():
            self.assertFalse(t.astype(str).isin(list(self.participants.index)).any().any(), name)


if __name__ == "__main__":
    unittest.main()
