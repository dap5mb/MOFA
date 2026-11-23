"""Issue 04: label-free preprocessing, the PPMI-only gene list, and the internal PPMI AUC."""

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ML_SCRIPTS = os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts"))
sys.path.insert(0, ML_SCRIPTS)
import shared
from normalize_cohort import normalize_cohort, read_participants
from prepare_features import prepare_features
from select_genes import read_gene_list, select_genes
from synthetic import make_synthetic_cohorts
from train import train_elastic_net
from train_ppmi import internal_auc_table

SMALL_GRID = {"C": [0.01, 0.1, 1.0], "l1_ratio": [0.5]}


def _limma_available():
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(["Rscript", "-e", "suppressPackageStartupMessages({library(limma); library(dplyr); library(ggplot2)})"],
                       capture_output=True)
    return r.returncode == 0


HAVE_LIMMA = _limma_available()
HAVE_PYDESEQ2 = importlib.util.find_spec("pydeseq2") is not None


def synthetic_counts(cohort, n=60, n_genes=40, n_unexpressed=5, seed=0):
    """Raw counts (genes x participants) with a plate shift and a PD block, and a participants file."""
    rng = np.random.default_rng(seed)
    ids = [f"{cohort}_{i:03d}" for i in range(n)]
    part = pd.DataFrame({
        "participant_id": ids, "cohort": cohort,
        "sex": rng.choice(["Female", "Male"], n),
        "diagnosis": rng.choice(["PD", "Control"], n),
        "plate": rng.choice(["P1", "P2", "P3"], n),
    })
    part["in_analysis_set"] = True
    part.loc[0, ["diagnosis", "in_analysis_set"]] = ["Other", False]
    mu = rng.uniform(50, 500, n_genes)[:, None] * rng.uniform(0.7, 1.3, n)[None, :]
    mu = mu * np.where(part["plate"] == "P2", 1.5, 1.0)[None, :]
    mu[:5] *= np.where(part["diagnosis"] == "PD", 1.8, 1.0)[None, :]
    counts = rng.poisson(mu)
    counts[-n_unexpressed:] = rng.poisson(0.2, (n_unexpressed, n))  # fail the low-count filter
    genes = [f"ENSG{j:011d}.{j % 3 + 1}" for j in range(n_genes)]
    return pd.DataFrame(counts, index=genes, columns=ids), part


class FakeNormalizer:
    """Stands in for MOFA's preprocess_rnaseq where pydeseq2 is missing: log1p CPM, centered."""

    def __init__(self):
        self.calls = []

    def __call__(self, counts_file, metadata_file, output_file, visit_filter, cohort_name, n_top_genes):
        meta = pd.read_csv(metadata_file)
        self.calls.append({"columns": list(meta.columns), "n_top_genes": n_top_genes})
        counts = pd.read_csv(counts_file, index_col=0)[list(meta["participant_id"])]
        x = np.log1p(counts / counts.sum() * 1e6)
        x.sub(x.mean(axis=1), axis=0).to_csv(output_file)


class RecordingPlateCorrection:
    def __init__(self):
        self.calls = []
        import plate_correction
        self.real = plate_correction.correct_plate

    def __call__(self, expr, metadata, include_diagnosis=True, **kw):
        self.calls.append({"columns": list(metadata.columns), "include_diagnosis": include_diagnosis})
        return self.real(expr, metadata, include_diagnosis=include_diagnosis, **kw)


@unittest.skipUnless(HAVE_LIMMA, "Rscript with limma, dplyr and ggplot2 not available")
class TestNormalizationIsLabelFree(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.counts, self.part = synthetic_counts("PPMI")
        self.counts_path = os.path.join(self.tmp, "counts.csv")
        self.counts.to_csv(self.counts_path)
        self.normalizer, self.plate = FakeNormalizer(), RecordingPlateCorrection()
        if not HAVE_PYDESEQ2:
            shared.preprocess_rnaseq = self.normalizer
        shared.correct_plate = self.plate

    def tearDown(self):
        for name in ("preprocess_rnaseq", "correct_plate"):
            shared.__dict__.pop(name, None)
        shutil.rmtree(self.tmp)

    def _run(self, part, cohort="PPMI"):
        path = os.path.join(self.tmp, "participants.csv")
        part.to_csv(path, index=False)
        return normalize_cohort(self.counts_path, read_participants(path), cohort)

    def test_permuting_diagnosis_leaves_the_matrix_unchanged(self):
        a = self._run(self.part)
        permuted = self.part.copy()
        rng = np.random.default_rng(1)
        permuted.loc[1:, "diagnosis"] = rng.permutation(permuted.loc[1:, "diagnosis"].to_numpy())
        b = self._run(permuted)
        pd.testing.assert_frame_equal(a, b, check_exact=True)

    def test_diagnosis_reaches_no_design(self):
        self._run(self.part)
        for call in self.plate.calls:
            self.assertFalse(call["include_diagnosis"])
            self.assertNotIn("diagnosis", call["columns"])
        for call in self.normalizer.calls:
            self.assertNotIn("diagnosis", call["columns"])

    def test_only_analysis_set_participants_are_kept(self):
        out = self._run(self.part)
        self.assertEqual(list(out.columns), list(self.part.loc[self.part["in_analysis_set"], "participant_id"]))

    def test_refuses_diagnosis_and_other_cohorts(self):
        with self.assertRaisesRegex(ValueError, "diagnosis"):
            normalize_cohort(self.counts_path, self.part.set_index("participant_id"), "PPMI")
        with self.assertRaisesRegex(ValueError, "PDBP normalization got"):
            self._run(self.part, cohort="PDBP")


@unittest.skipUnless(HAVE_LIMMA and HAVE_PYDESEQ2, "needs pydeseq2 and Rscript with limma (Terra)")
class TestNormalizationWithMofaCode(unittest.TestCase):
    """The real preprocess_rnaseq: PPMI keeps expressed genes only, PDBP keeps all genes."""

    def test_gene_sets_and_label_freedom(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = {}
            for cohort, seed in (("PPMI", 0), ("PDBP", 1)):
                counts, part = synthetic_counts(cohort, seed=seed)
                cpath, ppath = os.path.join(tmp, f"{cohort}_counts.csv"), os.path.join(tmp, f"{cohort}_p.csv")
                counts.to_csv(cpath)
                part.to_csv(ppath, index=False)
                out[cohort] = normalize_cohort(cpath, read_participants(ppath), cohort)
                part["diagnosis"] = np.random.default_rng(2).permutation(part["diagnosis"].to_numpy())
                part.to_csv(ppath, index=False)
                pd.testing.assert_frame_equal(out[cohort], normalize_cohort(cpath, read_participants(ppath), cohort))
            self.assertEqual(len(out["PPMI"]), 35)
            self.assertEqual(len(out["PDBP"]), 40)


class TestGeneSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ppmi_expr, cls.ppmi_meta, cls.pdbp_expr, cls.pdbp_meta = make_synthetic_cohorts("signal", seed=0)

    def test_list_comes_from_ppmi_only(self):
        genes, report = select_genes(self.ppmi_expr, self.ppmi_meta, 50, pdbp_ids=self.pdbp_meta.index)
        self.assertEqual(len(genes), 50)
        self.assertEqual(report["cohorts"], ["PPMI"])
        self.assertEqual(report["n_pdbp_participants_in_selection"], 0)
        self.assertFalse(report["labels_used"])

    def test_pdbp_participants_are_refused(self):
        mixed = pd.concat([self.ppmi_expr, self.pdbp_expr.iloc[:3]])
        meta = pd.concat([self.ppmi_meta, self.pdbp_meta])
        with self.assertRaisesRegex(ValueError, "PPMI participants only"):
            select_genes(mixed, meta, 50)
        # A PDBP sample mislabelled as PPMI is still caught by the ID check
        relabelled = meta.copy()
        relabelled.loc[self.pdbp_meta.index, "cohort"] = "PPMI"
        with self.assertRaisesRegex(ValueError, "PDBP participants entered"):
            select_genes(mixed, relabelled, 50, pdbp_ids=self.pdbp_meta.index)

    def test_selection_ignores_labels(self):
        meta = self.ppmi_meta.copy()
        meta["diagnosis"] = np.random.default_rng(3).permutation(meta["diagnosis"].to_numpy())
        self.assertEqual(select_genes(self.ppmi_expr, self.ppmi_meta, 50)[0],
                         select_genes(self.ppmi_expr, meta, 50)[0])

    def test_mofa_overlap_matches_across_id_styles(self):
        genes, _ = select_genes(self.ppmi_expr, self.ppmi_meta, 50)
        mofa = [f"{g}.7_RNA_seq" for g in genes[:40]] + ["OTHER1", "OTHER2"]
        _, report = select_genes(self.ppmi_expr, self.ppmi_meta, 50, mofa_genes=mofa)
        o = report["mofa_overlap"]
        self.assertEqual((o["n_shared"], o["n_reference"]), (40, 42))
        self.assertAlmostEqual(o["fraction_of_selected"], 0.8)


class TestPrepareFeatures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ppmi_expr, cls.ppmi_meta, cls.pdbp_expr, cls.pdbp_meta = make_synthetic_cohorts("signal", seed=0)
        cls.genes = list(cls.ppmi_expr.columns[[5, 1, 9]])

    def test_list_applied_unchanged_and_standardized_within_sex(self):
        out = prepare_features(self.pdbp_expr, self.pdbp_meta, self.genes)
        self.assertEqual(list(out.columns), self.genes)
        by_sex = out.groupby(self.pdbp_meta.loc[out.index, "sex"])
        np.testing.assert_allclose(by_sex.mean().to_numpy(), 0, atol=1e-12)
        np.testing.assert_allclose(by_sex.std().to_numpy(), 1, atol=1e-12)

    def test_missing_gene_stops_the_run(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            prepare_features(self.pdbp_expr.drop(columns=self.genes[0]), self.pdbp_meta, self.genes)

    def test_one_cohort_at_a_time(self):
        mixed = pd.concat([self.ppmi_expr, self.pdbp_expr])
        with self.assertRaisesRegex(ValueError, "one cohort"):
            prepare_features(mixed, pd.concat([self.ppmi_meta, self.pdbp_meta]), self.genes)


class TestInternalAucBySex(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        expr, cls.meta, _, _ = make_synthetic_cohorts("signal", seed=0)
        cls.genes = list(expr.columns[:100])
        cls.expr = expr
        cls.model, cls.oof = train_elastic_net(expr, cls.meta, cls.genes, grid=SMALL_GRID, return_oof=True)

    def test_every_participant_scored_once_out_of_fold(self):
        self.assertEqual(list(self.oof.index), list(self.expr.index))
        self.assertFalse(self.oof.isna().any())

    def test_auc_reported_overall_and_per_sex(self):
        cv = self.model.internal_cv_auc
        self.assertEqual(set(cv["by_sex"]), {"Female", "Male"})
        self.assertEqual(sum(r["n"] for r in cv["by_sex"].values()), self.model.training["n"])
        for r in cv["by_sex"].values():
            self.assertGreater(r["auc_pooled_oof"], 0.6)
        table = internal_auc_table(self.model)
        self.assertEqual(list(table["sex"]), ["all", "Female", "Male"])
        self.assertTrue((table["cohort"] == "PPMI").all())

    def test_same_model_as_before(self):
        # Out-of-fold bookkeeping changes neither the folds nor the refit model
        plain = train_elastic_net(self.expr, self.meta, self.genes, grid=SMALL_GRID)
        self.assertEqual(plain.to_dict(), self.model.to_dict())


class TestScriptsRunInOrder(unittest.TestCase):
    """select_genes -> prepare_features (both cohorts) -> train_ppmi, as the Terra notebook calls them."""

    def test_chain_writes_ppmi_results_and_no_pdbp_scores(self):
        ppmi_expr, ppmi_meta, pdbp_expr, pdbp_meta = make_synthetic_cohorts("signal", n_ppmi=200, n_genes=80, seed=0)
        with tempfile.TemporaryDirectory() as tmp:
            P = lambda name: os.path.join(tmp, name)  # noqa: E731
            for cohort, expr, meta in (("ppmi", ppmi_expr, ppmi_meta), ("pdbp", pdbp_expr, pdbp_meta)):
                expr.T.to_csv(P(f"{cohort}_expr.csv"))
                meta.assign(in_analysis_set=True).rename_axis("participant_id").to_csv(P(f"{cohort}_participants.csv"))

            def run(script, *args):
                subprocess.run([sys.executable, os.path.join(ML_SCRIPTS, script), *args], check=True,
                               capture_output=True, text=True)

            run("select_genes.py", "--expr", P("ppmi_expr.csv"), "--participants", P("ppmi_participants.csv"),
                "--pdbp_participants", P("pdbp_participants.csv"), "--n_top", "30",
                "--genes_out", P("genes.txt"), "--report", P("report.json"))
            for c in ("ppmi", "pdbp"):
                run("prepare_features.py", "--expr", P(f"{c}_expr.csv"), "--participants", P(f"{c}_participants.csv"),
                    "--genes", P("genes.txt"), "--output", P(f"{c}_features.csv"))
            run("train_ppmi.py", "--features", P("ppmi_features.csv"), "--participants", P("ppmi_participants.csv"),
                "--model", P("model.json"), "--auc", P("auc.tsv"), "--oof", P("oof.csv"))

            genes = read_gene_list(P("genes.txt"))
            self.assertEqual(len(genes), 30)
            with open(P("report.json")) as fh:
                self.assertEqual(json.load(fh)["n_pdbp_participants_in_selection"], 0)
            self.assertEqual(list(pd.read_csv(P("pdbp_features.csv"), index_col=0).columns), genes)
            auc = pd.read_csv(P("auc.tsv"), sep="\t")
            self.assertEqual(list(auc["sex"]), ["all", "Female", "Male"])
            self.assertTrue(pd.read_csv(P("oof.csv"), index_col=0).index.str.startswith("PPMI").all())
            self.assertFalse([f for f in os.listdir(tmp) if "score" in f])


if __name__ == "__main__":
    unittest.main()
