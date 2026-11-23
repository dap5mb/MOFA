import unittest
import numpy as np
import pandas as pd
import h5py
import tempfile
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))

from compare_models import match


def write_model(path, z, w, r2, groups=("all",)):
    """Minimal mofapy2-layout HDF5 with one view: Z per group (factors x samples), W (factors x features), R2."""
    n = z.shape[0]
    halves = np.array_split(np.arange(n), len(groups))
    with h5py.File(path, "w") as f:
        exp = f.create_group("expectations")
        zg, sg = exp.create_group("Z"), f.create_group("samples")
        r2g = f.create_group("variance_explained").create_group("r2_per_factor")
        for g, idx in zip(groups, halves):
            zg.create_dataset(g, data=z[idx].T)
            sg.create_dataset(g, data=np.array([f"S{i}".encode() for i in idx]))
            r2g.create_dataset(g, data=np.array([r2]))
        exp.create_group("W").create_dataset("RNA_seq", data=w.T)
        f.create_group("features").create_dataset("RNA_seq", data=np.array([f"g{j}".encode() for j in range(w.shape[0])]))


class TestCompareModels(unittest.TestCase):
    def test_matches_shared_factors_and_leaves_extra_factor_unmatched(self):
        rng = np.random.default_rng(0)
        n, d = 200, 60
        shared_z, shared_w = rng.normal(size=(n, 2)), rng.normal(size=(d, 2))
        sex_z, sex_w = rng.normal(size=(n, 1)), rng.normal(size=(d, 1))
        with tempfile.TemporaryDirectory() as tmp:
            a, b = os.path.join(tmp, "a.hdf5"), os.path.join(tmp, "b.hdf5")
            write_model(a, shared_z, shared_w, [30.0, 20.0], groups=("Female", "Male"))
            # Model B: factors reordered and sign-flipped, plus an extra "sex" factor
            write_model(b, np.column_stack([sex_z, -shared_z[:, 1], shared_z[:, 0]]),
                        np.column_stack([sex_w, -shared_w[:, 1], shared_w[:, 0]]), [15.0, 20.0, 30.0])
            a_to_b = match(a, b, 0.1).set_index("factor")
            b_to_a = match(b, a, 0.1).set_index("factor")
        self.assertEqual(a_to_b.loc["Factor1", "matched_factor"], "Factor3")
        self.assertEqual(a_to_b.loc["Factor2", "matched_factor"], "Factor2")
        self.assertGreater(a_to_b["z_abs_r"].min(), 0.99)
        self.assertGreater(a_to_b["w_abs_r_RNA_seq"].min(), 0.99)
        self.assertLess(b_to_a.loc["Factor1", "z_abs_r"], 0.3)  # the extra factor has no counterpart
        self.assertEqual(len(a_to_b), 2)
        self.assertAlmostEqual(b_to_a.loc["Factor1", "max_R2"], 15.0)


class TestSingleGroupProjection(unittest.TestCase):
    def test_single_group_model_ignores_sex_for_grouping(self):
        from project_pdbp import project_pdbp
        rng = np.random.default_rng(1)
        n, d = 80, 30
        w = rng.normal(size=(d, 2))
        ids = [f"S{i}" for i in range(n)]
        with tempfile.TemporaryDirectory() as tmp:
            model = os.path.join(tmp, "m.hdf5")
            with h5py.File(model, "w") as f:
                f.create_group("groups").create_dataset("groups", data=np.array([b"all"]))
                f.create_group("features").create_dataset("RNA_seq", data=np.array([f"g{j}".encode() for j in range(d)]))
                f.create_group("expectations").create_group("W").create_dataset("RNA_seq", data=w.T)
            os.makedirs(os.path.join(tmp, "in"))
            pd.DataFrame((rng.normal(size=(n, 2)) @ w.T).T, index=[f"g{j}" for j in range(d)], columns=ids).to_csv(
                os.path.join(tmp, "in", "pdbp_rna_seq_aligned.csv"))
            meta = os.path.join(tmp, "meta.csv")
            pd.DataFrame({"participant_id": ids, "sex": ["Female", "Male"] * (n // 2)}).to_csv(meta, index=False)
            project_pdbp(model, os.path.join(tmp, "in"), os.path.join(tmp, "out"), "PDBP", meta)
            z = pd.read_csv(os.path.join(tmp, "out", "pdbp_latent_factors_Z.csv"), index_col=0)
        self.assertEqual(len(z), n)
        self.assertEqual(list(z["group"].unique()), ["all"])


if __name__ == "__main__":
    unittest.main()
