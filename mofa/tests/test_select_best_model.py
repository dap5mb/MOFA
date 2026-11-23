import unittest
import numpy as np
import pandas as pd
import h5py
import tempfile
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))

from select_best_model import get_convergence, match_factors, load_z


def write_model(path, z, elbo, maxiter, groups=("Female", "Male")):
    """Minimal mofapy2-layout HDF5: Z per group as (factors, samples), ELBO trace, training_opts."""
    n = z.shape[0]
    halves = np.array_split(np.arange(n), len(groups))
    with h5py.File(path, "w") as f:
        zg = f.create_group("expectations").create_group("Z")
        sg = f.create_group("samples")
        for g, idx in zip(groups, halves):
            zg.create_dataset(g, data=z[idx].T)
            sg.create_dataset(g, data=np.array([f"S{i}".encode() for i in idx]))
        f.create_group("training_stats").create_dataset("elbo", data=np.asarray(elbo, dtype=float))
        opts = f.create_dataset("training_opts", data=np.array([maxiter, 5, 1, 0, 0, 0], dtype=float))
        opts.attrs["names"] = np.array([b"maxiter", b"freqELBO", b"start_elbo", b"gpu_mode", b"stochastic", b"seed"])


class TestSelectBestModel(unittest.TestCase):
    def test_convergence_detects_cap_and_early_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            capped, early = os.path.join(tmp, "capped.hdf5"), os.path.join(tmp, "early.hdf5")
            z = np.random.default_rng(0).normal(size=(20, 2))
            write_model(capped, z, [-100.0] + [np.nan] * 4 + [-50.0, -49.9], maxiter=6)  # 7 entries = cap
            write_model(early, z, [-100.0, np.nan, -50.0, -49.99], maxiter=10)           # truncated
            self.assertFalse(get_convergence(capped)["converged"])
            self.assertEqual(get_convergence(capped)["iterations"], 6)
            self.assertTrue(get_convergence(early)["converged"])

    def test_match_factors_handles_permutation_sign_and_constant_factor(self):
        rng = np.random.default_rng(1)
        z = rng.normal(size=(100, 3))
        z[:, 2] = 0.0  # inactive factor
        other = np.column_stack([-z[:, 1], z[:, 0] + 0.01 * rng.normal(size=100), rng.normal(size=100)])
        with tempfile.TemporaryDirectory() as tmp:
            write_model(os.path.join(tmp, "a.hdf5"), z, [-1.0, -0.5], 10)
            write_model(os.path.join(tmp, "b.hdf5"), other, [-1.0, -0.5], 10)
            r, idx = match_factors(load_z(os.path.join(tmp, "a.hdf5")), load_z(os.path.join(tmp, "b.hdf5")))
        self.assertGreater(r["Factor1"], 0.99)
        self.assertEqual(idx["Factor1"], "Factor2")
        self.assertGreater(r["Factor2"], 0.99)  # matched despite the sign flip
        self.assertEqual(idx["Factor2"], "Factor1")
        self.assertTrue(np.isnan(r["Factor3"]))
        self.assertTrue(pd.isna(idx["Factor3"]))


if __name__ == "__main__":
    unittest.main()
