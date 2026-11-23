import unittest
import pandas as pd
import numpy as np
import tempfile
import os
import sys
import h5py

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../scripts')))

class TestPDBPProjection(unittest.TestCase):
    def test_project_pdbp_mock(self):
        # We no longer mock entry_point, we test the new implementation with load_model_weights
        # This will be tricky because we need an HDF5 file. 
        # I'll update the test to create a dummy HDF5 file first.
        
        with tempfile.TemporaryDirectory() as tmpdir:
            model_path = os.path.join(tmpdir, "dummy_model.hdf5")
            with h5py.File(model_path, 'w') as f:
                exp = f.create_group('expectations')
                w = exp.create_group('W')
                view_w = w.create_group('RNA_seq')
                view_w.create_dataset('group1', data=np.array([[0.1, 0.2], [0.3, 0.4]]))
                feat = f.create_group('features')
                feat.create_dataset('RNA_seq', data=np.array([b'Gene1', b'Gene2']))
                
            input_dir = os.path.join(tmpdir, "input")
            output_dir = os.path.join(tmpdir, "output")
            os.makedirs(input_dir, exist_ok=True)
            
            # Create dummy aligned CSV
            df = pd.DataFrame({"sample1": [1.0, 2.0]}, index=["Gene1", "Gene2"])
            df.to_csv(os.path.join(input_dir, "pdbp_rna_seq_aligned.csv"))
            
            # Run project_pdbp
            from project_pdbp import project_pdbp
            project_pdbp(model_path, input_dir, output_dir, "PDBP")
            
            out_file = os.path.join(output_dir, "pdbp_latent_factors_Z.csv")
            self.assertTrue(os.path.exists(out_file))
            
            res = pd.read_csv(out_file, index_col=0)
            self.assertEqual(res.shape[0], 1)
            self.assertEqual([c for c in res.columns if c.startswith("Factor")], ["Factor1", "Factor2"])

    def test_projection_centers_within_group_and_uses_observed_features(self):
        """Mimics a mofapy2 multi-group HDF5: sex mean shifts must not leak into Z, missing views must not be zero-filled."""
        from project_pdbp import project_pdbp
        rng = np.random.default_rng(1)
        K, n = 2, 200
        dims = {"RNA_seq": 60, "CSF_Proteomics": 20}
        W = {v: rng.normal(size=(d, K)) for v, d in dims.items()}
        groups = np.array(["Female", "Male"] * (n // 2))
        Z_true = rng.normal(size=(n, K))
        ids = [f"S{i}" for i in range(n)]

        with tempfile.TemporaryDirectory() as tmpdir:
            model_path = os.path.join(tmpdir, "model.hdf5")
            with h5py.File(model_path, 'w') as f:
                f.create_group('views').create_dataset('views', data=np.array([b"CSF_Proteomics", b"RNA_seq"]))
                f.create_group('groups').create_dataset('groups', data=np.array([b"Female", b"Male"]))
                feat = f.create_group('features')
                w_grp = f.create_group('expectations').create_group('W')
                for v, d in dims.items():
                    feat.create_dataset(v, data=np.array([f"{v}_{j}".encode() for j in range(d)]))
                    w_grp.create_dataset(v, data=W[v].T)  # mofapy2 layout: factors x features

            input_dir = os.path.join(tmpdir, "input")
            os.makedirs(input_dir)
            for v, d in dims.items():
                X = Z_true @ W[v].T + 0.05 * rng.normal(size=(n, d))
                X[groups == "Male"] += 3.0  # sex mean shift that group centering must remove
                df = pd.DataFrame(X.T, index=[f"{v}_{j}" for j in range(d)], columns=ids)
                if v == "CSF_Proteomics":
                    df = df.iloc[:, : n // 2]  # half the samples have no CSF view
                df.to_csv(os.path.join(input_dir, f"pdbp_{v.lower()}_aligned.csv"))
            meta_path = os.path.join(tmpdir, "meta.csv")
            pd.DataFrame({"participant_id": ids, "sex": [g.lower() for g in groups]}).to_csv(meta_path, index=False)

            output_dir = os.path.join(tmpdir, "output")
            project_pdbp(model_path, input_dir, output_dir, "PDBP", meta_path)

            z = pd.read_csv(os.path.join(output_dir, "pdbp_latent_factors_Z.csv"), index_col=0).loc[ids]
            self.assertEqual(sorted(z["group"].unique()), ["Female", "Male"])
            for k in range(K):
                # Factors recovered up to the view-scaling constant, including samples missing CSF
                self.assertGreater(abs(np.corrcoef(z[f"Factor{k+1}"], Z_true[:, k])[0, 1]), 0.95)
                no_csf = ids[n // 2:]
                self.assertGreater(abs(np.corrcoef(z.loc[no_csf, f"Factor{k+1}"], Z_true[n // 2:, k])[0, 1]), 0.95)
                # Group means are ~0 in both sexes (sex shift removed)
                self.assertLess(z.groupby("group")[f"Factor{k+1}"].mean().abs().max(), 0.1)

            r2 = pd.read_csv(os.path.join(output_dir, "pdbp_variance_explained_per_factor.csv"))
            self.assertEqual(len(r2), len(dims) * 2)  # one row per view x group
            total = pd.read_csv(os.path.join(output_dir, "pdbp_variance_explained_total.csv"))
            self.assertTrue(((total["R2_total"] > 90) & (total["R2_total"] <= 100)).all())


class TestVarianceExplained(unittest.TestCase):
    def test_matches_mofa_definition(self):
        from project_pdbp import variance_explained
        rng = np.random.default_rng(2)
        Z, W = rng.normal(size=(50, 3)), rng.normal(size=(20, 3))
        Y = Z @ W.T
        Y[0, 0] = np.nan
        per_factor, total = variance_explained(Y, Z, W)
        self.assertAlmostEqual(total, 100.0, places=6)
        obs = ~np.isnan(Y)
        expected = (1 - np.sum((Y[obs] - np.outer(Z[:, 0], W[:, 0])[obs]) ** 2) / np.sum(Y[obs] ** 2)) * 100
        self.assertAlmostEqual(per_factor[0], expected, places=6)


class TestSavedTauAlphaZ(unittest.TestCase):
    def test_saved_alphaz_shrinks_factor_inactive_in_one_group(self):
        """A factor with huge AlphaZ in Male (inactive there) must project to ~0 for males only."""
        from project_pdbp import project_pdbp
        rng = np.random.default_rng(3)
        K, n, d = 2, 100, 40
        W = rng.normal(size=(d, K))
        groups = np.array(["Female", "Male"] * (n // 2))
        ids = [f"S{i}" for i in range(n)]
        X = rng.normal(size=(n, K)) @ W.T + 0.1 * rng.normal(size=(n, d))

        with tempfile.TemporaryDirectory() as tmpdir:
            model_path = os.path.join(tmpdir, "model.hdf5")
            with h5py.File(model_path, 'w') as f:
                f.create_group('groups').create_dataset('groups', data=np.array([b"Female", b"Male"]))
                f.create_group('features').create_dataset('RNA_seq', data=np.array([f"g{j}".encode() for j in range(d)]))
                exp = f.create_group('expectations')
                exp.create_group('W').create_dataset('RNA_seq', data=W.T)
                tau = exp.create_group('Tau').create_group('RNA_seq')
                alpha = exp.create_group('AlphaZ')
                for g in ("Female", "Male"):
                    tau.create_dataset(g, data=np.full((n // 2, d), 100.0))
                alpha.create_dataset("Female", data=np.array([1.0, 1.0]))
                alpha.create_dataset("Male", data=np.array([1.0, 1e8]))

            input_dir = os.path.join(tmpdir, "input")
            os.makedirs(input_dir)
            pd.DataFrame(X.T, index=[f"g{j}" for j in range(d)], columns=ids).to_csv(
                os.path.join(input_dir, "pdbp_rna_seq_aligned.csv"))
            meta_path = os.path.join(tmpdir, "meta.csv")
            pd.DataFrame({"participant_id": ids, "sex": groups}).to_csv(meta_path, index=False)

            output_dir = os.path.join(tmpdir, "output")
            project_pdbp(model_path, input_dir, output_dir, "PDBP", meta_path)
            z = pd.read_csv(os.path.join(output_dir, "pdbp_latent_factors_Z.csv"), index_col=0)

            self.assertLess(z.loc[z.group == "Male", "Factor2"].abs().max(), 1e-3)
            self.assertGreater(z.loc[z.group == "Female", "Factor2"].std(), 0.1)
            self.assertGreater(z.loc[z.group == "Male", "Factor1"].std(), 0.1)


if __name__ == "__main__":
    unittest.main()
