import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from factor_technical_check import dx_with_and_without_tech, library_stats, top_features  # noqa: E402


def test_library_stats_sums_across_chunks(tmp_path):
    counts = pd.DataFrame({"P1": [10, 0, 5], "P2": [100, 1, 0]}, index=["g1", "g2", "g3"])
    path = tmp_path / "counts.csv"
    counts.to_csv(path)
    out = library_stats(str(path), chunksize=2).set_index("participant_id")
    assert np.isclose(out.loc["P1", "log10_lib_size"], np.log10(15))
    assert np.isclose(out.loc["P2", "log10_lib_size"], np.log10(101))
    assert out.loc["P1", "genes_detected"] == 2 and out.loc["P2", "genes_detected"] == 2


def test_adjustment_removes_a_pure_library_size_effect():
    # PD depends on library size only; the factor is library size plus noise
    rng = np.random.default_rng(0)
    n = 2000
    lib = rng.normal(size=n)
    pd_status = (rng.random(n) < 1 / (1 + np.exp(-1.5 * lib))).astype(float)
    df = pd.DataFrame({"Factor1": lib + 0.3 * rng.normal(size=n), "pd": pd_status,
                       "age": rng.normal(60, 10, n), "male": rng.integers(0, 2, n).astype(float),
                       "log10_lib_size": 7 + 0.1 * lib, "genes_detected": 20000 + rng.normal(0, 100, n)})
    df["Factor1"] = (df["Factor1"] - df["Factor1"].mean()) / df["Factor1"].std()
    row = dx_with_and_without_tech(df, ["Factor1"]).iloc[0]
    assert row["p"] < 1e-10
    assert abs(row["estimate_adj_lib"]) < abs(row["estimate"]) / 3


def test_top_features_ignores_group_duplicates(tmp_path):
    w = pd.DataFrame({"Factor4": [0.9, -0.8, 0.1, 0.9], "view": ["RNA_seq"] * 4,
                      "group": ["Female", "Female", "Female", "Male"], "feature": ["A", "B", "C", "A"]})
    path = tmp_path / "w.csv"
    w.to_csv(path, index=False)
    top = top_features(str(path), "Factor4", n_top=30)
    pos = top[(top.sign == "positive")]
    assert list(pos.feature) == ["A", "C"]
    assert list(top[top.sign == "negative"].feature) == ["B"]
