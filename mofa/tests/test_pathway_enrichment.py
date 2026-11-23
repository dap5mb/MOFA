import os
import subprocess
import sys

import numpy as np
import pandas as pd

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, SCRIPTS)
from pathway_enrichment import load_view_weights, ranked_symbols  # noqa: E402


def test_ranked_symbols_drops_unmapped_and_keeps_largest_duplicate():
    w = pd.Series([0.5, -0.9, 0.1, 0.3], index=["ENSG1.2", "ENSG2.1", "ENSG3.1", "ENSG4.7"])
    symbols = pd.Series({"ENSG1": "A", "ENSG2": "A", "ENSG3": None, "ENSG4": "B"})
    r = ranked_symbols(w, symbols)
    assert r.to_dict() == {"B": 0.3, "A": -0.9}


def test_load_view_weights_dedupes_groups(tmp_path):
    df = pd.DataFrame({"Factor1": [1.0, 2.0, 1.0, 5.0], "view": ["RNA_seq", "RNA_seq", "RNA_seq", "CSF"],
                       "group": ["F", "F", "M", "F"], "feature": ["g1", "g2", "g1", "p1"]})
    path = tmp_path / "w.csv"
    df.to_csv(path, index=False)
    w = load_view_weights(str(path), "RNA_seq")
    assert list(w.index) == ["g1", "g2"]


def test_end_to_end_finds_planted_sets(tmp_path):
    rng = np.random.default_rng(1)
    n = 600
    ids = [f"ENSG{i:05d}.1" for i in range(n)]
    syms = [f"G{i}" for i in range(n)]
    f1 = rng.normal(size=n)
    up, down = syms[:25], syms[25:50]
    f1[:25] += 3  # planted positive set
    f1[25:50] -= 3  # planted negative set
    w = pd.DataFrame({"Factor1": f1, "Factor2": rng.normal(size=n), "view": "RNA_seq", "feature": ids})
    w.to_csv(tmp_path / "w.csv", index=False)
    pd.DataFrame({"ensembl_id": [i.split(".")[0] for i in ids], "symbol": syms}).to_csv(tmp_path / "map.csv", index=False)
    rand = list(rng.choice(syms[50:], 30, replace=False))
    with open(tmp_path / "toy.v2024.1.gmt", "w") as f:
        for name, genes in (("UP", up), ("DOWN", down), ("RANDOM", rand)):
            f.write("\t".join([name, "na"] + genes) + "\n")
    out = tmp_path / "out"
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "pathway_enrichment.py"), "--weights", str(tmp_path / "w.csv"),
                    "--factors", "Factor1", "--gene_sets", str(tmp_path / "toy.v2024.1.gmt"), "--gene_map",
                    str(tmp_path / "map.csv"), "--permutations", "200", "--processes", "1", "--output", str(out)],
                   check=True)
    res = pd.read_csv(out / "gsea_prerank_results.csv").set_index("term")
    assert res.loc["UP", "nes"] > 0 and res.loc["UP", "fdr"] < 0.05
    assert res.loc["DOWN", "nes"] < 0 and res.loc["DOWN", "fdr"] < 0.05
    assert res.loc["RANDOM", "fdr"] > 0.05
    assert set(res["library"]) == {"toy"}
