import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import string_networks as sn  # noqa: E402


def test_select_top_skips_unmapped_and_splits_directions():
    ranked = pd.Series({"A": 0.9, "LOC1": 0.8, "B": 0.5, "C": -0.2, "D": -0.7, "IGHG1": -0.95})
    mapped = pd.DataFrame({"symbol": ["A", "B", "C", "D"], "string_id": ["s.a", "s.b", "s.c", "s.d"],
                           "preferred_name": ["A", "B", "C", "D"]})
    top = sn.select_top(ranked, mapped, top_n=1)
    assert top.set_index("direction").symbol.to_dict() == {"positive": "A", "negative": "D"}
    assert list(top.string_id) == ["s.a", "s.d"]


def test_get_network_drops_reverse_duplicates(monkeypatch):
    tsv = ("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\n"
           "s.a\ts.b\tA\tB\t9606\t0.9\ns.b\ts.a\tB\tA\t9606\t0.9\ns.a\ts.d\tA\tD\t9606\t0.75\n")
    monkeypatch.setattr(sn, "post", lambda path, params, pause=1.0: tsv)
    edges = sn.get_network(["s.a", "s.b", "s.d"], 700)
    assert len(edges) == 2


def test_draw_networks_writes_files(tmp_path):
    nodes = pd.DataFrame({"symbol": ["A", "B", "C"], "weight": [0.5, -0.4, 0.1], "direction": ["positive", "negative", "positive"],
                          "string_id": ["s.a", "s.b", "s.c"], "preferred_name": ["A", "B", "C"]})
    edges = pd.DataFrame({"stringId_A": ["s.a"], "stringId_B": ["s.b"], "score": [0.8]})
    sn.draw_networks({"Toy": (nodes, edges), "Empty": (nodes, edges.iloc[0:0])}, str(tmp_path / "net"))
    assert (tmp_path / "net.png").exists() and (tmp_path / "net.pdf").exists()
