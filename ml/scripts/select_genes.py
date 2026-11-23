#!/usr/bin/env python3
"""
ml/scripts/select_genes.py

Feature selector (spec story 11; issue 04): the n most variable genes in label-free PPMI data
(normalize_cohort.py output), with no labels and no PDBP data. The gene list is written once
and applied unchanged to both cohorts by prepare_features.py.

Checks written to the report (aggregate only, may leave Terra):
  - every participant used is a PPMI analysis-set participant (select_top_variance refuses
    anything else), and with --pdbp_participants, that none of them is a PDBP participant;
  - with --mofa_genes, the overlap with MOFA's 5,000-gene RNA-seq list. MOFA chose its genes
    the same way (variance after size factors and log1p) but before plate correction, so the
    overlap is expected to be high, not complete.

    python3 ml/scripts/select_genes.py --expr results/ml/preprocessed/ppmi_expr.csv \\
        --participants results/ml/analysis_set/ppmi_participants.csv \\
        --pdbp_participants results/ml/analysis_set/pdbp_participants.csv \\
        --mofa_genes <mofa>/data/ppmi_rnaseq_adjusted.csv --n_top 5000 \\
        --genes_out results/ml/preprocessed/genes_top5000.txt \\
        --report results/ml/preprocessed/gene_selection_report.json
"""

import argparse
import json
import os
import re

import pandas as pd

from features import select_top_variance


def read_expr(path: str) -> pd.DataFrame:
    """A genes x participants CSV as participants x genes."""
    return pd.read_csv(path, index_col=0).T.rename_axis("participant_id")


def read_gene_list(path: str) -> list:
    """One gene per line (.txt), or the first column of a CSV (a gene list or a genes x samples matrix)."""
    if path.endswith(".txt"):
        with open(path) as fh:
            return [line.strip() for line in fh if line.strip()]
    return list(pd.read_csv(path, usecols=[0]).iloc[:, 0].astype(str))


def write_gene_list(genes: list, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(genes) + "\n")


def base_gene_id(g: str) -> str:
    """Ensembl ID without version or MOFA view suffix: ENSG00000000003.14_RNA_seq -> ENSG00000000003."""
    return re.sub(r"(_?RNA_seq)+$", "", str(g), flags=re.IGNORECASE).split(".")[0]


def overlap(genes: list, reference: list) -> dict:
    a, b = {base_gene_id(g) for g in genes}, {base_gene_id(g) for g in reference}
    shared_n = len(a & b)
    return {"n_selected": len(a), "n_reference": len(b), "n_shared": shared_n,
            "fraction_of_selected": shared_n / len(a) if a else float("nan"),
            "jaccard": shared_n / len(a | b) if a | b else float("nan")}


def select_genes(expr: pd.DataFrame, meta: pd.DataFrame, n_top: int, pdbp_ids=None,
                 mofa_genes=None):
    """Returns (gene list, report). expr: participants x genes; meta: cohort per participant."""
    genes = select_top_variance(expr, meta, n_top)
    used = expr.index
    report = {
        "method": "variance, label-free PPMI expression (normalize_cohort.py)",
        "n_top": n_top, "n_genes_available": int(expr.shape[1]), "n_selected": len(genes),
        "n_participants": int(len(used)),
        "cohorts": sorted(meta.loc[used, "cohort"].unique()),
        "labels_used": False,
    }
    if pdbp_ids is not None:
        report["n_pdbp_participants_checked"] = int(len(pdbp_ids))
        report["n_pdbp_participants_in_selection"] = int(len(used.intersection(pdbp_ids)))
        if report["n_pdbp_participants_in_selection"]:
            raise ValueError(f"{report['n_pdbp_participants_in_selection']} PDBP participants "
                             f"entered gene selection")
    if mofa_genes is not None:
        report["mofa_overlap"] = overlap(genes, mofa_genes)
    return genes, report


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--expr", required=True, help="PPMI normalize_cohort.py output (genes x participants)")
    p.add_argument("--participants", required=True, help="PPMI analysis_set.py output")
    p.add_argument("--pdbp_participants", help="PDBP analysis_set.py output, to check no PDBP sample is used")
    p.add_argument("--mofa_genes", help="MOFA's RNA-seq gene list (.txt) or a matrix whose rows are its genes")
    p.add_argument("--n_top", type=int, default=5000)
    p.add_argument("--genes_out", required=True, help="Gene list, one per line")
    p.add_argument("--report", required=True, help="Selection report (JSON, aggregate)")
    a = p.parse_args()

    expr = read_expr(a.expr)
    meta = pd.read_csv(a.participants, usecols=["participant_id", "cohort"], index_col="participant_id")
    pdbp_ids = (pd.read_csv(a.pdbp_participants, usecols=["participant_id"])["participant_id"]
                if a.pdbp_participants else None)
    mofa = read_gene_list(a.mofa_genes) if a.mofa_genes else None

    genes, report = select_genes(expr, meta, a.n_top, pdbp_ids=pdbp_ids, mofa_genes=mofa)
    report.update({"input": a.expr, "mofa_genes": a.mofa_genes})
    write_gene_list(genes, a.genes_out)
    with open(a.report, "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report, indent=2))
    print(f"Saved: {a.genes_out}, {a.report}")


if __name__ == "__main__":
    main()
