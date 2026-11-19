"""
Rank-based pathway enrichment (GSEA preranked) on MOFA factor weights.

For each factor, every feature of one view is ranked by its weight (the full weight vector,
not a top-k cut-off) and tested against gene-set libraries (GMT files with gene symbols, e.g.
MSigDB Hallmark, Reactome, GO Biological Process). The universe is the model's features
(5,000 variable genes for RNA-seq), so results describe which pathways a factor separates
among those genes, not genome-wide enrichment.

Sign: a positive NES means the gene set is concentrated among positive weights, i.e. higher
expression in participants with higher factor scores. Factor signs are arbitrary, so read NES
together with the factor's clinical association (e.g. Factor 7 is lower in PD).

Ensembl gene IDs are mapped to symbols with MyGene.info (needs internet); the mapping is
cached in --gene_map so later runs work offline.

Usage (from ~/mofa or locally):
  python3 scripts/pathway_enrichment.py --weights feature_weights_W_full.csv \
    --factors Factor7 Factor10 --gene_sets h.all.gmt c2.cp.reactome.gmt c5.go.bp.gmt \
    --gene_map gene_map.csv --output enrichment
"""

import argparse
import json
import os
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd


def load_view_weights(weights_path: str, view: str) -> pd.DataFrame:
    """Weights of one view, one row per feature (weights are shared across groups)."""
    w = pd.read_csv(weights_path)
    w = w[w["view"] == view].drop_duplicates(subset=["feature"])
    if w.empty:
        raise ValueError(f"No features for view '{view}' in {weights_path}")
    return w.set_index("feature")[[c for c in w.columns if c.startswith("Factor")]]


def strip_version(ensembl_id: str) -> str:
    return str(ensembl_id).split(".")[0]


def query_mygene(ids: list, batch: int = 1000) -> dict:
    out = {}
    for i in range(0, len(ids), batch):
        data = urllib.parse.urlencode({"q": ",".join(ids[i:i + batch]), "scopes": "ensembl.gene",
                                       "fields": "symbol", "species": "human"}).encode()
        req = urllib.request.Request("https://mygene.info/v3/query", data=data,
                                     headers={"User-Agent": "mofa-enrichment/1.0"})
        for hit in json.loads(urllib.request.urlopen(req, timeout=120).read()):
            if "symbol" in hit and hit["query"] not in out:
                out[hit["query"]] = hit["symbol"]
    return out


def map_to_symbols(features: list, cache_path: str = None) -> pd.Series:
    """Ensembl ID (version stripped) -> gene symbol, using and extending a CSV cache."""
    ids = sorted({strip_version(f) for f in features})
    cache = {}
    if cache_path and os.path.exists(cache_path):
        cache = pd.read_csv(cache_path).dropna().set_index("ensembl_id")["symbol"].to_dict()
    missing = [i for i in ids if i not in cache]
    if missing:
        print(f"Querying MyGene.info for {len(missing)} gene IDs ...")
        found = query_mygene(missing)
        cache.update(found)
        if cache_path:
            pd.DataFrame({"ensembl_id": list(cache), "symbol": list(cache.values())}).to_csv(cache_path, index=False)
    return pd.Series({i: cache.get(i) for i in ids})


def ranked_symbols(weights: pd.Series, symbols: pd.Series) -> pd.Series:
    """Weights indexed by symbol; unmapped IDs dropped, duplicate symbols keep the largest |weight|."""
    df = pd.DataFrame({"w": weights.values, "symbol": [symbols.get(strip_version(f)) for f in weights.index]})
    df = df.dropna(subset=["symbol"])
    df = df.loc[df["w"].abs().groupby(df["symbol"]).idxmax()]
    return df.set_index("symbol")["w"].sort_values(ascending=False)


def read_gmt(path: str) -> dict:
    sets = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = parts[2:]
    return sets


def normalize_result(res: pd.DataFrame) -> pd.DataFrame:
    """Common columns across gseapy versions: term, es, nes, pval, fdr, size, leading_edge."""
    r = res.reset_index() if "Term" not in res.columns else res.copy()
    r = r.rename(columns={"Term": "term", "ES": "es", "NES": "nes", "NOM p-val": "pval", "FDR q-val": "fdr",
                          "Lead_genes": "leading_edge", "ledge_genes": "leading_edge", "matched_size": "size",
                          "Tag %": "tag_pct"})
    if "size" not in r.columns and "tag_pct" in r.columns:
        r["size"] = r["tag_pct"].astype(str).str.split("/").str[1].astype(float)
    for c in ("es", "nes", "pval", "fdr", "size"):
        r[c] = pd.to_numeric(r[c], errors="coerce")
    return r[["term", "es", "nes", "pval", "fdr", "size", "leading_edge"]]


def run_prerank(rnk: pd.Series, gene_sets: dict, permutations: int, min_size: int, max_size: int,
                seed: int, processes: int) -> pd.DataFrame:
    import gseapy
    res = gseapy.prerank(rnk=rnk, gene_sets=gene_sets, permutation_num=permutations, min_size=min_size,
                         max_size=max_size, outdir=None, seed=seed, no_plot=True, processes=processes)
    return normalize_result(res.res2d)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--weights", required=True, help="feature_weights_W_full.csv from extract_factors.py")
    p.add_argument("--factors", nargs="+", default=["Factor7", "Factor10"])
    p.add_argument("--view", default="RNA_seq")
    p.add_argument("--gene_sets", nargs="+", required=True, help="GMT files (gene symbols)")
    p.add_argument("--gene_map", default=None, help="CSV cache of ensembl_id,symbol")
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--min_size", type=int, default=15)
    p.add_argument("--max_size", type=int, default=500)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--processes", type=int, default=4)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    os.makedirs(args.output, exist_ok=True)

    w = load_view_weights(args.weights, args.view)
    symbols = map_to_symbols(list(w.index), args.gene_map)
    print(f"{args.view}: {len(w)} features, {symbols.notna().sum()} mapped to gene symbols")

    results = []
    for factor in args.factors:
        rnk = ranked_symbols(w[factor], symbols)
        rnk.rename("weight").to_frame().to_csv(os.path.join(args.output, f"{factor}_ranked_genes.csv"))
        for gmt in args.gene_sets:
            library = os.path.basename(gmt).split(".v20")[0]
            res = run_prerank(rnk, read_gmt(gmt), args.permutations, args.min_size, args.max_size,
                              args.seed, args.processes)
            results.append(res.assign(factor=factor, library=library))
            n_sig = int((res["fdr"] < 0.05).sum())
            print(f"{factor} / {library}: {len(res)} gene sets tested, {n_sig} with FDR < 0.05")

    allres = pd.concat(results, ignore_index=True).sort_values(["factor", "library", "fdr", "pval"])
    allres = allres[["factor", "library", "term", "nes", "es", "pval", "fdr", "size", "leading_edge"]]
    allres.to_csv(os.path.join(args.output, "gsea_prerank_results.csv"), index=False)

    for factor in args.factors:
        top = allres[(allres.factor == factor) & (allres.fdr < 0.05)].sort_values("nes")
        print(f"\n{factor}: gene sets with FDR < 0.05 (most negative to most positive NES)")
        if top.empty:
            print("  none")
        for lib, g in top.groupby("library"):
            print(f"  {lib}:")
            for _, r in pd.concat([g.head(8), g.tail(8)]).drop_duplicates("term").iterrows():
                print(f"    NES {r.nes:+.2f}  FDR {r.fdr:.3f}  n={int(r['size'])}  {r.term}")
    print(f"\nSaved results to: {args.output}")


if __name__ == "__main__":
    main()
