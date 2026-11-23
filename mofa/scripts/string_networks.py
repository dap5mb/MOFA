"""
STRING protein-network views of a factor's top-weighted genes (secondary to GSEA preranked).

For each factor and each weight direction, the top genes that map to STRING proteins are taken
from <factor>_ranked_genes.csv (written by pathway_enrichment.py). The script then:
  - retrieves STRING interactions among them (default: combined score >= 0.7, high confidence);
  - runs STRING functional enrichment for each direction against the model's own genes as the
    background (all mapped genes of the ranked list), not the whole genome;
  - draws one network panel per factor, coloured by weight direction.

Top-weight cut-offs are used only here, because top-gene overlap across subsamples is moderate;
the primary enrichment is GSEA on the full weight vector. Genes without a STRING protein
(immunoglobulin segments, most lncRNAs and pseudogenes) cannot appear in the networks.

Needs internet (string-db.org API). Usage:
  python3 mofa/scripts/string_networks.py --enrichment_dir results/mofa/main_noplasma/enrichment \
    --factors Factor7 Factor10 --top_n 50 --output results/mofa/main_noplasma/string
"""

import argparse
import io
import json
import os
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

API = "https://string-db.org/api"
CALLER = "mofa_string_views"


def post(path: str, params: dict, pause: float = 1.0) -> str:
    """POST to the STRING API, pausing between calls as STRING asks."""
    data = urllib.parse.urlencode({**params, "caller_identity": CALLER}).encode()
    req = urllib.request.Request(API + path, data=data, headers={"User-Agent": "mofa-string/1.0"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                text = r.read().decode()
            time.sleep(pause)
            return text
        except urllib.error.HTTPError as e:
            if attempt == 3 or e.code < 500:
                raise
            time.sleep(5 * (attempt + 1))


def map_to_string(symbols: list, species: int = 9606, batch: int = 2000) -> pd.DataFrame:
    """Gene symbol -> STRING protein ID (best match). Unmapped symbols are absent from the result."""
    rows = []
    for i in range(0, len(symbols), batch):
        hits = json.loads(post("/json/get_string_ids", {"identifiers": "\r".join(symbols[i:i + batch]),
                                                        "species": species, "limit": 1, "echo_query": 1}))
        rows += [{"symbol": h["queryItem"], "string_id": h["stringId"], "preferred_name": h["preferredName"]}
                 for h in hits]
    return pd.DataFrame(rows, columns=["symbol", "string_id", "preferred_name"]).drop_duplicates("symbol")


def select_top(ranked: pd.Series, mapped: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Top `top_n` STRING-mapped genes by weight in each direction."""
    m = mapped.drop_duplicates("string_id").set_index("symbol")
    df = pd.DataFrame({"symbol": ranked.index, "weight": ranked.values})
    df = df[df.symbol.isin(m.index)]
    pos = df[df.weight > 0].nlargest(top_n, "weight").assign(direction="positive")
    neg = df[df.weight < 0].nsmallest(top_n, "weight").assign(direction="negative")
    out = pd.concat([pos, neg], ignore_index=True)
    out["string_id"] = out.symbol.map(m["string_id"])
    out["preferred_name"] = out.symbol.map(m["preferred_name"])
    return out


def get_network(string_ids: list, required_score: int, species: int = 9606) -> pd.DataFrame:
    text = post("/tsv/network", {"identifiers": "\r".join(string_ids), "species": species,
                                 "required_score": required_score})
    edges = pd.read_csv(io.StringIO(text), sep="\t")
    if edges.empty:
        return pd.DataFrame(columns=["stringId_A", "stringId_B", "score"])
    edges = edges[["stringId_A", "stringId_B", "preferredName_A", "preferredName_B", "score"]]
    key = edges[["stringId_A", "stringId_B"]].apply(lambda r: tuple(sorted(r)), axis=1)
    return edges.loc[~key.duplicated()].reset_index(drop=True)


def get_enrichment(string_ids: list, background_ids: list, species: int = 9606) -> pd.DataFrame:
    text = post("/json/enrichment", {"identifiers": "\r".join(string_ids), "species": species,
                                     "background_string_identifiers": "\r".join(background_ids)})
    res = pd.DataFrame(json.loads(text))
    if res.empty:
        return res
    keep = ["category", "term", "description", "number_of_genes", "number_of_genes_in_background", "fdr", "inputGenes"]
    return res[[c for c in keep if c in res.columns]]


def packed_layout(g, gap: float = 0.9) -> dict:
    """Lays out each connected component on its own (Kamada-Kawai), sized by its node count, and packs
    the components in rows, largest first, so one big cluster isn't crushed by distant small ones."""
    import networkx as nx

    comps = sorted(nx.connected_components(g), key=len, reverse=True)
    boxes = []
    for comp in comps:
        sub = g.subgraph(comp)
        pos = nx.kamada_kawai_layout(sub) if len(comp) > 2 else nx.circular_layout(sub)
        xy = np.array(list(pos.values()))
        xy = xy - xy.min(axis=0)
        span = xy.max(axis=0)
        span[span == 0] = 1
        side = 0.9 * np.sqrt(len(comp))
        xy = xy / span.max() * side
        boxes.append((list(pos), xy, xy.max(axis=0)))
    # The largest component on the left; the rest packed in rows in a column to its right
    out = {n: (px, -boxes[0][2][1] + py) for n, (px, py) in zip(boxes[0][0], boxes[0][1])}
    x0 = boxes[0][2][0] + 1.5 * gap
    width = max(boxes[0][2][0], 4.0)
    x, y, row_h = x0, 0.0, 0.0
    for nodes, xy, size in boxes[1:]:
        if x > x0 and x + size[0] > x0 + width:
            x, y, row_h = x0, y - row_h - gap, 0.0
        for n, (px, py) in zip(nodes, xy):
            out[n] = (x + px, y - size[1] + py)
        x += size[0] + gap
        row_h = max(row_h, size[1])
    return out


def draw_networks(panels: dict, out_path: str, pos_color="#eb6834", neg_color="#2a78d6"):
    """panels: {title: (nodes_df, edges_df)}; one network panel each, isolated nodes counted, not drawn."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "pdf.fonttype": 42, "savefig.dpi": 300})
    fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 3.3))
    for ax, (title, (nodes, edges)) in zip(np.atleast_1d(axes), panels.items()):
        g = nx.Graph()
        info = nodes.set_index("string_id")
        for _, e in edges.iterrows():
            g.add_edge(e.stringId_A, e.stringId_B, score=e.score)
        n_iso = len(nodes) - g.number_of_nodes()
        ax.set_title(title, loc="left", fontsize=9)
        ax.set_axis_off()
        if g.number_of_nodes() == 0:
            ax.text(0.5, 0.5, "No interactions at this confidence", ha="center", transform=ax.transAxes, fontsize=8)
            continue
        layout = packed_layout(g)
        w = info.weight.abs()
        size = 30 + 170 * (w - w.min()) / (w.max() - w.min() + 1e-12)
        nx.draw_networkx_edges(g, layout, ax=ax, edge_color="#c9c8c3",
                               width=[0.4 + 1.2 * (d["score"] - 0.7) / 0.3 for _, _, d in g.edges(data=True)])
        nlist = list(g.nodes)
        nx.draw_networkx_nodes(g, layout, nodelist=nlist, ax=ax, linewidths=0.8, edgecolors="white",
                               node_color=[pos_color if info.loc[n, "weight"] > 0 else neg_color for n in nlist],
                               node_size=[size.loc[n] for n in nlist])
        nx.draw_networkx_labels(g, {n: (x, y + 0.17) for n, (x, y) in layout.items()}, ax=ax, font_size=5.5,
                                labels={n: info.loc[n, "preferred_name"] for n in nlist}, font_color="#0b0b0b")
        ax.set_aspect("equal")
        ax.margins(0.08)
        n_pos = int((info.loc[nlist, "weight"] > 0).sum())
        ax.text(0.0, -0.02, f"{len(nlist)} of {len(nodes)} genes connected ({n_pos} positive, {len(nlist) - n_pos} "
                f"negative weight); {n_iso} without a high-confidence partner not shown.", transform=ax.transAxes,
                fontsize=6.5, color="#52514e", va="top")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", ms=6, mfc=c, mec="white", label=l)
               for c, l in ((pos_color, "Positive weight"), (neg_color, "Negative weight"))]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_path}.{ext}", bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--enrichment_dir", required=True, help="Folder with <factor>_ranked_genes.csv")
    p.add_argument("--factors", nargs="+", default=["Factor7", "Factor10"])
    p.add_argument("--top_n", type=int, default=50, help="STRING-mapped genes per direction")
    p.add_argument("--required_score", type=int, default=700)
    p.add_argument("--titles", nargs="*", default=None, help="Panel titles, one per factor")
    p.add_argument("--output", required=True)
    p.add_argument("--figure", default=None, help="Figure path without extension (default: <output>/string_networks)")
    args = p.parse_args()
    os.makedirs(args.output, exist_ok=True)

    ranked = {f: pd.read_csv(os.path.join(args.enrichment_dir, f"{f}_ranked_genes.csv"), index_col=0)["weight"]
              for f in args.factors}
    universe = sorted(set().union(*[set(r.index) for r in ranked.values()]))
    print(f"Mapping {len(universe)} gene symbols to STRING ...")
    mapped = map_to_string(universe)
    mapped.to_csv(os.path.join(args.output, "string_id_map.csv"), index=False)
    background = sorted(mapped.string_id.unique())
    print(f"{len(mapped)} of {len(universe)} symbols mapped; background = {len(background)} STRING proteins")

    panels, enrich_all = {}, []
    titles = args.titles or args.factors
    for f, title in zip(args.factors, titles):
        nodes = select_top(ranked[f], mapped, args.top_n)
        edges = get_network(list(nodes.string_id), args.required_score)
        nodes.to_csv(os.path.join(args.output, f"{f}_network_nodes.csv"), index=False)
        edges.to_csv(os.path.join(args.output, f"{f}_network_edges.csv"), index=False)
        print(f"\n{f}: {len(nodes)} genes, {len(edges)} interactions (score >= {args.required_score / 1000:.1f})")
        for direction, g in nodes.groupby("direction"):
            res = get_enrichment(list(g.string_id), background)
            n_sig = 0 if res.empty else int((res.fdr < 0.05).sum())
            print(f"  {direction} weights ({len(g)} genes): {n_sig} terms with FDR < 0.05")
            if n_sig:
                top = res[res.fdr < 0.05].sort_values("fdr").head(8)
                for _, r in top.iterrows():
                    print(f"    {r.category:10} FDR {r.fdr:.1e}  {r.number_of_genes:>3} genes  {r.description}")
                enrich_all.append(res.assign(factor=f, direction=direction))
        panels[title] = (nodes, edges)

    if enrich_all:
        pd.concat(enrich_all, ignore_index=True).to_csv(os.path.join(args.output, "string_enrichment.csv"), index=False)
    fig_path = args.figure or os.path.join(args.output, "string_networks")
    draw_networks(panels, fig_path)
    print(f"\nSaved tables to {args.output} and figure to {fig_path}.pdf/.png")


if __name__ == "__main__":
    main()
