"""
Checks MOFA factor scores against RNA-seq library size and lists one factor's top features.

For each cohort:
  - library size (total counts) and genes detected (count > 0) per participant, from the raw
    count matrix (genes x participants);
  - Spearman correlation of every factor with log10 library size and genes detected;
  - logit(PD) ~ factor + age + male, refitted with log10 library size and genes detected
    added, to see whether the factor's PD association survives adjustment.
Then the factor's top RNA-seq and CSF features by absolute weight.

Usage (from ~/mofa):
  python3 scripts/factor_technical_check.py \
    --ppmi_z ppmi_mofa_outputs_noplasma/latent_factors_Z_with_meta.csv --ppmi_metadata data/metadata_m0_PPMI.csv \
    --ppmi_counts <PPMI raw counts csv> \
    --pdbp_z pdbp_mofa_outputs_noplasma/pdbp_latent_factors_Z.csv --pdbp_metadata data/metadata_m0_PDBP.csv \
    --pdbp_counts data/pdbp/raw/pdbp_rnaseq.csv \
    --weights ppmi_mofa_outputs_noplasma/feature_weights_W_full.csv \
    --factor Factor4 --output ppmi_mofa_outputs_noplasma/technical_check
"""

import argparse
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import stats

from clinical_models import fit, load_cohort, prepare


def library_stats(counts_path: str, chunksize: int = 5000) -> pd.DataFrame:
    """Total counts and genes detected per participant, read in chunks of genes."""
    lib, detected = None, None
    for chunk in pd.read_csv(counts_path, index_col=0, chunksize=chunksize, low_memory=False):
        chunk = chunk.apply(pd.to_numeric, errors="coerce").fillna(0)
        s, d = chunk.sum(axis=0), (chunk > 0).sum(axis=0)
        lib = s if lib is None else lib.add(s, fill_value=0)
        detected = d if detected is None else detected.add(d, fill_value=0)
    out = pd.DataFrame({"log10_lib_size": np.log10(lib.clip(lower=1)), "genes_detected": detected})
    out.index = out.index.astype(str)
    out.index.name = "participant_id"
    return out.reset_index()


def technical_correlations(df: pd.DataFrame, factors: list, tech_cols: list) -> pd.DataFrame:
    rows = []
    for f in factors:
        for t in tech_cols:
            d = df[[f, t]].dropna()
            rho, p = stats.spearmanr(d[f], d[t])
            rows.append({"factor": f, "covariate": t, "spearman_rho": rho, "p": p, "n": len(d)})
    return pd.DataFrame(rows)


def dx_with_and_without_tech(df: pd.DataFrame, factors: list) -> pd.DataFrame:
    """logit(PD) ~ factor + age + male, without and with library-size covariates (log-odds per SD)."""
    d = df[df["pd"].notna()].copy()
    d["lib_z"] = (d["log10_lib_size"] - d["log10_lib_size"].mean()) / d["log10_lib_size"].std()
    d["det_z"] = (d["genes_detected"] - d["genes_detected"].mean()) / d["genes_detected"].std()
    rows = []
    for f in factors:
        base = fit(f"pd ~ {f} + age + male", d, f, [f, "pd", "age", "male", "lib_z", "det_z"], logistic=True)
        adj = fit(f"pd ~ {f} + age + male + lib_z + det_z", d, f, [f, "pd", "age", "male", "lib_z", "det_z"], logistic=True)
        rows.append({"factor": f, "estimate": base["estimate"], "p": base["p"],
                     "estimate_adj_lib": adj["estimate"], "p_adj_lib": adj["p"], "n": adj["n"]})
    return pd.DataFrame(rows)


def top_features(weights_path: str, factor: str, n_top: int) -> pd.DataFrame:
    w = pd.read_csv(weights_path)
    w = w.drop_duplicates(subset=["view", "feature"])  # weights are shared across groups
    out = []
    for view, g in w.groupby("view"):
        n = n_top if len(g) > 1000 else max(5, n_top // 3)
        g = g.assign(abs_w=g[factor].abs()).sort_values("abs_w", ascending=False)
        for sign, part in (("positive", g[g[factor] > 0]), ("negative", g[g[factor] < 0])):
            part = part.head(n)
            out.append(pd.DataFrame({"view": view, "sign": sign, "rank": range(1, len(part) + 1),
                                     "feature": part["feature"].values, "weight": part[factor].values}))
    return pd.concat(out, ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for c in ("ppmi", "pdbp"):
        p.add_argument(f"--{c}_z", required=True)
        p.add_argument(f"--{c}_metadata", required=True)
        p.add_argument(f"--{c}_counts", required=True, help="Raw counts, genes x participants")
    p.add_argument("--weights", required=True, help="feature_weights_W_full.csv from extract_factors.py")
    p.add_argument("--factor", default="Factor4")
    p.add_argument("--n_top", type=int, default=30)
    p.add_argument("--case_levels", nargs="+", default=["PD"])
    p.add_argument("--control_levels", nargs="+", default=["Ctrl"])
    p.add_argument("--output", required=True)
    args = p.parse_args()
    os.makedirs(args.output, exist_ok=True)

    prep_args = SimpleNamespace(diagnosis_col="diagnosis", case_levels=args.case_levels,
                                control_levels=args.control_levels, age_col="age_at_baseline", scores={})
    corr_all, dx_all = [], []
    for cohort in ("PPMI", "PDBP"):
        c = cohort.lower()
        df = load_cohort(getattr(args, f"{c}_z"), getattr(args, f"{c}_metadata"))
        factors = [col for col in df.columns if col.startswith("Factor")]
        print(f"\n{cohort}: reading library sizes from {getattr(args, f'{c}_counts')} ...")
        lib = library_stats(getattr(args, f"{c}_counts"))
        df = df.merge(lib, on="participant_id", how="left")
        print(f"{cohort}: library size for {df['log10_lib_size'].notna().sum()} of {len(df)} participants")
        df = prepare(df, factors, prep_args)

        corr = technical_correlations(df, factors, ["log10_lib_size", "genes_detected"]).assign(cohort=cohort)
        dx = dx_with_and_without_tech(df, factors).assign(cohort=cohort)
        corr_all.append(corr)
        dx_all.append(dx)

        print(f"\n{cohort}: Spearman rho of each factor with library size and genes detected")
        print(corr.pivot(index="factor", columns="covariate", values="spearman_rho").reindex(factors).round(3).to_string())
        print(f"\n{cohort}: PD log-odds per SD, without and with library-size adjustment")
        print(dx.set_index("factor")[["estimate", "p", "estimate_adj_lib", "p_adj_lib", "n"]].round(4).to_string())

    pd.concat(corr_all).to_csv(os.path.join(args.output, "factor_library_size_correlations.csv"), index=False)
    pd.concat(dx_all).to_csv(os.path.join(args.output, "pd_association_library_adjusted.csv"), index=False)

    top = top_features(args.weights, args.factor, args.n_top)
    top.to_csv(os.path.join(args.output, f"{args.factor}_top_features.csv"), index=False)
    print(f"\nTop {args.factor} features by absolute weight:")
    for (view, sign), g in top.groupby(["view", "sign"]):
        print(f"  {view} ({sign}): " + ", ".join(g["feature"].astype(str)))
    print(f"\nSaved results to: {args.output}")


if __name__ == "__main__":
    main()
