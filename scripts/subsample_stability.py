#!/usr/bin/env python3
"""
mofa/scripts/subsample_stability.py

Factor stability by subsampling. mofapy2 initializes factors from a PCA of the data, so seed
restarts are near-identical (see select_best_model.py); retraining on random subsets of
participants is the meaningful robustness check.

For each subset:
  1. Draw `fraction` of participants within each MOFA group (sex), without replacement.
  2. Train MOFA on that subset with the same settings as the full model (train_mofa.build_and_train).
  3. Match every factor of the full (reference) model to its closest subset factor by |Pearson r|
     of factor scores Z on the shared participants (sign and order invariant).
  4. For the matched pair, compare weights per view: |r| of the weight vectors and the Jaccard
     overlap of the top `top_pct` features by |weight| (the features pathway enrichment would use).

Subset models are saved as subset_<i>.hdf5 and reused if present, so an interrupted run resumes.

Outputs in --output:
  subset_<i>_participants.txt           participants in each subset
  subset_<i>.hdf5                       trained subset models
  subset_factor_matching_long.csv       one row per subset x reference factor x view
  subset_stability_summary.csv          one row per reference factor
"""

import os
import argparse
import numpy as np
import pandas as pd

from select_best_model import load_z, match_factors, get_max_r2_per_factor
from project_pdbp import load_model_weights


def draw_subset(data: pd.DataFrame, fraction: float, rng: np.random.Generator) -> list:
    """Participants sampled without replacement within each group, keeping the sex balance."""
    pairs = data[["sample", "group"]].drop_duplicates()
    chosen = []
    for _, sub in pairs.groupby("group"):
        n = int(round(fraction * len(sub)))
        chosen.extend(rng.choice(sub["sample"].values, size=n, replace=False).tolist())
    return sorted(chosen)


def top_jaccard(w_ref: pd.Series, w_sub: pd.Series, top_pct: float) -> float:
    n = max(1, int(len(w_ref) * top_pct))
    a = set(w_ref.abs().nlargest(n).index)
    b = set(w_sub.abs().nlargest(n).index)
    return len(a & b) / len(a | b)


def compare_to_reference(ref_path: str, sub_path: str, subset_id: int, top_pct: float) -> pd.DataFrame:
    z_r, z_idx = match_factors(load_z(ref_path), load_z(sub_path))
    w_ref, w_sub = load_model_weights(ref_path), load_model_weights(sub_path)
    rows = []
    for factor in z_r.index:
        matched = z_idx[factor]
        for view in w_ref:
            row = {"subset": subset_id, "factor": factor, "matched_factor": matched,
                   "z_abs_r": z_r[factor], "view": view, "w_abs_r": np.nan, "top_jaccard": np.nan}
            if isinstance(matched, str) and view in w_sub:
                a = w_ref[view][factor]
                b = w_sub[view].reindex(a.index)[matched]
                if a.std() > 0 and b.std() > 0:
                    row["w_abs_r"] = abs(np.corrcoef(a, b)[0, 1])
                    row["top_jaccard"] = top_jaccard(a, b, top_pct)
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Assess MOFA factor stability by retraining on participant subsets.")
    parser.add_argument("--input", required=True, help="Long-format MOFA input used for the full model (sample, feature, value, view, group).")
    parser.add_argument("--reference_model", required=True, help="Full-data model HDF5 to compare against (e.g. the best seed).")
    parser.add_argument("--output", default="ppmi_mofa_outputs/subsample_stability", help="Output directory.")
    parser.add_argument("--n_subsets", type=int, default=10, help="Number of subsets (default: 10).")
    parser.add_argument("--fraction", type=float, default=0.8, help="Fraction of participants per subset, within each group (default: 0.8).")
    parser.add_argument("--n_iter", type=int, default=1000, help="Maximum training iterations (default: 1000, as for the full model).")
    parser.add_argument("--seed", type=int, default=0, help="Seed for drawing subsets (default: 0).")
    parser.add_argument("--top_pct", type=float, default=0.05, help="Top fraction of features by |weight| for overlap (default: 0.05).")
    parser.add_argument("--match_threshold", type=float, default=0.8, help="|r| at or above which a factor counts as reproduced (default: 0.8).")
    parser.add_argument("--active_r2", type=float, default=1.0, help="Max per-view R2 (%%) for a reference factor to count as active (default: 1.0).")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    data = pd.read_csv(args.input)
    ref_weights = load_model_weights(args.reference_model)
    n_factors = next(iter(ref_weights.values())).shape[1]
    print(f"Reference model: {args.reference_model} ({n_factors} factors)")
    print(f"Participants: {data['sample'].nunique()} | groups: {data.groupby('group')['sample'].nunique().to_dict()}")

    rng = np.random.default_rng(args.seed)
    results = []
    for i in range(args.n_subsets):
        subset = draw_subset(data, args.fraction, rng)  # drawn every time so subsets are identical on resume
        with open(os.path.join(args.output, f"subset_{i}_participants.txt"), "w") as f:
            f.write("\n".join(subset))

        model_path = os.path.join(args.output, f"subset_{i}.hdf5")
        if os.path.exists(model_path):
            print(f"\n--- Subset {i}: reusing {model_path} ---")
        else:
            print(f"\n--- Subset {i}: training on {len(subset)} participants ---")
            from train_mofa import build_and_train, SAVED_EXPECTATIONS
            ent = build_and_train(data[data["sample"].isin(subset)], n_factors=n_factors, n_iter=args.n_iter, seed=i)
            ent.save(model_path, expectations=SAVED_EXPECTATIONS)

        res = compare_to_reference(args.reference_model, model_path, i, args.top_pct)
        results.append(res)
        z_summary = res.drop_duplicates("factor").set_index("factor")["z_abs_r"]
        print(f"Subset {i} |r| of factor scores: " + ", ".join(f"{k}={v:.2f}" for k, v in z_summary.items()))

    long = pd.concat(results, ignore_index=True)
    long.to_csv(os.path.join(args.output, "subset_factor_matching_long.csv"), index=False)

    # One row per reference factor: Z agreement across subsets, mean weight agreement over views
    z = long.drop_duplicates(["subset", "factor"]).pivot(index="factor", columns="subset", values="z_abs_r")
    factors = list(next(iter(ref_weights.values())).columns)
    z = z.reindex(factors)
    w = long.groupby("factor")[["w_abs_r", "top_jaccard"]].median().reindex(factors)
    summary = pd.DataFrame({
        "ref_max_R2": get_max_r2_per_factor(args.reference_model).reindex(factors),
        "z_abs_r_median": z.median(axis=1),
        "z_abs_r_min": z.min(axis=1),
        "frac_subsets_reproduced": (z >= args.match_threshold).mean(axis=1),
        "w_abs_r_median": w["w_abs_r"],
        f"top{int(args.top_pct * 100)}pct_jaccard_median": w["top_jaccard"],
    })
    summary["stable"] = (summary["frac_subsets_reproduced"] >= 0.9) & (summary["ref_max_R2"] >= args.active_r2)
    summary.to_csv(os.path.join(args.output, "subset_stability_summary.csv"))

    print(f"\n=== Subsample stability ({args.n_subsets} subsets of {args.fraction:.0%} participants) ===")
    print(summary.round(2).to_string())
    print(f"\nA factor is marked stable if it is active (max per-view R2 >= {args.active_r2}%) and matched with "
          f"|r| >= {args.match_threshold} in at least 90% of subsets.")
    print(f"Saved results to: {args.output}")


if __name__ == "__main__":
    main()


###### Example use (from ~/mofa on Terra, after selecting the best seed) ##########
#python3 scripts/subsample_stability.py \
#  --input ppmi_mofa_inputs/ppmi_mofa_input_long.csv \
#  --reference_model $(cat ppmi_mofa_outputs/best_model_path.txt) \
#  --output ppmi_mofa_outputs/subsample_stability \
#  --n_subsets 10
