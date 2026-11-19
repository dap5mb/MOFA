#!/usr/bin/env python3
"""
mofa/scripts/compare_models.py

Matches the factors of two MOFA models trained on the same participants, e.g. the sex-grouped
main model and the single-group sensitivity model (pdbp-projection issue 06).

For each factor of one model, finds the factor of the other model with the highest |Pearson r|
of factor scores Z on the shared participants (sign and order invariant), then compares the
matched weights per view: |r| of the weight vectors and the Jaccard overlap of the top
`top_pct` features by |weight|. Matching is run in both directions, because a factor can have
no counterpart (e.g. a sex factor that only the single-group model can form).

Outputs in --output: factor_matching_<a>_to_<b>.csv and factor_matching_<b>_to_<a>.csv.
"""

import os
import argparse
import pandas as pd

from subsample_stability import compare_to_reference
from select_best_model import get_max_r2_per_factor


def match(ref_path: str, other_path: str, top_pct: float) -> pd.DataFrame:
    long = compare_to_reference(ref_path, other_path, 0, top_pct).drop(columns="subset")
    # One row per factor: the match, then weight agreement for each view as its own columns
    wide = long.drop_duplicates("factor")[["factor", "matched_factor", "z_abs_r"]].reset_index(drop=True)
    for view, sub in long.groupby("view"):
        wide = wide.merge(sub[["factor", "w_abs_r", "top_jaccard"]].rename(
            columns={"w_abs_r": f"w_abs_r_{view}", "top_jaccard": f"top_jaccard_{view}"}), on="factor", how="left")
    wide.insert(1, "max_R2", wide["factor"].map(get_max_r2_per_factor(ref_path)))
    order = {f: int(f.replace("Factor", "")) for f in wide["factor"]}
    return wide.sort_values("factor", key=lambda s: s.map(order)).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description="Match factors between two MOFA models trained on the same participants.")
    parser.add_argument("--model_a", required=True, help="First model HDF5 (e.g. the sex-grouped main model).")
    parser.add_argument("--model_b", required=True, help="Second model HDF5 (e.g. the single-group model).")
    parser.add_argument("--label_a", default="grouped", help="Short name for model A (default: grouped).")
    parser.add_argument("--label_b", default="single", help="Short name for model B (default: single).")
    parser.add_argument("--output", default="model_comparison", help="Output directory.")
    parser.add_argument("--top_pct", type=float, default=0.05, help="Top fraction of features by |weight| for overlap (default: 0.05).")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    for (ref, ref_label), (other, other_label) in (((args.model_a, args.label_a), (args.model_b, args.label_b)),
                                                   ((args.model_b, args.label_b), (args.model_a, args.label_a))):
        res = match(ref, other, args.top_pct)
        out = os.path.join(args.output, f"factor_matching_{ref_label}_to_{other_label}.csv")
        res.to_csv(out, index=False)
        print(f"\n=== Each {ref_label} factor's best match in the {other_label} model ===")
        print(res.round(2).to_string(index=False))
        print(f"Saved to: {out}")


if __name__ == "__main__":
    main()


###### Example use (from ~/mofa on Terra) ##########
#python3 scripts/compare_models.py \
#  --model_a ppmi_mofa_outputs_noplasma/model_seed_0.hdf5 \
#  --model_b ppmi_mofa_outputs_singlegroup/model_seed_0.hdf5 \
#  --output ppmi_mofa_outputs_singlegroup/model_comparison
