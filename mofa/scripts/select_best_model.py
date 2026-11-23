#!/usr/bin/env python3
"""
mofa/scripts/select_best_model.py

Identifies the best MOFA seed by final ELBO, and checks two things ELBO ranking alone does not:
  1. Convergence: whether each seed converged or stopped at the iteration cap (mofapy2 truncates
     the stored ELBO trace on convergence, so a full-length trace means the cap was reached), and
     the relative ELBO change over the last recorded checkpoints.
  2. Seed stability: for every factor of the best seed, the best-matching factor in each other seed
     by absolute Pearson correlation of factor scores Z (all seeds share the same samples; this is
     the comparison MOFA2's compare_factors() makes). Factors that are active in the best seed but
     have no close match elsewhere are not reproducible.

Note: mofapy2 (0.7.5) initializes Z from a PCA of the data and W at zero, so the seed only matters
through scikit-learn's randomized PCA solver on large matrices. Seeds are therefore near-identical
restarts, and agreement across them is weak evidence of stability; a subsampling check (retraining
on random subsets of participants) is the stronger test.

Outputs in --model_dir: best_model_path.txt, seed_convergence_summary.csv,
seed_factor_matching.csv (|r| of best match) and seed_factor_matching_index.csv (matched factor).
"""
import os
import argparse
import sys
import pandas as pd
import h5py

import numpy as np

def get_final_elbo(model_path: str) -> float:
    """Extracts final non-NaN ELBO value from the saved MOFA HDF5 file."""
    with h5py.File(model_path, 'r') as f:
        if 'training_stats' in f and 'elbo' in f['training_stats']:
            elbo = f['training_stats']['elbo'][:]
            # Filter out NaNs before getting the final value
            valid_elbo = elbo[~np.isnan(elbo)]
            if len(valid_elbo) > 0:
                return float(valid_elbo[-1])
    return -float('inf')

def get_convergence(model_path: str, default_maxiter: int = None) -> dict:
    """
    Reports iterations run, the iteration cap, whether training converged, and the relative
    ELBO change (%) between the last two recorded ELBO checkpoints.
    """
    with h5py.File(model_path, 'r') as f:
        elbo = f['training_stats']['elbo'][:]
        maxiter = default_maxiter
        if 'training_opts' in f:
            opts = f['training_opts']
            names = [n.decode() if isinstance(n, bytes) else str(n) for n in opts.attrs.get('names', [])]
            if 'maxiter' in names:
                maxiter = int(opts[:][names.index('maxiter')])
    # The trace holds the pre-training ELBO plus one entry per iteration
    n_iter = len(elbo) - 1
    valid = elbo[~np.isnan(elbo)]
    rel_change = np.nan
    if len(valid) >= 2 and valid[-2] != 0:
        rel_change = 100 * abs((valid[-1] - valid[-2]) / valid[-2])
    converged = None if maxiter is None else n_iter < maxiter
    return {"iterations": n_iter, "maxiter": maxiter, "converged": converged,
            "last_rel_elbo_change_pct": rel_change}

def load_z(model_path: str) -> pd.DataFrame:
    """Factor scores Z (samples x factors) across all groups, indexed by sample name."""
    with h5py.File(model_path, 'r') as f:
        frames = []
        for group in f['expectations']['Z'].keys():
            z = f['expectations']['Z'][group][:]
            samples = [s.decode() if isinstance(s, bytes) else str(s) for s in f['samples'][group][:]]
            if z.shape[0] != len(samples):
                z = z.T  # mofapy2 stores Z as (factors, samples)
            frames.append(pd.DataFrame(z, index=samples, columns=[f"Factor{i+1}" for i in range(z.shape[1])]))
    return pd.concat(frames)

def get_max_r2_per_factor(model_path: str) -> pd.Series:
    """Largest per-factor R2 (%) over views and groups, to flag which factors are active."""
    with h5py.File(model_path, 'r') as f:
        if 'variance_explained' not in f or 'r2_per_factor' not in f['variance_explained']:
            return pd.Series(dtype=float)
        r2 = np.max([f['variance_explained']['r2_per_factor'][g][:] for g in f['variance_explained']['r2_per_factor']], axis=(0, 1))
    return pd.Series(r2, index=[f"Factor{i+1}" for i in range(len(r2))])

def match_factors(z_ref: pd.DataFrame, z_other: pd.DataFrame):
    """For each reference factor: best |Pearson r| with any factor of the other seed, and which one."""
    common = z_ref.index.intersection(z_other.index)
    def standardize(m):
        # Inactive factors can be constant (sd 0): leave them as NaN so they match nothing
        sd = m.std(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(sd > 1e-12, (m - m.mean(0)) / sd, np.nan)

    a = standardize(z_ref.loc[common].values)
    b = standardize(z_other.loc[common].values)
    corr = np.abs(a.T @ b / len(common))  # ref factors x other factors (NaN for constant factors)
    best_r = np.nanmax(np.where(np.isnan(corr), -1, corr), axis=1)
    best_idx = np.where(np.isnan(corr), -1, corr).argmax(axis=1)
    matched = [z_other.columns[i] if r >= 0 else None for i, r in zip(best_idx, best_r)]
    return (pd.Series(np.where(best_r >= 0, best_r, np.nan), index=z_ref.columns),
            pd.Series(matched, index=z_ref.columns))

def main():
    parser = argparse.ArgumentParser(description="Identify best MOFA model seed based on final ELBO.")
    parser.add_argument("--model_dir", required=True, help="Directory containing model_seed_*.hdf5 files.")
    parser.add_argument("--maxiter", type=int, default=None, help="Iteration cap, if not stored in the model files.")
    parser.add_argument("--match_threshold", type=float, default=0.8, help="|r| at or above which a factor counts as reproduced (default: 0.8).")
    parser.add_argument("--active_r2", type=float, default=1.0, help="Max per-view R2 (%%) for a factor to count as active (default: 1.0).")
    args = parser.parse_args()

    if not os.path.exists(args.model_dir):
        print(f"Error: Directory '{args.model_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    all_files = os.listdir(args.model_dir)
    print(f"Found {len(all_files)} total file(s) in '{args.model_dir}'")

    matching_files = sorted(f for f in all_files if f.startswith("model_seed_") and f.endswith(".hdf5"))
    print(f"Found {len(matching_files)} matching 'model_seed_*.hdf5' file(s)")

    best_elbo = -float('inf')
    best_model = None
    records = []

    for filename in matching_files:
        path = os.path.join(args.model_dir, filename)
        try:
            elbo = get_final_elbo(path)
            conv = get_convergence(path, args.maxiter)
            status = {True: "converged", False: "stopped at iteration cap", None: "unknown (pass --maxiter)"}[conv["converged"]]
            print(f"{filename}: Final ELBO = {elbo} | {conv['iterations']} iterations, {status} | "
                  f"last ELBO change = {conv['last_rel_elbo_change_pct']:.2e}%")
            records.append({"model": filename, "final_elbo": elbo, **conv})
            if elbo > best_elbo:
                best_elbo = elbo
                best_model = path
        except Exception as e:
            print(f"Could not read {filename}: {e}")

    if not best_model:
        print("No valid models found with non-inf ELBO.")
        return

    print(f"\nBest model identified: {best_model} (ELBO: {best_elbo})")
    with open(os.path.join(args.model_dir, "best_model_path.txt"), "w") as f:
        f.write(best_model)

    summary = pd.DataFrame(records).sort_values("final_elbo", ascending=False)
    summary["elbo_gap_to_best"] = best_elbo - summary["final_elbo"]
    summary.to_csv(os.path.join(args.model_dir, "seed_convergence_summary.csv"), index=False)
    if (summary["converged"] == False).any():
        print("Warning: some seeds stopped at the iteration cap. Check that the best seed's ELBO curve has "
              "flattened, or that its ELBO lead is larger than its last-checkpoint change.")

    # Seed stability: match every factor of the best seed against each other seed
    others = [r["model"] for r in records if os.path.join(args.model_dir, r["model"]) != best_model]
    if not others:
        return
    z_best = load_z(best_model)
    active = get_max_r2_per_factor(best_model)
    match_r, match_idx = {}, {}
    for filename in others:
        r, idx = match_factors(z_best, load_z(os.path.join(args.model_dir, filename)))
        seed = filename.replace("model_", "").replace(".hdf5", "")
        match_r[seed], match_idx[seed] = r, idx
    match_r = pd.DataFrame(match_r)
    match_idx = pd.DataFrame(match_idx)
    if not active.empty:
        match_r.insert(0, "best_model_max_R2", active.reindex(match_r.index))
    match_r.to_csv(os.path.join(args.model_dir, "seed_factor_matching.csv"))
    match_idx.to_csv(os.path.join(args.model_dir, "seed_factor_matching_index.csv"))

    seed_cols = [c for c in match_r.columns if c != "best_model_max_R2"]
    print(f"\nSeed stability: best |r| of each factor of {os.path.basename(best_model)} with any factor of the other seeds")
    print(match_r.round(2).to_string())
    for factor, row in match_r.iterrows():
        is_active = "best_model_max_R2" not in row or row["best_model_max_R2"] >= args.active_r2
        n_rep = int((row[seed_cols] >= args.match_threshold).sum())
        if is_active and n_rep < len(seed_cols):
            print(f"  {factor}: active but reproduced (|r| >= {args.match_threshold}) in only {n_rep} / {len(seed_cols)} other seeds")
    print(f"\nSaved {os.path.join(args.model_dir, 'seed_convergence_summary.csv')}, seed_factor_matching.csv and seed_factor_matching_index.csv")

if __name__ == "__main__":
    main()
