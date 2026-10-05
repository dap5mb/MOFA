#!/usr/bin/env python3
"""
mofa/scripts/project_pdbp.py

Projects preprocessed PDBP validation omic data into the trained PPMI MOFA factor space
using weights directly loaded from the HDF5 model file.

The projection mirrors what mofapy2 does to the training data:
  1. Center every feature within each group (sex)            -> center_groups=True
  2. Scale each view to unit variance                         -> scale_views=True
  3. Per sample, use only the features that were observed (no zero-filling) and weight
     each feature by its noise precision (tau), estimated from the PPMI training residuals:
         z = (W_o^T T_o W_o + I)^-1 W_o^T T_o y_o
     This is the posterior mean of Z under MOFA's N(0, I) factor prior.

Variance explained uses the same definition as mofapy2 (BayesNet.calculate_variance_explained),
per view and per group, on the centered and scaled data:
    R2_k = 1 - sum((Y - z_k w_k^T)^2) / sum(Y^2)        (observed entries only, in %)
"""

import os
import argparse
import h5py
import pandas as pd
import numpy as np

VIEWS = ["RNA_seq", "CSF_Proteomics", "Plasma_Proteomics"]


def _decode(values) -> list:
    return [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in values]


def _first_dataset(node):
    """Returns node if it is a dataset, else the first dataset inside the group."""
    if isinstance(node, h5py.Group):
        return node[list(node.keys())[0]]
    return node


def load_model_weights(model_path: str) -> dict:
    """
    Extracts weight expectations W (SW or W) for each view from the trained MOFA HDF5 file.
    Returns dict: {view_name: dataframe (features x factors)}
    """
    weights = {}
    with h5py.File(model_path, 'r') as f:
        exp_group = f['expectations']
        w_key = 'SW' if 'SW' in exp_group else ('W' if 'W' in exp_group else None)

        if w_key is None:
            raise KeyError(f"Could not find weight expectations (W/SW) in {model_path}")

        w_group = exp_group[w_key]

        for view_name in w_group.keys():
            w_matrix = _first_dataset(w_group[view_name])[:]

            if 'features' in f and view_name in f['features']:
                feature_names = _decode(_first_dataset(f['features'][view_name])[:])
            else:
                feature_names = [f"feat_{i}" for i in range(w_matrix.shape[0])]

            # mofapy2 stores W as (factors, features)
            if w_matrix.shape[0] != len(feature_names) and w_matrix.shape[1] == len(feature_names):
                w_matrix = w_matrix.T

            factor_names = [f"Factor{i+1}" for i in range(w_matrix.shape[1])]
            weights[view_name] = pd.DataFrame(w_matrix, index=feature_names, columns=factor_names)

    return weights


def load_model_groups(model_path: str) -> list:
    """Returns the group names the model was trained with (e.g. ['Female', 'Male'])."""
    with h5py.File(model_path, 'r') as f:
        if 'groups' in f:
            g_node = f['groups']['groups'] if isinstance(f['groups'], h5py.Group) else f['groups']
            return _decode(g_node[:])
        return list(f['expectations']['Z'].keys()) if 'Z' in f['expectations'] else []


def estimate_noise_precision(model_path: str, model_weights: dict) -> dict:
    """
    Estimates per-feature noise precision tau = 1 / mean(residual^2) for each view and group
    from the training data and expectations stored in the HDF5 (mofapy2 saves data, W and Z
    by default but not Tau). Returns {view: {group: pd.Series(tau, index=features)}}, or
    an empty dict when the model file holds no training data.
    """
    tau = {}
    with h5py.File(model_path, 'r') as f:
        if 'data' not in f or 'Z' not in f['expectations']:
            return tau
        for view, df_w in model_weights.items():
            if view not in f['data']:
                continue
            tau[view] = {}
            for group in f['data'][view].keys():
                Y = f['data'][view][group][:]                   # Samples x Features
                Z = f['expectations']['Z'][group][:]            # Factors x Samples
                if Z.shape[0] != Y.shape[0]:
                    Z = Z.T
                resid = Y - Z @ df_w.values.T
                mse = np.nanmean(resid ** 2, axis=0)
                mse = np.where(np.isfinite(mse) & (mse > 0), mse, np.nanmedian(mse))
                tau[view][group] = pd.Series(1.0 / mse, index=df_w.index)
    return tau


def center_scale_like_mofa(df_x: pd.DataFrame, sample_groups: pd.Series) -> pd.DataFrame:
    """
    Applies mofapy2's training-time transformation to a features x samples matrix:
    center each feature within each group, then scale the whole view to unit variance.
    """
    out = df_x.copy()
    for group in sample_groups.unique():
        cols = sample_groups.index[sample_groups == group]
        out[cols] = out[cols].sub(out[cols].mean(axis=1, skipna=True), axis=0)
    view_sd = np.nanstd(out.values)
    if view_sd > 0:
        out = out / view_sd
    return out


def project_samples(Y: np.ndarray, W: np.ndarray, tau: np.ndarray) -> np.ndarray:
    """
    Posterior mean of Z for each sample given fixed W, using only observed features.
    Y: Samples x Features (NaN = missing), W: Features x Factors, tau: Samples x Features.
    """
    n_samples, n_factors = Y.shape[0], W.shape[1]
    Z = np.full((n_samples, n_factors), np.nan)
    identity = np.eye(n_factors)
    for i in range(n_samples):
        obs = ~np.isnan(Y[i])
        if not obs.any():
            continue
        W_o = W[obs]
        t_o = tau[i, obs]
        A = W_o.T @ (t_o[:, None] * W_o) + identity
        Z[i] = np.linalg.solve(A, W_o.T @ (t_o * Y[i, obs]))
    return Z


def variance_explained(Y: np.ndarray, Z: np.ndarray, W: np.ndarray):
    """
    mofapy2's variance explained (%) on observed entries of Y (Samples x Features).
    Returns (per_factor array of length K, total).
    """
    obs = ~np.isnan(Y)
    Y0 = np.where(obs, Y, 0.0)
    ss = np.sum(Y0 ** 2)
    if ss == 0:
        return np.zeros(W.shape[1]), 0.0
    per_factor = np.empty(W.shape[1])
    for k in range(W.shape[1]):
        pred = np.where(obs, np.outer(Z[:, k], W[:, k]), 0.0)
        per_factor[k] = (1.0 - np.sum((Y0 - pred) ** 2) / ss) * 100.0
    pred = np.where(obs, Z @ W.T, 0.0)
    total = (1.0 - np.sum((Y0 - pred) ** 2) / ss) * 100.0
    return per_factor, total


def load_sample_groups(metadata_path: str, model_groups: list, group_col: str = "sex") -> pd.Series:
    """Maps participant_id -> model group name using the metadata grouping column."""
    meta = pd.read_csv(metadata_path)
    id_col = "participant_id" if "participant_id" in meta.columns else meta.columns[0]
    meta = meta.drop_duplicates(subset=[id_col]).set_index(id_col)
    if group_col not in meta.columns:
        raise ValueError(f"Metadata has no '{group_col}' column to assign MOFA groups.")
    lookup = {g.lower(): g for g in model_groups}
    raw = meta[group_col].dropna().astype(str)
    groups = raw.str.lower().map(lookup)
    unmatched = sorted(set(raw[groups.isna()]))
    if unmatched:
        print(f"Warning: metadata '{group_col}' values {unmatched} do not match model groups {model_groups}; those samples are skipped.")
    return groups.dropna()


def project_pdbp(model_path: str, input_dir: str, output_dir: str, cohort: str, metadata_path: str = None):
    """
    Loads model weights W, applies MOFA's group centering and view scaling to the cohort
    matrices, projects samples into factor space and computes per-view, per-group R2.
    """
    os.makedirs(output_dir, exist_ok=True)

    # 1. Load trained weights, groups and noise precision from HDF5
    print(f"Loading trained weights from: {model_path}")
    model_weights = load_model_weights(model_path)
    model_groups = load_model_groups(model_path)
    print(f"Model groups: {model_groups}")
    tau = estimate_noise_precision(model_path, model_weights)
    if not tau:
        print("Warning: no training data in model file; using equal noise precision for all features.")

    # 2. Assign samples to model groups
    if metadata_path:
        sample_groups = load_sample_groups(metadata_path, model_groups)
    elif len(model_groups) <= 1:
        sample_groups = None
    else:
        raise ValueError(f"Model has groups {model_groups}; pass --metadata so samples can be centered within each group.")

    # 3. Load, align, center and scale each view
    view_data = {}
    for view in VIEWS:
        file_path = os.path.join(input_dir, f"{cohort.lower()}_{view.lower()}_aligned.csv")
        if not os.path.exists(file_path):
            print(f"Warning: Missing view {view} at {file_path}")
            continue
        if view not in model_weights:
            print(f"Warning: View {view} not found in model weights. Skipping.")
            continue

        df_x = pd.read_csv(file_path, index_col=0)  # Features x Samples
        df_w = model_weights[view]
        df_x = df_x.reindex(index=df_w.index)

        observed = df_x.notna().any(axis=1)
        print(f"View {view}: {observed.sum()} / {len(df_w)} model features observed in {cohort}.")

        if sample_groups is None:
            groups_v = pd.Series(model_groups[0] if model_groups else "group0", index=df_x.columns)
        else:
            keep = df_x.columns.intersection(sample_groups.index)
            dropped = len(df_x.columns) - len(keep)
            if dropped:
                print(f"View {view}: skipping {dropped} samples without a group assignment.")
            df_x = df_x[keep]
            groups_v = sample_groups.loc[keep]

        view_data[view] = (center_scale_like_mofa(df_x, groups_v), groups_v)

    if not view_data:
        raise ValueError("No matching aligned feature data found for projection.")

    # 4. Build one samples x features matrix across views (NaN where a view or feature is missing)
    all_samples = sorted(set().union(*[set(x.columns) for x, _ in view_data.values()]))
    sample_group = pd.concat([g for _, g in view_data.values()])
    sample_group = sample_group[~sample_group.index.duplicated()].reindex(all_samples)

    Y_blocks, W_blocks, tau_blocks = [], [], []
    for view, (df_x, _) in view_data.items():
        df_w = model_weights[view]
        Y_blocks.append(df_x.reindex(columns=all_samples).T.values)
        W_blocks.append(df_w.values)
        if view in tau:
            t = np.vstack([
                tau[view][g].values if g in tau[view] else np.mean(list(tau[view].values()), axis=0)
                for g in sample_group.values
            ])
        else:
            t = np.ones((len(all_samples), len(df_w)))
        tau_blocks.append(t)

    Y_full = np.hstack(Y_blocks)
    W_full = np.vstack(W_blocks)
    tau_full = np.hstack(tau_blocks)
    factor_names = list(next(iter(model_weights.values())).columns)

    Z_proj = project_samples(Y_full, W_full, tau_full)
    df_z = pd.DataFrame(Z_proj, index=all_samples, columns=factor_names)
    df_z = df_z.dropna(how="all")
    df_z["group"] = sample_group.reindex(df_z.index)

    # 5. Variance explained per view and group (same definition as mofapy2)
    r2_records, total_records = [], []
    for view, (df_x, groups_v) in view_data.items():
        W_v = model_weights[view].values
        for group in sorted(groups_v.unique()):
            samples = groups_v.index[groups_v == group].intersection(df_z.index)
            if len(samples) == 0:
                continue
            Y_v = df_x[samples].T.values
            Z_g = df_z.loc[samples, factor_names].values
            per_factor, total = variance_explained(Y_v, Z_g, W_v)
            row = dict(zip(factor_names, per_factor))
            row["group"] = group
            row["view"] = view
            r2_records.append(row)
            total_records.append({"view": view, "group": group, "R2_total": total, "n_samples": len(samples)})

    if r2_records:
        r2_out_path = os.path.join(output_dir, f"{cohort.lower()}_variance_explained_per_factor.csv")
        pd.DataFrame(r2_records).to_csv(r2_out_path, index=False)
        print(f"Saved projected variance explained R2 to: {r2_out_path}")
        tot_out_path = os.path.join(output_dir, f"{cohort.lower()}_variance_explained_total.csv")
        pd.DataFrame(total_records).to_csv(tot_out_path, index=False)
        print(f"Saved projected total variance explained to: {tot_out_path}")

    # 6. Save projected factors Z
    out_path = os.path.join(output_dir, f"{cohort.lower()}_latent_factors_Z.csv")
    df_z.to_csv(out_path)

    print(f"\nSuccessfully projected {df_z.shape[0]} samples into {len(factor_names)} MOFA factors.")
    print(f"Saved projected factors Z to: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Project PDBP data into PPMI MOFA space using HDF5 model.")
    parser.add_argument("--model_path", required=True, help="Path to trained PPMI model.")
    parser.add_argument("--input_dir", required=True, help="Directory with _aligned.csv files.")
    parser.add_argument("--output_dir", required=True, help="Output directory.")
    parser.add_argument("--cohort", default="PDBP", help="Cohort name.")
    parser.add_argument("--metadata", default=None, help="Metadata CSV with participant_id and sex (required for multi-group models).")

    args = parser.parse_args()
    project_pdbp(args.model_path, args.input_dir, args.output_dir, args.cohort, args.metadata)

if __name__ == "__main__":
    main()


###### Example use ##########
#python scripts/project_pdbp.py \
#  --model_path ppmi_mofa_outputs/model_seed_3.hdf5 \
#  --input_dir pdbp_mofa_inputs/ \
#  --output_dir pdbp_mofa_outputs/ \
#  --metadata data/metadata_m0_PDBP.csv \
#  --cohort PDBP
