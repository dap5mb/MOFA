#!/usr/bin/env python3
"""
mofa/scripts/evaluate_factors.py

Reusable evaluation script for MOFA latent factors (Z) across discovery (PPMI) 
and validation (PDBP) cohorts. Evaluates statistical associations with clinical 
and technical covariates (Pearson, Spearman, ANOVA) and compares variance 
explained replication ($R^2$).
"""

import os
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

def dedupe_metadata(meta_df: pd.DataFrame, id_col: str = "participant_id") -> pd.DataFrame:
    """
    Keeps one row per participant (the first). Clinical joins can repeat a participant, for example
    several baseline exam records, which inflates n and significance downstream. Columns whose values
    differ between a participant's rows are reported so the first-row choice can be reviewed.
    """
    dup_mask = meta_df[id_col].duplicated(keep=False)
    if not dup_mask.any():
        return meta_df
    dups = meta_df[dup_mask]
    n_people = dups[id_col].nunique()
    print(f"Metadata: {len(dups)} rows for {n_people} participants with more than one row; keeping the first row per participant.")
    conflicts = {}
    for col in meta_df.columns:
        if col == id_col:
            continue
        n_conflict = (dups.groupby(id_col)[col].nunique(dropna=True) > 1).sum()
        if n_conflict:
            conflicts[col] = int(n_conflict)
    if conflicts:
        print("  Columns whose values differ between a participant's rows (participants affected): " +
              ", ".join(f"{c} ({n})" for c, n in sorted(conflicts.items(), key=lambda kv: -kv[1])))
    return meta_df.drop_duplicates(subset=[id_col], keep="first")

def load_and_merge_factors_metadata(z_path: str, metadata_path: str) -> pd.DataFrame:
    """Loads latent factor Z CSV and merges with sample clinical metadata with ID namespace diagnostics."""
    z_df = pd.read_csv(z_path)
    meta_df = pd.read_csv(metadata_path)
    
    # Identify sample / participant ID column in Z
    id_col = None
    for col in ["participant_id", "patient_id", "sample", "Unnamed: 0"]:
        if col in z_df.columns:
            id_col = col
            break
    if id_col is None:
        id_col = z_df.columns[0]
        
    meta_id_col = "participant_id" if "participant_id" in meta_df.columns else meta_df.columns[0]
    
    # Rename id col in z_df to participant_id if needed
    if id_col != "participant_id":
        z_df = z_df.rename(columns={id_col: "participant_id"})
    if meta_id_col != "participant_id":
        meta_df = meta_df.rename(columns={meta_id_col: "participant_id"})
    meta_df = dedupe_metadata(meta_df)
    z_df = dedupe_metadata(z_df)

    z_ids = set(z_df["participant_id"].astype(str))
    meta_ids = set(meta_df["participant_id"].astype(str))
    overlap = z_ids.intersection(meta_ids)
    
    print(f"ID Namespace Check: {len(z_ids)} IDs in Z matrix, {len(meta_ids)} IDs in metadata. Overlap: {len(overlap)}")
    if len(overlap) == 0:
        print(f"WARNING: Zero ID overlap between Z matrix (e.g., {list(z_ids)[:3]}) and metadata (e.g., {list(meta_ids)[:3]}). Please verify participant ID formatting/namespaces.")

    # Merge
    merged = pd.merge(z_df, meta_df, on="participant_id", how="inner", suffixes=("", "_meta"))
    print(f"Merged latent factors with metadata: {len(merged)} samples matched.")
    return merged

def evaluate_factor_associations(df: pd.DataFrame, output_dir: str, covariates: list):
    """
    Correlates MOFA latent factors (Z) with clinical and technical metadata 
    (Pearson correlation for numeric, ANOVA for categorical).
    """
    os.makedirs(output_dir, exist_ok=True)
    factor_cols = [col for col in df.columns if col.startswith("Factor")]
    print(f"Found latent factors: {factor_cols}")
    
    available_covs = [c for c in covariates if c in df.columns]
    missing_covs = [c for c in covariates if c not in df.columns]
    if missing_covs:
        print(f"Warning: The following requested covariates were not found: {missing_covs}")
    
    if not available_covs:
        raise ValueError(f"None of the requested covariates {covariates} were found in dataframe columns.")    
    
    print(f"Evaluating associations against: {available_covs}")
    results = []
    
    for factor in factor_cols:
        for cov in available_covs:
            sub = df[[factor, cov]].dropna()
            if len(sub) < 10:
                continue
                
            y = sub[factor]
            x = sub[cov]
            
            if pd.api.types.is_numeric_dtype(x):
                corr, pval = stats.pearsonr(x, y)
                test_type = "Pearson"
                stat_val = corr
            else:
                categories = x.unique()
                groups = [y[x == cat] for cat in categories if len(y[x == cat]) > 1]
                if len(groups) > 1:
                    f_val, pval = stats.f_oneway(*groups)
                    test_type = "ANOVA"
                    stat_val = f_val
                else:
                    continue
                    
            results.append({
                "Factor": factor,
                "Metadata": cov,
                "Test": test_type,
                "Statistic": stat_val,
                "P-value": pval
            })
            
    res_df = pd.DataFrame(results)
    out_csv = os.path.join(output_dir, "factor_metadata_associations.csv")
    res_df.to_csv(out_csv, index=False)
    print(f"Saved association statistics to: {out_csv}")
    
    # Generate Heatmap of -log10(p-value)
    if not res_df.empty and HAS_PLOT:
        pivot_p = res_df.pivot(index="Factor", columns="Metadata", values="P-value")
        log_p = -np.log10(pivot_p.clip(lower=1e-300))
        
        plt.figure(figsize=(10, max(4, len(factor_cols) * 0.5)))
        sns.heatmap(log_p, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={'label': '-log10(p-value)'})
        plt.title("Latent Factor vs. Metadata Associations (-log10 p-value)", fontsize=14, fontweight="bold")
        plt.tight_layout()
        
        out_plot = os.path.join(output_dir, "factor_metadata_associations_heatmap.png")
        plt.savefig(out_plot, dpi=300)
        plt.close()
        print(f"Saved association heatmap to: {out_plot}")

def compare_variance_replication(ppmi_r2_path: str, pdbp_r2_path: str, output_dir: str):
    """
    Compares per-factor variance explained ($R^2$) between PPMI and PDBP for each view and
    group (sex). R2 values from different views or groups have different denominators, so
    they are compared side by side and never summed.
    """
    if not os.path.exists(ppmi_r2_path) or not os.path.exists(pdbp_r2_path):
        print("Skipping R2 replication comparison (missing PPMI or PDBP R2 file).")
        return

    ppmi_r2 = pd.read_csv(ppmi_r2_path)
    pdbp_r2 = pd.read_csv(pdbp_r2_path)

    print("Comparing variance explained ($R^2$) replication across cohorts...")

    for df in (ppmi_r2, pdbp_r2):
        for col in ("view", "group"):
            if col not in df.columns:
                raise ValueError(f"R2 file is missing a '{col}' column; regenerate it with extract_factors.py / project_pdbp.py.")
            df[col] = df[col].astype(str)

    factor_cols = [c for c in ppmi_r2.columns if c.startswith("Factor")]
    keys = ["view", "group"]
    ppmi_long = ppmi_r2.melt(id_vars=keys, value_vars=factor_cols, var_name="Factor", value_name="PPMI_R2")
    pdbp_long = pdbp_r2.melt(id_vars=keys, value_vars=factor_cols, var_name="Factor", value_name="PDBP_R2")

    # Join on lower-cased view and group names
    for df in (ppmi_long, pdbp_long):
        df["view_key"] = df["view"].str.lower()
        df["group_key"] = df["group"].str.lower()
    comp = pd.merge(ppmi_long, pdbp_long.drop(columns=keys), on=["view_key", "group_key", "Factor"], how="inner")

    if comp.empty:
        raise ValueError(
            f"No overlapping (view, group) pairs between PPMI {sorted(set(zip(ppmi_r2.view, ppmi_r2.group)))} "
            f"and PDBP {sorted(set(zip(pdbp_r2.view, pdbp_r2.group)))}."
        )
    comp = comp.drop(columns=["view_key", "group_key"])

    out_csv = os.path.join(output_dir, "cohort_variance_replication_comparison.csv")
    comp.to_csv(out_csv, index=False)
    print(f"Saved replication comparison table to: {out_csv}")

    # Agreement of the per-factor R2 profile within each view and group
    for (view, group), sub in comp.groupby(keys):
        if len(sub) > 2:
            rho, _ = stats.spearmanr(sub["PPMI_R2"], sub["PDBP_R2"])
            print(f"  {view} / {group}: Spearman rho of per-factor R2 = {rho:.2f}")

    if HAS_PLOT:
        views = sorted(comp["view"].unique())
        groups = sorted(comp["group"].unique())
        fig, axes = plt.subplots(len(views), len(groups), figsize=(5 * len(groups), 3.5 * len(views)),
                                 squeeze=False, sharex=True)
        for i, view in enumerate(views):
            for j, group in enumerate(groups):
                ax = axes[i][j]
                sub = comp[(comp["view"] == view) & (comp["group"] == group)].set_index("Factor").reindex(factor_cols)
                sub[["PPMI_R2", "PDBP_R2"]].plot(kind="bar", ax=ax, colormap="tab10", width=0.7, legend=(i == 0 and j == 0))
                ax.set_title(f"{view} / {group}", fontsize=11)
                ax.set_xlabel("")
                ax.set_ylabel("Variance Explained (%)" if j == 0 else "")
                ax.axhline(0, color="black", linewidth=0.5)
                ax.tick_params(axis="x", rotation=90)
        fig.suptitle("Variance Explained ($R^2$) Replication per View and Group: PPMI vs. PDBP", fontsize=13, fontweight="bold")
        plt.tight_layout()

        out_path = os.path.join(output_dir, "cohort_variance_replication_comparison.png")
        plt.savefig(out_path, dpi=300)
        plt.close()
        print(f"Saved replication comparison plot to: {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Evaluate MOFA latent factors and metadata associations.")
    parser.add_argument("--input_z", required=True, help="Path to latent factors Z CSV or pre-merged CSV")
    parser.add_argument("--metadata", default=None, help="Path to master metadata CSV (optional if already merged)")
    parser.add_argument("--output", default="../results", help="Output directory")
    parser.add_argument("--covariates", nargs="+", default=["diagnosis", "sex", "age_at_baseline", "mds_updrs_3", "moca"], help="Covariates to evaluate")
    parser.add_argument("--ppmi_r2", default=None, help="Path to PPMI R2 CSV for replication comparison (optional)")
    parser.add_argument("--pdbp_r2", default=None, help="Path to PDBP R2 CSV for replication comparison (optional)")
    
    args = parser.parse_args()
    
    if args.metadata and os.path.exists(args.metadata):
        merged_df = load_and_merge_factors_metadata(args.input_z, args.metadata)
    else:
        merged_df = pd.read_csv(args.input_z)
        
    evaluate_factor_associations(merged_df, args.output, args.covariates)
    
    if args.ppmi_r2 and args.pdbp_r2:
        compare_variance_replication(args.ppmi_r2, args.pdbp_r2, args.output)
        
    print("\n=== Evaluation Complete ===")

if __name__ == "__main__":
    main()
