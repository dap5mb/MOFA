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
    """Compares view-wise variance explained ($R^2$) between PPMI and PDBP."""
    if not os.path.exists(ppmi_r2_path) or not os.path.exists(pdbp_r2_path):
        print("Skipping R2 replication comparison (missing PPMI or PDBP R2 file).")
        return
        
    ppmi_r2 = pd.read_csv(ppmi_r2_path)
    pdbp_r2 = pd.read_csv(pdbp_r2_path)
    
    print("Comparing variance explained ($R^2$) replication across cohorts...")
    if HAS_PLOT:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Aggregate total R2 per view for comparison if view column exists
        if "view" in ppmi_r2.columns and "view" in pdbp_r2.columns:
            factor_cols_ppmi = [c for c in ppmi_r2.columns if c.startswith("Factor")]
            factor_cols_pdbp = [c for c in pdbp_r2.columns if c.startswith("Factor")]
            
            ppmi_summary = ppmi_r2.groupby("view")[factor_cols_ppmi].sum().sum(axis=1).reset_index(name="PPMI_R2")
            pdbp_summary = pdbp_r2.groupby("view")[factor_cols_pdbp].sum().sum(axis=1).reset_index(name="PDBP_R2")
            
            comp_df = pd.merge(ppmi_summary, pdbp_summary, on="view", how="inner")
            if not comp_df.empty:
                comp_df.set_index("view")[["PPMI_R2", "PDBP_R2"]].plot(kind="bar", ax=ax, colormap="tab10")
                ax.set_title("Variance Explained ($R^2$) Replication: PPMI vs. PDBP", fontsize=14, fontweight="bold")
                ax.set_xlabel("Omic View", fontsize=12)
                ax.set_ylabel("Total Variance Explained (%)", fontsize=12)
                plt.xticks(rotation=0)
                plt.legend(title="Cohort", frameon=True)
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
