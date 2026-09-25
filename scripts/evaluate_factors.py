import os
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

def evaluate_factor_associations(z_meta_path: str, output_dir: str, covariates: list):
    """
    Correlates MOFA latent factors (Z) with clinical and technical metadata 
    to check for biological signal vs. technical batch effects.
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"Loading Z matrix and metadata from: {z_meta_path}")
    df = pd.read_csv(z_meta_path)
    
    # Identify factor columns
    factor_cols = [col for col in df.columns if col.startswith("Factor")]
    print(f"Found latent factors: {factor_cols}")
    
    # Filter covariates to those present in the dataframe
    available_covs = [c for c in covariates if c in df.columns]
    missing_covs = [c for c in covariates if c not in df.columns]
    if missing_covs:
        print(f"Warning: The following requested covariates were not found in metadata: {missing_covs}")
    
    if not available_covs:
        raise ValueError(f"None of the requested covariates {covariates} were found in the dataframe columns: {list(df.columns)}")    
    
    print(f"Evaluating associations against: {available_covs}")
    
    results = []
    
    for factor in factor_cols:
        for cov in available_covs:
            sub = df[[factor, cov]].dropna()
            if len(sub) < 10:
                continue
                
            y = sub[factor]
            x = sub[cov]
            
            # Determine data type of covariate
            if pd.api.types.is_numeric_dtype(x):
                # Pearson correlation for numeric variables (e.g., age)
                corr, pval = stats.pearsonr(x, y)
                test_type = "Pearson"
                stat_val = corr
            else:
                # ANOVA (F-test) for categorical variables (e.g., sex, diagnosis, plate)
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
    
    # Save association table
    out_csv = os.path.join(output_dir, "factor_metadata_associations.csv")
    res_df.to_csv(out_csv, index=False)
    print(f"Saved association statistics to: {out_csv}")
    
    # Generate Heatmap of -log10(p-value)
    if not res_df.empty:
        pivot_p = res_df.pivot(index="Factor", columns="Metadata", values="P-value")
        # Convert p-values to -log10(p) for visualization
        log_p = -np.log10(pivot_p.clip(lower=1e-300))
        
        plt.figure(figsize=(10, len(factor_cols) * 0.5 + 2))
        sns.heatmap(log_p, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={'label': '-log10(p-value)'})
        plt.title("Latent Factor vs. Metadata Associations (-log10 p-value)")
        plt.tight_layout()
        
        out_plot = os.path.join(output_dir, "factor_metadata_associations_heatmap.png")
        plt.savefig(out_plot, dpi=300)
        plt.close()
        print(f"Saved association heatmap to: {out_plot}")
        
    print("\n=== Evaluation Complete ===")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate MOFA factor correlations with metadata.")
    parser.add_argument("--input", required=True, help="Path to latent_factors_Z_with_meta.csv")
    parser.add_argument("--output", default="../ppmi_mofa_outputs", help="Output directory")
    parser.add_argument("--covariates", nargs="+", default=["diagnosis", "sex", "age", "plate", "batch"], help="List of metadata columns to test")
    
    args = parser.parse_args()
    evaluate_factor_associations(args.input, args.output, args.covariates)
