import os
import argparse
import pandas as pd
import numpy as np
from mofapy2.run.load_model import load_model

def extract_mofa_results(model_path: str, metadata_path: str, output_dir: str, top_pct: float = 0.05):
    """
    Loads a trained MOFA model and extracts variance explained (R2), 
    latent factor coordinates (Z), and feature weights (W).
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"Loading trained model from: {model_path}")
    model = load_model(model_path)
    
    print("\n=== 1. Extracting Variance Explained (R2) == @==")
    # mofapy2 provides get_variance_explained()
    r2_dict = model.get_variance_explained()
    
    # r2_dict typically contains 'r2_per_factor' and 'r2_total'
    if "r2_per_factor" in r2_dict:
        r2_df = r2_dict["r2_per_factor"]
        # Save variance explained summary
        r2_out_path = os.path.join(output_dir, "variance_explained_per_factor.csv")
        r2_df.to_csv(r2_out_path)
        print(f"Saved variance explained per factor to: {r2_out_path}")
    
    print("\n=== 2. Extracting Latent Factors (Z Matrix) ===")
    # Get expectations for Z (latent factors) across groups
    # model.get_expectations('Z') returns a dict/list keyed by group
    Z_dict = model.get_expectations('Z')
    groups = list(Z_dict.keys())
    
    z_dfs = []
    for group in groups:
        # Z_dict[group] is shape (samples x factors)
        group_z = Z_dict[group]
        # Get sample names for this group
        sample_names = model.get_samples_names()[groups.index(group)]
        factor_names = [f"Factor{i+1}" for i in range(group_z.shape[1])]
        
        df_z = pd.DataFrame(group_z, index=sample_names, columns=factor_names)
        df_z["group"] = group
        df_z["participant_id"] = df_z.index
        z_dfs.append(df_z)
        
    combined_z = pd.concat(z_dfs, axis=0)
    
    # Merge with master metadata if provided
    if os.path.exists(metadata_path):
        print(f"Merging latent factors with metadata from: {metadata_path}")
        meta = pd.read_csv(metadata_path)
        if "participant_id" in meta.columns:
            meta = meta.set_index("participant_id")
        
        # Join metadata columns
        combined_z = combined_z.set_index("participant_id").join(meta, how="left").reset_index()
    
    z_out_path = os.path.join(output_dir, "latent_factors_Z_with_meta.csv")
    combined_z.to_csv(z_out_path, index=False)
    print(f"Saved latent factors Z to: {z_out_path}")
    
    print("\n=== 3. Extracting Feature Weights (W Matrix) ===")
    # model.get_expectations('W') returns weights per view
    W_dict = model.get_expectations('W')
    view_names = model.get_views_names()
    
    weight_dfs = []
    top_weight_dfs = []
    
    for view_idx, view_name in enumerate(view_names):
        # W_dict[view_idx] is shape (features x factors x groups) or similar depending on mofapy2 version
        # Let's inspect shape and structure safely
        view_w = W_dict[view_name] # view_w is usually a dict of groups or a numpy array
        
        # If view_w is a dict per group:
        if isinstance(view_w, dict):
            for group, g_weights in view_w.items():
                feature_names = model.get_features_names()[view_idx]
                factor_names = [f"Factor{i+1}" for i in range(g_weights.shape[1])]
                df_w = pd.DataFrame(g_weights, index=feature_names, columns=factor_names)
                df_w["view"] = view_name
                df_w["group"] = group
                df_w["feature"] = df_w.index
                weight_dfs.append(df_w)
                
                # Extract top X% drivers per factor (by absolute weight)
                for col in factor_names:
                    top_n = int(len(df_w) * top_pct)
                    top_features = df_w.reindex(df_w[col].abs().sort_values(ascending=False).index).head(top_n)
                    top_weight_dfs.append(top_features)
        else:
            # Fallback if array structure differs
            print(f"Note: View {view_name} weights have non-dict structure. Skipping standard extraction.")
            
    if weight_dfs:
        full_w_df = pd.concat(weight_dfs, axis=0)
        full_w_path = os.path.join(output_dir, "feature_weights_W_full.csv")
        full_w_df.to_csv(full_w_path, index=False)
        print(f"Saved full feature weights to: {full_w_path}")
        
    if top_weight_dfs:
        top_w_df = pd.concat(top_weight_dfs, axis=0)
        top_w_path = os.path.join(output_dir, f"top_{int(top_pct*100)}pct_feature_weights.csv")
        top_w_df.to_csv(top_w_path, index=False)
        print(f"Saved top {int(top_pct*100)}% feature drivers to: {top_w_path}")
        
    print("\n=== Extraction Complete ===")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract MOFA factors and weights.")
    parser.add_argument("--model", required=True, help="Path to trained HDF5 model file")
    parser.add_argument("--metadata", required=True, help="Path to master metadata CSV")
    parser.add_argument("--output", default="../ppmi_mofa_outputs", help="Output directory")
    parser.add_argument("--top_pct", type=float, default=0.05, help="Top fraction of features to extract (default: 0.05)")
    
    args = parser.parse_args()
    extract_mofa_results(args.model, args.metadata, args.output, args.top_pct)
