# Use:
#   python3 extract_factors.py \  
#    --model ../ppmi_mofa_outputs/model_seed_0.hdf5 \                                                                                 --metadata ../ppmi_mofa_inputs/ppmi_master_metadata_mofa_ready.csv \                                                             --output ../ppmi_mofa_outputs \                                                                                                   --top_pct 0.05   

#######################
import os
import argparse
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_variance_explained(r2_df: pd.DataFrame, output_dir: str):
    """
    Generates summary plots (stacked bar plot and heatmap) of variance 
    explained per factor across omic views.
    """
    try:
        # Prepare pivot table for plotting (Views on columns, Factors on index)
        # Average across groups if multiple groups exist
        r2_pivot = r2_df.pivot_table(index="view", columns=None, values=[c for c in r2_df.columns if c.startswith("Factor")], aggfunc="mean")
        r2_pivot = r2_pivot.T  # Index: Factors, Columns: Views
        
        # 1. Stacked Bar Chart
        fig, ax = plt.subplots(figsize=(10, 6))
        r2_pivot.plot(kind="bar", stacked=True, ax=ax, colormap="viridis")
        
        ax.set_title("Variance Explained per Factor by Omic View", fontsize=14, fontweight="bold")
        ax.set_xlabel("Latent Factors", fontsize=12)
        ax.set_ylabel("Variance Explained (%)", fontsize=12)
        ax.legend(title="Views", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=True)
        plt.xticks(rotation=0)
        plt.tight_layout()
        
        bar_plot_path = os.path.join(output_dir, "variance_explained_bar_plot.png")
        plt.savefig(bar_plot_path, dpi=300)
        plt.close()
        print(f"Saved variance explained bar plot to: {bar_plot_path}")

        # 2. Heatmap
        plt.figure(figsize=(10, 6))
        sns.heatmap(r2_pivot.T, annot=True, fmt=".2f", cmap="Blues", cbar_kws={'label': 'Variance Explained (%)'})
        plt.title("Variance Explained Heatmap (Views x Factors)", fontsize=14, fontweight="bold")
        plt.xlabel("Latent Factors", fontsize=12)
        plt.ylabel("Omic View", fontsize=12)
        plt.tight_layout()
        
        heatmap_path = os.path.join(output_dir, "variance_explained_heatmap.png")
        plt.savefig(heatmap_path, dpi=300)
        plt.close()
        print(f"Saved variance explained heatmap to: {heatmap_path}")

    except Exception as e:
        print(f"Warning: Could not auto-generate variance plot due to data shape error: {e}")

def extract_mofa_results(model_path: str, metadata_path: str, output_dir: str, top_pct: float = 0.05):
    """
    Loads a trained MOFA HDF5 model directly and extracts variance explained (R2), 
    latent factor coordinates (Z), and feature weights (W).
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"Loading trained model from: {model_path}")
    
    with h5py.File(model_path, 'r') as f:
        
        # mofapy2 stores view names as a dataset at views/views, in training order
        # (which is also the row order of variance_explained/r2_per_factor)
        if 'views' in f:
            v_node = f['views']['views'] if isinstance(f['views'], h5py.Group) else f['views']
            view_names = [v.decode('utf-8') if isinstance(v, bytes) else str(v) for v in v_node[:]]
        else:
            view_names = list(f['features'].keys())
        print(f"Views: {view_names}")

        print("\n=== 1. Extracting Variance Explained (R2) ===")
        if 'variance_explained' in f:
            v_grp = f['variance_explained']
            if 'r2_per_factor' in v_grp:
                r2_data = v_grp['r2_per_factor']
                r2_dfs = []
                
                if isinstance(r2_data, h5py.Group):
                    for group_name in r2_data.keys():
                        g_r2 = r2_data[group_name][:] # Shape: (views, factors)
                        factor_names = [f"Factor{i+1}" for i in range(g_r2.shape[1])]
                        
                        if len(view_names) != g_r2.shape[0]:
                            raise ValueError(f"Found {len(view_names)} view names but R2 has {g_r2.shape[0]} view rows")

                        df_r2 = pd.DataFrame(g_r2, index=view_names, columns=factor_names)
                        df_r2["group"] = group_name
                        df_r2["view"] = df_r2.index
                        r2_dfs.append(df_r2)
                    
                    if r2_dfs:
                        full_r2_df = pd.concat(r2_dfs, axis=0)
                        r2_out_path = os.path.join(output_dir, "variance_explained_per_factor.csv")
                        full_r2_df.to_csv(r2_out_path, index=False)
                        print(f"Saved variance explained per factor to: {r2_out_path}")
                        
                        # Generate visualization plots
                        plot_variance_explained(full_r2_df, output_dir)

            if 'r2_total' in v_grp:
                total_records = []
                for group_name in v_grp['r2_total'].keys():
                    g_tot = v_grp['r2_total'][group_name][:]  # Shape: (views,)
                    for view_name, val in zip(view_names, g_tot):
                        total_records.append({"view": view_name, "group": group_name, "R2_total": val})
                tot_out_path = os.path.join(output_dir, "variance_explained_total.csv")
                pd.DataFrame(total_records).to_csv(tot_out_path, index=False)
                print(f"Saved total variance explained per view and group to: {tot_out_path}")

        print("\n=== 2. Extracting Latent Factors (Z Matrix) ===")
        z_dfs = []
        z_grp = f['expectations']['Z']
        groups = list(z_grp.keys())
        
        for group in groups:
            raw_z = z_grp[group][:] 
            sample_names = [s.decode('utf-8') if isinstance(s, bytes) else str(s) for s in f['samples'][group][:]]
            
            # Transpose if shape is (factors, samples)
            if raw_z.shape[0] != len(sample_names) and raw_z.shape[1] == len(sample_names):
                group_z = raw_z.T
            else:
                group_z = raw_z
                
            factor_names = [f"Factor{i+1}" for i in range(group_z.shape[1])]
            
            df_z = pd.DataFrame(group_z, index=sample_names, columns=factor_names)
            df_z["group"] = group
            df_z["participant_id"] = df_z.index
            z_dfs.append(df_z)
            
        combined_z = pd.concat(z_dfs, axis=0)
        
        if os.path.exists(metadata_path):
            print(f"Merging latent factors with metadata from: {metadata_path}")
            meta = pd.read_csv(metadata_path)
            if "participant_id" in meta.columns:
                meta = meta.set_index("participant_id")
            combined_z = combined_z.set_index("participant_id").join(meta, how="left").reset_index()
        
        z_out_path = os.path.join(output_dir, "latent_factors_Z_with_meta.csv")
        combined_z.to_csv(z_out_path, index=False)
        print(f"Saved latent factors Z to: {z_out_path}")
        
        print("\n=== 3. Extracting Feature Weights (W Matrix) ===")
        w_grp = f['expectations']['W']
        views = list(w_grp.keys())
        
        weight_dfs = []
        top_weight_dfs = []
        
        for view_name in views:
            view_w_grp = w_grp[view_name]
            feature_names = [feat.decode('utf-8') if isinstance(feat, bytes) else str(feat) for feat in f['features'][view_name][:]]
            
            if isinstance(view_w_grp, h5py.Group):
                for group in view_w_grp.keys():
                    raw_w = view_w_grp[group][:]
                    
                    g_weights = raw_w.T if raw_w.shape[0] != len(feature_names) and raw_w.shape[1] == len(feature_names) else raw_w
                    factor_names = [f"Factor{i+1}" for i in range(g_weights.shape[1])]
                    
                    df_w = pd.DataFrame(g_weights, index=feature_names, columns=factor_names)
                    df_w["view"] = view_name
                    df_w["group"] = group
                    df_w["feature"] = df_w.index
                    weight_dfs.append(df_w)
                    
                    for col in factor_names:
                        top_n = max(1, int(len(df_w) * top_pct))
                        top_features = df_w.reindex(df_w[col].abs().sort_values(ascending=False).index).head(top_n)
                        top_weight_dfs.append(top_features)
            else:
                raw_w = view_w_grp[:]
                g_weights = raw_w.T if raw_w.shape[0] != len(feature_names) and raw_w.shape[1] == len(feature_names) else raw_w
                factor_names = [f"Factor{i+1}" for i in range(g_weights.shape[1])]
                
                df_w = pd.DataFrame(g_weights, index=feature_names, columns=factor_names)
                df_w["view"] = view_name
                df_w["feature"] = df_w.index
                weight_dfs.append(df_w)
                
                for col in factor_names:
                    top_n = max(1, int(len(df_w) * top_pct))
                    top_features = df_w.reindex(df_w[col].abs().sort_values(ascending=False).index).head(top_n)
                    top_weight_dfs.append(top_features)
                
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