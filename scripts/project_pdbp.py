#!/usr/bin/env python3
"""
mofa/scripts/project_pdbp.py

Projects preprocessed PDBP validation omic data into the trained PPMI MOFA factor space
using weights directly loaded from the HDF5 model file.
"""

import os
import argparse
import h5py
import pandas as pd
import numpy as np

def load_model_weights(model_path: str) -> dict:
    """
    Extracts weight expectations W (SW or W) for each view from the trained MOFA HDF5 file.
    Returns dict: {view_name: dataframe (features x factors)}
    """
    weights = {}
    with h5py.File(model_path, 'r') as f:
        # Check where expectations/W or expectations/SW resides
        exp_group = f['expectations']
        w_key = 'SW' if 'SW' in exp_group else ('W' if 'W' in exp_group else None)
        
        if w_key is None:
            raise KeyError(f"Could not find weight expectations (W/SW) in {model_path}")
            
        w_group = exp_group[w_key]
        
        # Load weights and feature names per view
        for view_name in w_group.keys():
            w_matrix = w_group[view_name][:]  # Shape: (features, factors) or (factors, features)
            
            # Extract feature names for alignment
            if 'features' in f and view_name in f['features']:
                feat_ds = f['features'][view_name]
                if isinstance(feat_ds, h5py.Group):
                    feat_ds = feat_ds[list(feat_ds.keys())[0]]
                feature_names = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in feat_ds[:]]
            else:
                feature_names = [f"feat_{i}" for i in range(w_matrix.shape[0])]
            
            # Ensure shape is (features, factors)
            if w_matrix.shape[0] != len(feature_names) and w_matrix.shape[1] == len(feature_names):
                w_matrix = w_matrix.T
                
            num_factors = w_matrix.shape[1]
            factor_names = [f"Factor{i+1}" for i in range(num_factors)]
            
            weights[view_name] = pd.DataFrame(w_matrix, index=feature_names, columns=factor_names)
            
    return weights

def project_pdbp(model_path: str, input_dir: str, output_dir: str, cohort: str):
    """
    Loads model weights W and aligns PDBP matrices to compute projected latent factors Z.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Load trained weight matrices W from HDF5
    print(f"Loading trained weights from: {model_path}")
    model_weights = load_model_weights(model_path)
    
    views = ["RNA_seq", "CSF_Proteomics", "Plasma_Proteomics"]
    concat_X = []
    concat_W = []
    
    # 2. Align views and aggregate across available omics
    for view in views:
        file_path = os.path.join(input_dir, f"{cohort.lower()}_{view.lower()}_aligned.csv")
        
        if not os.path.exists(file_path):
            print(f"Warning: Missing view {view} at {file_path}")
            continue
            
        df_x = pd.read_csv(file_path, index_col=0) # Shape: (features, samples)
        
        if view not in model_weights:
            print(f"Warning: View {view} not found in model weights. Skipping.")
            continue
            
        df_w = model_weights[view] # Shape: (features, factors)
        
        # Align features
        common_features = df_x.index.intersection(df_w.index)
        print(f"View {view}: Matched {len(common_features)} / {len(df_w)} features.")
        
        df_x_aligned = df_x.loc[common_features]
        df_w_aligned = df_w.loc[common_features]
        
        concat_X.append(df_x_aligned)
        concat_W.append(df_w_aligned)
        
    if not concat_X:
        raise ValueError("No matching aligned feature data found for projection.")
        
    # Concatenate features across all available views
    X_full = pd.concat(concat_X, axis=0) # Features x Samples
    W_full = pd.concat(concat_W, axis=0) # Features x Factors
    
    # Fill missing values with 0 (mean-centered)
    X_mat = X_full.fillna(0).values.T # Shape: (Samples, Features)
    W_mat = W_full.values              # Shape: (Features, Factors)
    
    # 3. Calculate latent factor projection: Z = X * W * (W^T * W)^(-1)
    # Ridge regularization term to ensure stable inversion
    lambda_reg = 1e-3
    wtw = np.dot(W_mat.T, W_mat) + lambda_reg * np.eye(W_mat.shape[1])
    inv_wtw = np.linalg.inv(wtw)
    
    # Z_proj shape: (Samples, Factors)
    Z_proj = np.dot(X_mat, np.dot(W_mat, inv_wtw))
    
    # 4. Save projected factors Z
    df_z = pd.DataFrame(Z_proj, index=X_full.columns, columns=W_full.columns)
    out_path = os.path.join(output_dir, f"{cohort.lower()}_latent_factors_Z.csv")
    df_z.to_csv(out_path)
    
    print(f"\nSuccessfully projected {df_z.shape[0]} samples into {df_z.shape[1]} MOFA factors.")
    print(f"Saved projected factors Z to: {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Project PDBP data into PPMI MOFA space using HDF5 model.")
    parser.add_argument("--model_path", required=True, help="Path to trained PPMI model.")
    parser.add_argument("--input_dir", required=True, help="Directory with _aligned.csv files.")
    parser.add_argument("--output_dir", required=True, help="Output directory.")
    parser.add_argument("--cohort", default="PDBP", help="Cohort name.")
    
    args = parser.parse_args()
    project_pdbp(args.model_path, args.input_dir, args.output_dir, args.cohort)

if __name__ == "__main__":
    main()
    
    
###### Example use ##########
#python scripts/project_pdbp.py \
#  --model_path ppmi_mofa_outputs/model_seed_1.hdf5 \
#  --input_dir pdbp_mofa_inputs/ \
#  --output_dir pdbp_mofa_outputs/ \
#  --cohort PDBP