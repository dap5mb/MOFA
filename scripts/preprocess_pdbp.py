#!/usr/bin/env python3
"""
mofa/scripts/preprocess_pdbp.py

Preprocesses PDBP validation omic matrices (RNA-seq and Proteomics) and aligns them
to the exact feature set of the trained PPMI MOFA model.
"""

import os
import sys
import h5py
import pandas as pd
import numpy as np

def align_features_to_model(df: pd.DataFrame, model_features: list) -> pd.DataFrame:
    """
    Subsets and reindexes a dataframe (features x samples) to match the exact feature set
    from the trained MOFA model weights. Missing features are filled with NaN.
    """
    # Reindex to model features, filling missing features with NaN
    aligned = df.reindex(index=model_features)
    return aligned

def get_model_features(model_path: str, view_name: str) -> list:
    """
    Extracts feature names for a given view from a trained MOFA HDF5 model.
    """
    with h5py.File(model_path, 'r') as f:
        # mofapy2 hdf5 structure for features: /features/<view_name>/features or similar
        # Let's inspect potential keys
        if 'features' in f:
            if view_name in f['features']:
                feat_data = f['features'][view_name]
                if isinstance(feat_data, h5py.Dataset):
                    features = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in feat_data[:]]
                    return features
                elif isinstance(feat_data, h5py.Group):
                    # dataset inside group
                    for k in feat_data.keys():
                        feat_ds = feat_data[k]
                        if isinstance(feat_ds, h5py.Dataset):
                            features = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in feat_ds[:]]
                            return features
        raise KeyError(f"View '{view_name}' features not found in model {model_path}")

def preprocess_pdbp_view(raw_df: pd.DataFrame, view_type: str, model_features: list) -> pd.DataFrame:
    """
    Applies normalization/transformation and aligns features to model.
    """
    if view_type == "RNA_seq":
        # Apply log1p normalization if raw counts
        # (Assuming raw or logCPM input, applying log1p as standard stable fallback)
        norm_df = np.log1p(raw_df)
        # Center features
        norm_df = norm_df.sub(norm_df.mean(axis=1, skipna=True), axis=0)
    elif "Proteomics" in view_type:
        # Apply log2 + median centering if not already log-transformed
        if raw_df.max().max() > 50:
            norm_df = np.log2(raw_df + 1)
        else:
            norm_df = raw_df.copy()
        sample_medians = norm_df.median(axis=0, skipna=True)
        global_median = sample_medians.median(skipna=True)
        norm_df = norm_df.sub(sample_medians - global_median, axis=1)
        norm_df = norm_df.sub(norm_df.mean(axis=1, skipna=True), axis=0)
    else:
        norm_df = raw_df.copy()
        
    # Align features to trained model
    aligned_df = align_features_to_model(norm_df, model_features)
    return aligned_df

def main():
    print("PDBP Preprocessing and Feature Alignment script ready.")

if __name__ == "__main__":
    main()
