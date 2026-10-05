#!/usr/bin/env python3
"""
mofa/scripts/preprocess_pdbp.py

Preprocesses PDBP validation omic matrices (RNA-seq and Proteomics), aligns them
to the exact feature set of the trained PPMI MOFA model, and exports both wide
and long-format data structures for downstream analysis/mofapy2.
"""

import os
import sys
import argparse
import h5py
import pandas as pd
import numpy as np
import re


def align_features_to_model(df: pd.DataFrame, model_features: list, view_type: str) -> pd.DataFrame:
    """
    Subsets and reindexes a dataframe (features x samples) to match the exact feature set
    from the trained MOFA model weights by stripping all occurrences of view suffixes.
    """
    # 1. Dynamically strip all repeating occurrences of view_type suffixes (e.g. O00391CSF_ProteomicsCSF_Proteomics -> O00391)
    suffix_pattern = r'(_?' + re.escape(view_type) + r')+$'
    clean_model_features = [re.sub(suffix_pattern, '', f, flags=re.IGNORECASE) for f in model_features]
    
    # 2. Handle Ensembl version stripping and duplicate resolution for RNA-seq
    if view_type == "RNA_seq":
        df.index = df.index.str.split('.').str[0]
        clean_model_features = [f.split('.')[0] for f in clean_model_features]
        if df.index.has_duplicates:
            df = df.groupby(df.index).mean()
    elif df.index.has_duplicates:
        df = df.groupby(df.index).mean()

    # 3. Reindex raw matrix using cleaned feature names
    aligned = df.reindex(index=clean_model_features)
    
    # 4. Restore original model feature names (with view suffixes) for mofapy2 model compatibility
    aligned.index = model_features
    return aligned

def get_model_features(model_path: str, view_name: str) -> list:
    """Extracts feature names for a given view from a trained MOFA HDF5 model."""
    with h5py.File(model_path, 'r') as f:
        if 'features' in f and view_name in f['features']:
            feat_data = f['features'][view_name]
            if isinstance(feat_data, h5py.Dataset):
                return [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in feat_data[:]]
            elif isinstance(feat_data, h5py.Group):
                for k in feat_data.keys():
                    feat_ds = feat_data[k]
                    if isinstance(feat_ds, h5py.Dataset):
                        return [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in feat_ds[:]]
        raise KeyError(f"View '{view_name}' features not found in model {model_path}")

def preprocess_pdbp_view(raw_df: pd.DataFrame, view_type: str, model_features: list) -> pd.DataFrame:
    """Applies normalization/transformation and aligns features to model."""
    if view_type == "RNA_seq":
        # Check if already log-transformed (e.g. from adjust_batch_effects.R)
        flat_vals = raw_df.values.flatten()
        flat_vals = flat_vals[~np.isnan(flat_vals)]
        max_val = np.max(flat_vals) if len(flat_vals) > 0 else 0
        
        if max_val > 50:
            print("PDBP RNA_seq: Raw count range detected. Applying log1p transformation...")
            norm_df = np.log1p(raw_df)
        else:
            print("PDBP RNA_seq: Log-transformed/normalized range detected. Skipping log1p transformation...")
            norm_df = raw_df.copy()
            
        norm_df = norm_df.sub(norm_df.mean(axis=1, skipna=True), axis=0)
    elif "Proteomics" in view_type:
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
        
    aligned_df = align_features_to_model(norm_df, model_features, view_type)
    return aligned_df

def create_long_format(aligned_views: dict, metadata_df: pd.DataFrame, output_dir: str):
    """Converts aligned wide matrices to long-format for mofapy2."""
    mofa_long_list = []
    for view_name, df in aligned_views.items():
        # Melt wide matrix (features x samples) into long format (feature, sample, value)
        long_df = df.reset_index().melt(id_vars="index", var_name="sample", value_name="value")
        long_df.rename(columns={"index": "feature"}, inplace=True)
        long_df["view"] = view_name
        
        # Map sample to group (e.g., sex or single group default)
        if metadata_df is not None and "sex" in metadata_df.columns:
            long_df["group"] = long_df["sample"].map(metadata_df["sex"]).fillna("single_group")
        else:
            long_df["group"] = "single_group"
            
        mofa_long_list.append(long_df)
        
    final_mofa_df = pd.concat(mofa_long_list, ignore_index=True)
    out_file = os.path.join(output_dir, "pdbp_mofa_input_long.csv")
    final_mofa_df.to_csv(out_file, index=False)
    print(f"Saved long-format matrix to {out_file} (Shape: {final_mofa_df.shape})")

def main():
    parser = argparse.ArgumentParser(description="Preprocess PDBP validation data aligned to trained MOFA model.")
    parser.add_argument("--rna_seq", default=None, help="Path to RNA-seq CSV file.")
    parser.add_argument("--csf_proteomics", default=None, help="Path to CSF Proteomics CSV file.")
    parser.add_argument("--plasma_proteomics", default=None, help="Path to Plasma Proteomics CSV file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save preprocessed matrices.")
    parser.add_argument("--cohort", default="PDBP", help="Cohort name.")
    parser.add_argument("--model_path", required=True, help="Path to trained MOFA model HDF5.")
    parser.add_argument("--metadata_path", default=None, help="Path to master metadata CSV (optional).")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Starting PDBP preprocessing for cohort: {args.cohort}")

    # Load metadata if provided
    metadata_df = None
    if args.metadata_path and os.path.exists(args.metadata_path):
        print(f"Loading metadata from {args.metadata_path}...")
        metadata_df = pd.read_csv(args.metadata_path, index_col=0)

    view_paths = {
        "RNA_seq": args.rna_seq,
        "CSF_Proteomics": args.csf_proteomics,
        "Plasma_Proteomics": args.plasma_proteomics
    }

    aligned_views = {}

    for view, input_file in view_paths.items():
        if not input_file:
            print(f"Skipping {view}: No input file path provided.")
            continue

        if not os.path.exists(input_file):
            print(f"Skipping {view}: File not found at '{input_file}'")
            continue

        print(f"Processing {view} from {input_file}...")
        raw_df = pd.read_csv(input_file, index_col=0)

        try:
            model_features = get_model_features(args.model_path, view)
        except KeyError as e:
            print(f"Error fetching features for {view}: {e}")
            continue

        aligned_df = preprocess_pdbp_view(raw_df, view, model_features)

        # Save aligned wide matrix
        output_file = os.path.join(args.output_dir, f"{args.cohort.lower()}_{view.lower()}_aligned.csv")
        aligned_df.to_csv(output_file)
        print(f"Saved aligned {view} matrix to {output_file} (Shape: {aligned_df.shape})")

        # Store in dict for long format generation
        aligned_views[view] = aligned_df

    # Generate long-format output if any views were processed
    if aligned_views:
        create_long_format(aligned_views, metadata_df, args.output_dir)

    print("\nPreprocessing complete.")

if __name__ == "__main__":
    main()