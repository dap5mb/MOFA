import os
import numpy as np
import pandas as pd

def run_preflight_check(rna_df, csf_df, plasma_df, metadata_df):
    """
    Verifies that participant IDs in metadata are present and aligned across all omic matrices.
    Returns a filtered metadata dataframe and cleaned views.
    """
    print("--- Running Pre-flight Integrity Check ---")
    
    # 1. Check for missing sex
    meta = metadata_df.copy()
    if 'sex' not in meta.columns:
        raise ValueError("Metadata must contain a 'sex' column.")
    
    initial_count = len(meta)
    meta = meta.dropna(subset=['sex'])
    if len(meta) < initial_count:
        print(f"Dropped {initial_count - len(meta)} samples with missing sex.")

    # 2. Align IDs across all sources
    meta_ids = set(meta.index)
    rna_ids = set(rna_df.columns)
    csf_ids = set(csf_df.columns)
    pla_ids = set(plasma_df.columns)
    
    # Ensure every participant in metadata exists in at least one view
    all_view_ids = rna_ids | csf_ids | pla_ids
    common_ids = meta_ids.intersection(all_view_ids)
    
    print(f"Metadata IDs: {len(meta_ids)}, IDs in views: {len(all_view_ids)}")
    print(f"Final aligned participant set: {len(common_ids)}")
    
    return meta.loc[list(common_ids)]

def validate_and_intersect_views(
    rna_df: pd.DataFrame,
    csf_df: pd.DataFrame,
    plasma_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    output_dir: str = "./mofa_inputs"
):
    """
    Task 1.2.4: Final Matrix Validation & Sample Intersection
    """
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== Step 1: Validating Input Matrices & Metadata ===")
    
    # Pre-flight check and metadata cleanup
    if "participant_id" in metadata_df.columns:
        metadata_df = metadata_df.set_index("participant_id")
    
    metadata_df = run_preflight_check(rna_df, csf_df, plasma_df, metadata_df)
    
    views = {
        "RNA_seq": rna_df,
        "CSF_Proteomics": csf_df,
        "Plasma_Proteomics": plasma_df
    }

    # --------------------------------------------------------------------------
    # 1. Inspect Sample Overlap Across Views
    # --------------------------------------------------------------------------
    rna_samples = set(rna_df.columns)
    csf_samples = set(csf_df.columns)
    plasma_samples = set(plasma_df.columns)
    meta_samples = set(metadata_df.index)
    
    # Find complete case intersection (present in ALL 3 views + metadata)
    complete_cases = rna_samples.intersection(csf_samples).intersection(plasma_samples).intersection(meta_samples)
    
    # Find union of all samples present in metadata and at least ONE omic view
    all_omic_samples = (rna_samples | csf_samples | plasma_samples).intersection(meta_samples)
    
    print(f"Total Unique Participants in Metadata with >= 1 view: {len(all_omic_samples)}")
    print(f"Complete Case Participants (Present in all 3 views): {len(complete_cases)}")
    print(f"  - RNA-seq samples aligned: {len(rna_samples.intersection(meta_samples))}")
    print(f"  - CSF Proteomics samples aligned: {len(csf_samples.intersection(meta_samples))}")
    print(f"  - Plasma Proteomics samples aligned: {len(plasma_samples.intersection(meta_samples))}\n")
    
    # --------------------------------------------------------------------------
    # 2. Re-align and Format Matrices
    # --------------------------------------------------------------------------
    # Master column ordering based on union of participants passing QC
    master_participant_id = sorted(list(all_omic_samples))
    aligned_metadata = metadata_df.loc[master_participant_id].copy()
    
    aligned_views = {}
    for view_name, df in views.items():
        # Reindex columns to master_participant_id (places NaN for missing participant assays)
        reindexed_df = df.reindex(columns=master_participant_id)
        
        # Verify row names (features) are unique strings
        assert reindexed_df.index.is_unique, f"Duplicate feature IDs found in {view_name}!"
        
        # Verify no entire row/feature is completely NaN
        all_na_features = reindexed_df.isna().all(axis=1)
        if all_na_features.any():
            print(f"Warning: Dropping {all_na_features.sum()} features with 100% missing values in {view_name}.")
            reindexed_df = reindexed_df.loc[~all_na_features]
            
        aligned_views[view_name] = reindexed_df
        print(f"Formated {view_name}: {reindexed_df.shape[0]} features x {reindexed_df.shape[1]} participants")

    # --------------------------------------------------------------------------
    # 3. Format into Long DataFrame for pyMOFA / mofapy2
    # --------------------------------------------------------------------------
    # mofapy2 accepts a dictionary of dataframes (features x samples) or long format
    mofa_long_list = []
    
    for view_name, df in aligned_views.items():
        # Melt matrix to long format (sample, feature, value, view, group)
        long_df = df.reset_index().rename(columns={df.index.name or "index": "feature"})
        long_df = long_df.melt(id_vars="feature", var_name="sample", value_name="value")
        long_df["view"] = view_name
        
        # Add primary grouping variable ('sex') for Task 1.3 Multi-Group setup
        long_df["group"] = long_df["sample"].map(aligned_metadata["sex"])
        
        # Drop missing values to keep the input file compact
        long_df = long_df.dropna(subset=["value"])
        mofa_long_list.append(long_df)

    final_mofa_df = pd.concat(mofa_long_list, ignore_index=True)
    
    # --------------------------------------------------------------------------
    # 4. Export Final Matrices & Metadata
    # --------------------------------------------------------------------------
    print("\n=== Step 4: Exporting Final Files ===")
    
    # Export wide matrices (features x samples)
    for view_name, df in aligned_views.items():
        out_path = os.path.join(output_dir, f"ppmi_{view_name.lower()}_mofa_ready.csv")
        df.to_csv(out_path)
        print(f"Saved: {out_path}")
        
    # Export Master Metadata
    meta_out_path = os.path.join(output_dir, "ppmi_master_metadata_mofa_ready.csv")
    aligned_metadata.to_csv(meta_out_path)
    print(f"Saved: {meta_out_path}")
    
    # Export Long DataFrame (Direct mofapy2 entry point)
    long_out_path = os.path.join(output_dir, "ppmi_mofa_input_long.csv")
    final_mofa_df.to_csv(long_out_path, index=False)
    print(f"Saved Long-Format MOFA Table: {long_out_path}")
    
    return aligned_views, aligned_metadata, final_mofa_df

# ------------------------------------------------------------------------------
# Example Execution
# ------------------------------------------------------------------------------
# rna_corrected: from Task 1.2.2 (genes x participant_id)
# csfp_clean: log2 + normalized from Task 1.2.3 (proteins x participant_id)
# plap_clean: log2 + normalized from Task 1.2.3 (proteins x participant_id)
# ppmi_metadata: master sample manifest (indexed by participant_id)

# aligned_views, aligned_meta, mofa_entry_df = validate_and_intersect_views(
#    rna_df=rna_corrected,
#    csf_df=csfp_clean,
#    plasma_df=plap_clean,
#    metadata_df=ppmi_metadata,
#    output_dir="./ppmi_mofa_inputs"
# )
