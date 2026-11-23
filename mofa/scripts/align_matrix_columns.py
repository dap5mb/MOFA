import pandas as pd

def load_and_align_matrix(matrix_path: str, 
                          master_meta: pd.DataFrame, 
                          sample_id_col: str, 
                          feature_id_col: str | int = 0,
                          output_path: str = None) -> pd.DataFrame:
    """
    Reads an expression or proteomics matrix, retains the feature ID column,
    filters sample columns to match baseline samples, and maps assay-specific
    sample IDs to universal participant_id.
    
    Parameters:
    -----------
    matrix_path : str
        Path to the raw expression/proteomics matrix file (CSV/TSV).
    master_meta : pd.DataFrame
        The master baseline metadata dataframe containing the sample mapping.
    sample_id_col : str
        The column name in master_meta corresponding to this view's sample IDs 
        (e.g., 'rnaseq_sample_id', 'csf_pro_sample_id', 'plasma_pro_sample_id').
    feature_id_col : str or int, default=0
        Column name or integer position for the feature IDs (e.g., gene or protein identifiers).
        
    Returns:
    --------
    pd.DataFrame
        The aligned matrix with feature IDs retained and sample column names converted to participant_id.
    """
    # Load data without setting an index
    matrix_df = pd.read_csv(matrix_path)
    
    # Resolve feature ID column name if integer position was passed
    if isinstance(feature_id_col, int):
        feature_col_name = matrix_df.columns[feature_id_col]
    else:
        feature_col_name = feature_id_col
    
    # Create mapping dictionary: sample_id -> participant_id
    id_map = master_meta.dropna(subset=[sample_id_col]).set_index(sample_id_col)['participant_id'].to_dict()
    
    # Identify sample columns present in the baseline mapping dictionary
    valid_sample_cols = [c for c in matrix_df.columns if c in id_map]
    
    # Subset matrix: keep feature ID column + valid sample columns
    aligned_matrix = matrix_df[[feature_col_name] + valid_sample_cols].copy()
    
    # Rename sample columns from assay sample IDs to participant_id
    aligned_matrix = aligned_matrix.rename(columns=id_map)
    
    # Deduplicate participant ID columns (keeping feature column intact)
    sample_cols = aligned_matrix.columns[1:]
    dedup_sample_cols = sample_cols[~sample_cols.duplicated()]
    aligned_matrix = aligned_matrix[[feature_col_name] + list(dedup_sample_cols)]
    
    if output_path:
        aligned_matrix.to_csv(output_path, index=False)
        print(f"Aligned matrix successfully saved to {output_path}")
    
    print(f"Loaded matrix from {matrix_path}")
    print(f"  - Original sample columns: {matrix_df.shape[1] - 1}")
    print(f"  - Aligned baseline sample columns (mapped to participant_id): {aligned_matrix.shape[1] - 1}")
    
    return aligned_matrix