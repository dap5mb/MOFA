import pandas as pd

def filter_by_metadata(
    data_path: str,
    metadata_df: pd.DataFrame,
    sample_col: str = "sample_id",
    index_col: int = 0
) -> pd.DataFrame:
    """
    Filters a gene/protein expression matrix by sample IDs from filtered metadata file.
    
    Parameters:
    - data_path: Path to the main expression CSV (genes/proteins x sample_ids).
    - metadata_df: Pandas DataFrame containing the metadata.
    - sample_col: Name of the column in the metadata containing sample IDs.
    - output_path: Path to save the filtered CSV file (set to None if you don't want to save).
    - index_col: Column index or name to use as the row index (genes/proteins).
    
    Returns:
    - filtered_df: The filtered pandas DataFrame.
    """
    # 1. Load data and metadata
    data_df = pd.read_csv(data_path, index_col=index_col)
    
    # 2. Validate that the sample column exists in the metadata
    if sample_col not in metadata_df.columns:
        raise ValueError(
            f"Column '{sample_col}' not found in metadata. "
            f"Available columns: {list(metadata_df.columns)}"
        )
    
    # 3. Extract unique target sample IDs
    target_samples = metadata_df[sample_col].dropna().unique()
    
    # 4. Find the intersection to prevent key errors
    samples_to_keep = data_df.columns.intersection(target_samples)
    
    print(f"Total unique samples in metadata: {len(target_samples)}")
    print(f"Total sample columns in expression data: {len(data_df.columns)}")
    print(f"Samples matched and kept: {len(samples_to_keep)}")
    
    if len(samples_to_keep) == 0:
        print("Warning: No matching sample IDs found between the two files!")
        
    # 5. Filter the dataframe
    filtered_df = data_df[samples_to_keep]
        
    return filtered_df