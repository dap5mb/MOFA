import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet

def preprocess_rnaseq(counts_file=str, 
                           metadata_file=str, 
                           output_file=str, 
                           visit_filter=0.0, 
                           cohort_name=str, 
                           min_count=10, 
                           min_prop=0.20):
    print("Loading data...")
    
    # 1. Load data
    counts_raw = pd.read_csv(counts_file, index_col=0, low_memory=False)
    col_data = pd.read_csv(metadata_file, low_memory=False)
    
    # Handle duplicate gene rows in counts if any exist
    if counts_raw.index.duplicated().any():
        print("Warning: Duplicate feature names found in counts file. Aggregating duplicates by summing...")
        counts_raw = counts_raw.groupby(counts_raw.index).sum()
        
    # Ensure metadata uses participant_id as index to match the updated count columns
    id_col = 'participant_id' if 'participant_id' in col_data.columns else col_data.columns[0]
    
    if col_data[id_col].duplicated().any():
        print(f"Warning: Duplicate entries found for '{id_col}' in metadata. Keeping the first occurrence...")
        col_data = col_data.drop_duplicates(subset=[id_col], keep='first')
        
    col_data = col_data.set_index(id_col)
    
    # 2. Filter metadata strictly to the target study/cohort
    if 'study' in col_data.columns:
        col_data = col_data[col_data['study'] == cohort_name]
        print(f"Filtered metadata to {len(col_data)} samples belonging to study: {cohort_name}")
    else:
        print("Warning: 'study' column not found in metadata. Proceeding with all rows.")

    # 3. Filter for baseline visit month
    if 'visit_month' in col_data.columns:
        col_data = col_data[col_data['visit_month'] == visit_filter]
        print(f"Filtered to {len(col_data)} samples where visit_month == {visit_filter}")

    # 4. Intersect sample IDs between counts columns and metadata index
    common_samples = counts_raw.columns.intersection(col_data.index)
    counts = counts_raw[common_samples].T # PyDESeq2 expects samples as rows, genes as columns
    col_data = col_data.loc[common_samples]
    
    print(f"Matrix aligned. Total samples for analysis: {counts.shape[0]}")

    # 5. Filter Lowly Expressed Genes
    min_samples = int(np.ceil(min_prop * counts.shape[0]))
    keep_genes = (counts >= min_count).sum(axis=0) >= min_samples
    counts_filtered = counts.loc[:, keep_genes]
    print(f"Retained {counts_filtered.shape[1]} genes out of {counts.shape[1]} total genes after low-count filtering.")

    # PyDESeq2 requires integer counts
    counts_filtered = counts_filtered.round().astype(int)

    # 6. Apply Size-Factor Normalization + log1p Transformation (Robust Fallback)
    print("Applying size-factor normalization and log1p transformation...")
    dds = DeseqDataSet(
        counts=counts_filtered,
        metadata=col_data,
        design_factors="1",
        ref_level=None
    )
    
    # Run fitting to compute reliable size factors
    dds.deseq2()
    
    # Extract size factors and normalize counts (Genes x Samples format)
    # counts_filtered has samples as rows, genes as columns in PyDESeq2
    size_factors = dds.obs['size_factors']
    
    # Normalize counts by sequencing depth / size factors, then apply log1p
    normed_counts = counts_filtered.div(size_factors, axis=0)
    expr_matrix_vst = np.log1p(normed_counts).T  # Transpose back to Genes x Samples
    
    # 7. Select Top Highly Variable Genes (HVGs) for MOFA
    print(f"Selecting top 5000 highly variable genes for MOFA...")
    
    # Calculate Variance or Median Absolute Deviation (MAD) across samples
    gene_variances = expr_matrix_vst.var(axis=1) # or scipy.stats.median_abs_deviation
    top_genes = gene_variances.nlargest(5000).index
    
    expr_matrix_mofa = expr_matrix_vst.loc[top_genes]
    
    # 8. Feature-wise Standardization (Crucial for MOFA)
    # MOFA performs best when features in each view are zero-centered
    expr_matrix_mofa = expr_matrix_mofa.sub(expr_matrix_mofa.mean(axis=1), axis=0)
    
    # Save output
    expr_matrix_mofa.to_csv(output_file)
    print(f"Preprocessing complete! Output saved to {output_file}")
    
    return expr_matrix_mofa