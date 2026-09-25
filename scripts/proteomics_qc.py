import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def verify_proteomics_scaling(df: pd.DataFrame, view_name: str):
    """
    Checks distribution, scale, and feature variances of pre-processed proteomics matrix.
    df: shape (proteins x participant_id)
    """
    # Extract flattened array of non-NA values
    vals = df.values.flatten()
    vals = vals[~np.isnan(vals)]
    
    val_min, val_max = vals.min(), vals.max()
    val_median = np.median(vals)
    val_mean = np.mean(vals)
    
    print(f"=== {view_name} Distribution Summary ===")
    print(f"Min: {val_min:.4f} | Max: {val_max:.4f}")
    print(f"Median: {val_median:.4f} | Mean: {val_mean:.4f}")
    
    # Infer transformation status
    if val_median > 50:
        print("--> Status: Appears to be on RAW intensity scale. Log2 transformation recommended.")
    else:
        print("--> Status: Appears to be LOG-TRANSFORMED / CENTERED. No further log transform needed.")
        
    # Check feature variances across proteins
    gene_vars = df.var(axis=1, skipna=True)
    print(f"Median Feature Variance: {gene_vars.median():.4f}")
    print(f"Variance IQR: {gene_vars.quantile(0.25):.4f} - {gene_vars.quantile(0.75):.4f}\n")
    
    return gene_vars

def normalize_proteomics(raw_df: pd.DataFrame, log2_transform=True) -> pd.DataFrame:
    # Step 1: Optional Log2 Transformation
    if log2_transform:
        norm_df = np.log2(raw_df + 1)
    else:
        norm_df = raw_df.copy()
    
    # Step 2: Sample-Level Median Centering
    sample_medians = norm_df.median(axis=0, skipna=True)
    global_median = sample_medians.median(skipna=True)
    
    # Subtract sample offset so every sample shares the exact same median
    norm_df = norm_df.sub(sample_medians - global_median, axis=1)
    
    return norm_df

# Load batch-corrected protein matrices
#csf_df = pd.read_csv("PLA-PDIA_protein_batch-corrected_matrix_csf.csv", index_col=0)
#plasma_df = pd.read_csv("PLA-PDIA_protein_batch-corrected_matrix_plasma.csv", index_col=0)

# Run verification checks
#csf_vars = verify_proteomics_scaling(csf_df, "CSF Proteomics")
#plasma_vars = verify_proteomics_scaling(plasma_df, "Plasma Proteomics")