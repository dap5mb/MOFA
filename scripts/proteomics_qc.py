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

def filter_proteomics(raw_df: pd.DataFrame, max_sample_zero_frac=0.50, max_protein_zero_frac=0.30,
                      view_name: str = "Proteomics") -> pd.DataFrame:
    """
    QC filter on a raw intensity matrix (proteins x participant_id) where exact zeros are
    non-detects. Drops samples with more than max_sample_zero_frac zeros first, then proteins
    with more than max_protein_zero_frac zeros among the remaining samples.
    Pass None for either threshold to skip that filter.
    """
    df = raw_df
    if max_sample_zero_frac is not None:
        sample_zero = (df == 0).mean(axis=0)
        bad_samples = sample_zero.index[sample_zero > max_sample_zero_frac]
        print(f"{view_name}: dropping {len(bad_samples)} / {df.shape[1]} samples with > {max_sample_zero_frac:.0%} zeros")
        df = df.drop(columns=bad_samples)
    if max_protein_zero_frac is not None:
        protein_zero = (df == 0).mean(axis=1)
        bad_proteins = protein_zero.index[protein_zero > max_protein_zero_frac]
        print(f"{view_name}: dropping {len(bad_proteins)} / {df.shape[0]} proteins with > {max_protein_zero_frac:.0%} zeros")
        df = df.drop(index=bad_proteins)
    return df

def normalize_proteomics(raw_df: pd.DataFrame, log2_transform=True, zeros_as_missing=True) -> pd.DataFrame:
    # Step 0: Exact zeros in raw intensities are non-detects, not measurements. Kept as values they
    # become log2(0 + 1) = 0, an extreme outlier (~17 log2 units below a typical protein) that
    # dominated the plasma view's variance. Mark them missing; MOFA handles NaN natively.
    if zeros_as_missing:
        raw_df = raw_df.mask(raw_df == 0)

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