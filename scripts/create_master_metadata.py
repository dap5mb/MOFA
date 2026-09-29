import os
import pandas as pd

def create_master_metadata(rnaseq_meta_path: str, proteomics_meta_path: str, cohort: str, output_path: str = None) -> pd.DataFrame:
    """
    Harmonizes and merges RNA-seq and proteomics metadata for baseline visits (visit_month == 0.0)
    to create a master sample manifest for MOFA analysis for a specific cohort.
    
    Parameters:
    -----------
    rnaseq_meta_path : str
        Path to the RNA-seq metadata CSV file.
    proteomics_meta_path : str
        Path to the proteomics metadata CSV file.
    output_path : str, optional
        Path to save the resulting master CSV file.
    cohort : str, optional
        Specific cohort to filter by ('PPMI' or 'PDBP'). 
        
    Returns:
    --------
    pd.DataFrame
        The unified master metadata dataframe for the specified cohort.
    """
    # 1. Load the metadata files
    rnaseq_df = pd.read_csv(rnaseq_meta_path)
    prot_df = pd.read_csv(proteomics_meta_path)
    
    # 2. Filter for baseline (visit_month == 0.0) and the specific study/cohort
    rnaseq_baseline = rnaseq_df[(rnaseq_df['visit_month'] == 0.0) & (rnaseq_df['study'] == cohort)].copy()
    prot_baseline = prot_df[(prot_df['visit_month'] == 0.0) & (prot_df['cohort'] == cohort)].copy()

    
    # 3. Rename RNA-seq sample ID to clearly indicate its view
    rnaseq_baseline = rnaseq_baseline.rename(columns={'sample_id': 'rnaseq_sample_id'})
    
    # 4. Pivot proteomics metadata to separate CSF and Plasma sample IDs by participant
    # Ensure tissue_type values match your data (e.g., 'CSF', 'Plasma')
    prot_pivot = prot_baseline.pivot_table(
        index='participant_id',
        columns='tissue_type',
        values='sample_id',
        aggfunc='first'
    ).reset_index()
    
    # Rename pivoted columns to standard identifiers
    rename_cols = {}
    if 'CSF' in prot_pivot.columns:
        rename_cols['CSF'] = 'csf_pro_sample_id'
    if 'Plasma' in prot_pivot.columns:
        rename_cols['Plasma'] = 'plasma_pro_sample_id'
    prot_pivot = prot_pivot.rename(columns=rename_cols)
    
    # Extract static clinical/demographic metadata attributes from proteomics to merge back
    # (Taking the first instance per participant to avoid duplication)
    prot_covariates = prot_baseline.drop_duplicates(subset=['participant_id'])[[
        'participant_id', 'APOE_genotype', 'APOE_E4_status', 
        'has_known_GBA_mutation_in_WGS', 'has_known_LRRK2_mutation_in_WGS',
        'has_known_SNCA_mutation_in_WGS', 'education_level_years'
    ]]
    
    # 5. Merge RNA-seq metadata with the pivoted proteomics IDs and covariates
    # Using an outer join to retain single-omic participants (handled natively by MOFA)
    master_df = pd.merge(
        rnaseq_baseline, 
        prot_pivot[['participant_id', 'csf_pro_sample_id', 'plasma_pro_sample_id']], 
        on='participant_id', 
        how='outer'
    )
    
    master_df = pd.merge(
        master_df,
        prot_covariates,
        on='participant_id',
        how='left'
    )
    
    # Clean up column overlap if study/cohort or diagnoses labels exist in both
    if 'cohort' in master_df.columns and 'study' in master_df.columns:
        master_df['study'] = master_df['study'].fillna(master_df['cohort'])
        master_df = master_df.drop(columns=['cohort'])

    # 6. Optional cohort filtering (already handled in step 2)
    master_df = master_df[master_df['study'] == cohort].copy()

    # Drop missing sex records if present
    if 'sex' in master_df.columns:
        master_df = master_df.dropna(subset=['sex'])

    # 7. Optional: Save to file if output path is provided
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        master_df.to_csv(output_path, index=False)
        print(f"Master metadata successfully saved to {output_path}")
        
    print(f"Total unique baseline participants ({cohort}): {master_df['participant_id'].nunique()}")
    print(f"Participants with RNA-seq: {master_df['rnaseq_sample_id'].notna().sum()}")
    print(f"Participants with CSF Proteomics: {master_df['csf_pro_sample_id'].notna().sum()}")
    print(f"Participants with Plasma Proteomics: {master_df['plasma_pro_sample_id'].notna().sum()}")
    
    return master_df

# --- Example Usage ---
# master_metadata = create_master_metadata(
#     rnaseq_meta_path='path/to/rnaseq_meta.csv',
#     proteomics_meta_path='path/to/proteomics_meta.csv',
#     output_path='path/to/master_baseline_metadata.csv',
#     cohort = "PPMI"
# )