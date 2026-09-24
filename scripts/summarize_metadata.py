import pandas as pd
import sys
from scipy.stats import mannwhitneyu

def summarize_metadata(input_file):
    """
    Reads metadata, filters by study, visit_month, and diagnosis,
    produces a summary table, performs a statistical test on age,
    and returns the filtered DataFrame.
    """
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return None

    # Ensure required columns exist
    required_columns = {'study', 'visit_month', 'diagnosis', 'participant_id', 'age_at_baseline', 'sex'}
    if not required_columns.issubset(df.columns):
        print(f"Error: Missing required columns. Need: {required_columns}")
        return None

    # Filter
    # Assuming visit_month is numeric, 0 matches both 0 and 0.0.
    filtered_df = df[
        (df['study'].isin(['PDBP', 'PPMI'])) &
        (df['visit_month'] == 0) &
        (df['diagnosis'].isin(['Ctrl', 'PD']))
    ].copy()

    # Grouping and Aggregation
    summary = filtered_df.groupby(['study', 'diagnosis', 'sex']).agg(
        sample_count=('participant_id', 'count'),
        age=('age_at_baseline', 'median')
    ).reset_index()

    print("--- Summary Table ---")
    print(summary)

    # Statistical Test (Mann-Whitney U for age between studies)
    ppmi_ages = filtered_df[filtered_df['study'] == 'PPMI']['age_at_baseline'].dropna()
    pdbp_ages = filtered_df[filtered_df['study'] == 'PDBP']['age_at_baseline'].dropna()
    
    if len(ppmi_ages) > 0 and len(pdbp_ages) > 0:
        stat, p_value = mannwhitneyu(ppmi_ages, pdbp_ages)
        print("\n--- Statistical Test (Mann-Whitney U) ---")
        print("Comparison: PPMI vs PDBP cohorts (Age)")
        print(f"Statistic: {stat:.4f}")
        print(f"P-value: {p_value:.4e}")
    else:
        print("\n--- Statistical Test ---")
        print("Insufficient data for Mann-Whitney U test between PPMI and PDBP cohorts.")

    return filtered_df
