#!/usr/bin/env python3
"""
mofa/scripts/prepare_pdbp_metadata.py

Harmonizes and preprocesses clinical and demographic metadata for the PDBP baseline validation cohort (visit_month == 0.0).
"""

import os
import pandas as pd

def process_pdbp_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters for PDBP participants at baseline (visit_month == 0.0), cleans missing sex,
    and harmonizes demographic and clinical features.
    """
    # 1. Filter for PDBP and baseline visit
    if "study" in df.columns:
        filtered = df[(df["study"] == "PDBP") & (df["visit_month"] == 0.0)].copy()
    elif "cohort" in df.columns:
        filtered = df[(df["cohort"] == "PDBP") & (df["visit_month"] == 0.0)].copy()
    else:
        raise ValueError("Metadata must contain 'study' or 'cohort' column.")

    # 2. Drop rows with missing sex or participant_id
    initial_count = len(filtered)
    if "sex" in filtered.columns:
        filtered = filtered.dropna(subset=["sex"])
    if "participant_id" in filtered.columns:
        filtered = filtered.dropna(subset=["participant_id"])
        
    dropped = initial_count - len(filtered)
    if dropped > 0:
        print(f"Dropped {dropped} PDBP baseline records with missing sex or participant_id.")

    # 3. Ensure required harmonization columns exist
    expected_cols = ["participant_id", "sex", "age_at_baseline", "diagnosis"]
    for col in expected_cols:
        if col not in filtered.columns:
            print(f"Warning: Expected column '{col}' not found in PDBP metadata.")

    # 4. Audit availability of secondary clinical scores (MDS-UPDRS, MoCA, etc.)
    clinical_scores = ["mds_updrs_3", "moca", "updrs_iii"]
    available_scores = [c for c in clinical_scores if c in filtered.columns]
    print(f"Audited PDBP clinical scores available: {available_scores}")

    return filtered

def main():
    import sys
    input_path = sys.argv[1] if len(sys.argv) > 1 else "sample_summary.csv"
    output_path = sys.argv[2] if len(sys.argv) > 2 else "mofa/data/pdbp_baseline_metadata.csv"
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    df = pd.read_csv(input_path)
    processed = process_pdbp_metadata(df)
    processed.to_csv(output_path, index=False)
    print(f"PDBP baseline metadata processed: {len(processed)} records saved to {output_path}")

if __name__ == "__main__":
    main()
