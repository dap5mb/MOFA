#!/usr/bin/env python3
"""
scripts/01_prepare_metadata.py

Filter clinical summary data for PPMI and PDBP participants at baseline (visit_month = 0.0)
and output baseline clinical metadata table.
"""

import os
import sys
import pandas as pd

def filter_baseline_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters clinical metadata for PPMI and PDBP participants at visit_month == 0.0.
    """
    valid_studies = {"PPMI", "PDBP"}
    filtered = df[(df["study"].isin(valid_studies)) & (df["visit_month"] == 0.0)].copy()
    return filtered

def main():
    # Handle snakemake or standalone execution
    if "snakemake" in globals():
        summary_path = snakemake.input.summary
        output_path = snakemake.output.metadata
    else:
        if len(sys.argv) > 2:
            summary_path = sys.argv[1]
            output_path = sys.argv[2]
        else:
            summary_path = "sample_summary.csv"
            output_path = "results/data/baseline_metadata.tsv"

    if not os.path.exists(summary_path):
        raise FileNotFoundError(f"Input file non-existent: {summary_path}")

    df = pd.read_csv(summary_path)
    baseline_df = filter_baseline_metadata(df)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    baseline_df.to_csv(output_path, sep="\t", index=False)
    print(f"Extracted {len(baseline_df)} baseline PPMI/PDBP cohort records -> {output_path}")

if __name__ == "__main__":
    main()
