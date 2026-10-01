#!/usr/bin/env python3
"""
mofa/scripts/query_comprehensive_metadata.py

Converted and expanded Python script to query Google BigQuery directly from a Terra workspace,
joining demographics, case-control status, WGS mutations, MDS-UPDRS Part III motor scores, 
and MoCA cognitive scores at baseline (visit_month == 0.0).
"""

import os
import time
import pandas as pd
from google.cloud import bigquery
import firecloud.api as fapi

def get_workspace_attributes():
    """Fetches workspace attributes from FireCloud API in a Terra environment."""
    namespace = os.getenv('WORKSPACE_NAMESPACE')
    workspace = os.getenv('WORKSPACE_NAME')
    
    if not namespace or not workspace:
        raise ValueError("WORKSPACE_NAMESPACE and WORKSPACE_NAME environment variables must be set.")
        
    resp = fapi.get_workspace(namespace, workspace)
    if resp.status_code != 200:
        raise RuntimeError(f"Failed to fetch workspace attributes: {resp.text}")
        
    return resp.json()['workspace']['attributes']

def query_master_clinical_metadata(cohort: str = "PDBP") -> pd.DataFrame:
    """Queries BigQuery for comprehensive clinical, demographic, and omic metadata."""
    start_time = time.time()
    
    # Initialize BigQuery client
    project_id = os.getenv('GOOGLE_PROJECT')
    client = bigquery.Client(project=project_id)
    
    # Get dataset IDs from workspace attributes
    attrs = get_workspace_attributes()
    
    tier1_ds = attrs.get('bq_ds_amp_pd_tier1_release', '').replace(':', '.')
    tier2_ds = attrs.get('bq_ds_amp_pd_tier2_release', '').replace(':', '.')
    rnaseq_ds = attrs.get('bq_ds_amp_pd_transcriptomics_rnaseq_release', '').replace(':', '.')
    
    print(f"Tier 1 Dataset: {tier1_ds}")
    print(f"Tier 2 Dataset: {tier2_ds}")
    print(f"RNA-seq Dataset: {rnaseq_ds}")
    
    query = f"""
    WITH baseline_clinical AS (
        SELECT
            p.participant_id,
            p.study,
            d.sex,
            d.age_at_baseline,
            d.ethnicity,
            d.race,
            d.education_level_years,
            cc.diagnosis_at_baseline,
            m.has_known_PD_Mutation_in_WGS,
            u3.mds_updrs_part_iii_summary_score,
            moca.moca_total_score
        FROM `{tier1_ds}.amp_pd_participants` p
        LEFT JOIN `{tier1_ds}.Demographics` d ON p.participant_id = d.participant_id AND d.visit_month = 0.0
        LEFT JOIN `{tier1_ds}.amp_pd_case_control` cc ON p.participant_id = cc.participant_id
        LEFT JOIN `{tier2_ds}.amp_pd_participant_mutations` m ON p.participant_id = m.participant_id
        LEFT JOIN `{tier1_ds}.MDS_UPDRS_Part_III` u3 ON p.participant_id = u3.participant_id AND u3.visit_month = 0.0
        LEFT JOIN `{tier1_ds}.MOCA` moca ON p.participant_id = moca.participant_id AND moca.visit_month = 0.0
        WHERE p.study = @cohort
    )
    SELECT * FROM baseline_clinical
    ORDER BY participant_id
    """
    
    job_config = bigquery.QueryJobConfig(
        query_parameters=[
            bigquery.ScalarQueryParameter("cohort", "STRING", cohort)
        ]
    )
    
    print(f"Executing BigQuery query for cohort {cohort}...")
    df = client.query(query, job_config=job_config).to_dataframe()
    
    elapsed = time.time() - start_time
    print(f"Query completed in {elapsed:.2f} seconds. Retrieved {len(df)} records.")
    return df

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Query comprehensive clinical and omic metadata from BigQuery on Terra.")
    parser.add_argument("--cohort", default="PDBP", help="Cohort name (e.g. PDBP or PPMI)")
    parser.add_argument("--output", default="comprehensive_baseline_clinical.csv", help="Output CSV path")
    
    args = parser.parse_args()
    
    df = query_master_clinical_metadata(cohort=args.cohort)
    df.to_csv(args.output, index=False)
    print(f"Saved comprehensive metadata to: {args.output}")

if __name__ == "__main__":
    main()
