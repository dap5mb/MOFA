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
    
    # Some participants have more than one baseline MDS-UPDRS Part III record (repeated exams,
    # e.g. off and on medication; visit_name "#2"/"#3" marks repeats). Keep one per participant,
    # preferring the off-medication exam, then untreated, then state not recorded, then on
    # medication; ties go to the first record by visit_name.
    query = f"""
    WITH u3_baseline AS (
        SELECT
            participant_id,
            mds_updrs_part_iii_summary_score,
            upd23a_medication_for_pd AS updrs3_medication_for_pd,
            upd23b_clinical_state_on_medication AS updrs3_clinical_state,
            visit_name AS updrs3_visit_name,
            COUNT(*) OVER (PARTITION BY participant_id) AS updrs3_n_baseline_records
        FROM `{tier1_ds}.MDS_UPDRS_Part_III`
        WHERE visit_month = 0.0
        QUALIFY ROW_NUMBER() OVER (
            PARTITION BY participant_id
            ORDER BY
                CASE
                    WHEN LOWER(upd23b_clinical_state_on_medication) LIKE '%off%' THEN 0
                    WHEN LOWER(upd23a_medication_for_pd) IN ('no', '0') THEN 1
                    WHEN upd23b_clinical_state_on_medication IS NULL THEN 2
                    ELSE 3
                END,
                visit_name
        ) = 1
    ),
    baseline_clinical AS (
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
            u3.updrs3_medication_for_pd,
            u3.updrs3_clinical_state,
            u3.updrs3_visit_name,
            u3.updrs3_n_baseline_records,
            moca.moca_total_score
        FROM `{tier1_ds}.amp_pd_participants` p
        LEFT JOIN `{tier1_ds}.Demographics` d ON p.participant_id = d.participant_id AND d.visit_month = 0.0
        LEFT JOIN `{tier1_ds}.amp_pd_case_control` cc ON p.participant_id = cc.participant_id
        LEFT JOIN `{tier2_ds}.amp_pd_participant_mutations` m ON p.participant_id = m.participant_id
        LEFT JOIN u3_baseline u3 ON p.participant_id = u3.participant_id
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
    print(f"Query completed in {elapsed:.2f} seconds. Retrieved {len(df)} records for {df['participant_id'].nunique()} participants.")

    # Show which MDS-UPDRS III exam was kept, so the selection rule can be checked against the codes
    repeated = df.drop_duplicates("participant_id")
    repeated = repeated[repeated["updrs3_n_baseline_records"] > 1]
    print(f"Participants with more than one baseline MDS-UPDRS III record: {len(repeated)}")
    if len(repeated):
        print("  Kept exam, clinical state on medication (3b): " +
              ", ".join(f"{k} ({v})" for k, v in repeated["updrs3_clinical_state"].value_counts(dropna=False).items()))
        print("  Kept exam, medication for PD (3a): " +
              ", ".join(f"{k} ({v})" for k, v in repeated["updrs3_medication_for_pd"].value_counts(dropna=False).items()))
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
