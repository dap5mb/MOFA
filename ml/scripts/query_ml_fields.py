#!/usr/bin/env python3
"""
ml/scripts/query_ml_fields.py

BigQuery (Terra only): the baseline fields the cross-cohort ML analysis set needs on top of
MOFA's metadata, one row per participant of the cohort. MOFA's own query
(query_comprehensive_metadata.py, via shared) already supplies sex, age at baseline,
diagnosis_at_baseline, has_known_PD_Mutation_in_WGS, MDS-UPDRS III (off-medication exam rule)
and MoCA; this query adds only:

  PD_Medical_History (tier 1)
    age_at_diagnosis   static, so taken from any record (baseline first, then nearest visit)
    use_of_pd_medication, on_levodopa, on_dopamine_agonist, on_other_pd_medications
                       from the latest record at or before baseline (visit_month <= 0:
                       baseline, else screening); LOG records without a visit month are not used
    medication_visit_name  which record the medication fields came from
  Clinically_Reported_Genetic_Status (tier 2)
    genetic_status_enrollment, genetic_status_wgs

Raw values are written unchanged; analysis_set.py interprets them and stops on any value it
does not recognize. The output is participant-level and stays on Terra.

    python3 ml/scripts/query_ml_fields.py --cohort PPMI --output results/ml/analysis_set/ppmi_ml_fields.csv
    python3 ml/scripts/query_ml_fields.py --cohort PDBP --output results/ml/analysis_set/pdbp_ml_fields.csv
"""

import argparse
import os

import pandas as pd

import shared

MEDICATION_FIELDS = ["use_of_pd_medication", "on_levodopa", "on_dopamine_agonist", "on_other_pd_medications"]
GENETIC_FIELDS = ["genetic_status_enrollment", "genetic_status_wgs"]


def build_query(tier1_ds: str, tier2_ds: str) -> str:
    med_cols = ", ".join(MEDICATION_FIELDS)
    return f"""
    WITH dx_age AS (
        SELECT participant_id, age_at_diagnosis
        FROM `{tier1_ds}.PD_Medical_History`
        WHERE age_at_diagnosis IS NOT NULL
        QUALIFY ROW_NUMBER() OVER (
            PARTITION BY participant_id ORDER BY ABS(IFNULL(visit_month, 999)), visit_name) = 1
    ),
    med AS (
        SELECT participant_id, {med_cols}, visit_name AS medication_visit_name
        FROM `{tier1_ds}.PD_Medical_History`
        WHERE visit_month <= 0 AND COALESCE({med_cols}) IS NOT NULL
        QUALIFY ROW_NUMBER() OVER (
            PARTITION BY participant_id ORDER BY visit_month DESC, visit_name) = 1
    ),
    gs AS (
        SELECT participant_id, {", ".join(GENETIC_FIELDS)}
        FROM `{tier2_ds}.Clinically_Reported_Genetic_Status`
        WHERE COALESCE({", ".join(GENETIC_FIELDS)}) IS NOT NULL
        QUALIFY ROW_NUMBER() OVER (
            PARTITION BY participant_id ORDER BY ABS(IFNULL(visit_month, 999)), visit_name) = 1
    )
    SELECT p.participant_id, dx_age.age_at_diagnosis,
           {", ".join("med." + c for c in MEDICATION_FIELDS)}, med.medication_visit_name,
           {", ".join("gs." + c for c in GENETIC_FIELDS)}
    FROM `{tier1_ds}.amp_pd_participants` p
    LEFT JOIN dx_age ON p.participant_id = dx_age.participant_id
    LEFT JOIN med ON p.participant_id = med.participant_id
    LEFT JOIN gs ON p.participant_id = gs.participant_id
    WHERE p.study = @cohort
    ORDER BY p.participant_id
    """


def query_ml_fields(cohort: str) -> pd.DataFrame:
    from google.cloud import bigquery

    attrs = shared.get_workspace_attributes()  # loaded on use: needs firecloud (Terra)
    tier1 = attrs.get("bq_ds_amp_pd_tier1_release", "").replace(":", ".")
    tier2 = attrs.get("bq_ds_amp_pd_tier2_release", "").replace(":", ".")
    client = bigquery.Client(project=os.getenv("GOOGLE_PROJECT"))
    config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ScalarQueryParameter("cohort", "STRING", cohort)])
    df = client.query(build_query(tier1, tier2), job_config=config).to_dataframe()
    if df["participant_id"].duplicated().any():
        raise RuntimeError("query returned more than one row for some participants")

    # Aggregate value counts only, so the vocabulary can be checked against analysis_set.py
    print(f"{cohort}: {len(df)} participants")
    for col in MEDICATION_FIELDS + GENETIC_FIELDS + ["medication_visit_name"]:
        counts = df[col].value_counts(dropna=False)
        print(f"  {col}: " + ", ".join(f"{k} ({v})" for k, v in counts.items()))
    print(f"  age_at_diagnosis recorded: {df['age_at_diagnosis'].notna().sum()}")
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cohort", required=True, choices=["PPMI", "PDBP"])
    p.add_argument("--output", required=True, help="Participant-level CSV (stays on Terra)")
    a = p.parse_args()
    df = query_ml_fields(a.cohort)
    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    df.to_csv(a.output, index=False)
    print(f"Saved: {a.output}")


if __name__ == "__main__":
    main()
