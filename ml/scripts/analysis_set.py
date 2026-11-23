#!/usr/bin/env python3
"""
ml/scripts/analysis_set.py

Baseline analysis set for one cohort (spec step 2, story 1-7; issue 03). Participants are the
MOFA participant definition restricted to whole-blood RNA-seq: baseline (visit_month = 0),
recorded sex, one row per participant (first record kept, as MOFA does), and present in the
cohort's MOFA RNA-seq matrix. Diagnosis is harmonized to PD / Control / Other; Other is kept
in the output with in_analysis_set = False so the coverage report can count it.

Inputs
  MOFA metadata  data/metadata_m0_<COHORT>.csv on Terra: MOFA's master metadata already merged
                 with query_comprehensive_metadata.py (sex, age_at_baseline, diagnosis columns,
                 has_known_PD_Mutation_in_WGS, MDS-UPDRS III by the off-medication exam rule, MoCA)
  ML fields      query_ml_fields.py output (age at diagnosis, baseline medication, genetic status)
  RNA-seq matrix optional; its header (participant IDs) defines who has RNA-seq. Without it,
                 a non-missing rnaseq_sample_id does.

Derived fields
  diagnosis      PPMI from 'diagnosis', PDBP from 'diagnosis_at_baseline' (as in the MOFA
                 write-up); labels outside PD_LABELS / CONTROL_LABELS become Other.
  carrier        True if has_known_PD_Mutation_in_WGS is yes OR genetic_status_enrollment names a
                 PD mutation; False if neither is positive and at least one is negative; missing
                 if both are missing. Both inputs are kept (carrier_wgs, carrier_enrollment).
  disease_duration_years
                 age_at_baseline - age_at_diagnosis, PD only. Ages are whole years, so -1 <= d < 0
                 is read as rounding and set to 0; d < -1 is set to missing (duration_status
                 records which).
  pd_medication  True if any PD_Medical_History medication field (use_of_pd_medication,
                 on_levodopa, on_dopamine_agonist, on_other_pd_medications) is yes; False if
                 use_of_pd_medication is no, or all three drug-class fields are no; otherwise the
                 MDS-UPDRS 3a answer (medication for PD) of the selected baseline exam; otherwise
                 missing. pd_medication_source records which.
  moca           baseline MoCA, else the latest screening score (query_ml_fields.py); this
                 replaces the baseline-only score in MOFA's metadata. moca_visit_month < 0
                 marks a screening score.

Any value outside the known vocabularies stops the run with the list of values, so a new
coding is mapped on purpose rather than silently.

Duration bands and the PPMI-like subset were fixed from the coverage report in the step-2 spec
addendum of 2026-11-23 (FIXED_DEFINITIONS; duration_band() and ppmi_like() apply them).

    python3 ml/scripts/analysis_set.py --cohort PDBP --mofa_metadata data/metadata_m0_PDBP.csv \\
        --ml_fields results/ml/analysis_set/pdbp_ml_fields.csv \\
        --rnaseq_matrix data/pdbp_rnaseq_adjusted.csv \\
        --output results/ml/analysis_set/pdbp_participants.csv
"""

import argparse
import os

import numpy as np
import pandas as pd

COHORTS = ("PPMI", "PDBP")
SEXES = ("Female", "Male")
DIAGNOSIS_SOURCE = {"PPMI": "diagnosis", "PDBP": "diagnosis_at_baseline"}
PD_LABELS = {"pd", "parkinson's disease"}
CONTROL_LABELS = {"ctrl", "control", "healthy control", "no pd nor other neurological disorder"}

_YES = {"yes", "y", "1", "1.0", "true"}
_NO = {"no", "n", "0", "0.0", "false"}
_MISSING = {"", "nan", "unknown", "not tested", "na", "n/a", "not available", "not done",
            "not applicable", "unknown/not collected as enrollment criterion"}
_GENETIC_NEGATIVE = {"no mutation", "no mutations", "negative", "no", "non-carrier", "noncarrier",
                     "not a carrier", "wild type", "wildtype", "wt", "none", "no pd mutation"}
_GENETIC_POSITIVE_TOKENS = ("lrrk2", "gba", "snca", "+", "positive", "carrier", "mutation")

MEDICATION_HISTORY_FIELDS = ["use_of_pd_medication", "on_levodopa", "on_dopamine_agonist",
                             "on_other_pd_medications"]
UPDRS3_MEDICATION = "updrs3_medication_for_pd"

# Step-2 spec addendum 2026-11-23, fixed before any PDBP classifier score. Durations are whole
# years: bands [0, 2), [2, 5), [5, 10), [10, inf) are 0-1, 2-4, 5-9 and >= 10 years. The PPMI-like
# subset is PD with duration <= 2 years (no medication criterion) plus all controls.
FIXED_DEFINITIONS = {
    "addendum": "2026-11-23",
    "duration_band_edges": [0, 2, 5, 10, np.inf],
    "duration_band_labels": ["0-1", "2-4", "5-9", ">=10"],
    "ppmi_like_max_duration": 2,
    "ppmi_like_unmedicated": False,
}


REQUIRED_METADATA = ["sex", "age_at_baseline", "has_known_PD_Mutation_in_WGS",
                     "mds_updrs_part_iii_summary_score"]
REQUIRED_ML_FIELDS = (["participant_id", "age_at_diagnosis", "genetic_status_enrollment", "moca_total_score",
                       "moca_visit_month"] + MEDICATION_HISTORY_FIELDS)


def _require(df, cols, what):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{what} is missing columns {missing}")


def require_fixed_definitions() -> dict:
    """Cohort-shift steps call this first: no bands or subset before the addendum fixes them."""
    if FIXED_DEFINITIONS is None:
        raise RuntimeError("Duration bands and the PPMI-like subset are not fixed yet: run the "
                           "coverage report, write the step-2 spec addendum, then set "
                           "analysis_set.FIXED_DEFINITIONS")
    return FIXED_DEFINITIONS


def duration_band(participants: pd.DataFrame) -> pd.Series:
    """Fixed duration band per PD case; missing for controls and PD without duration."""
    d = require_fixed_definitions()
    band = pd.cut(participants["disease_duration_years"], d["duration_band_edges"], right=False,
                  labels=d["duration_band_labels"])
    return band.where(participants["diagnosis"] == "PD")


def ppmi_like(participants: pd.DataFrame) -> pd.Series:
    """True for controls and for PD cases meeting the fixed PPMI-like definition."""
    d = require_fixed_definitions()
    is_pd = participants["diagnosis"] == "PD"
    keep = is_pd & (participants["disease_duration_years"] <= d["ppmi_like_max_duration"])
    if d["ppmi_like_unmedicated"]:
        keep &= participants["pd_medication"].eq(False).fillna(False)
    return (participants["diagnosis"] == "Control") | keep


def _norm(v) -> str:
    return "" if pd.isna(v) else str(v).strip().lower()


def _unrecognized(values, field):
    raise ValueError(f"{field}: unrecognized values {sorted(set(map(str, values)))}; "
                     f"add them to the vocabulary in analysis_set.py")


def parse_yes_no(series: pd.Series, field: str) -> pd.Series:
    """True / False / <NA>; stops on any other value."""
    s = series.map(_norm)
    bad = s[~s.isin(_YES | _NO | _MISSING)]
    if len(bad):
        _unrecognized(series[bad.index].unique(), field)
    return s.map(lambda v: True if v in _YES else (False if v in _NO else pd.NA)).astype("boolean")


def parse_genetic_status(series: pd.Series, field: str) -> pd.Series:
    """True if the status names a PD mutation, False if it says none; <NA> if missing."""
    def one(v):
        n = _norm(v)
        if n in _MISSING:
            return pd.NA
        if n in _GENETIC_NEGATIVE or n.startswith(("no ", "non-", "not ")):
            return False
        if any(t in n for t in _GENETIC_POSITIVE_TOKENS):
            return True
        return "?"
    out = series.map(one)
    bad = out[out.astype(str) == "?"]
    if len(bad):
        _unrecognized(series[bad.index].unique(), field)
    return out.astype("boolean")


def harmonize_diagnosis(raw: pd.Series) -> pd.Series:
    n = raw.map(_norm)
    return pd.Series(np.select([n.isin(PD_LABELS), n.isin(CONTROL_LABELS)], ["PD", "Control"], "Other"),
                     index=raw.index)


def derive_carrier(wgs: pd.Series, enrollment: pd.Series) -> pd.Series:
    either_pos = (wgs == True).fillna(False) | (enrollment == True).fillna(False)  # noqa: E712
    any_neg = (wgs == False).fillna(False) | (enrollment == False).fillna(False)  # noqa: E712
    out = pd.Series(pd.NA, index=wgs.index, dtype="boolean")
    out[any_neg] = False
    out[either_pos] = True
    return out


def derive_duration(age_at_baseline, age_at_diagnosis, is_pd):
    d = pd.to_numeric(age_at_baseline, errors="coerce") - pd.to_numeric(age_at_diagnosis, errors="coerce")
    status = pd.Series(np.select(
        [~is_pd, d.isna(), d < -1, d < 0], ["not_pd", "missing", "negative_set_missing", "rounded_to_zero"],
        "ok"), index=d.index)
    d = d.where(is_pd & (d >= -1)).clip(lower=0)
    return d, status


def derive_medication(df: pd.DataFrame):
    hist = pd.DataFrame({f: parse_yes_no(df[f], f) if f in df else pd.Series(pd.NA, index=df.index, dtype="boolean")
                         for f in MEDICATION_HISTORY_FIELDS})
    any_yes = hist.eq(True).fillna(False).any(axis=1)
    classes = hist[MEDICATION_HISTORY_FIELDS[1:]]
    hist_no = hist["use_of_pd_medication"].eq(False).fillna(False) | classes.eq(False).fillna(False).all(axis=1)

    u3 = (parse_yes_no(df[UPDRS3_MEDICATION], UPDRS3_MEDICATION) if UPDRS3_MEDICATION in df
          else pd.Series(pd.NA, index=df.index, dtype="boolean"))

    med = pd.Series(pd.NA, index=df.index, dtype="boolean")
    src = pd.Series(pd.NA, index=df.index, dtype="object")
    use_u3 = ~any_yes & ~hist_no & u3.notna()
    med[use_u3] = u3[use_u3]
    src[use_u3] = "updrs3_3a"
    med[hist_no] = False
    med[any_yes] = True
    src[hist_no | any_yes] = "medical_history"
    return med, src


def build_participants(mofa_meta: pd.DataFrame, ml_fields: pd.DataFrame, cohort: str,
                       rnaseq_ids=None) -> pd.DataFrame:
    """One row per MOFA RNA-seq participant of the cohort, indexed by participant_id."""
    if cohort not in COHORTS:
        raise ValueError(f"cohort must be one of {COHORTS}")
    meta = mofa_meta.reset_index() if "participant_id" not in mofa_meta.columns else mofa_meta.copy()
    if "study" in meta:
        meta = meta[meta["study"] == cohort]
    if "visit_month" in meta:
        meta = meta[meta["visit_month"] == 0.0]
    meta = meta.dropna(subset=["participant_id", "sex"])
    meta = meta.drop_duplicates("participant_id", keep="first").set_index("participant_id")

    if rnaseq_ids is not None:
        ids = pd.Index(rnaseq_ids).drop_duplicates()
        missing = ids.difference(meta.index)
        if len(missing):
            raise ValueError(f"{len(missing)} RNA-seq participants have no baseline MOFA metadata")
        meta = meta.loc[ids]
    else:
        meta = meta[meta["rnaseq_sample_id"].notna()]

    bad_sex = set(meta["sex"]) - set(SEXES)
    if bad_sex:
        raise ValueError(f"sex must be one of {SEXES}; found {sorted(map(str, bad_sex))}")

    _require(meta, REQUIRED_METADATA + [DIAGNOSIS_SOURCE[cohort]], f"{cohort} MOFA metadata")
    _require(ml_fields, REQUIRED_ML_FIELDS, "ML fields")
    fields = ml_fields.drop_duplicates("participant_id").set_index("participant_id")
    meta = meta.drop(columns=[c for c in fields.columns if c in meta.columns]).join(fields, how="left")

    dx_col = DIAGNOSIS_SOURCE[cohort]
    out = pd.DataFrame(index=meta.index)
    out["cohort"] = cohort
    out["sex"] = meta["sex"]
    out["diagnosis_raw"] = meta[dx_col]
    out["diagnosis"] = harmonize_diagnosis(meta[dx_col])
    out["in_analysis_set"] = out["diagnosis"].isin(["PD", "Control"])
    out["age_at_baseline"] = pd.to_numeric(meta["age_at_baseline"], errors="coerce")

    out["carrier_wgs"] = parse_yes_no(meta["has_known_PD_Mutation_in_WGS"], "has_known_PD_Mutation_in_WGS")
    out["carrier_enrollment"] = parse_genetic_status(meta["genetic_status_enrollment"], "genetic_status_enrollment")
    out["carrier"] = derive_carrier(out["carrier_wgs"], out["carrier_enrollment"])

    is_pd = out["diagnosis"] == "PD"
    out["age_at_diagnosis"] = pd.to_numeric(meta["age_at_diagnosis"], errors="coerce")
    out["disease_duration_years"], out["duration_status"] = derive_duration(
        out["age_at_baseline"], out["age_at_diagnosis"], is_pd)

    out["pd_medication"], out["pd_medication_source"] = derive_medication(meta)
    out["updrs3"] = pd.to_numeric(meta["mds_updrs_part_iii_summary_score"], errors="coerce")
    # MoCA from query_ml_fields.py (baseline, else screening), which replaces MOFA's baseline-only score
    out["moca"] = pd.to_numeric(meta["moca_total_score"], errors="coerce")
    out["moca_visit_month"] = pd.to_numeric(meta["moca_visit_month"], errors="coerce")
    for col in ("medication_visit_name", "updrs3_clinical_state", "plate", "rnaseq_sample_id"):
        if col in meta:
            out[col] = meta[col]
    out.index.name = "participant_id"
    return out


def analysis_set(participants: pd.DataFrame) -> pd.DataFrame:
    """PD and control participants only: what the classifiers use."""
    return participants[participants["in_analysis_set"]].copy()


def read_matrix_ids(path: str) -> list:
    """Participant IDs from a genes x participants matrix header (first column = gene IDs)."""
    return list(pd.read_csv(path, nrows=0).columns[1:])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cohort", required=True, choices=COHORTS)
    p.add_argument("--mofa_metadata", required=True)
    p.add_argument("--ml_fields", required=True, help="query_ml_fields.py output")
    p.add_argument("--rnaseq_matrix", help="MOFA RNA-seq matrix (genes x participants) defining who has RNA-seq")
    p.add_argument("--output", required=True, help="Participant-level CSV (stays on Terra)")
    a = p.parse_args()

    ids = read_matrix_ids(a.rnaseq_matrix) if a.rnaseq_matrix else None
    df = build_participants(pd.read_csv(a.mofa_metadata), pd.read_csv(a.ml_fields), a.cohort, ids)
    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    df.to_csv(a.output)
    print(f"{a.cohort}: {len(df)} RNA-seq participants, {int(df['in_analysis_set'].sum())} PD or control")
    print(pd.crosstab(df["diagnosis"], df["sex"], margins=True).to_string())
    print(f"Saved: {a.output}")


if __name__ == "__main__":
    main()
