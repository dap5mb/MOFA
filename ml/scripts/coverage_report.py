#!/usr/bin/env python3
"""
ml/scripts/coverage_report.py

Aggregate report on the analysis set (issue 03): participant counts against the MOFA write-up,
carrier flags and how the two sources agree, how complete disease duration and baseline PD
medication are among PD cases by cohort and sex, the duration distribution, and the size of
the PPMI-like subset under candidate definitions. Run before any PDBP classifier score exists;
the duration bands and PPMI-like definition are then fixed in the step-2 spec addendum.

Every output is a count, a percentage or a percentile of a group: nothing participant-level,
so the report directory is what may leave Terra.

    python3 ml/scripts/coverage_report.py \\
        --participants results/ml/analysis_set/ppmi_participants.csv \\
                       results/ml/analysis_set/pdbp_participants.csv \\
        --output_dir results/ml/analysis_set/coverage_report
"""

import argparse
import os

import numpy as np
import pandas as pd

from analysis_set import SEXES

# MOFA write-up (Methods, "Participants"). PDBP: 1,248 RNA-seq participants. PPMI: 1,308 with
# RNA-seq; the diagnosis and sex counts there describe all 1,329 participants with any view.
MOFA_REFERENCE = [
    ("PDBP", "RNA-seq participants", 1248),
    ("PDBP", "PD", 767), ("PDBP", "Control", 481), ("PDBP", "Other", 0),
    ("PDBP", "Female", 554), ("PDBP", "Male", 694),
    ("PPMI", "RNA-seq participants", 1308),
    ("PPMI", "PD (any view, n = 1,329)", 753), ("PPMI", "Control (any view, n = 1,329)", 555),
    ("PPMI", "Other (any view, n = 1,329)", 21),
]
CANDIDATE_MAX_DURATION = [1, 2, 3, 5]
CANDIDATE_BAND_EDGES = [0, 2, 5, 10, np.inf]
PERCENTILES = [10, 25, 50, 75, 90]


def _pct(num, den):
    return round(100 * num / den, 1) if den else np.nan


def _label(v):
    return "missing" if pd.isna(v) else ("yes" if v is True or v == "True" else "no")


def _bool(s: pd.Series) -> pd.Series:
    """Booleans read back from CSV come as True/False/NaN objects or strings."""
    return s.map(lambda v: pd.NA if pd.isna(v) else str(v) == "True").astype("boolean")


def load(paths) -> pd.DataFrame:
    df = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    for col in ("in_analysis_set", "carrier", "carrier_wgs", "carrier_enrollment", "pd_medication"):
        df[col] = _bool(df[col])
    return df


def reference_check(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, what, expected in MOFA_REFERENCE:
        c = df[df["cohort"] == cohort]
        key = what.split(" (")[0]
        if key == "RNA-seq participants":
            observed = len(c)
        elif key in SEXES:
            observed = int((c["sex"] == key).sum())
        else:
            observed = int((c["diagnosis"] == key).sum())
        rows.append({"cohort": cohort, "count": what, "mofa_writeup": expected, "observed": observed,
                     "difference": observed - expected})
    return pd.DataFrame(rows)


def counts(df: pd.DataFrame) -> pd.DataFrame:
    t = df.groupby(["cohort", "sex", "diagnosis"]).size().unstack("diagnosis", fill_value=0)
    return t.reset_index()


def excluded_diagnoses(df: pd.DataFrame) -> pd.DataFrame:
    other = df[df["diagnosis"] == "Other"]
    return other.groupby(["cohort", other["diagnosis_raw"].fillna("missing")]).size().rename("n").reset_index()


def carriers(df: pd.DataFrame) -> pd.DataFrame:
    a = df[df["in_analysis_set"]].assign(carrier=lambda d: d["carrier"].map(_label))
    return a.groupby(["cohort", "sex", "diagnosis", "carrier"]).size().unstack("carrier", fill_value=0).reset_index()


def carrier_sources(df: pd.DataFrame) -> pd.DataFrame:
    a = df[df["in_analysis_set"]]
    t = a.assign(wgs=a["carrier_wgs"].map(_label), enrollment=a["carrier_enrollment"].map(_label))
    return t.groupby(["cohort", "wgs", "enrollment"]).size().rename("n").reset_index()


def coverage(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    pd_cases = df[df["diagnosis"] == "PD"]
    for (cohort, sex), g in pd_cases.groupby(["cohort", "sex"]):
        n = len(g)
        dur = g["disease_duration_years"].notna()
        med = g["pd_medication"].notna()
        src = g["pd_medication_source"]
        rows.append({
            "cohort": cohort, "sex": sex, "n_pd": n,
            "duration_n": int(dur.sum()), "duration_pct": _pct(dur.sum(), n),
            "duration_negative_set_missing": int((g["duration_status"] == "negative_set_missing").sum()),
            "duration_rounded_to_zero": int((g["duration_status"] == "rounded_to_zero").sum()),
            "medication_n": int(med.sum()), "medication_pct": _pct(med.sum(), n),
            "medication_from_history": int((src == "medical_history").sum()),
            "medication_from_updrs3_3a": int((src == "updrs3_3a").sum()),
            "medicated_n": int(g["pd_medication"].eq(True).sum()),
            "unmedicated_n": int(g["pd_medication"].eq(False).sum()),
            "both_n": int((dur & med).sum()), "both_pct": _pct((dur & med).sum(), n),
            "updrs3_n": int(g["updrs3"].notna().sum()), "moca_n": int(g["moca"].notna().sum()),
        })
    return pd.DataFrame(rows)


def duration_distribution(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    pd_cases = df[(df["diagnosis"] == "PD") & df["disease_duration_years"].notna()]
    for (cohort, sex), g in pd_cases.groupby(["cohort", "sex"]):
        d = g["disease_duration_years"]
        row = {"cohort": cohort, "sex": sex, "n": len(d), "mean": round(d.mean(), 2)}
        row.update({f"p{q}": round(float(np.percentile(d, q)), 1) for q in PERCENTILES})
        bands = pd.cut(d, CANDIDATE_BAND_EDGES, right=False)
        row.update({f"band_{iv.left:g}-{iv.right:g}": int(n) for iv, n in bands.value_counts(sort=False).items()})
        rows.append(row)
    return pd.DataFrame(rows)


def ppmi_like_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """PDBP PD cases meeting each candidate definition, by sex, with the PDBP controls they'd be scored against."""
    pdbp = df[(df["cohort"] == "PDBP") & df["in_analysis_set"]]
    rows = []
    for max_d in CANDIDATE_MAX_DURATION:
        for unmedicated in (True, False):
            row = {"max_duration_years": max_d, "unmedicated_required": unmedicated}
            for sex in SEXES:
                g = pdbp[(pdbp["sex"] == sex)]
                pd_g = g[g["diagnosis"] == "PD"]
                keep = pd_g["disease_duration_years"] <= max_d
                if unmedicated:
                    keep &= pd_g["pd_medication"].eq(False).fillna(False)
                row[f"{sex}_pd"] = int(keep.sum())
                row[f"{sex}_controls"] = int((g["diagnosis"] == "Control").sum())
            rows.append(row)
    return pd.DataFrame(rows)


def coverage_report(df: pd.DataFrame) -> dict:
    return {
        "reference_check": reference_check(df),
        "counts": counts(df),
        "excluded_diagnoses": excluded_diagnoses(df),
        "carriers": carriers(df),
        "carrier_sources": carrier_sources(df),
        "coverage_pd": coverage(df),
        "duration_distribution": duration_distribution(df),
        "ppmi_like_candidates": ppmi_like_candidates(df),
    }


SECTIONS = {
    "reference_check": "Counts against the MOFA write-up",
    "counts": "Analysis set by cohort, sex and diagnosis (RNA-seq participants)",
    "excluded_diagnoses": "Excluded: diagnosis other than PD or control (raw label)",
    "carriers": "Carrier status (analysis set)",
    "carrier_sources": "Carrier sources: WGS-confirmed mutation vs enrollment genetic status",
    "coverage_pd": "Coverage of duration, medication and clinical scores among PD cases",
    "duration_distribution": "Disease duration (years) among PD cases with duration",
    "ppmi_like_candidates": "PDBP PPMI-like subset under candidate definitions (PD cases; all PDBP controls)",
}


def _markdown(t: pd.DataFrame) -> str:
    """Pipe table without the optional tabulate dependency."""
    rows = [list(map(str, t.columns))] + [["" if pd.isna(v) else str(v) for v in r] for r in t.itertuples(index=False)]
    out = ["| " + " | ".join(rows[0]) + " |", "|" + "---|" * len(rows[0])]
    return "\n".join(out + ["| " + " | ".join(r) + " |" for r in rows[1:]])


def write_report(tables: dict, out_dir: str) -> str:
    os.makedirs(out_dir, exist_ok=True)
    lines = ["# Analysis set and field-coverage report (issue 03)", "",
             "Aggregate counts only. Written before any PDBP classifier score exists.", ""]
    for key, title in SECTIONS.items():
        t = tables[key]
        t.to_csv(os.path.join(out_dir, f"{key}.tsv"), sep="\t", index=False)
        lines += [f"## {title}", "", _markdown(t) if len(t) else "_none_", ""]
    path = os.path.join(out_dir, "coverage_report.md")
    with open(path, "w") as f:
        f.write("\n".join(lines))
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--participants", nargs="+", required=True, help="analysis_set.py outputs")
    p.add_argument("--output_dir", required=True)
    a = p.parse_args()
    path = write_report(coverage_report(load(a.participants)), a.output_dir)
    print(open(path).read())
    print(f"Saved: {path}")


if __name__ == "__main__":
    main()
