#!/usr/bin/env python3
"""
mofa/scripts/clinical_models.py

Regression models linking MOFA factors to clinical variables in PPMI (discovery) and PDBP
(validation), including the sex x diagnosis interaction tests (pdbp-projection issue 04).

Factors are z-scored within each cohort, and continuous clinical scores are z-scored within the
PD cases of each cohort, so estimates are per SD and comparable across cohorts. Metadata is
reduced to one row per participant before modelling (evaluate_factors.dedupe_metadata).

Analyses (per factor):
  interaction_dx       factor ~ pd * male + age                         (cases + controls)
  within_sex_dx        factor ~ pd + age, fitted separately in each sex  (pd effect per sex)
  interaction_<score>  factor ~ score * male + age                       (PD cases only)
  within_sex_<score>   factor ~ score + age, fitted separately in each sex (PD cases only)
  main_dx              logit(pd) ~ factor + age + male                   (log-odds per SD; OR = exp(estimate))
  main_<score>         score ~ factor + age + male                       (PD cases only)

Severity models use PD cases only so that disease severity isn't confounded with diagnosis.
p-values are Benjamini-Hochberg corrected across factors within each cohort, analysis and term.
A result counts as replicated when PPMI q < 0.05, PDBP p < 0.05 and the estimates share a sign.

The diagnosis column's labels differ between data sources, so --case_levels and --control_levels
are required. Run without them to print the labels found in each cohort.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

from evaluate_factors import dedupe_metadata

try:
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests
except ImportError:
    sys.exit("clinical_models.py needs statsmodels: pip install statsmodels")

try:
    import matplotlib.pyplot as plt
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False


def load_cohort(z_path: str, metadata_path: str) -> pd.DataFrame:
    """Factor scores (Factor columns only) merged with one metadata row per participant."""
    z = pd.read_csv(z_path)
    id_col = next((c for c in ["participant_id", "patient_id", "sample", "Unnamed: 0"] if c in z.columns), z.columns[0])
    z = z.rename(columns={id_col: "participant_id"})
    z = z[["participant_id"] + [c for c in z.columns if c.startswith("Factor")]]
    z["participant_id"] = z["participant_id"].astype(str)
    z = dedupe_metadata(z)

    meta = pd.read_csv(metadata_path)
    meta_id = "participant_id" if "participant_id" in meta.columns else meta.columns[0]
    meta = meta.rename(columns={meta_id: "participant_id"})
    meta["participant_id"] = meta["participant_id"].astype(str)
    meta = dedupe_metadata(meta)
    meta = meta[[c for c in meta.columns if not c.startswith("Factor")]]
    return z.merge(meta, on="participant_id", how="inner")


def prepare(df: pd.DataFrame, factors: list, args) -> pd.DataFrame:
    """Adds pd (1 case / 0 control / NaN other), male (1/0), z-scored factors and scores."""
    out = df.copy()
    dx = out[args.diagnosis_col].astype(str)
    out["pd"] = np.where(dx.isin(args.case_levels), 1.0, np.where(dx.isin(args.control_levels), 0.0, np.nan))
    sex = out["sex"].astype(str).str.lower().str[0]
    out["male"] = np.where(sex == "m", 1.0, np.where(sex == "f", 0.0, np.nan))
    out["age"] = pd.to_numeric(out[args.age_col], errors="coerce")
    for f in factors:
        out[f] = (out[f] - out[f].mean()) / out[f].std()
    cases = out["pd"] == 1
    for name, col in args.scores.items():
        if col in out.columns:
            v = pd.to_numeric(out[col], errors="coerce")
            out[name] = (v - v[cases].mean()) / v[cases].std()
    return out


def fit(formula: str, data: pd.DataFrame, term: str, variables: list, logistic: bool = False) -> dict:
    """Fits on complete cases for the given variables and returns the estimate for one term."""
    d = data.dropna(subset=variables)
    if len(d) < 20:
        return {"estimate": np.nan, "se": np.nan, "p": np.nan, "n": len(d)}
    try:
        model = (smf.logit(formula, d).fit(disp=0) if logistic else smf.ols(formula, d).fit())
        return {"estimate": model.params[term], "se": model.bse[term], "p": model.pvalues[term], "n": int(model.nobs)}
    except Exception as e:  # e.g. perfect separation or a constant term
        print(f"  Could not fit '{formula}' ({term}): {e}")
        return {"estimate": np.nan, "se": np.nan, "p": np.nan, "n": len(d)}


def run_cohort(df: pd.DataFrame, cohort: str, factors: list, args) -> pd.DataFrame:
    rows = []
    def add(analysis, factor, term, res, sex="all"):
        rows.append({"cohort": cohort, "analysis": analysis, "factor": factor, "term": term, "sex": sex, **res})

    dx_data = df[df["pd"].notna()]
    pd_cases = df[df["pd"] == 1]
    sexes = {"Female": 0.0, "Male": 1.0}
    for f in factors:
        add("interaction_dx", f, "pd:male", fit(f"{f} ~ pd * male + age", dx_data, "pd:male", [f, "pd", "male", "age"]))
        for sex, code in sexes.items():
            add("within_sex_dx", f, "pd", fit(f"{f} ~ pd + age", dx_data[dx_data["male"] == code], "pd", [f, "pd", "age"]), sex)
        add("main_dx", f, f, fit(f"pd ~ {f} + age + male", dx_data, f, [f, "pd", "male", "age"], logistic=True))
        for score in args.scores:
            if score not in df.columns:
                continue
            add(f"interaction_{score}", f, f"{score}:male", fit(f"{f} ~ {score} * male + age", pd_cases, f"{score}:male", [f, score, "male", "age"]))
            for sex, code in sexes.items():
                add(f"within_sex_{score}", f, score, fit(f"{f} ~ {score} + age", pd_cases[pd_cases["male"] == code], score, [f, score, "age"]), sex)
            add(f"main_{score}", f, f, fit(f"{score} ~ {f} + age + male", pd_cases, f, [f, score, "male", "age"]))

    res = pd.DataFrame(rows)
    res["q"] = np.nan
    for _, idx in res.groupby(["analysis", "term", "sex"]).groups.items():
        p = res.loc[idx, "p"]
        ok = p.notna()
        if ok.any():
            res.loc[p[ok].index, "q"] = multipletests(p[ok], method="fdr_bh")[1]
    return res


def replication_table(res: pd.DataFrame) -> pd.DataFrame:
    keys = ["analysis", "factor", "term", "sex"]
    ppmi = res[res.cohort == "PPMI"].set_index(keys)[["estimate", "p", "q", "n"]].add_prefix("PPMI_")
    pdbp = res[res.cohort == "PDBP"].set_index(keys)[["estimate", "p", "q", "n"]].add_prefix("PDBP_")
    rep = ppmi.join(pdbp, how="outer").reset_index()
    rep["replicated"] = ((rep["PPMI_q"] < 0.05) & (rep["PDBP_p"] < 0.05) &
                         (np.sign(rep["PPMI_estimate"]) == np.sign(rep["PDBP_estimate"])))
    return rep


def plot_within_sex_dx(res: pd.DataFrame, factors: list, out_path: str):
    sub = res[res.analysis == "within_sex_dx"]
    cohorts = [c for c in ["PPMI", "PDBP"] if c in sub.cohort.unique()]
    fig, axes = plt.subplots(1, len(cohorts), figsize=(5 * len(cohorts), 0.45 * len(factors) + 1.5), sharey=True, squeeze=False)
    y = np.arange(len(factors))
    for ax, cohort in zip(axes[0], cohorts):
        for offset, (sex, color) in zip((-0.15, 0.15), (("Female", "tab:red"), ("Male", "tab:blue"))):
            d = sub[(sub.cohort == cohort) & (sub.sex == sex)].set_index("factor").reindex(factors)
            ax.errorbar(d["estimate"], y + offset, xerr=1.96 * d["se"], fmt="o", color=color, label=sex, capsize=2)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_title(f"{cohort}: PD vs control, within sex")
        ax.set_xlabel("PD effect on factor (SD, age-adjusted, 95% CI)")
    axes[0][0].set_yticks(y)
    axes[0][0].set_yticklabels(factors)
    axes[0][0].invert_yaxis()
    axes[0][0].legend()
    plt.tight_layout(w_pad=3)
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f"Saved within-sex diagnosis plot to: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Clinical regression and sex x diagnosis interaction models for MOFA factors.")
    parser.add_argument("--ppmi_z", required=True, help="PPMI factor scores CSV (e.g. latent_factors_Z_with_meta.csv).")
    parser.add_argument("--ppmi_metadata", required=True, help="PPMI clinical metadata CSV.")
    parser.add_argument("--pdbp_z", default=None, help="PDBP projected factor scores CSV (optional).")
    parser.add_argument("--pdbp_metadata", default=None, help="PDBP clinical metadata CSV (required with --pdbp_z).")
    parser.add_argument("--output", default="clinical_models", help="Output directory.")
    parser.add_argument("--factors", nargs="+", default=None, help="Factors to model (default: all).")
    parser.add_argument("--diagnosis_col", default="diagnosis", help="Diagnosis column (default: diagnosis).")
    parser.add_argument("--case_levels", nargs="+", default=None, help="Diagnosis labels counted as PD.")
    parser.add_argument("--control_levels", nargs="+", default=None, help="Diagnosis labels counted as controls.")
    parser.add_argument("--age_col", default="age_at_baseline")
    parser.add_argument("--updrs_col", default="mds_updrs_part_iii_summary_score")
    parser.add_argument("--moca_col", default="moca_total_score")
    args = parser.parse_args()
    args.scores = {"updrs3": args.updrs_col, "moca": args.moca_col}

    cohorts = {"PPMI": load_cohort(args.ppmi_z, args.ppmi_metadata)}
    if args.pdbp_z:
        if not args.pdbp_metadata:
            sys.exit("--pdbp_metadata is required with --pdbp_z")
        cohorts["PDBP"] = load_cohort(args.pdbp_z, args.pdbp_metadata)

    for name, df in cohorts.items():
        print(f"\n{name}: {len(df)} participants with factor scores and metadata")
        print(f"  {args.diagnosis_col} labels: " + ", ".join(f"{k} ({v})" for k, v in df[args.diagnosis_col].value_counts(dropna=False).items()))
        print(f"  sex: " + ", ".join(f"{k} ({v})" for k, v in df["sex"].value_counts(dropna=False).items()))

    if not args.case_levels or not args.control_levels:
        print("\nPass --case_levels and --control_levels (exact labels from the lists above) to run the models.")
        return

    os.makedirs(args.output, exist_ok=True)
    all_res = []
    for name, df in cohorts.items():
        factors = args.factors or [c for c in df.columns if c.startswith("Factor")]
        data = prepare(df, factors, args)
        n_case, n_ctrl = int((data.pd == 1).sum()), int((data.pd == 0).sum())
        dropped = int(data.pd.isna().sum())
        print(f"\n{name}: {n_case} PD, {n_ctrl} controls ({dropped} with other or missing diagnosis excluded from diagnosis models)")
        print("  PD / control by sex: " + ", ".join(
            f"{s} {int(((data.pd == 1) & (data.male == c)).sum())} / {int(((data.pd == 0) & (data.male == c)).sum())}"
            for s, c in (("F", 0.0), ("M", 1.0))))
        all_res.append(run_cohort(data, name, factors, args))

    res = pd.concat(all_res, ignore_index=True)
    res.to_csv(os.path.join(args.output, "clinical_models_long.csv"), index=False)

    if "PDBP" in cohorts:
        rep = replication_table(res)
        rep.to_csv(os.path.join(args.output, "clinical_models_replication.csv"), index=False)
        inter = rep[rep.analysis.str.startswith("interaction")]
        print("\n=== Sex x clinical interaction terms (PPMI / PDBP) ===")
        cols = ["analysis", "factor", "PPMI_estimate", "PPMI_q", "PDBP_estimate", "PDBP_p", "replicated"]
        print(inter[cols].round(4).to_string(index=False))
        hits = rep[rep.replicated]
        print(f"\n=== Replicated results (PPMI q < 0.05, PDBP p < 0.05, same sign): {len(hits)} ===")
        if len(hits):
            print(hits[["analysis", "factor", "term", "sex", "PPMI_estimate", "PPMI_q", "PDBP_estimate", "PDBP_p"]].round(4).to_string(index=False))

    if HAS_PLOT:
        factors = args.factors or sorted({f for f in res.factor}, key=lambda f: int(f.replace("Factor", "")))
        plot_within_sex_dx(res, factors, os.path.join(args.output, "within_sex_diagnosis_effects.png"))
    print(f"\nSaved results to: {args.output}")


if __name__ == "__main__":
    main()


###### Example use (from ~/mofa on Terra) ##########
# 1. List the diagnosis labels in each cohort:
#python3 scripts/clinical_models.py \
#  --ppmi_z ppmi_mofa_outputs_noplasma/latent_factors_Z_with_meta.csv --ppmi_metadata data/metadata_m0_PPMI.csv \
#  --pdbp_z pdbp_mofa_outputs_noplasma/pdbp_latent_factors_Z.csv --pdbp_metadata data/metadata_m0_PDBP.csv
# 2. Run the models with the labels chosen from that list:
#python3 scripts/clinical_models.py ... --case_levels "<PD label>" --control_levels "<control label>" \
#  --output ppmi_mofa_outputs_noplasma/clinical_models
