#!/usr/bin/env python3
"""
mofa/scripts/plasma_diagnostics.py

Diagnoses why plasma proteomics is barely explained by the PPMI MOFA model (~2.5% R2)
and does not replicate in PDBP. CSF proteomics is run alongside as a reference view.

Main hypothesis: exact zeros in the raw intensity matrices are non-detects. The current
pipeline turns them into log2(0 + 1) = 0, i.e. ~17 log2 units below a typical protein,
so after centering they become extreme outliers that can dominate a view's variance.

For each view (CSF, Plasma) the script reports:
  1. Raw scale: min / median / max, % exact zeros, % NaN, % negative (per cohort)
  2. Zero (non-detect) pattern: zeros per protein and per sample, and whether the
     per-sample zero fraction differs by plate
  3. Share of the processed view's total sum of squares that comes from zero entries
  4. PCA structure in PPMI under two treatments:
        "current"     : log2(x + 1) with zeros kept (as in the pipeline)
        "zero_as_nan" : zeros set to missing before log2 (mean-imputed for PCA only)
  5. Associations of the top PPMI PCs with plate, sex, age, diagnosis and MOFA factors
  6. Replication: how much PDBP variance the top-k PPMI PCs explain, compared with
     PDBP's own top-k PCs (upper bound)

Inputs are proteins x participant_id matrices before any log transform (the files
produced by align_matrix_columns.load_and_align_matrix).
"""

import os
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats
try:
    import matplotlib.pyplot as plt
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

TREATMENTS = ["current", "zero_as_nan"]
N_PCS = 10
N_ASSOC_PCS = 5


def load_metadata(path: str) -> pd.DataFrame:
    meta = pd.read_csv(path)
    id_col = "participant_id" if "participant_id" in meta.columns else meta.columns[0]
    return meta.drop_duplicates(subset=[id_col]).set_index(id_col)


def raw_scale_summary(df: pd.DataFrame) -> dict:
    vals = df.values.ravel()
    finite = vals[~np.isnan(vals)]
    return {
        "n_proteins": df.shape[0],
        "n_samples": df.shape[1],
        "min": np.min(finite),
        "median": np.median(finite),
        "max": np.max(finite),
        "pct_zero": 100 * np.mean(finite == 0),
        "pct_nan": 100 * np.mean(np.isnan(vals)),
        "pct_negative": 100 * np.mean(finite < 0),
    }


def process(df: pd.DataFrame, treatment: str) -> pd.DataFrame:
    """
    Mirrors proteomics_qc.normalize_proteomics (log2(x + 1), sample-median centering)
    followed by feature centering. With 'zero_as_nan', zeros are masked first.
    """
    x = df.astype(float)
    if treatment == "zero_as_nan":
        x = x.mask(x == 0)
    x = np.log2(x + 1)
    sample_medians = x.median(axis=0, skipna=True)
    x = x.sub(sample_medians - sample_medians.median(skipna=True), axis=1)
    return x.sub(x.mean(axis=1, skipna=True), axis=0)


def zero_share_of_variance(raw: pd.DataFrame, processed_current: pd.DataFrame) -> float:
    """Percent of the processed view's total sum of squares contributed by raw zero entries."""
    sq = processed_current.values ** 2
    total = np.nansum(sq)
    return 100 * np.nansum(np.where(raw.values == 0, sq, 0.0)) / total if total > 0 else np.nan


def pca(x: pd.DataFrame, n_pcs: int):
    """
    PCA on a centered proteins x samples matrix (NaN mean-imputed, i.e. set to 0).
    Returns (variance explained % per PC, loadings proteins x PCs, scores samples x PCs).
    """
    X = np.nan_to_num(x.values.T)  # samples x proteins
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    var_pct = 100 * s ** 2 / np.sum(s ** 2)
    k = min(n_pcs, len(s))
    pcs = [f"PC{i+1}" for i in range(k)]
    loadings = pd.DataFrame(Vt[:k].T, index=x.index, columns=pcs)
    scores = pd.DataFrame(U[:, :k] * s[:k], index=x.columns, columns=pcs)
    return var_pct[:k], loadings, scores


def variance_captured(x: pd.DataFrame, loadings: pd.DataFrame) -> float:
    """Percent of x's total variance reconstructed by projecting onto the given (orthonormal) loadings."""
    X = np.nan_to_num(x.loc[loadings.index].values.T)
    V = loadings.values
    resid = X - X @ V @ V.T
    total = np.sum(X ** 2)
    return 100 * (1 - np.sum(resid ** 2) / total) if total > 0 else np.nan


def associate(scores: pd.DataFrame, covariates: pd.DataFrame) -> pd.DataFrame:
    """-log10 p of each PC against each covariate (Pearson for numeric, ANOVA for categorical)."""
    records = []
    for pc in scores.columns:
        for cov in covariates.columns:
            sub = pd.concat([scores[pc], covariates[cov]], axis=1, join="inner").dropna()
            if len(sub) < 10:
                continue
            y, x = sub.iloc[:, 0], sub.iloc[:, 1]
            if pd.api.types.is_numeric_dtype(x):
                _, p = stats.pearsonr(x, y)
            else:
                groups = [y[x == c] for c in x.unique() if (x == c).sum() > 1]
                if len(groups) < 2:
                    continue
                _, p = stats.f_oneway(*groups)
            records.append({"PC": pc, "covariate": cov, "neg_log10_p": -np.log10(max(p, 1e-300))})
    return pd.DataFrame(records)


def run_view(view: str, ppmi_raw: pd.DataFrame, pdbp_raw: pd.DataFrame, ppmi_meta: pd.DataFrame,
             ppmi_z: pd.DataFrame, output_dir: str) -> dict:
    print(f"\n{'=' * 70}\n{view}\n{'=' * 70}")
    common = ppmi_raw.index.intersection(pdbp_raw.index)
    print(f"Proteins: PPMI {len(ppmi_raw)}, PDBP {len(pdbp_raw)}, shared {len(common)}")
    ppmi_raw, pdbp_raw = ppmi_raw.loc[common], pdbp_raw.loc[common]

    # 1. Raw scale
    print("\n[1] Raw scale (before log2)")
    scale = pd.DataFrame({"PPMI": raw_scale_summary(ppmi_raw), "PDBP": raw_scale_summary(pdbp_raw)}).T
    print(scale.round(3).to_string())

    # 2. Zero pattern
    print("\n[2] Zero (non-detect) pattern")
    zero_rows = []
    for cohort, raw in (("PPMI", ppmi_raw), ("PDBP", pdbp_raw)):
        per_protein = 100 * (raw == 0).mean(axis=1)
        per_sample = 100 * (raw == 0).mean(axis=0)
        zero_rows.append({
            "cohort": cohort,
            "proteins_with_any_zero": int((per_protein > 0).sum()),
            "proteins_>30pct_zero": int((per_protein > 30).sum()),
            "proteins_>50pct_zero": int((per_protein > 50).sum()),
            "median_pct_zero_per_sample": per_sample.median(),
            "max_pct_zero_per_sample": per_sample.max(),
        })
        per_protein.rename("pct_zero").to_csv(os.path.join(output_dir, f"{view.lower()}_{cohort.lower()}_zeros_per_protein.csv"))
    print(pd.DataFrame(zero_rows).set_index("cohort").round(2).to_string())

    if "plate" in ppmi_meta.columns:
        per_sample = (100 * (ppmi_raw == 0).mean(axis=0)).rename("pct_zero")
        assoc = associate(per_sample.to_frame(), ppmi_meta[["plate"]])
        if not assoc.empty:
            print(f"PPMI per-sample zero fraction vs plate: -log10 p = {assoc['neg_log10_p'].iloc[0]:.1f}")

    # 3. Zero share of variance under the current pipeline
    ppmi_current = process(ppmi_raw, "current")
    pdbp_current = process(pdbp_raw, "current")
    print("\n[3] Share of processed-view sum of squares coming from zero entries (current pipeline)")
    zshare = {"PPMI": zero_share_of_variance(ppmi_raw, ppmi_current), "PDBP": zero_share_of_variance(pdbp_raw, pdbp_current)}
    for cohort, val in zshare.items():
        print(f"  {cohort}: {val:.1f}%")

    # 4-6. PCA structure, associations, replication per treatment
    summary = {"view": view, "zero_share_PPMI": zshare["PPMI"], "zero_share_PDBP": zshare["PDBP"]}
    scree = {}
    for treatment in TREATMENTS:
        ppmi_x = ppmi_current if treatment == "current" else process(ppmi_raw, treatment)
        pdbp_x = pdbp_current if treatment == "current" else process(pdbp_raw, treatment)
        var_ppmi, load_ppmi, scores_ppmi = pca(ppmi_x, N_PCS)
        var_pdbp, _, _ = pca(pdbp_x, N_PCS)
        scree[treatment] = (var_ppmi, var_pdbp)

        print(f"\n[4] PCA, treatment = {treatment}")
        print("  PPMI variance % PC1-5: " + ", ".join(f"{v:.1f}" for v in var_ppmi[:5]))
        print("  PDBP variance % PC1-5: " + ", ".join(f"{v:.1f}" for v in var_pdbp[:5]))

        top_load = load_ppmi["PC1"].abs().sort_values(ascending=False).head(5)
        zero_pct = 100 * (ppmi_raw.loc[top_load.index] == 0).mean(axis=1)
        print("  Top PC1 proteins (PPMI) and their % zeros: " +
              ", ".join(f"{p} ({z:.0f}%)" for p, z in zero_pct.items()))

        covs = [c for c in ["plate", "sex", "age_at_baseline", "diagnosis"] if c in ppmi_meta.columns]
        cov_df = ppmi_meta[covs].copy()
        if ppmi_z is not None:
            cov_df = cov_df.join(ppmi_z, how="left")
        assoc = associate(scores_ppmi.iloc[:, :N_ASSOC_PCS], cov_df)
        if not assoc.empty:
            pivot = assoc.pivot(index="PC", columns="covariate", values="neg_log10_p").reindex(scores_ppmi.columns[:N_ASSOC_PCS])
            pivot.to_csv(os.path.join(output_dir, f"{view.lower()}_{treatment}_ppmi_pc_associations.csv"))
            print(f"\n[5] PPMI PC associations (-log10 p), treatment = {treatment}")
            print(pivot.round(1).to_string())

        print(f"\n[6] Replication, treatment = {treatment}")
        for k in (1, 3, 5):
            ppmi_in_pdbp = variance_captured(pdbp_x, load_ppmi.iloc[:, :k])
            pdbp_own = float(np.sum(var_pdbp[:k]))
            ratio = ppmi_in_pdbp / pdbp_own if pdbp_own > 0 else np.nan
            print(f"  top {k} PCs: PPMI PCs explain {ppmi_in_pdbp:.1f}% of PDBP variance "
                  f"(PDBP's own top {k}: {pdbp_own:.1f}%; ratio {ratio:.2f})")
            summary[f"{treatment}_PPMI_PC1to{k}_var"] = float(np.sum(var_ppmi[:k]))
            summary[f"{treatment}_replication_ratio_k{k}"] = ratio

    if HAS_PLOT:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        for ax, treatment in zip(axes, TREATMENTS):
            var_ppmi, var_pdbp = scree[treatment]
            pcs = np.arange(1, len(var_ppmi) + 1)
            ax.plot(pcs, var_ppmi, marker="o", label="PPMI")
            ax.plot(np.arange(1, len(var_pdbp) + 1), var_pdbp, marker="s", label="PDBP")
            ax.set_title(f"{view}: {treatment}")
            ax.set_xlabel("Principal component")
            ax.set_xticks(pcs)
        axes[0].set_ylabel("Variance explained (%)")
        axes[0].legend()
        plt.tight_layout()
        out_plot = os.path.join(output_dir, f"{view.lower()}_pca_scree.png")
        plt.savefig(out_plot, dpi=200)
        plt.close()
        print(f"\nSaved scree plot to: {out_plot}")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Diagnose plasma (and CSF) proteomics signal and replication.")
    parser.add_argument("--ppmi_plasma", required=True, help="PPMI raw plasma matrix (proteins x participant_id).")
    parser.add_argument("--pdbp_plasma", required=True, help="PDBP raw plasma matrix (proteins x participant_id).")
    parser.add_argument("--ppmi_csf", default=None, help="PPMI raw CSF matrix (optional reference view).")
    parser.add_argument("--pdbp_csf", default=None, help="PDBP raw CSF matrix (optional reference view).")
    parser.add_argument("--ppmi_metadata", required=True, help="PPMI metadata CSV (participant_id, plate, sex, age_at_baseline, diagnosis).")
    parser.add_argument("--ppmi_z", default=None, help="PPMI latent factors CSV (participant_id + Factor columns) to correlate with PCs.")
    parser.add_argument("--output", default="plasma_diagnostics", help="Output directory.")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    ppmi_meta = load_metadata(args.ppmi_metadata)

    ppmi_z = None
    if args.ppmi_z and os.path.exists(args.ppmi_z):
        z = pd.read_csv(args.ppmi_z)
        id_col = "participant_id" if "participant_id" in z.columns else z.columns[0]
        z = z.drop_duplicates(subset=[id_col]).set_index(id_col)
        ppmi_z = z[[c for c in z.columns if c.startswith("Factor")]]

    views = {"Plasma": (args.ppmi_plasma, args.pdbp_plasma)}
    if args.ppmi_csf and args.pdbp_csf:
        views["CSF"] = (args.ppmi_csf, args.pdbp_csf)

    summaries = []
    for view, (ppmi_path, pdbp_path) in views.items():
        ppmi_raw = pd.read_csv(ppmi_path, index_col=0)
        pdbp_raw = pd.read_csv(pdbp_path, index_col=0)
        summaries.append(run_view(view, ppmi_raw, pdbp_raw, ppmi_meta, ppmi_z, args.output))

    out_csv = os.path.join(args.output, "proteomics_diagnostics_summary.csv")
    pd.DataFrame(summaries).to_csv(out_csv, index=False)
    print(f"\nSaved summary to: {out_csv}")
    print("\n=== Diagnostics Complete ===")


if __name__ == "__main__":
    main()


###### Example use (from ~/mofa on Terra) ##########
#python3 scripts/plasma_diagnostics.py \
#  --ppmi_plasma data/ppmi/raw/ppmi_pla_prot.csv \
#  --pdbp_plasma data/pdbp/raw/plasma_proteomics.csv \
#  --ppmi_csf data/ppmi/raw/ppmi_csf_prot.csv \
#  --pdbp_csf data/pdbp/raw/csf_proteomics.csv \
#  --ppmi_metadata data/metadata_m0_PPMI.csv \
#  --ppmi_z ppmi_mofa_outputs/latent_factors_Z_with_meta.csv \
#  --output plasma_diagnostics
