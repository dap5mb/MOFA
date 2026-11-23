"""
Figures for the MOFA sex-differences write-up (docs/mofa_writeup.md), drawn from the CSVs in
the local results folder. Each figure is skipped, with a message, if its inputs are missing.

  fig2_variance_explained   R² per factor, view, sex and cohort (main model)
  fig3_associations         unadjusted factor-covariate associations, PPMI and PDBP
  fig4_within_sex_pd        PD effect within each sex (age-adjusted), with the interaction p
  fig5_enrichment           GSEA preranked results for the PD and aging factors
  figS1_stability           subsampling stability of the main model's factors

Usage:
  python3 mofa/scripts/make_figures.py --results results/mofa --output figures/mofa/writeup
"""

import argparse
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# Reference palette (light mode): categorical slots 1-2, text and surface tokens
SEX_COLORS = {"Female": "#eb6834", "Male": "#2a78d6"}
POS, NEG = "#eb6834", "#2a78d6"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f4f7fc", "#a9c7ee", "#2a78d6", "#0d3b75"])
FACTORS = [f"Factor{i}" for i in range(1, 11)]
FLABEL = {f: f"F{f[6:]}" for f in FACTORS}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42,
    "savefig.dpi": 300, "figure.facecolor": "white",
})


def save(fig, out_dir, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {name}.pdf / .png")


def need(*paths):
    missing = [p for p in paths if not os.path.exists(p)]
    for p in missing:
        print(f"  missing input: {p}")
    return not missing


def fig_design(out_dir):
    """Study design schematic. Counts are those in docs/mofa_writeup.md (Methods: Cohorts and data)."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7.2)
    ax.set_axis_off()
    fill, model_fill, edge = "#f2f1ed", "#e6eefa", "#b9b8b2"

    def box(x, y, w, h, title, body, face=fill):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=face, ec=edge, lw=0.6))
        ax.text(x + 0.15, y + h - 0.14, title, fontsize=7.8, fontweight="bold", va="top", color=INK)
        ax.text(x + 0.15, y + h - 0.44, body, fontsize=6.6, va="top", color=INK2, linespacing=1.35)

    def arrow(x1, y1, x2, y2, label=None, lx=0, ly=0):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9, color=INK2, lw=0.8))
        if label:
            ax.text((x1 + x2) / 2 + lx, (y1 + y2) / 2 + ly, label, fontsize=6.4, color=INK2, ha="center", va="center")

    # Cohorts
    box(0.1, 5.75, 4.3, 1.35, "PPMI: discovery cohort (baseline)",
        "1,329 participants: 613 women, 716 men\n753 PD, 555 controls, 21 other diagnoses\n"
        "Whole-blood RNA-seq 1,308; CSF proteomics 354")
    box(5.6, 5.75, 4.3, 1.35, "PDBP: validation cohort (baseline)",
        "1,248 participants: 554 women, 694 men\n767 PD, 481 controls\nWhole-blood RNA-seq 1,248; CSF proteomics 107")
    # Processing
    box(0.1, 4.05, 4.3, 1.35, "Processing (within PPMI)",
        "RNA-seq: DESeq2 size factors, log, 5,000 most\nvariable genes; sequencing plate removed\n"
        "CSF: non-detects set to missing, log2, median-\ncentered; 249 proteins (plasma excluded)")
    box(5.6, 4.05, 4.3, 1.35, "Processing (within PDBP)",
        "Same steps, estimated within PDBP\nGenes and proteins fixed to the PPMI model\nNo PDBP data used for feature selection")
    # Model and projection
    box(0.1, 2.35, 4.3, 1.35, "Multi-group MOFA (mofapy2)",
        "2 views, 10 factors, groups = sex\nFeatures centered within each sex; views scaled\nShared weights; per-sex factor activity",
        face=model_fill)
    box(5.6, 2.35, 4.3, 1.35, "Projection into PDBP",
        "PPMI weights held fixed\nCentered within sex, scaled as in training\nObserved features only, noise-weighted",
        face=model_fill)
    # Analyses
    box(0.1, 0.1, 9.8, 1.9, "Analyses (both cohorts unless stated)",
        "Variance explained per factor, view and sex (same formula in both cohorts)\n"
        "Stability: retraining on 10 subsets of 80% of PPMI (|r| ≥ 0.8 in ≥ 9 of 10)\n"
        "Clinical models: factor ~ PD × sex + age; MDS-UPDRS III and MoCA × sex in PD; PD within each sex\n"
        "Replicated if PPMI q < 0.05 (BH across factors), PDBP p < 0.05 and same sign\n"
        "Sensitivity: single-group model (no sex grouping); three-view model with plasma; technical check (library size)\n"
        "Enrichment (PPMI weights): GSEA preranked and STRING networks for the PD and age factors")
    for x in (2.25, 7.75):
        for y in (5.75, 4.05, 2.35):
            arrow(x, y, x, y - 0.33)
    arrow(4.44, 3.0, 5.56, 3.0, "weights,\nnoise, prior", ly=0.3)
    arrow(4.44, 4.7, 5.56, 4.7, "features", ly=0.18)
    save(fig, out_dir, "fig1_study_design")


def fig_variance(results, out_dir):
    path = os.path.join(results, "main_noplasma", "cohort_variance_replication_comparison.csv")
    if not need(path):
        return print("Skipped fig2_variance_explained")
    df = pd.read_csv(path)
    floor = 0.01
    views = [("RNA_seq", "Whole-blood RNA-seq"), ("CSF_Proteomics", "CSF proteomics")]
    offsets = {("PPMI", "Female"): -0.24, ("PDBP", "Female"): -0.08, ("PPMI", "Male"): 0.08, ("PDBP", "Male"): 0.24}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
    for ax, (view, title) in zip(axes, views):
        d = df[df.view == view]
        for (cohort, sex), off in offsets.items():
            s = d[d.group == sex].set_index("Factor").reindex(FACTORS)[f"{cohort}_R2"].clip(lower=floor)
            x = np.arange(len(FACTORS)) + off
            ax.scatter(x, s, s=22, marker="o", linewidths=1.2, zorder=3, edgecolors=SEX_COLORS[sex],
                       facecolors=SEX_COLORS[sex] if cohort == "PPMI" else "white")
        ax.set_yscale("log")
        ax.set_ylim(floor * 0.7, 100)
        ax.set_xticks(range(len(FACTORS)), [FLABEL[f] for f in FACTORS])
        ax.set_title(title, loc="left")
        ax.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
        ax.set_xlim(-0.6, len(FACTORS) - 0.4)
    axes[0].set_ylabel("Variance explained (%)")
    axes[0].set_yticks([0.01, 0.1, 1, 10, 100], ["≤0.01", "0.1", "1", "10", "100"])
    handles = [Line2D([], [], marker="o", ls="", ms=5, mfc=SEX_COLORS[s] if c == "PPMI" else "white",
                      mec=SEX_COLORS[s], mew=1.2, label=f"{c}, {'women' if s == 'Female' else 'men'}")
               for c, s in offsets]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.06))
    fig.tight_layout()
    save(fig, out_dir, "fig2_variance_explained")


def fig_associations(results, out_dir):
    paths = {c: os.path.join(results, "main_noplasma", f"factor_metadata_associations_{c.lower()}.csv") for c in ("PPMI", "PDBP")}
    if not need(*paths.values()):
        return print("Skipped fig3_associations")
    covs = [("age_at_baseline", "Age"), ("diagnosis", "Diagnosis"), ("mds_updrs_part_iii_summary_score", "MDS-UPDRS III"),
            ("moca_total_score", "MoCA"), ("sex", "Sex"), ("plate", "Plate")]
    cap = 30
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.4))
    for ax, (cohort, path) in zip(axes, paths.items()):
        a = pd.read_csv(path)
        a["mlp"] = -np.log10(a["P-value"].clip(lower=1e-300))
        m = a.pivot_table(index="Factor", columns="Metadata", values="mlp").reindex(index=FACTORS, columns=[c for c, _ in covs])
        im = ax.imshow(m.clip(upper=cap).values, cmap=SEQ, vmin=0, vmax=cap, aspect="auto")
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                v = m.values[i, j]
                if np.isfinite(v) and v >= 2:
                    ax.text(j, i, f"{v:.0f}" if v >= 10 else f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                            color="white" if min(v, cap) > cap * 0.45 else INK)
        ax.set_xticks(range(len(covs)), [l for _, l in covs], rotation=40, ha="right")
        ax.set_yticks(range(len(FACTORS)), [FLABEL[f] for f in FACTORS])
        ax.set_title(cohort, loc="left")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    cb = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.02)
    cb.set_label(f"−log10 p (capped at {cap})")
    cb.outline.set_visible(False)
    save(fig, out_dir, "fig3_associations")


def fig_within_sex(results, out_dir):
    path = os.path.join(results, "main_noplasma", "clinical_models_long.csv")
    if not need(path):
        return print("Skipped fig4_within_sex_pd")
    df = pd.read_csv(path)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True)
    for ax, cohort in zip(axes, ("PPMI", "PDBP")):
        d = df[df.cohort == cohort]
        for k, f in enumerate(FACTORS):
            y = len(FACTORS) - 1 - k
            for sex, off in (("Female", 0.15), ("Male", -0.15)):
                r = d[(d.analysis == "within_sex_dx") & (d.factor == f) & (d.sex == sex)]
                if r.empty:
                    continue
                est, se = r.estimate.iloc[0], r.se.iloc[0]
                ax.plot([est - 1.96 * se, est + 1.96 * se], [y + off] * 2, color=SEX_COLORS[sex], lw=1.4,
                        solid_capstyle="round", zorder=2)
                ax.scatter(est, y + off, s=20, color=SEX_COLORS[sex], edgecolors="white", linewidths=1, zorder=3)
            ia = d[(d.analysis == "interaction_dx") & (d.factor == f)]
            if not ia.empty:
                ax.text(1.02, y, f"p = {ia.p.iloc[0]:.2f}", transform=ax.get_yaxis_transform(), va="center",
                        fontsize=6.5, color=INK2)
        ax.axvline(0, color=INK2, lw=0.6, zorder=1)
        ax.set_yticks(range(len(FACTORS)), [FLABEL[f] for f in reversed(FACTORS)])
        ax.set_xlabel("PD − control (SD of factor score)")
        ax.set_title(cohort, loc="left")
        ax.text(1.02, len(FACTORS) - 0.3, "PD × sex", transform=ax.get_yaxis_transform(), fontsize=6.5, color=INK2)
        ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    handles = [Line2D([], [], color=c, marker="o", ms=4, lw=1.4, label=("Women" if s == "Female" else "Men") + " (95% CI)")
               for s, c in SEX_COLORS.items()]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.04))
    fig.tight_layout(w_pad=4)
    save(fig, out_dir, "fig4_within_sex_pd")


ACRONYMS = {"BCR", "FCGR", "FCGR3A", "FCERI", "IL10", "IL1", "IL2", "IL6", "STAT5", "STAT3", "JAK", "NF", "KB", "CD22",
             "LAT2", "NTAL", "LAB", "NMD", "SRP", "P53", "KRAS", "TNFA", "MTORC1", "MYC", "E2F", "G2M", "UV", "DNA", "RNA",
             "MHC", "II", "TCR", "NK", "ER", "ATP", "GTP", "TGF", "WNT"}


def clean_term(term, width=50):
    lib = term.split("_")[0]
    words = re.sub(r"^(HALLMARK|REACTOME|GOBP)_", "", term).split("_")
    words = [w if w in ACRONYMS else w.lower() for w in words]
    words = ["B" if w == "b" and i + 1 < len(words) and words[i + 1] == "cell" else w for i, w in enumerate(words)]
    t = " ".join(words)
    t = t[0].upper() + t[1:]
    t = t.replace("NF KB", "NF-κB").replace("Gram negative", "Gram-negative").replace("gram negative", "Gram-negative")
    tag = {"HALLMARK": "H", "REACTOME": "R", "GOBP": "GO"}.get(lib, "")
    t = t if len(t) <= width else t[: width - 1] + "…"
    return f"{t} [{tag}]"


def pick_terms(d, per_library=2, per_direction=6):
    """Up to `per_library` terms per library and direction (largest |NES|), at most `per_direction` per direction."""
    picks = []
    for sign in (1, -1):
        s = d[np.sign(d.nes) == sign].assign(abs_nes=lambda x: x.nes.abs())
        top = s.sort_values("abs_nes", ascending=False).groupby("library").head(per_library)
        picks.append(top.nlargest(per_direction, "abs_nes"))
    return pd.concat(picks).sort_values("nes")


def fig_enrichment(results, out_dir, factors=("Factor7", "Factor10"), fdr=0.05):
    path = os.path.join(results, "main_noplasma", "enrichment", "gsea_prerank_results.csv")
    if not need(path):
        return print("Skipped fig5_enrichment")
    df = pd.read_csv(path)
    titles = {"Factor7": "Factor 7 (lower in PD)", "Factor10": "Factor 10 (higher with age)"}
    chosen = {f: pick_terms(df[(df.factor == f) & (df.fdr < fdr)].drop_duplicates("term")) for f in factors}
    heights = [max(len(c), 1) for c in chosen.values()]
    fig, axes = plt.subplots(len(factors), 1, figsize=(5.6, 0.22 * sum(heights) + 1.4), sharex=True,
                             gridspec_kw={"height_ratios": heights, "hspace": 0.35})
    for ax, f in zip(np.atleast_1d(axes), factors):
        top = chosen[f]
        ax.set_title(titles.get(f, f), loc="left")
        ax.axvline(0, color=INK2, lw=0.6)
        ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
        if top.empty:
            ax.text(0.5, 0.5, f"No gene set with FDR < {fdr}", ha="center", va="center", transform=ax.transAxes, color=INK2)
            ax.set_yticks([])
            continue
        y = np.arange(len(top))
        ax.barh(y, top.nes, color=[POS if v > 0 else NEG for v in top.nes], height=0.62, zorder=2)
        ax.set_yticks(y, [clean_term(t) for t in top.term], fontsize=6.5)
        ax.set_ylim(-0.6, len(top) - 0.4)
        ax.tick_params(axis="y", length=0)
    np.atleast_1d(axes)[-1].set_xlabel("Normalized enrichment score (NES)")
    fig.text(0.0, -0.01, f"GSEA preranked on all RNA-seq weights with a gene symbol. Gene sets with FDR < {fdr}; up to 2 per library and "
             "direction. NES > 0: higher\nin participants with higher factor scores. H MSigDB Hallmark, R Reactome, "
             "GO GO biological process.", fontsize=6.2, color=INK2, va="top")
    save(fig, out_dir, "fig5_enrichment")


def fig_stability(results, out_dir):
    path = os.path.join(results, "main_noplasma", "subset_stability_summary.csv")
    if not need(path):
        return print("Skipped figS1_stability")
    s = pd.read_csv(path, index_col=0).reindex(FACTORS)
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    y = np.arange(len(FACTORS))[::-1]
    ax.hlines(y, s.z_abs_r_min, 1, color=GRID, lw=2, zorder=1)
    ax.scatter(s.z_abs_r_min, y, s=18, color=INK2, zorder=3, label="Worst subset")
    ax.scatter(s.z_abs_r_median, y, s=24, color=NEG, zorder=4, label="Median of 10 subsets")
    for yi, (_, r) in zip(y, s.iterrows()):
        ax.text(1.03, yi, f"{int(round(r.frac_subsets_reproduced * 10))}/10", va="center", fontsize=6.5, color=INK2,
                transform=ax.get_yaxis_transform())
    ax.axvline(0.8, color=INK2, lw=0.6, ls=(0, (3, 2)))
    ax.set_yticks(y, [FLABEL[f] for f in FACTORS])
    left = min(0.5, np.floor((s.z_abs_r_min.min() - 0.05) * 10) / 10)
    ax.set_xlim(max(left, 0), 1.0)
    ax.set_ylim(-0.6, len(FACTORS) - 0.2)
    ax.set_xlabel("|r| of factor scores, subset vs full model")
    ax.legend(frameon=False, loc="lower left", fontsize=6.5)
    ax.text(0.795, len(FACTORS) - 0.45, "stable if ≥ 0.8\nin ≥ 9/10 subsets", ha="right", va="center", fontsize=6.5,
            color=INK2)
    ax.text(1.03, len(FACTORS) - 0.45, "Reproduced", transform=ax.get_yaxis_transform(), va="center", fontsize=6.5,
            color=INK2)
    save(fig, out_dir, "figS1_stability")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    os.makedirs(args.output, exist_ok=True)
    fig_design(args.output)
    for fn in (fig_variance, fig_associations, fig_within_sex, fig_enrichment, fig_stability):
        fn(args.results, args.output)


if __name__ == "__main__":
    main()
