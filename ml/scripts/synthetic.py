"""
ml/scripts/synthetic.py

Two simulated cohorts (PPMI-like and PDBP-like) on the same genes, for testing the pipeline
end to end without AMP-PD data. Each cohort gets its own per-gene offset (cohort shift) and a
sex mean difference on some genes; the 'signal' scenario shifts a block of genes in PD in
both cohorts, the 'null' scenario plants nothing.

Output matches what the trainer and scorer take: expression participants x genes (log scale),
metadata indexed by participant_id with cohort, sex and diagnosis ("PD" / "Control").
"""

import numpy as np
import pandas as pd

SCENARIOS = ("signal", "null")


def _cohort(cohort, n, gene_means, gene_sds, signal_genes, effect, sex_genes, pd_frac, offset_sd, rng):
    n_genes = len(gene_means)
    ids = pd.Index([f"{cohort}_{i:04d}" for i in range(n)], name="participant_id")
    sex = rng.choice(["Female", "Male"], n)
    pd_case = rng.random(n) < pd_frac

    expr = gene_means + rng.normal(0, offset_sd, n_genes) + rng.normal(size=(n, n_genes)) * gene_sds
    expr[:, sex_genes] += np.where(sex == "Male", 1.0, 0.0)[:, None] * gene_sds[sex_genes]
    expr[:, signal_genes] += pd_case[:, None] * effect * gene_sds[signal_genes]

    genes = [f"G{j:04d}" for j in range(n_genes)]
    meta = pd.DataFrame({"cohort": cohort, "sex": sex,
                         "diagnosis": np.where(pd_case, "PD", "Control")}, index=ids)
    return pd.DataFrame(expr, index=ids, columns=genes), meta


def make_synthetic_cohorts(scenario: str = "signal", n_ppmi: int = 300, n_pdbp: int = 300,
                           n_genes: int = 300, n_signal: int = 20, effect: float = 0.5,
                           pd_frac: float = 0.55, seed: int = 0):
    """Returns (ppmi_expr, ppmi_meta, pdbp_expr, pdbp_meta). effect is in gene SDs."""
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}")
    rng = np.random.default_rng(seed)
    gene_means = rng.uniform(4, 10, n_genes)
    gene_sds = rng.lognormal(np.log(0.5), 0.3, n_genes)
    order = rng.permutation(n_genes)
    signal_genes = order[:n_signal]
    sex_genes = order[n_signal:n_signal + 10]
    eff = effect if scenario == "signal" else 0.0

    ppmi = _cohort("PPMI", n_ppmi, gene_means, gene_sds, signal_genes, eff, sex_genes, pd_frac, 0.2, rng)
    pdbp = _cohort("PDBP", n_pdbp, gene_means, gene_sds, signal_genes, eff, sex_genes, pd_frac, 0.2, rng)
    return (*ppmi, *pdbp)
