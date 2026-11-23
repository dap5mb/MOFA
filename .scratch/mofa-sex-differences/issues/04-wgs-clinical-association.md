# 04: Post-hoc WGS variant association and clinical correlation

**What to build:** Run linear models associating learned factor scores with WGS variants (factor-GWAS / QTL) to detect genetic drivers without matrix sparsity issues during factor training. Compute correlation matrices between latent factors and baseline clinical endpoints (MDS-UPDRS, MoCA, DaTSCAN).

**Blocked by:** `pdbp-projection` 04 (sex × diagnosis interaction test), which settles which factors to carry forward

**Status:** needs-triage (queued after `pdbp-projection` 04–06)

- [ ] Execute linear models testing factor scores against WGS variant matrices adjusting for age, sex, and ancestry PCs ($F_k \sim \text{Genotype} + \text{Covariates}$).
- [ ] Compute Pearson and Spearman correlations between factor scores and MDS-UPDRS subscores, MoCA scores, and DaTSCAN SBR measures.
- [ ] Identify clinical and genetic associations that differ significantly by sex.
- [ ] Export association summary tables and p-value distributions.

**Note (2026-10-07):** the factors to carry forward are the stable, replicated ones in the main model (RNA-seq + CSF): Factor 7 (PD axis) and Factor 10 (aging axis). The clinical-correlation item overlaps with `pdbp-projection` 04; this issue should keep the WGS work plus anything 04 doesn't cover (MDS-UPDRS subscores, DaTSCAN SBR). DaTSCAN isn't yet in the metadata pulled for the evaluation runs.
