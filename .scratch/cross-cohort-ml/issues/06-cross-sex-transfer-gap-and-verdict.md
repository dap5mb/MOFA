# 06: Cross-sex transfer gap and verdict

**What to build:** The primary result. Sex-stratified frozen models (genes plus covariates) trained on PPMI women and on PPMI men, each scored in PDBP women and PDBP men, giving the **cross-sex transfer gap** in both directions. Training is matched on N (50 subsamples of the larger sex at the smaller sex's size, by diagnosis; PDBP scores averaged across subsamples). Intervals combine 2,000 bootstrap resamples of PDBP participants with subsample variation. Verdict: **shared** if the 90% interval lies within ±0.05 AUC, **sex-specific** if the 95% interval excludes 0 in both directions, otherwise **inconclusive**. A pooled model with sex as a covariate is compared with the stratified models within each PDBP sex as a secondary result.

**Blocked by:** 05

**Status:** ready-for-agent

- [ ] Seam 1 synthetic tests: a signal planted in both sexes → shared; in women only → sex-specific; weak → inconclusive; a gap created only by unequal sex sizes disappears with matched-N training.
- [ ] Real-data table: within-sex and cross-sex AUC per training sex, gap, both intervals, verdict.
- [ ] Pooled-vs-stratified comparison reported as secondary.
