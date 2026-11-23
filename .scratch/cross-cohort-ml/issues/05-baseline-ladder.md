# 05: Baseline ladder

**What to build:** The **baseline ladder** on real data: a covariates-only model (age, library size, genes detected, cell-composition scores; sex in pooled models only; no plate), Factor 7 plus covariates, all 10 projected MOFA factors plus covariates, and genes plus unpenalized covariates. Composition scores (neutrophils, B cells, T/NK cells, monocytes) use a published marker set retrieved from PubMed and passed through the project's reference checker; markers taken from MOFA results are not allowed. Each model is frozen, reported internally in PPMI and scored once in PDBP, as AUC and as gain over covariates only.

**Blocked by:** 04

**Status:** ready-for-agent

- [ ] The marker-set reference passes the reference checker and is added to the verified references.
- [ ] Gene-model covariates are unpenalized.
- [ ] Synthetic test: a signal driven only by composition is captured by the covariates-only model and gives little or no gene gain.
- [ ] Table: each rung's PPMI internal AUC, PDBP AUC with interval, and gain over covariates only.
