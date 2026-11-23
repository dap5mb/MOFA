# AMP-PD Multi-Omics Sex-Specific Factor Analysis

This document defines terms and domain concepts for analyzing sex-specific differences in Parkinson's Disease using AMP-PD data.

## Cohort definitions

- **PPMI**: Parkinson's Progression Markers Initiative cohort. Serves as the primary discovery cohort for multi-omics factor modeling.
- **PDBP**: Parkinson's Disease Biomarkers Program cohort. Serves as the primary validation cohort for factor projection and replication.
- **Baseline visit (`M0`)**: Clinical and biospecimen data collected at `visit_month = 0.0`.

## Omic views

- **Whole blood RNAseq (`rnaseq_WB_RWTS`)**: Gene expression view that captures systemic transcriptomic variation.
- **Proteomics (`proteomics_PLA_PDIA` / `proteomics_CSF_PDIA`)**: Plasma and cerebrospinal fluid protein abundance views.
- **WGS view**: Whole-genome sequencing data tested post hoc for factor associations (GWAS, QTL, or gene burden) rather than included directly in the MOFA feature matrix.

## Factor modeling strategy

- **Group-MOFA**: Multi-group multi-omics factor analysis (`MOFA2` with `groups = sex`). Decomposes variation into factors shared between sexes and factors active in one sex.
- **Anchored factor training**: Trains MOFA on a complete-case subset ($N \approx 300$) with concurrent RNAseq and proteomics data to prevent RNAseq sample counts from dominating factor priors. The model then projects factors onto the full cohort ($N \approx 1,500$).

## Cross-cohort predictive analysis

A gene-level companion to the factor analysis: does the blood RNA-seq signal that separates PD from controls differ between women and men? Classifiers are built in PPMI and tested in PDBP.

- **Frozen model**: A classifier whose features, scaling and coefficients were fixed after training in PPMI. Only frozen models are scored in PDBP.
  _Avoid_: final model, locked model
- **Cross-sex transfer gap**: Within-sex AUC minus cross-sex AUC in PDBP (for example, a model trained on PPMI women scored on PDBP women vs PDBP men), with training sets matched on size.
- **Shared predictive signal**: A PD-vs-control signal whose cross-sex transfer gap lies within the pre-specified equivalence margin.
- **Sex-specific predictive signal**: A PD-vs-control signal whose within-sex performance exceeds its cross-sex performance in PDBP.
- **Inconclusive**: The verdict when the cross-sex transfer gap neither falls within the equivalence margin nor excludes zero.
- **PPMI-like subset**: PDBP participants with short disease duration and no PD medication at baseline, resembling PPMI's de novo enrollment.
- **Baseline ladder**: The ordered comparison of covariates-only, MOFA-factor, and gene-plus-covariate models; a gene model's value is its gain over the covariates-only model.

### Ambiguous terms

- **Sex-specific** means different things at different levels. For a MOFA factor it means *active in one sex*; for the predictive analysis it means a *sex-specific predictive signal* as defined above. The earlier PPMI IR3 gene categories (sex-neutral, sex-specific, sex-dimorphic, sex-modulated) classify single genes and are not used in this project.
