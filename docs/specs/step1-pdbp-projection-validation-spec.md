# Spec: PDBP Validation Cohort Projection and Multi-Omics Factor Replication

## Problem Statement

Researchers investigating sex-specific molecular differences in Parkinson's Disease across AMP-PD cohorts have successfully trained and validated a multi-group MOFA model (`groups = sex`) on the PPMI discovery cohort. To establish scientific rigor and generalizability, these discovered latent factors and sex-specific axes must be validated on an independent cohort. The user needs an automated, robust pipeline to preprocess the independent PDBP baseline validation cohort, project PDBP samples into the trained PPMI latent factor space, evaluate structural and clinical replication, and run downstream pathway enrichment.

## Solution

Implement an automated PDBP validation and projection pipeline consisting of:
1. Independent preprocessing and batch normalization of PDBP whole blood RNA-seq and CSF/plasma proteomics using protocols identical to PPMI.
2. Strict feature alignment of PDBP matrices against trained PPMI model weights ($W$).
3. Latent factor projection of PDBP baseline samples into the PPMI MOFA factor space utilizing `mofapy2`'s native missing data handling.
4. Multi-faceted evaluation comparing view-wise variance explained ($R^2$), clinical score correlations (MDS-UPDRS, MoCA, diagnosis), and sex-specific factor distributions.
5. Automated pathway and network enrichment on top feature weights of replicated factors.
6. Modular workflow integration of projection (`project_pdbp.py`) and evaluation (`evaluate_pdbp_factors.py`) scripts.

## User Stories

1. As a computational biologist, I want to preprocess PDBP baseline whole blood RNA-seq and CSF/plasma proteomics using the identical normalization and filtering pipeline as PPMI, so that technical artifacts do not bias cohort comparison.
2. As a researcher, I want to align PDBP omic matrices to the exact feature set of the trained PPMI MOFA model, so that matrix projection shapes match without dimension errors.
3. As a researcher, I want to project preprocessed PDBP baseline samples into the PPMI MOFA factor space using native missing data handling, so that single-view and multi-view validation subjects are fully scored.
4. As a downstream analyst, I want to evaluate view-wise variance explained ($R^2$) and factor replication across PPMI and PDBP, so that I can confirm structural generalization of discovered latent axes.
5. As a clinical researcher, I want to automatically test projected PDBP latent factors against diagnosis (`Ctrl` vs. `PD`) and clinical scores (MDS-UPDRS, MoCA), so that disease associations are independently replicated.
6. As a molecular biologist, I want to run automated pathway and network enrichment (GSEA, Reactome, STRING PPI) on top feature weights of replicated factors, so that the biological processes driving sex-specific axes are annotated.
7. As a workflow maintainer, I want PDBP projection and evaluation scripts (`project_pdbp.py`, `evaluate_pdbp_factors.py`) integrated into the project's modular workflow architecture, so that the entire discovery-to-validation pipeline runs reproducibly.

## Implementation Decisions

- **Feature Alignment:** PDBP omic matrices are strictly subset to the exact feature set (rows / genes and proteins) present in the trained PPMI model weights ($W$) to ensure mathematical consistency during projection.
- **Batch Correction & Normalization Strategy:** Technical batch correction and normalization (VST for RNA-seq, log2 + plate median centering for proteomics) are performed independently within PDBP using identical protocols to avoid data leakage or artificial batch alignment.
- **Projection Execution:** Implemented via a dedicated Python script (`project_pdbp.py`) within `mofa/scripts/` that loads the trained PPMI HDF5 model, loads preprocessed PDBP data views and master metadata, invokes `mofapy2`'s projection/prediction methods, and saves PDBP latent factors ($Z_{\text{PDBP}}$) and replication metrics.
- **Missing Data Handling:** Relies on `mofapy2`'s native missing data handling during projection, allowing single-view or partial-view PDBP subjects to be projected into the latent space using available views without artificial imputation noise.
- **Replication & Evaluation Metrics:** Multi-faceted evaluation including (1) view-wise variance explained ($R^2$) correlation between PPMI and PDBP, (2) association significance of projected factors $Z_{\text{PDBP}}$ with diagnosis (`Ctrl` vs. `PD`) and clinical scores (MDS-UPDRS), and (3) conservation of sex-specific factor loading patterns.
- **Clinical & Sex-Stratified Association Pipeline:** Implemented in `evaluate_pdbp_factors.py` to perform linear/logistic regression between projected PDBP factors and clinical metadata, generating comparative replication plots.
- **Pathway Enrichment:** Automated functional annotation on top feature weights (top 1%–5%) of replicated factors to link statistical factors back to biological processes.
- **Modular Integration:** Integrated into the project's modular script and workflow structure alongside existing PPMI scripts.

## Testing Decisions

- **Testing Approach:** Test at the highest integration seam where preprocessed PDBP matrices and trained PPMI model weights meet. Tests will verify external behavior (that output latent factor matrices $Z_{\text{PDBP}}$ have expected dimensions, match sample metadata, and successfully output validation metrics such as $R^2$ and clinical regression $p$-values) rather than internal implementation details.
- **Prior Art:** Aligns with prior pre-flight checks and validation reporting established in PPMI discovery workflows (`prep_utils.py`, `evaluate_factors.py`, and `mofa-step1-validation-report.md`).

## Out of Scope

- Modifying the core PPMI model hyperparameters or re-training the PPMI discovery MOFA model.
- Longitudinal trajectory projection (focus is strictly on baseline visit `visit_month = 0.0`).
- Inclusion of WGS variants in the primary MOFA matrix (retained for post-hoc GWAS/QTL factor association).

## Further Notes

- Adheres to ADR-0001 ("Multi-omics sex-specific factor analysis strategy for AMP-PD").
- Since GitHub CLI (`gh`) is not installed in the current environment, this spec is published as a local version-controlled Markdown artifact in `docs/specs/step1-pdbp-projection-validation-spec.md`.
