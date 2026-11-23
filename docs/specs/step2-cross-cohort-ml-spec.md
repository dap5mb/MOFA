# Spec: Cross-Cohort Predictive Analysis of Sex Differences in the PD Blood Transcriptome

Status: `ready-for-agent` · Date: 2026-11-23 · Decisions: ADR-0002 · Glossary: `CONTEXT.md` ("Cross-cohort predictive analysis")

## Problem Statement

The multi-group MOFA analysis found no large sex-specific multi-omic axis in PD, but it can only see effects big enough to form a factor (the replicated PD factor explains about 1% of RNA-seq variance) and its interaction tests had power only for differences of about 0.3 SD. The earlier PPMI-only classifiers (PPMI IR3 counts, one 70/30 split, sex-stratified LASSO/RF/XGBoost/SVM) reported high AUCs (LASSO 0.90 in women, 0.87 in men) but were never tested in an independent cohort, were not compared with technical or cell-composition baselines, and cannot say whether the predictive signal differs between women and men. The researcher needs a gene-level predictive analysis that asks the sex question with external validation in PDBP, defends against leakage and technical confounding, and can return a "shared" verdict rather than only "no significant difference".

## Solution

A separate analysis module, alongside but not inside the MOFA code, that:

1. Builds a baseline analysis set from AMP-PD (2023_v4) whole-blood RNA-seq for PPMI and PDBP, using the MOFA participant definition (baseline, PD vs control), with genetic-carrier, disease-duration and medication flags.
2. Normalizes each cohort on its own and applies leakage guards: diagnosis never enters preprocessing, and nothing from PDBP influences feature selection, scaling parameters or hyperparameters.
3. Trains elastic-net (primary) and XGBoost (comparison) classifiers in PPMI with nested cross-validation, refits on all PPMI, and saves **frozen models**.
4. Scores PDBP once per pre-specified model and evaluates a **baseline ladder** (covariates only → MOFA factors → genes plus covariates).
5. Estimates the **cross-sex transfer gap** with training matched on N, and returns a three-way verdict: **shared predictive signal**, **sex-specific predictive signal**, or **inconclusive**.
6. Runs the pre-specified sensitivity analyses, a secondary PD-only regression of MDS-UPDRS III and MoCA, and a descriptive stability-selection layer.
7. Feeds a separate write-up that will become a section of the MOFA paper.

## User Stories

### Analysis set
1. As a researcher, I want both cohorts taken from the AMP-PD release with the same RNA-seq processing, so that a PPMI→PDBP performance drop cannot be blamed on different pipelines.
2. As a researcher, I want the analysis set to match the MOFA participant definition (baseline visit, PD or control, recorded sex, whole-blood RNA-seq), so that the two analyses describe the same people.
3. As a researcher, I want PPMI participants with a diagnosis other than PD or control excluded, so that the classification target is unambiguous.
4. As a researcher, I want each participant flagged as a genetic carrier or not (WGS-confirmed PD mutation and enrollment genetic status), so that I can run a sensitivity analysis without carriers.
5. As a researcher, I want disease duration (age at baseline minus age at diagnosis) and PD medication use at baseline derived for every PD participant, so that I can describe cohort shift and define the PPMI-like subset.
6. As a researcher, I want a coverage report of duration and medication fields in PDBP, so that I know before scoring whether the PPMI-like subset is large enough to use.
7. As a researcher, I want MDS-UPDRS III selected with the same exam-selection rule as MOFA (off-medication exam first), so that clinical targets are comparable across analyses.

### Preprocessing and leakage guards
8. As a researcher, I want counts normalized within each cohort (median-of-ratios size factors, log transform), so that cohorts are not normalized against each other.
9. As a researcher, I want plate corrected within each cohort without diagnosis in the design, so that the label never enters the features.
10. As a researcher, I want a no-plate-correction version of every primary model, so that I can see how much of the signal follows plate.
11. As a researcher, I want the 5,000-gene feature list recomputed from PPMI only, by variance, after label-free preprocessing, so that feature selection uses no labels and no PDBP data.
12. As a researcher, I want an all-expressed-genes feature set as a sensitivity analysis, so that conclusions do not depend on the variance filter.
13. As a researcher, I want genes standardized within cohort and sex (as in the MOFA projection), so that expression offsets between cohorts and sex mean differences do not drive predictions.
14. As a researcher, I want a sensitivity analysis that scales PDBP with PPMI's means and SDs, so that I can see whether cohort offsets matter.

### Covariates and baselines
15. As a researcher, I want cell-composition scores (neutrophils, B cells, T/NK cells, monocytes) computed from a published, PubMed-verified marker set, so that composition can be controlled for without blood counts and without reusing MOFA-derived markers.
16. As a researcher, I want a covariates-only model (age, library size, genes detected, composition scores; sex in pooled models only; no plate), so that I know how far technical and composition variables alone predict PD.
17. As a researcher, I want MOFA-factor models (Factor 7 alone, then all 10 projected factors, each with covariates), so that the gene models are compared against the factor-level result.
18. As a researcher, I want gene models to carry the covariates unpenalized, so that any gain from genes is measured on top of covariates that are always included.
19. As a reader, I want each gene model reported as its gain in AUC over the covariates-only model, so that "high AUC" cannot hide a composition or technical signal.

### Training and freezing
20. As a researcher, I want elastic-net logistic regression as the primary model and XGBoost as a nonlinear comparison, so that the analysis is interpretable and still checks for nonlinear signal.
21. As a researcher, I want nested cross-validation in PPMI (outer 5-fold for the internal estimate, inner 10-fold for hyperparameters), so that the internal AUC is not optimistically biased.
22. As a researcher, I want class weights instead of down-sampling, so that no participants are discarded.
23. As a researcher, I want each model refit on all PPMI with the chosen hyperparameters and saved as a frozen model (features, scaling rule, hyperparameters, coefficients), so that PDBP is scored by exactly one fixed model.
24. As a reader, I want the internal PPMI AUC reported next to the PDBP AUC, so that the drop between cohorts is visible.

### PDBP scoring and the sex question
25. As a researcher, I want PDBP scored only from frozen model files, once per model named in this spec, so that the external validation is defensible.
26. As a researcher, I want any model added after PDBP is scored labelled exploratory, so that the reported results keep to the committed plan.
27. As a researcher, I want sex-stratified models trained on PPMI women and on PPMI men and each scored in both PDBP women and PDBP men, so that I can estimate the cross-sex transfer gap in both directions.
28. As a researcher, I want training sets matched on N by drawing 50 subsamples of the larger sex at the smaller sex's size (and per diagnosis), averaging their PDBP scores, so that a gap cannot come from sample size.
29. As a researcher, I want confidence intervals that combine a 2,000-resample bootstrap over PDBP participants with variation across training subsamples, so that both sources of uncertainty are included.
30. As a reader, I want the verdict "shared" when the 90% interval of the gap lies inside ±0.05 AUC, "sex-specific" when the 95% interval excludes 0 in both directions, and "inconclusive" otherwise, so that a null result is not reported as "no difference".
31. As a researcher, I want a pooled model (sex as a covariate) compared with the sex-stratified models within each PDBP sex, so that I have the secondary stratified-vs-pooled comparison.
32. As a researcher, I want the reverse direction (train in PDBP, score PPMI) run with the same pipeline as a secondary analysis, so that conclusions do not depend on which cohort trained the model.

### Cohort shift and sensitivity analyses
33. As a researcher, I want PDBP AUC reported within disease-duration bands and by medication use, so that I can see where performance falls.
34. As a researcher, I want the primary models also scored in the PPMI-like subset of PDBP, so that a PPMI→PDBP drop is read as overfitting only if it persists there.
35. As a researcher, I want the carrier-excluded analysis, the no-plate-correction analysis, the all-genes feature set and PPMI-scale standardization each run as named sensitivity analyses, so that each robustness claim has its own result.

### Clinical scores
36. As a researcher, I want MDS-UPDRS III and MoCA regressed on expression in PD cases only (train in PPMI, test in PDBP, same leakage guards), so that I can see whether the signal tracks severity or only diagnosis.
37. As a reader, I want clinical results reported as Spearman ρ with a confidence interval, so that weak associations are not overstated.

### Description and reporting
38. As a researcher, I want stability selection (100 × 80% subsamples, elastic net) run within each sex in PPMI and, separately, in PDBP, so that I can list genes that are stable in both cohorts.
39. As a reader, I want only genes stable in both cohorts interpreted biologically, so that gene lists reflect replicated signal.
40. As a reader, I want ROC curves, a transfer-gap plot with the equivalence margin, and a baseline-ladder figure, so that the main results can be read at a glance.
41. As a writer, I want a separate write-up for this analysis, written to slot into the MOFA paper as a section, so that the two analyses stay physically separate but tell one story.

### Workflow
42. As a workflow maintainer, I want the module in its own folder next to the MOFA code, as Snakemake rules with pytest tests, running on Terra, so that the AMP-PD data stays where the data use agreement requires and the pipeline reruns end to end. *(Amended 2026-11-23: Snakemake cannot run on Terra; see Addendum.)*
43. As a workflow maintainer, I want the module to consume MOFA outputs (projected factor scores, preprocessing utilities, baseline metadata) rather than copy them, so that there is one source of truth.

## Implementation Decisions

- **Primary question:** Sex-specificity of the PD-vs-control predictive signal; classifier performance is reported but is not the headline (ADR-0002).
- **Location:** A new module beside the MOFA code, not inside it. It reuses MOFA preprocessing, metadata and PDBP-projection utilities and reads MOFA factor scores as inputs. Python (scikit-learn, xgboost), Snakemake, pytest, run on Terra. *(Amended 2026-11-23: see Addendum.)*
- **Data:** AMP-PD 2023_v4, baseline (`visit_month = 0`), whole-blood RNA-seq only. CSF proteomics is out (about 41 PDBP women with CSF).
- **Analysis set:** The MOFA participant definition, restricted to PD and control. Carriers stay in the primary analysis; a sensitivity analysis excludes them.
- **Validation direction:** PPMI trains and PDBP tests (primary); PDBP→PPMI is secondary. Pooled training is not used.
- **Preprocessing modules:**
  - *Cohort normalizer:* size factors and log transform within a cohort; plate correction without diagnosis (primary) or none (sensitivity).
  - *Feature selector:* top 5,000 genes by variance in label-free PPMI data (primary); all expressed genes (sensitivity). Output is a gene list consumed unchanged by PDBP.
  - *Scaler:* standardizes within cohort and sex (primary) or with PPMI parameters (sensitivity).
- **Covariate builder:** age, library size, genes detected, marker-based composition scores (published marker set, verified through the project's reference checker); sex only in pooled models; plate never used as a feature.
- **Model trainer:** nested CV (outer 5, inner 10); elastic net tunes alpha and lambda, XGBoost tunes depth, learning rate and rounds; class weights; covariates unpenalized in elastic net. Output is a frozen model artifact holding the gene list, scaling rule, covariate definition, hyperparameters, coefficients/trees and the internal CV AUC.
- **Scorer:** the only entry point that touches PDBP for classification. It takes a frozen model artifact and a PDBP matrix and returns per-participant scores. It cannot fit anything.
- **Transfer-gap evaluator:** takes PDBP scores, labels and sex; computes within-sex and cross-sex AUC per training sex, the gap, bootstrap-plus-subsample intervals and the verdict.
  - Matched-N training: 50 subsamples, stratified by diagnosis.
  - Bootstrap: 2,000 resamples over PDBP participants.
  - Equivalence margin: ±0.05 AUC.
  - Verdict rule:

    | Verdict | Rule |
    |---|---|
    | Shared | 90% interval inside the margin |
    | Sex-specific | 95% interval excludes 0, in both directions |
    | Inconclusive | Anything else |
- **Baseline ladder:** covariates only; Factor 7 + covariates; all 10 factors + covariates; genes + covariates. Each is reported as AUC and as gain over covariates only.
- **Cohort-shift reporting:** PDBP AUC by disease-duration band and by medication use; separate scoring of the PPMI-like subset (bands and thresholds fixed before scoring, once the field-coverage check is in).
- **Clinical regression:** elastic-net regression of MDS-UPDRS III and MoCA in PD cases, same pipeline and guards, Spearman ρ with bootstrap interval.
- **Stability description:** 100 × 80% subsamples within each sex and cohort; selection frequency ≥ 0.8 counts as stable; genes stable in both cohorts are reported. No test depends on it.
- **Locked plan:** this spec and ADR-0002 are fixed before any PDBP classifier score exists. Additions after scoring are labelled exploratory.
- **Write-up:** a separate document for this analysis, intended as a section of the MOFA paper. References are PubMed-retrieved and pass the project's reference checker.

## Testing Decisions

- **What makes a good test:** tests drive the module through its public entry points and check outputs (verdicts, AUCs, frozen-model contents, scores). They do not inspect internal helpers.
- **Seam 1, whole pipeline on synthetic cohorts (primary):** two simulated cohorts (PPMI-like and PDBP-like) with known structure, run from processed matrices to verdict:
  - a signal planted in both sexes yields **shared**;
  - a signal planted in women only yields **sex-specific**;
  - a weak or noisy signal yields **inconclusive**;
  - a signal tied to plate or composition is captured by the covariates-only model, or appears as a gap between plate correction and no correction, rather than as a gene-model gain;
  - matched-N training removes a gap that is caused only by unequal sex sizes.
- **Seam 2, the train/score boundary (leakage guard):**
  - changing PDBP values or labels leaves every frozen parameter (gene list, scaling rule, hyperparameters, coefficients) unchanged;
  - the scorer cannot fit;
  - passing PDBP participants to the trainer or feature selector raises an error;
  - diagnosis never appears in the plate-correction design.
- **No separate unit tests** for deconvolution, subsampling or the bootstrap; Seam 1 covers them.
- **Prior art:** `test_pdbp_projection.py` (synthetic model artifact plus synthetic PDBP matrix through the public projection entry point); the toy check for model comparison, which planted a sex-gene block and confirmed it formed its own factor; `test_clinical_models.py` for clinical-model outputs.

## Out of Scope

- CSF and plasma proteomics as classifier features.
- Longitudinal visits and progression.
- Pooled PPMI + PDBP training as a primary analysis.
- Gene × sex interaction terms as a primary test.
- Gene-level classification categories from the IR3 work (sex-neutral, sex-specific, sex-dimorphic, sex-modulated).
- Retraining or changing the MOFA models.
- WGS associations.
- Random forest and SVM.

## Further Notes

- Follows ADR-0001 for the discovery/validation split and ADR-0002 for the locked validation and estimand.
- Supersedes the PPMI IR3 ML results, which are retired.
- Published as a local Markdown spec, like the step-1 spec; GitHub Issues were not used.
- Fields to confirm on Terra at build time: how carriers are flagged (`has_known_PD_Mutation_in_WGS`, `genetic_status_enrollment`), `age_at_diagnosis`, `use_of_pd_medication`, `on_levodopa`. These are in the AMP-PD schema; their coverage in PDBP is not yet known.
- Known glossary inconsistency (not fixed here): "Anchored factor training" in `CONTEXT.md` and ADR-0001 describe complete-case training, but the MOFA write-up trains on all participants with at least one view.

## Addendum 2026-11-23: how the workflow runs on Terra

Snakemake is not allowed on Terra, and pytest is not installed there. This changes how the pipeline is driven, not what it computes; no analysis decision, threshold or test in this spec changes, and no PDBP classifier score existed when it was made.

- **Steps are scripts.** Every step is a standalone script with explicit input and output paths, which reads files and writes files. The analysis logic lives in the scripts and the modules they import, never in the driver.
- **On Terra, a notebook drives the run**, calling the scripts in order, as `mofa_wf.ipynb` does for MOFA. This is the run of record for real data.
- **`ml/Snakefile` is optional and local only**, for rerunning the synthetic tracer and other synthetic checks. Nothing on Terra depends on it, and each of its rules only calls a script.
- **Tests are `unittest` test cases**, so they run on Terra with `python -m unittest <file>` and locally with either `unittest` or pytest.
- Story 42 and the "Location" decision above are read with these changes.
