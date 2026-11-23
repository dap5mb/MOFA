# 04: Label-free preprocessing and the first real-data run (PPMI internal AUC only)

**What to build:** Real RNA-seq for both cohorts, prepared with the leakage guards: within-cohort size-factor normalization and log transform, plate correction without diagnosis, the 5,000 most variable genes chosen in PPMI only and applied unchanged to PDBP, and standardization within cohort and sex. The tracer model from 02 is run on real PPMI and reports its internal nested-CV AUC. PDBP is prepared but **not scored** in this ticket.

**Blocked by:** 01, 02, 03

**Status:** built 2026-11-23; waiting for the Terra run

- [ ] The gene list is derived from PPMI alone; a check confirms no PDBP sample entered its selection. *(Built: `select_genes.py` refuses non-PPMI rows and, given the PDBP participants file, fails if any PDBP ID was used; the report records `n_pdbp_participants_in_selection`. To tick after the Terra run.)*
- [x] Diagnosis does not appear in any preprocessing design. *(`normalize_cohort.py` reads the participants file without its diagnosis columns and passes only sex and plate on; tested by permuting diagnosis and by recording what preprocessing and plate correction receive.)*
- [ ] The overlap with MOFA's 5,000-gene list is reported (expected to be high). *(`gene_selection_report.json`, `mofa_overlap`; needs the Terra run.)*
- [ ] Internal PPMI nested-CV AUC is reported, overall and per sex; no PDBP classifier scores exist yet. *(`internal_auc.tsv`; needs the Terra run. Nothing in this issue calls the scorer.)*
- [x] Each step is a standalone script with explicit input and output paths; on Terra a notebook calls them in order (no Snakemake on Terra; spec addendum 2026-11-23). *(`ml/ml_wf.ipynb`.)*

## What was built (2026-11-23)

| Step | Script | Output |
|---|---|---|
| Normalize each cohort | `normalize_cohort.py` | `results/ml/preprocessed/{ppmi,pdbp}_expr.csv` (genes x participants) |
| Choose genes in PPMI | `select_genes.py` | `genes_top5000.txt`, `gene_selection_report.json` (aggregate) |
| Apply list, standardize within cohort and sex | `prepare_features.py` | `{ppmi,pdbp}_features_top5000.csv` |
| Train the tracer model on PPMI | `train_ppmi.py` | `results/ml/ppmi_elastic_net/frozen_model.json`, `internal_auc.tsv` (aggregate), `ppmi_oof_scores.csv` |

- **Normalization** reuses MOFA's `preprocess_rnaseq` and `correct_plate` through `shared.py`; MOFA's code is unchanged. PPMI keeps every gene passing MOFA's low-count filter (counts ≥ 10 in ≥ 20% of participants): this is the "all expressed genes" set for issue 09. It is obtained by setting `n_top_genes` to the number of genes, so MOFA's filter is used, not copied. PDBP keeps all genes. Both are restricted to the analysis set (PD or control) before size factors are estimated.
- **Plate correction** uses the label-free design `~ sex`. `--no_plate_correction` gives the issue 09 sensitivity matrix.
- **Gene selection** differs from MOFA's in one respect: MOFA chose its 5,000 genes before plate correction, and these are chosen after label-free plate correction. That is why the overlap is expected to be high rather than complete.
- **Internal AUC per sex:** the trainer now runs the outer CV loop itself (same splits as before, so the tracer's models are unchanged; tested) and keeps the out-of-fold scores. Each sex's AUC is the AUC of the pooled out-of-fold scores within that sex. The model is the pooled tracer model (both sexes, no covariates), so these are per-sex readouts of one model, not sex-stratified models (issue 06).
- **Tests:** `ml/tests/test_preprocessing.py` (16 tests). One test, which runs the real `preprocess_rnaseq`, needs pydeseq2 and is skipped locally; it runs on Terra. All ML test files and MOFA's `test_plate_correction` pass locally, run one file at a time.

## Terra run

All commands for this issue, with their inputs and outputs, are in `ml/STEPS.md` (section "Issue 04").

Run `ml/ml_wf.ipynb` from the repository root after setting `MOFA` in its first cell. Before running, check two input paths against the MOFA notebook: raw counts at `<mofa>/data/ppmi_rnaseq.csv` and `<mofa>/data/pdbp/raw/pdbp_rnaseq.csv`, and MOFA's 5,000 genes as the rows of `<mofa>/data/ppmi_rnaseq_adjusted.csv`. Nested CV on about 1,300 participants and 5,000 genes means 1,050 inner fits plus refits, so use `--n_jobs -1`.

Also run `python -m unittest test_preprocessing` in `ml/tests` on Terra, so the pydeseq2 test runs.

Bring back only `gene_selection_report.json` and `internal_auc.tsv`.
