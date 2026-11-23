# Cross-cohort ML: run steps

Every step of the cross-cohort ML analysis (spec: `docs/specs/step2-cross-cohort-ml-spec.md`), in run order, with one section per issue (`.scratch/cross-cohort-ml/issues/`). When an issue adds or changes a step, update its own ticket, this file, and `ml/ml_wf.ipynb`.

## Conventions

- **Where:** real data runs only on Terra, where `ml/ml_wf.ipynb` calls these scripts in order. Snakemake is not allowed on Terra; `ml/Snakefile` is for local synthetic runs only (spec addendum 2026-11-23).
- **From:** the repository root. `<mofa>` is the MOFA working folder that holds `data/`.
- **Folders:** `AS = results/ml/analysis_set`, `PRE = results/ml/preprocessed`, `EN = results/ml/ppmi_elastic_net`.
- **What may leave Terra:** only outputs marked **aggregate** (counts, percentages, AUCs, gene-selection summaries). Everything else is participant-level and stays under `results/`.
- **PDBP scoring:** PDBP is scored only from frozen model files, once per model named in the spec. Any model added after PDBP is scored is labelled exploratory.

## Tests

Run each test file on its own (`python -m unittest <file>` in `ml/tests`; pytest is not on Terra):

```bash
cd ml/tests && for f in test_shared test_tracer test_train_score_boundary test_analysis_set test_preprocessing; do python -m unittest $f; done
```

Some tests need packages that are on Terra but may be missing locally (pydeseq2 in `test_preprocessing`); they are skipped when the package is missing.

## Issue 01: shared preprocessing (done 2026-11-23)

No run step. The ML module imports MOFA's metadata query, normalization, plate correction and PDBP preparation through `ml/scripts/shared.py`. Plate correction takes `include_diagnosis=False` for the label-free design.

## Issue 02: tracer on synthetic cohorts (done 2026-11-23)

Synthetic data only; no AMP-PD data. Optional check that the frozen-model path works end to end.

```bash
python3 ml/scripts/run_tracer.py --scenario signal --model results/ml/synthetic/signal/frozen_model.json --scores results/ml/synthetic/signal/pdbp_scores.csv --auc results/ml/synthetic/signal/auc.tsv
```

```bash
python3 ml/scripts/run_tracer.py --scenario null --model results/ml/synthetic/null/frozen_model.json --scores results/ml/synthetic/null/pdbp_scores.csv --auc results/ml/synthetic/null/auc.tsv
```

Locally, `snakemake -s ml/Snakefile --cores 1` runs both.

## Issue 03: analysis set and coverage report (done 2026-11-23)

| Step | Script | Inputs | Outputs |
|---|---|---|---|
| 3.1 | `query_ml_fields.py` | AMP-PD 2023_v4 (BigQuery) | `AS/{ppmi,pdbp}_ml_fields.csv` |
| 3.2 | `analysis_set.py` | MOFA metadata, ML fields, MOFA RNA-seq matrix (header only) | `AS/{ppmi,pdbp}_participants.csv` |
| 3.3 | `coverage_report.py` | both participants files | `AS/coverage_report/` (**aggregate**) |

```bash
python3 ml/scripts/query_ml_fields.py --cohort PPMI --output results/ml/analysis_set/ppmi_ml_fields.csv
```

```bash
python3 ml/scripts/query_ml_fields.py --cohort PDBP --output results/ml/analysis_set/pdbp_ml_fields.csv
```

```bash
python3 ml/scripts/analysis_set.py --cohort PPMI --mofa_metadata <mofa>/data/metadata_m0_PPMI.csv --ml_fields results/ml/analysis_set/ppmi_ml_fields.csv --rnaseq_matrix <mofa>/data/ppmi_rnaseq_adjusted.csv --output results/ml/analysis_set/ppmi_participants.csv
```

```bash
python3 ml/scripts/analysis_set.py --cohort PDBP --mofa_metadata <mofa>/data/metadata_m0_PDBP.csv --ml_fields results/ml/analysis_set/pdbp_ml_fields.csv --rnaseq_matrix <mofa>/data/pdbp_rnaseq_adjusted.csv --output results/ml/analysis_set/pdbp_participants.csv
```

```bash
python3 ml/scripts/coverage_report.py --participants results/ml/analysis_set/ppmi_participants.csv results/ml/analysis_set/pdbp_participants.csv --output_dir results/ml/analysis_set/coverage_report
```

The duration bands and the PPMI-like subset were fixed from this report (spec addendum 2026-11-23) and are coded in `analysis_set.FIXED_DEFINITIONS`.

## Issue 04: label-free preprocessing and the internal PPMI AUC (built 2026-11-23; Terra run pending)

| Step | Script | Inputs | Outputs |
|---|---|---|---|
| 4.1 | `normalize_cohort.py` | raw counts, participants file | `PRE/{ppmi,pdbp}_expr.csv` |
| 4.2 | `select_genes.py` | PPMI `expr`, both participants files, MOFA's gene list | `PRE/genes_top5000.txt`, `PRE/gene_selection_report.json` (**aggregate**) |
| 4.3 | `prepare_features.py` | each cohort's `expr`, participants file, gene list | `PRE/{ppmi,pdbp}_features_top5000.csv` |
| 4.4 | `train_ppmi.py` | PPMI features, PPMI participants | `EN/frozen_model.json`, `EN/internal_auc.tsv` (**aggregate**), `EN/ppmi_oof_scores.csv` |

Check the input paths against the MOFA notebook before the first run: raw counts at `<mofa>/data/ppmi_rnaseq.csv` and `<mofa>/data/pdbp/raw/pdbp_rnaseq.csv`; MOFA's 5,000 genes are the rows of `<mofa>/data/ppmi_rnaseq_adjusted.csv`.

**4.1 Normalize each cohort on its own.** Size factors and log1p, then plate correction with `~ sex` (no diagnosis). PPMI keeps the genes passing MOFA's low-count filter; PDBP keeps all genes. Add `--no_plate_correction` for the issue 09 sensitivity matrix.

```bash
python3 ml/scripts/normalize_cohort.py --cohort PPMI --counts <mofa>/data/ppmi_rnaseq.csv --participants results/ml/analysis_set/ppmi_participants.csv --output results/ml/preprocessed/ppmi_expr.csv
```

```bash
python3 ml/scripts/normalize_cohort.py --cohort PDBP --counts <mofa>/data/pdbp/raw/pdbp_rnaseq.csv --participants results/ml/analysis_set/pdbp_participants.csv --output results/ml/preprocessed/pdbp_expr.csv
```

**4.2 Choose the 5,000 most variable genes in PPMI only.** The step fails if any PDBP participant was used, and it reports the overlap with MOFA's genes.

```bash
python3 ml/scripts/select_genes.py --expr results/ml/preprocessed/ppmi_expr.csv --participants results/ml/analysis_set/ppmi_participants.csv --pdbp_participants results/ml/analysis_set/pdbp_participants.csv --mofa_genes <mofa>/data/ppmi_rnaseq_adjusted.csv --n_top 5000 --genes_out results/ml/preprocessed/genes_top5000.txt --report results/ml/preprocessed/gene_selection_report.json
```

**4.3 Apply the list unchanged to both cohorts** and standardize within cohort and sex. PDBP is prepared here but not scored.

```bash
python3 ml/scripts/prepare_features.py --expr results/ml/preprocessed/ppmi_expr.csv --participants results/ml/analysis_set/ppmi_participants.csv --genes results/ml/preprocessed/genes_top5000.txt --output results/ml/preprocessed/ppmi_features_top5000.csv
```

```bash
python3 ml/scripts/prepare_features.py --expr results/ml/preprocessed/pdbp_expr.csv --participants results/ml/analysis_set/pdbp_participants.csv --genes results/ml/preprocessed/genes_top5000.txt --output results/ml/preprocessed/pdbp_features_top5000.csv
```

**4.4 Train the tracer elastic net on PPMI.** Nested CV (outer 5, inner 10, 21 grid points) on about 1,300 × 5,000; `--n_jobs -1` uses all cores. The internal AUC is reported overall and per sex.

```bash
python3 ml/scripts/train_ppmi.py --features results/ml/preprocessed/ppmi_features_top5000.csv --participants results/ml/analysis_set/ppmi_participants.csv --model results/ml/ppmi_elastic_net/frozen_model.json --auc results/ml/ppmi_elastic_net/internal_auc.tsv --oof results/ml/ppmi_elastic_net/ppmi_oof_scores.csv --n_jobs -1
```

Bring back: `gene_selection_report.json` and `internal_auc.tsv`.

## Issues 05–13

Not built yet. Each adds its section here when it is built.
