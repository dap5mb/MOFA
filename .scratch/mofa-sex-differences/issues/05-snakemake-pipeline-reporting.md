# 05: End-to-end Snakemake pipeline orchestration and reporting

**What to build:** Integrate all data preparation, model training, projection, and association steps into an automated Snakemake pipeline. Generate a summary HTML report displaying variance explained, sex-stratified factor loadings, validation metrics, and test coverage using synthetic fixture data.

**Blocked by:** 04: Post-hoc WGS variant association and clinical correlation

**Status:** needs-triage (last in queue)

- [ ] Connect Snakemake rules across all execution stages from harmonization to reporting.
- [ ] Implement an RMarkdown report template rendering multi-group factor plots, variance explained heatmaps, and validation replication metrics.
- [ ] Build a synthetic downsampled multi-omics dataset for CI test suites.
- [ ] Run end-to-end integration tests verifying execution from clean environment to final HTML report output.

**Note (2026-10-07):** the existing `Snakefile` and `config.yaml` predate the current scripts. `config.yaml` says 15 factors and the old raw paths, while the pipeline runs from `mofa_wf.ipynb` on Terra. The 20-test suite in `tests/` already uses synthetic fixtures and can seed the CI item.
