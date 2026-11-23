# 05: Automated Pathway & Network Enrichment on Top Feature Weights

**What to build:** An automated enrichment pipeline for stable, replicated factors. Use rank-based GSEA on the full weight vector as the primary method, with Reactome and STRING PPI on top-weighted features as secondary views, to annotate the PD and aging axes and any sex-differential factors.

**Blocked by:** 04: PDBP Factor Replication & Clinical Association Evaluation

**Status:** closed (2026-11-19)

- [x] (2026-11-19) Extract full weight vectors for the stable, replicated factors (Factors 7 and 10 in the main model). There are no sex-differential factors to add (04, 06), and Factor 4 is excluded as probably technical. Use top-weight cut-offs only for secondary network views, because top-5% overlap across subsets is only moderate (median 0.37 for Factor 7; validation report, Section 2).
- [x] (2026-11-19) GSEA preranked on the full RNA-seq weight vector (`scripts/pathway_enrichment.py`; MSigDB 2024.1 Hallmark, Reactome, GO BP; run locally). Factor 7: 156 gene sets with FDR < 0.05; neutrophil degranulation and inflammatory response on the PD side (lower Factor 7) vs B-cell receptor and immunoglobulin sets on the other, which reads as blood cell composition. Factor 10: 3 sets only; its top genes suggest cytotoxic vs naive lymphocytes.
- [x] (2026-11-19) STRING PPI network views on top-weighted features (`scripts/string_networks.py`; top 50 STRING-mapped genes per direction, score ≥ 0.7, enrichment against the model genes). Factor 7: an 18-gene neutrophil-granule cluster on the PD side (Reactome neutrophil degranulation, FDR < 1e-32); the B-cell side is mostly immunoglobulin genes, which STRING lacks. Factor 10: one cytotoxic cluster (ZNF683, GZMH, LAG3, IFNG, B3GAT1, VCAM1); no pathway terms. Tables in `results/main_noplasma/string/`; figure `figures/writeup/figS2_string_networks`.
- [x] (2026-11-19) Summary table (`results/main_noplasma/enrichment/gsea_prerank_results.csv`) and figure (`figures/writeup/fig5_enrichment`, from `scripts/make_figures.py`).
- [x] Archive enrichment reports alongside validation outputs (`results/main_noplasma/enrichment/`; write-up Results section).
