# 1. Multi-omics sex-specific factor analysis strategy for AMP-PD

* Status: Accepted
* Date: 2025-05-18

## Context

We need an automated workflow to identify sex-specific molecular differences in Parkinson's disease using baseline data (`visit_month = 0.0`) from the PPMI and PDBP cohorts in AMP-PD. Available data includes whole blood RNAseq ($N \approx 1,500$ per cohort), plasma and CSF proteomics ($N \approx 100-300$ per cohort), WGS, and clinical assessments (MDS-UPDRS, MoCA, DaTSCAN).

This study addresses three analytical constraints:
1. Isolating shared molecular variation from sex-specific variation.
2. Handling sample size differences across modalities (RNAseq $N \gg$ proteomics $N$).
3. Managing the high dimensionality and sparsity of WGS variants.

## Decision

We will implement a Snakemake workflow using multi-group MOFA (`MOFA2` / `mofapy2`) with this structure:

1. **Discovery and validation split**: Train MOFA on PPMI baseline data and project learned factors onto PDBP baseline data.
2. **Multi-group MOFA (`groups = sex`)**: Set sex as the grouping variable in MOFA to separate shared factors from sex-specific factors.
3. **Omic views and feature selection**:
   - View 1: Whole blood RNAseq (variance-stabilizing transformation applied to top variable genes).
   - View 2: Plasma and CSF proteomics abundance.
   - Post-hoc WGS: Exclude WGS from the initial MOFA matrix to avoid sparsity issues. Test variant associations against learned factor weights post hoc using factor-GWAS or QTL models ($F_k \sim \text{Genotype}$).
4. **Anchored training strategy**: Train initial MOFA weights on the complete-case subset ($N \approx 300$ with both RNAseq and proteomics) to balance view contributions. Project the remaining single-view subjects ($N \approx 1,200$) into this factor space.
5. **Workflow orchestration**: Implement the pipeline in Snakemake using modular R and Python scripts.

## Consequences

### Positive
- Keeps RNAseq sample counts ($N \approx 1,500$) from dominating proteomics priors ($N \approx 300$).
- Isolates sex-specific factors using a multi-group Bayesian framework.
- Evaluates factor replication on an independent clinical cohort (PDBP).

### Trade-offs
- The complete-case subset ($N \approx 300$) limits the initial training sample size, requiring convergence checks across multiple random seeds.
- Excludes WGS from direct multi-omic factor construction.
