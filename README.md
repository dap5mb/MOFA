# MOFA Multi-Omics Workflow

This directory contains the pipeline for multi-group (sex-specific) factor analysis on AMP-PD baseline data.

## Workflow Execution

1. **Environment**: Ensure `mofapy2`, `pandas`, `numpy`, and `matplotlib` are installed.
2. **Preprocessing**:
   - Run `mofa_wf.ipynb` to process raw RNA-seq and Proteomics matrices.
   - This generates the long-format CSV required by `mofapy2`.
3. **Training**:
   - Execute `train_mofa.py` to run the multi-group MOFA model.
   - Outputs (HDF5 model files and ELBO diagnostic plots) will be saved to the configured output directory.

## File Structure

- `prep_utils.py`: Contains `validate_and_intersect_views` for sample alignment (Task 1.2.4).
- `proteomics_qc.py`: Handles median-centering and normalization.
- `train_mofa.py`: Handles multi-seed training and stability diagnostics.
- `mofa_wf.ipynb`: Main execution notebook.
