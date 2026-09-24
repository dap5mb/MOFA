# scripts/03_prepare_proteomics.R

library(dplyr)
library(readr)

# Load inputs
metadata <- read_tsv(snakemake@input[["metadata"]])
proteomics <- read_csv(snakemake@input[["proteomics"]])

# Join and filter for baseline participants
proteomics_filtered <- proteomics %>%
  filter(participant_id %in% metadata$participant_id)

# Dummy normalization (e.g., Z-score)
proteomics_norm <- proteomics_filtered %>%
  mutate(across(-participant_id, ~ (. - mean(.)) / sd(.)))

# Save output
write_tsv(proteomics_norm, snakemake@output[["proteomics_norm"]])
