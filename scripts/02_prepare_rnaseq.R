# scripts/02_prepare_rnaseq.R

library(dplyr)
library(readr)

# Load inputs
metadata <- read_tsv(snakemake@input[["metadata"]])
rnaseq <- read_csv(snakemake@input[["rnaseq"]])

# Join and filter for baseline participants
# Assuming participant_id is the key
rnaseq_filtered <- rnaseq %>%
  filter(participant_id %in% metadata$participant_id)

# Dummy VST and selection (in a real pipeline we'd use DESeq2/vst)
rnaseq_vst <- rnaseq_filtered %>%
  mutate(across(-participant_id, ~ log1p(.))) # Mock VST

# Save output
write_tsv(rnaseq_vst, snakemake@output[["rnaseq_vst"]])
