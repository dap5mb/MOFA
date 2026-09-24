preprocess_rnaseq <- function(counts_file, metadata_file, visit_filter = 0.0, cohort_name = "PPMI", min_count = 10, min_prop = 0.20) {
  suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
  })
  
  cat("Loading data...\n")
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  colData <- read.csv(metadata_file, row.names = 1)
  
  # 1. Filter metadata strictly to the target cohort (e.g., PPMI)
  if ("study" %in% colnames(colData)) {
    colData <- colData %>% filter(study == cohort_name)
    cat(sprintf("Filtered metadata to %d samples belonging to study: %s\n", nrow(colData), cohort_name))
  }
  
  # 2. Intersect sample IDs between counts and the filtered metadata
  common_samples <- intersect(colnames(counts), rownames(colData))
  counts <- counts[, common_samples]
  colData <- colData[common_samples, ]
  
  # 3. Filter for baseline/specified visit month if available
  if ("visit_month" %in% colnames(colData)) {
    colData <- colData %>% filter(visit_month == visit_filter)
    counts <- counts[, rownames(colData)]
    cat(sprintf("Filtered to %d samples where visit_month == %.1f\n", ncol(counts), visit_filter))
  }
  
  # 4. Construct DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData = colData,
                                design = ~ 1)
  
  # 5. Filter Lowly Expressed Genes
  min_samples <- ceiling(min_prop * ncol(dds))
  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  dds <- dds[keep, ]
  cat(sprintf("Retained %d genes out of %d total genes after low-count filtering.\n", 
              nrow(dds), length(keep)))
  
  # 6. Apply Variance-Stabilizing Transformation (VST)
  cat("Applying Variance-Stabilizing Transformation (VST)...\n")
  vsd <- vst(dds, blind = TRUE)
  
  expr_matrix_vst <- assay(vsd)
  cat("Preprocessing complete!\n")
  
  write.csv(expr_matrix_vst, file = "ppmi_rnaseq_vst_filtered.csv")
  cat("Output saved to ppmi_rnaseq_vst_filtered.csv\n")
}

# --- Command Line Argument Handling ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  counts_arg <- args[1]
  metadata_arg <- args[2]
  cohort_arg <- if (length(args) >= 3) args[3] else "PPMI"
  
  preprocess_rnaseq(counts_file = counts_arg, 
                         metadata_file = metadata_arg, 
                         cohort_name = cohort_arg)
} else {
  cat("Error: Missing required arguments.\n")
  cat("Usage: Rscript preprocess_rnaseq.R <counts_file> <metadata_file> [cohort_name]\n")
}