#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Parse Command-Line Options (Base R)
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(expression = NULL, metadata = NULL, output = "adjusted_rna_matrix.csv", plot_prefix = NULL)
  for (i in seq_along(args)) {
    if (args[i] %in% c("-e", "--expression"))  opts$expression  <- args[i + 1]
    if (args[i] %in% c("-m", "--metadata"))    opts$metadata    <- args[i + 1]
    if (args[i] %in% c("-o", "--output"))      opts$output      <- args[i + 1]
    if (args[i] %in% c("-p", "--plot_prefix"))  opts$plot_prefix <- args[i + 1]
  }
  return(opts)
}

opt <- parse_args(args)

if (is.null(opt$expression) || is.null(opt$metadata)) {
  stop("Error: Both --expression (-e) and --metadata (-m) are required.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 2. Load & Align Data
# ------------------------------------------------------------------------------
logcpm_matrix <- read.csv(opt$expression, row.names = 1, check.names = FALSE)
metadata      <- read.csv(opt$metadata, row.names = 1, check.names = FALSE)

if (!"participant_id" %in% colnames(metadata)) metadata$participant_id <- rownames(metadata)

common_samples <- intersect(colnames(logcpm_matrix), metadata$participant_id)
if (length(common_samples) == 0) stop("Error: No matching sample IDs found.", call. = FALSE)

metadata <- metadata %>% 
  filter(participant_id %in% common_samples) %>% 
  slice(match(common_samples, participant_id))

logcpm_matrix <- as.matrix(logcpm_matrix[, metadata$participant_id])

metadata$plate     <- as.factor(metadata$plate)
metadata$sex       <- as.factor(metadata$sex)
metadata$diagnosis <- as.factor(metadata$diagnosis)

# ------------------------------------------------------------------------------
# 3. Batch Adjustment
# ------------------------------------------------------------------------------
design_bio <- model.matrix(~ sex + diagnosis, data = metadata)

adjusted_rna_matrix <- removeBatchEffect(
  x      = logcpm_matrix,
  batch  = metadata$plate,
  design = design_bio
)

write.csv(adjusted_rna_matrix, opt$output, quote = FALSE)
cat(sprintf("[SUCCESS] Saved adjusted matrix to %s\n", opt$output))

# ------------------------------------------------------------------------------
# 4. Optional PCA Plots
# ------------------------------------------------------------------------------
if (!is.null(opt$plot_prefix)) {
  get_pca_df <- function(mat, label) {
    pca <- prcomp(t(mat), scale. = TRUE)
    data.frame(
      PC1 = pca$x[, 1], PC2 = pca$x[, 2],
      Plate = metadata$plate, Diagnosis = metadata$diagnosis,
      Stage = factor(label, levels = c("Before Adjustment", "After Adjustment"))
    )
  }

  pca_df <- rbind(get_pca_df(logcpm_matrix, "Before Adjustment"), 
                  get_pca_df(adjusted_rna_matrix, "After Adjustment"))

  pdf(paste0(opt$plot_prefix, "_pca.pdf"), width = 10, height = 4.5)
  print(ggplot(pca_df, aes(x = PC1, y = PC2, color = Plate)) + geom_point(size = 3, alpha = 0.8) + facet_wrap(~Stage) + theme_bw() + labs(title = "PCA colored by Plate"))
  print(ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis)) + geom_point(size = 3, alpha = 0.8) + facet_wrap(~Stage) + theme_bw() + labs(title = "PCA colored by Diagnosis"))
  dev.off()
  
  cat(sprintf("[SUCCESS] Saved PCA plots to %s_pca.pdf\n", opt$plot_prefix))
}