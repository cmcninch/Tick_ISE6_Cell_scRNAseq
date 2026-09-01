# ==============================================================================
# scripts/pipeline/01_decontX.R
# Estimate ambient RNA contamination per sample using decontX.
# Input : CellRanger filtered + raw feature-barcode matrices
# Output: per-sample decontX estimate TSVs
# ==============================================================================

library(Seurat)
library(tidyverse)
library(celda)
library(SingleCellExperiment)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
cell_ranger_out_paths <- here("results", "CellRanger", "filtered_feature_bc_matrices")
out_dir <- here("results", "decontX")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# decontX function: model ambient RNA using raw matrix as background
# ------------------------------------------------------------------------------
get_decontX_seurat_obj <- function(filtered_h5_file, sample, output_path) {
  counts     <- Read10X_h5(filtered_h5_file)
  counts.raw <- Read10X_h5(gsub("filtered", "raw", filtered_h5_file))

  sce     <- SingleCellExperiment(list(counts = counts))
  sce.raw <- SingleCellExperiment(list(counts = counts.raw))

  sce <- decontX(sce, background = sce.raw)

  colData(sce) %>%
    as.data.frame() %>%
    mutate(Sample = sample) %>%
    rownames_to_column(var = "BarcodeID") %>%
    write.table(file = output_path, quote = FALSE, row.names = FALSE,
                col.names = TRUE, sep = "\t")
}

# ------------------------------------------------------------------------------
# Run across all filtered matrices
# ------------------------------------------------------------------------------
filtered_h5_files <- list.files(cell_ranger_out_paths, full.names = TRUE)

for (i in seq_along(filtered_h5_files)) {
  sample <- basename(gsub("_filtered_feature_bc_matrix.h5", "", filtered_h5_files[i]))
  message(sample)
  get_decontX_seurat_obj(
    filtered_h5_file = filtered_h5_files[i],
    sample           = sample,
    output_path      = file.path(out_dir, paste0(sample, "_decontX_estimates.tsv"))
  )
}
