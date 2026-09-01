# ==============================================================================
# scripts/pipeline/02_merge_qc_seurat.R
# Build merged (unfiltered) Seurat object with %mito, decontX, sample metadata.
# Also writes aggregated per-sample counts + QC metadata table.
# Output: unfiltered merged Seurat object (.rds)
# ==============================================================================

library(Seurat)
library(tidyverse)
library(rtracklayer)
library(readxl)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
h5_file_dir          <- here("results", "CellRanger", "filtered_feature_bc_matrices")
gtf_file_path        <- here("data", "genome_files", "genes.gtf.gz")
decontX_results_path <- here("results", "decontX")
sample_metadata_path <- here("data", "metadata", "Samples.xlsx")
out_dir              <- here("data", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Merge all sample matrices into one Seurat object
# ------------------------------------------------------------------------------
h5_files <- list.files(h5_file_dir, pattern = "*h5", full.names = TRUE)
samples  <- gsub("_filtered_feature_bc_matrix.h5", "", basename(h5_files))

h5_read   <- lapply(h5_files, Read10X_h5)
h5_seurat <- lapply(h5_read, CreateSeuratObject)
merged_seurat_object <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)],
                              add.cell.ids = samples)

merged_seurat_object$Sample <- str_extract(colnames(merged_seurat_object),
                                            "^spl[0-9]+_[A-Za-z0-9]+")
merged_seurat_object$BarcodeID <- str_replace(colnames(merged_seurat_object),
                                              "^spl[0-9]+_[A-Za-z0-9]+_", "")

# ------------------------------------------------------------------------------
# % mitochondrial gene expression (stringent + lenient definitions)
# ------------------------------------------------------------------------------
mitochondrial_genes_stringent <- rtracklayer::import(gtf_file_path) %>%
  as.data.frame() %>%
  distinct(gene_id, gene, product) %>%
  filter(str_detect(product, ", mitochondrial")) %>%
  filter(!str_detect(product, ", mitochondrial-like")) %>%
  pull(gene) %>% unique()

mitochondrial_genes_lenient <- rtracklayer::import(gtf_file_path) %>%
  as.data.frame() %>%
  distinct(gene_id, gene, product) %>%
  filter(str_detect(product, ", mitochondrial")) %>%
  pull(gene) %>% unique()

merged_seurat_object[["percent.mt.stringent"]] <-
  PercentageFeatureSet(merged_seurat_object, features = mitochondrial_genes_stringent)
merged_seurat_object[["percent.mt.lenient"]] <-
  PercentageFeatureSet(merged_seurat_object, features = mitochondrial_genes_lenient)

# ------------------------------------------------------------------------------
# Join decontX contamination estimates
# ------------------------------------------------------------------------------
combined_decontX_estimates <- do.call(
  bind_rows,
  lapply(list.files(decontX_results_path, pattern = "*tsv", full.names = TRUE),
         read_tsv, show_col_types = FALSE)
) %>%
  unite("BarcodeID", Sample, BarcodeID, sep = "_")

merged_seurat_object@meta.data <- merged_seurat_object@meta.data %>%
  rownames_to_column(var = "BarcodeID_Full") %>%
  left_join(combined_decontX_estimates, by = c("BarcodeID_Full" = "BarcodeID")) %>%
  column_to_rownames(var = "BarcodeID_Full")

# ------------------------------------------------------------------------------
# Join remaining sample metadata
# ------------------------------------------------------------------------------
merged_seurat_object@meta.data <- merged_seurat_object@meta.data %>%
  rownames_to_column(var = "BarcodeID_Full") %>%
  left_join(
    read_xlsx(sample_metadata_path) %>% dplyr::select(Sample_Name, Control:Replicate),
    by = c("Sample" = "Sample_Name")
  ) %>%
  column_to_rownames(var = "BarcodeID_Full")

# ------------------------------------------------------------------------------
# Aggregated per-sample counts (for downstream QC / PCA) + QC metadata
# ------------------------------------------------------------------------------
AggregateExpression(
  merged_seurat_object, assays = "RNA", slot = "counts",
  return.seurat = FALSE, group.by = c("Sample")
)$RNA %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene_ID") %>%
  write_tsv(file.path(out_dir, "aggregated_gene_counts_per_sample.tsv"))

merged_seurat_object@meta.data %>%
  write_tsv(file.path(out_dir, "iSE6_tick_cell_line_Muris_Eauclairensis_infected_metadata.tsv"))

# ------------------------------------------------------------------------------
# Save unfiltered merged object
# ------------------------------------------------------------------------------
saveRDS(
  merged_seurat_object,
  file.path(out_dir,
            "iSE6_tick_cell_line_Muris_Eauclairensis_infected_unfiltered_merged_seurat_object.rds")
)
