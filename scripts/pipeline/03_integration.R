# ==============================================================================
# scripts/pipeline/03_integration.R
# Filter cells, SCTransform, RPCA integration, clustering, UMAP.
# QC cutoffs: %mito < 5, 500 < nFeature_RNA < 5000, decontX < 0.5,
#             drop sample 'spl4_Sept22'
# Output: integrated Seurat object (.rds) consumed by all figure scripts
# ==============================================================================

library(Seurat)
library(tidyverse)
library(glmGamPoi)
library(here)

# ------------------------------------------------------------------------------
# File paths
#   NOTE: The output .rds filename below MUST match the string every figure
#   script reads. Replace the abbreviated placeholder with one literal name
#   used consistently across the whole repo (see README / confidence report).
# ------------------------------------------------------------------------------
merged_seurat_object_path <- here(
  "data", "processed",
  "iSE6_tick_cell_line_Muris_Eauclairensis_infected_unfiltered_merged_seurat_object.rds"
)
out_dir <- here("data", "processed")

# ------------------------------------------------------------------------------
# Load and filter for good-quality cells
# ------------------------------------------------------------------------------
merged_seurat_object <- readRDS(merged_seurat_object_path) %>%
  subset(subset = `percent.mt.stringent` < 5 &
                  nFeature_RNA > 500 &
                  nFeature_RNA < 5000 &
                  decontX_contamination < 0.5 &
                  Sample != "spl4_Sept22")

# ------------------------------------------------------------------------------
# SCTransform + PCA
# ------------------------------------------------------------------------------
options(future.globals.maxSize = 150000 * 1024^2)

merged_seurat_object <- SCTransform(merged_seurat_object, assay = "RNA", verbose = FALSE)
merged_seurat_object <- RunPCA(merged_seurat_object, assay = "SCT", npcs = 30, verbose = FALSE)

# ------------------------------------------------------------------------------
# RPCA integration (top 30 PCs)
# ------------------------------------------------------------------------------
integrated_seurat_object <- IntegrateLayers(
  object = merged_seurat_object, method = RPCAIntegration,
  normalization.method = "SCT", orig.reduction = "pca",
  new.reduction = "integrated.rpca", dims = 1:30, verbose = FALSE
)

# ------------------------------------------------------------------------------
# Neighbors, clustering across resolutions, UMAP
# ------------------------------------------------------------------------------
integrated_seurat_object <- FindNeighbors(integrated_seurat_object,
                                          reduction = "integrated.rpca", dims = 1:30)
integrated_seurat_object <- FindClusters(integrated_seurat_object,
                                         resolution = seq(0.05, 0.25, 0.05))
integrated_seurat_object <- RunUMAP(integrated_seurat_object,
                                    reduction = "integrated.rpca",
                                    reduction.name = "umap.rpca",
                                    dims = 1:30, return.model = TRUE, verbose = FALSE)

# ------------------------------------------------------------------------------
# Save integrated object (canonical name used by every figure script)
# ------------------------------------------------------------------------------
saveRDS(
  integrated_seurat_object,
  file.path(out_dir,
            "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds")
)
