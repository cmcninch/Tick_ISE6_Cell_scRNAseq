# ==============================================================================
# Figure_03.R
# Defining transcriptional heterogeneity of ISE6 cells (UMAPs)
# Panels: 3A (UMAP by sample group), 3B (UMAP by 15 clusters, res 0.2)
# ==============================================================================

library(scCustomize)
library(Seurat)
library(tidyverse)
library(pals)
library(readxl)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
integrated_seurat_object_path <- here(
  "data", "processed",
  "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds"
)

# ------------------------------------------------------------------------------
# Load integrated Seurat object and derive infection-status metadata
# ------------------------------------------------------------------------------
integrated_seurat_object <- readRDS(integrated_seurat_object_path)

integrated_seurat_object@meta.data <- integrated_seurat_object@meta.data %>%
  rownames_to_column(var = "BarcodeID_Long") %>%
  mutate(
    Infection_Status = ifelse(
      str_detect(integrated_seurat_object$Sample_Group, "EME"), "EME", "UI"
    ),
    Sample_Group_UI_Together = case_when(
      str_detect(Sample_Group, "EME") ~ Sample_Group,
      str_detect(Sample_Group, "UI") ~ "UI"
    )
  ) %>%
  column_to_rownames(var = "BarcodeID_Long")

# ------------------------------------------------------------------------------
## Panel A: UMAP of unsupervised clustering (by sample group)
# ------------------------------------------------------------------------------
panel_A <- DimPlot(
  integrated_seurat_object,
  group.by  = "Sample_Group",
  reduction = "umap.rpca",
  pt.size   = 0.0001,
  alpha     = 0.5
) +
  labs(title = NULL)

print(panel_A)

# ------------------------------------------------------------------------------
## Panel B: UMAP of the 15 transcriptional clusters (resolution 0.2)
# ------------------------------------------------------------------------------
panel_B <- DimPlot(
  integrated_seurat_object,
  group.by   = "SCT_snn_res.0.2",
  split.by   = "Infection_Status",
  reduction  = "umap.rpca",
  pt.size    = 0.0001,
  cols       = cols25(16),
  alpha      = 0.75,
  label      = TRUE,
  label.size = 5
) +
  labs(title = NULL)

print(panel_B)
