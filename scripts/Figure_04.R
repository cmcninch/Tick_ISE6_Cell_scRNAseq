# ==============================================================================
# Figure_04.R
# Dot plots of putative transcriptional markers across clusters
# Panels: 4A (uninfected only), 4B (EME-infected only)
# Markers computed from uninfected cells, resolution 0.2
# ==============================================================================

library(Seurat)
library(scCustomize)
library(tidyverse)
library(readxl)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
integrated_seurat_object_path <- here(
  "data", "processed",
  "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds"
)
marker_genes_path <- here(
  "data", "metadata", "Genes_to_Use_for_Marker_Genes_Plot_2025_12_15.xlsx"
)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
integrated_seurat_object <- readRDS(integrated_seurat_object_path)

`%nin%` <- Negate(`%in%`)

marker_genes <- read_xlsx(marker_genes_path) %>%
  dplyr::select(gene:product) %>%
  mutate(cluster = factor(cluster, levels = 0:14)) %>%
  pull(gene) %>%
  unique() %>%
  as.character()

# ------------------------------------------------------------------------------
# Define sample groupings
# ------------------------------------------------------------------------------
uninfected_samples <- c(
  "spl1_Sept22", "spl2_Sept22", "spl1_Nov22",
  "spl5_Sept22", "spl6_Sept22", "spl5_Nov22"
)

uninfected_integrated_seurat_object <- subset(
  integrated_seurat_object, subset = Sample %in% uninfected_samples
)
infected_integrated_seurat_object <- subset(
  integrated_seurat_object, subset = Sample %nin% uninfected_samples
)

# ------------------------------------------------------------------------------
## Panel A: Marker dot plot -- uninfected samples only
# ------------------------------------------------------------------------------
Idents(uninfected_integrated_seurat_object) <- "SCT_snn_res.0.2"

panel_A <- DotPlot(
  uninfected_integrated_seurat_object,
  features = marker_genes,
  cols     = c("lightgrey", "firebrick1"),
  col.min  = 0,
  col.max  = 2
) +
  RotatedAxis() +
  labs(x = NULL, y = "Cluster Identity", title = "UI Samples Only") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 16, hjust = 0.5))

print(panel_A)

# ------------------------------------------------------------------------------
## Panel B: Marker dot plot -- EME-infected samples only
# ------------------------------------------------------------------------------
Idents(infected_integrated_seurat_object) <- "SCT_snn_res.0.2"

panel_B <- DotPlot(
  infected_integrated_seurat_object,
  features = marker_genes,
  cols     = c("lightgrey", "firebrick1"),
  col.min  = 0,
  col.max  = 2
) +
  RotatedAxis() +
  labs(x = NULL, y = "Cluster Identity", title = "EME Samples Only") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 16, hjust = 0.5))

print(panel_B)
