# ==============================================================================
# Figure_09.R
# Transcription factor expression across clusters (UI vs EME)
# Panels: 9A (TF heatmap), 9B (UMAP of composite transcription score),
#         9C (violin plot of transcription score by cluster)
# ==============================================================================

library(scCustomize)
library(Seurat)
library(tidyverse)
library(readxl)
library(pals)
library(RColorBrewer)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
integrated_seurat_object_path <- here(
  "data", "processed",
  "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds"
)
transcriptional_genes_path <- here("data", "metadata", "Transcriptional markers.xlsx")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
integrated_seurat_object <- readRDS(integrated_seurat_object_path)
transcriptional_genes <- read_xlsx(transcriptional_genes_path) %>% pull(Symbol)

# ------------------------------------------------------------------------------
## Panel A: Heatmap of transcription-factor expression (150 random cells)
# ------------------------------------------------------------------------------
Idents(integrated_seurat_object) <- "SCT_snn_res.0.2"

random_cells_150 <- Random_Cells_Downsample(
  seurat_object = integrated_seurat_object, num_cells = 150, allow_lower = TRUE
)

scaled_genes  <- rownames(GetAssayData(integrated_seurat_object, layer = "scale.data"))
genes_to_plot <- intersect(transcriptional_genes, scaled_genes)

mat <- GetAssayData(integrated_seurat_object, layer = "scale.data")[genes_to_plot,
                                                                    random_cells_150]

# Order TF rows by expression similarity
gene_clust    <- hclust(dist(mat))
ordered_genes <- gene_clust$labels[gene_clust$order]

panel_A <- DoHeatmap(
  object = integrated_seurat_object, features = ordered_genes,
  cells = random_cells_150, angle = 0, hjust = 0.5
) +
  ggtitle("Cluster Identity") + ylab("Transcription Factor") +
  theme(axis.text.y = element_text(size = 0),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  guides(color = "none")

print(panel_A)

# ------------------------------------------------------------------------------
# Composite transcription score (used for Panels B & C)
# ------------------------------------------------------------------------------
transcriptional_genes_final <- intersect(
  transcriptional_genes, rownames(integrated_seurat_object@assays$SCT@counts)
)

integrated_seurat_object <- AddModuleScore(
  object = integrated_seurat_object,
  features = list(transcriptional_genes_final),
  name = "Transcription_Score", nbin = 5
)

integrated_seurat_object@meta.data <- integrated_seurat_object@meta.data %>%
  rownames_to_column(var = "BarcodeID_Long") %>%
  mutate(
    Infection_Status = ifelse(str_detect(integrated_seurat_object$Sample_Group, "EME"), "EME", "UI"),
    Sample_Group_UI_Together = case_when(
      str_detect(Sample_Group, "EME") ~ Sample_Group,
      str_detect(Sample_Group, "UI") ~ "UI")
  ) %>%
  column_to_rownames(var = "BarcodeID_Long")

# ------------------------------------------------------------------------------
## Panel B: UMAP of composite transcription score (UI / EME_d2 / EME_d4)
# ------------------------------------------------------------------------------
integrated_seurat_object$TEMP_Sample_Group_UI_Together <-
  factor(integrated_seurat_object$Sample_Group_UI_Together,
         levels = c("UI", "EME_d2", "EME_d4"))

panel_B <- FeaturePlot_scCustom(
  integrated_seurat_object, features = "Transcription_Score1",
  reduction = "umap.rpca", split.by = "TEMP_Sample_Group_UI_Together",
  pt.size = 0.0001, colors_use = viridis_plasma_dark_high,
  na_color = "lightgray", combine = TRUE
)

print(panel_B)

# ------------------------------------------------------------------------------
## Panel C: Violin plot of transcription score by cluster (split by group)
# ------------------------------------------------------------------------------
original_colors <- scCustomize_Palette(num_groups = 3)
new_colors <- original_colors[c(3, 1, 2)]   # UI, EME_d2, EME_d4

integrated_seurat_object$Sample_Group_UI_Together <- factor(
  integrated_seurat_object$Sample_Group_UI_Together,
  levels = c("UI", "EME_d2", "EME_d4")
)

panel_C <- VlnPlot_scCustom(
  seurat_object = integrated_seurat_object, features = "Transcription_Score1",
  group.by = "SCT_snn_res.0.2", split.by = "Sample_Group_UI_Together",
  colors_use = new_colors
)

print(panel_C)
