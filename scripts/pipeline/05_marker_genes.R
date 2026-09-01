# ==============================================================================
# scripts/pipeline/05_marker_genes.R
# Cluster marker discovery on UNINFECTED cells only, resolution 0.20.
# Thresholds: adjusted p-value < 0.01 & |avg_log2FC| > 0.05
# Output: marker gene TSV + annotated XLSX (feeds Fig 4 & 10 gene lists)
# ==============================================================================

library(Seurat)
library(tidyverse)
library(readxl)
library(openxlsx)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
integrated_seurat_object_path <- here(
  "data", "processed", "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds"
)
clustering_res_path   <- here("results", "markers")
annotation_table_path <- here("data", "metadata", "consolidated_annotations.tsv")
dir.create(clustering_res_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load object and subset to uninfected samples
# ------------------------------------------------------------------------------
integrated_seurat_object <- readRDS(integrated_seurat_object_path)

uninfected_samples <- c("spl1_Sept22", "spl2_Sept22", "spl1_Nov22",
                        "spl5_Sept22", "spl6_Sept22", "spl5_Nov22")

uninfected_integrated_seurat_object <- subset(
  integrated_seurat_object, subset = Sample %in% uninfected_samples
)

# ------------------------------------------------------------------------------
# Find markers (Wilcoxon), resolution 0.20
# ------------------------------------------------------------------------------
fdr_thresh    <- 0.01
log2fc_thresh <- 0.05
resolutions   <- 0.20

uninfected_integrated_seurat_object <- PrepSCTFindMarkers(uninfected_integrated_seurat_object)

for (i in seq_along(resolutions)) {
  Idents(uninfected_integrated_seurat_object) <- paste0("SCT_snn_res.", resolutions[i])

  markers_all_clusters <- FindAllMarkers(
    uninfected_integrated_seurat_object,
    only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5
  ) %>%
    filter(p_val_adj < fdr_thresh & abs(avg_log2FC) > log2fc_thresh) %>%
    dplyr::select(gene, cluster, avg_log2FC, p_val_adj, `pct.1`, `pct.2`)

  write_tsv(markers_all_clusters,
            file.path(clustering_res_path,
                      paste0("marker_genes_uninfected_only_resolution_",
                             resolutions[i], ".tsv")))
}

# ------------------------------------------------------------------------------
# Annotate and export XLSX
# ------------------------------------------------------------------------------
markers_all_clusters <- markers_all_clusters %>%
  mutate(Resolution = 0.20) %>%
  filter(p_val_adj < fdr_thresh & abs(avg_log2FC) > log2fc_thresh) %>%
  left_join(read_tsv(annotation_table_path, show_col_types = FALSE)) %>%
  distinct(gene, cluster, .keep_all = TRUE)

write.xlsx(
  list("Res0.20_UninfectedOnly" = markers_all_clusters),
  firstRow = TRUE, colWidths = "auto",
  file.path(clustering_res_path,
            "marker_genes_uninfected_only_resolution_0.2.xlsx")
)
