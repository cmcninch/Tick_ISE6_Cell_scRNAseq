# ==============================================================================
# Figure_05.R
# Global transcriptional response of ISE6 cells across EME infection
# Panels: 5A (UpSet of shared DEGs), 5B (EME_d2 vs UI_d2 heatmap),
#         5C (EME_d4 vs UI_d4 heatmap), 5D (EME_d2 vs EME_d4 heatmap)
# Diff Exp: FDR-adj p < 0.01 & |log2FC| > 1  (pseudobulk, DESeq2)
# ==============================================================================

library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(edgeR)
library(circlize)
library(scales)
library(UpSetR)
library(here)

# ------------------------------------------------------------------------------
# File paths (pre-computed pseudobulk DEG tables, res 0.20)
# ------------------------------------------------------------------------------
integrated_seurat_object_path <- here(
  "data", "processed",
  "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds"
)
EMEd2_vs_UId2_path    <- here("results", "degs", "EMEd2_vs_UId2_pseudobulk_res0.20.tsv")
EMEd4_vs_UId4_path    <- here("results", "degs", "EMEd4_vs_UId4_pseudobulk_res0.20.tsv")
UI_d4_vs_UI_d2_path   <- here("results", "degs", "UI_d4_vs_UI_d2_pseudobulk_res0.20.tsv")
EME_d4_vs_EME_d2_path <- here("results", "degs", "EME_d4_vs_EME_d2_pseudobulk_res0.20.tsv")

# ------------------------------------------------------------------------------
# Load integrated object
# ------------------------------------------------------------------------------
integrated_seurat_object <- readRDS(integrated_seurat_object_path)

# ------------------------------------------------------------------------------
## Panel A: UpSet plot of shared DEGs ("All Cells" comparisons)
# ------------------------------------------------------------------------------
consolidated_degs <- bind_rows(
  read_tsv(EMEd2_vs_UId2_path, show_col_types = FALSE) %>%
    dplyr::rename("Avg_Log2FC" = `Avg_Log2FC (pos=higher in EME-d2)`) %>%
    dplyr::select(Cluster, Gene, P_Value_Adj, Avg_Log2FC) %>%
    mutate(Comparison = "EME_d2_vs_UI_d2"),
  read_tsv(EMEd4_vs_UId4_path, show_col_types = FALSE) %>%
    dplyr::rename("Avg_Log2FC" = `Avg_Log2FC (pos=higher in EME-d4)`) %>%
    dplyr::select(Cluster, Gene, P_Value_Adj, Avg_Log2FC) %>%
    mutate(Comparison = "EME_d4_vs_UI_d4"),
  read_tsv(UI_d4_vs_UI_d2_path, show_col_types = FALSE) %>%
    dplyr::rename("Avg_Log2FC" = `Avg_Log2FC (pos=higher in UI-d4)`) %>%
    dplyr::select(Cluster, Gene, P_Value_Adj, Avg_Log2FC) %>%
    mutate(Comparison = "UI_d4_vs_UI_d2"),
  read_tsv(EME_d4_vs_EME_d2_path, show_col_types = FALSE) %>%
    dplyr::rename("Avg_Log2FC" = `Avg_Log2FC (pos=higher in EME-d4)`) %>%
    dplyr::select(Cluster, Gene, P_Value_Adj, Avg_Log2FC) %>%
    mutate(Comparison = "EME_d4_vs_EME_d2")
)

upset_matrix <- consolidated_degs %>%
  filter(Cluster == "All Cells") %>%
  mutate(DEG = 1) %>%
  dplyr::select(Gene, Comparison, DEG) %>%
  pivot_wider(names_from = Comparison, values_from = DEG, values_fill = 0) %>%
  as.data.frame()

upset(upset_matrix)

# ------------------------------------------------------------------------------
# Recompute pseudobulk DEG gene sets + TMM matrix for heatmap panels
# ------------------------------------------------------------------------------
pseudobulk_integrated_seurat_object <- AggregateExpression(
  integrated_seurat_object, assays = "RNA", return.seurat = TRUE,
  group.by = c("Sample", "Sample_Group")
)
Idents(pseudobulk_integrated_seurat_object) <- "Sample_Group"

get_deg_set <- function(id1, id2) {
  FindMarkers(pseudobulk_integrated_seurat_object,
              ident.1 = id1, ident.2 = id2, test.use = "DESeq2") %>%
    filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.01) %>%
    rownames()
}

EMEd2_vs_UId2degs  <- get_deg_set("EME-d2", "UI-d2")
EMEd4_vs_UId4degs  <- get_deg_set("EME-d4", "UI-d4")
EMEd4_vs_EMEd2degs <- get_deg_set("EME-d4", "EME-d2")

# TMM-normalized counts
counts_for_de <- pseudobulk_integrated_seurat_object[["RNA"]]$counts
count_matrix  <- counts_for_de[rowSums(counts_for_de) > 0, ]
tmms <- cpm(calcNormFactors(DGEList(counts = count_matrix), method = "TMM")) %>%
  as.data.frame()

# ------------------------------------------------------------------------------
# Shared heatmap helper (preserves original styling & Z-scoring)
# ------------------------------------------------------------------------------
draw_deg_heatmap <- function(genes_to_plot, group_1, group_2,
                             palette_idx_1, palette_idx_2) {
  metadata <- pseudobulk_integrated_seurat_object@meta.data %>%
    filter(Sample_Group %in% c(group_1, group_2)) %>%
    mutate(Sample_Group = gsub("-", "_", Sample_Group))

  scaled_numeric_matrix <- tmms %>%
    rownames_to_column(var = "gene_id") %>%
    dplyr::select(gene_id, all_of(rownames(metadata))) %>%
    filter(gene_id %in% genes_to_plot) %>%
    column_to_rownames(var = "gene_id") %>%
    as.matrix() %>% t() %>% scale() %>% t() %>%
    as.data.frame() %>%
    dplyr::select(all_of(contains(group_1)), all_of(contains(group_2))) %>%
    as.matrix()

  annotation <- metadata %>% dplyr::select(Sample_Group)
  annotation <- annotation[colnames(scaled_numeric_matrix), , drop = FALSE]

  col_fun <- colorRamp2(c(-2, 0, 2), c("mediumblue", "white", "red3"))
  grp1 <- gsub("-", "_", group_1); grp2 <- gsub("-", "_", group_2)

  ha_cols <- setNames(list(
    setNames(c(hue_pal()(4)[palette_idx_1], hue_pal()(4)[palette_idx_2]),
             c(grp1, grp2))
  ), "Sample_Group")

  heatmap <- Heatmap(
    scaled_numeric_matrix,
    col = col_fun(seq(floor(min(scaled_numeric_matrix)),
                      ceiling(max(scaled_numeric_matrix)), 1)),
    use_raster = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
    show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = TRUE,
    show_column_dend = FALSE,
    heatmap_legend_param = list(
      title = "Z-score", title_gp = gpar(fontsize = 16, fontface = "bold"),
      at = c(-2, -1, 0, 1, 2), labels_gp = gpar(fontsize = 12),
      title_position = "leftcenter-rot"),
    top_annotation = HeatmapAnnotation(
      df = annotation, col = ha_cols, border = TRUE,
      simple_anno_size = unit(1.5, "cm"),
      annotation_name_gp = gpar(fontsize = 0),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 16),
        at = c(grp1, grp2), labels = c(grp1, grp2)))
  )
  draw(heatmap, heatmap_legend_side = "right",
       annotation_legend_side = "right", merge_legend = TRUE)
}

# ------------------------------------------------------------------------------
## Panel B: EME_d2 vs UI_d2 DEG heatmap
# ------------------------------------------------------------------------------
draw_deg_heatmap(EMEd2_vs_UId2degs, "UI-d2", "EME-d2",
                 palette_idx_1 = 3, palette_idx_2 = 1)

# ------------------------------------------------------------------------------
## Panel C: EME_d4 vs UI_d4 DEG heatmap
# ------------------------------------------------------------------------------
draw_deg_heatmap(EMEd4_vs_UId4degs, "UI-d4", "EME-d4",
                 palette_idx_1 = 4, palette_idx_2 = 2)

# ------------------------------------------------------------------------------
## Panel D: EME_d2 vs EME_d4 DEG heatmap
# ------------------------------------------------------------------------------
draw_deg_heatmap(EMEd4_vs_EMEd2degs, "EME-d2", "EME-d4",
                 palette_idx_1 = 1, palette_idx_2 = 2)
