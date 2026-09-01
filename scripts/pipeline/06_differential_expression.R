# ==============================================================================
# scripts/pipeline/06_differential_expression.R
# Pseudobulk DESeq2 DEGs for all 4 comparisons, per-cluster (res 0.20) +
# "All Cells". Writes unfiltered tables; figure scripts apply display filters.
# Diff Exp cutoffs (for figures): FDR-adj p < 0.01 & |log2FC| > 1
# Output: results/degs/*_pseudobulk_res0.20.tsv
# ==============================================================================

library(Seurat)
library(tidyverse)
library(readxl)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
integrated_seurat_object_path <- here(
  "data", "processed", "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds"
)
out_dir <- here("results", "degs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

EMEd2_vs_UId2_path    <- file.path(out_dir, "EMEd2_vs_UId2_pseudobulk_res0.20.tsv")
EMEd4_vs_UId4_path    <- file.path(out_dir, "EMEd4_vs_UId4_pseudobulk_res0.20.tsv")
UI_d4_vs_UI_d2_path   <- file.path(out_dir, "UI_d4_vs_UI_d2_pseudobulk_res0.20.tsv")
EME_d4_vs_EME_d2_path <- file.path(out_dir, "EME_d4_vs_EME_d2_pseudobulk_res0.20.tsv")

integrated_seurat_object <- readRDS(integrated_seurat_object_path)

# ------------------------------------------------------------------------------
# Pseudobulk DEG function: "All Cells" + per-cluster (res 0.20)
# NOTE: fdr_thresh / log2fc_thresh set permissively here to export full tables;
#       display thresholds (0.01 / 1) are applied downstream in figure scripts.
# ------------------------------------------------------------------------------
compute_degs_func <- function(ident_1, ident_2, fdr_thresh, log2fc_thresh, out_path) {

  # All-cells pseudobulk
  pb_all <- AggregateExpression(
    integrated_seurat_object, assays = "RNA", return.seurat = TRUE,
    group.by = c("Sample", "Sample_Group")
  )
  Idents(pb_all) <- "Sample_Group"

  degs <- FindMarkers(pb_all, ident.1 = ident_1, ident.2 = ident_2,
                      test.use = "DESeq2") %>%
    rownames_to_column(var = "Gene") %>%
    mutate(Cluster = "All Cells") %>%
    dplyr::select(Cluster, Gene, p_val_adj, avg_log2FC, pct.1, pct.2)

  # Per-cluster pseudobulk
  pb_cluster <- AggregateExpression(
    integrated_seurat_object, assays = "RNA", return.seurat = TRUE,
    group.by = c("Sample", "Sample_Group", "SCT_snn_res.0.2")
  )
  Idents(pb_cluster) <- "Sample_Group"

  clusters_to_compare <- 0:max(as.numeric(pb_cluster$SCT_snn_res.0.2))

  for (i in seq_along(clusters_to_compare)) {
    temp_cluster <- clusters_to_compare[i]
    tryCatch({
      temp_degs <- FindMarkers(
        object = subset(pb_cluster, subset = `SCT_snn_res.0.2` == temp_cluster),
        ident.1 = ident_1, ident.2 = ident_2, test.use = "DESeq2"
      ) %>%
        rownames_to_column(var = "Gene") %>%
        mutate(Cluster = as.character(temp_cluster)) %>%
        dplyr::select(Cluster, Gene, p_val_adj, avg_log2FC, pct.1, pct.2)
    }, error = function(e) {})

    degs <- bind_rows(degs, temp_degs) %>%
      filter(p_val_adj < fdr_thresh & abs(avg_log2FC) > log2fc_thresh)
  }

  colnames(degs) <- c(
    "Cluster", "Gene", "P_Value_Adj",
    paste0("Avg_Log2FC (pos=higher in ", ident_1, ")"),
    paste0("Pct_Cells_Exp_in_", ident_1),
    paste0("Pct_Cells_Exp_in_", ident_2)
  )

  write_tsv(degs, out_path)
}

# ------------------------------------------------------------------------------
# Run all four comparisons (fdr_thresh = 1, log2fc_thresh = 0 -> full tables)
# ------------------------------------------------------------------------------
fdr_thresh    <- 1
log2fc_thresh <- 0

compute_degs_func("EME-d2", "UI-d2",  fdr_thresh, log2fc_thresh, EMEd2_vs_UId2_path)
compute_degs_func("EME-d4", "UI-d4",  fdr_thresh, log2fc_thresh, EMEd4_vs_UId4_path)
compute_degs_func("UI-d4",  "UI-d2",  fdr_thresh, log2fc_thresh, UI_d4_vs_UI_d2_path)
compute_degs_func("EME-d4", "EME-d2", fdr_thresh, log2fc_thresh, EME_d4_vs_EME_d2_path)
