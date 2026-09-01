# ==============================================================================
# Figure_10.R
# Comparison of tissue-specific transcriptional signatures: bulk atlas vs ISE6
# Panels: 10A (tissue-specific marker heatmap across bulk tissues),
#         10B (tissue marker heatmap across 15 ISE6 clusters)
# Tissue-specific: Tau == 1, then DESeq2 filter FDR < 0.01 & log2FC > 1
# ==============================================================================

library(tispec)
library(tidyverse)
library(openxlsx)
library(rtracklayer)
library(readxl)
library(pals)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(DESeq2)
library(scCustomize)
library(Seurat)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
tpms_path                   <- here("data", "bulk_rnaseq", "TPMCountFile_rsemgenes.txt")
tpms_rolandelli_path        <- here("data", "bulk_rnaseq", "RSEM.genes.TPM.all_samples.txt")
gene_counts_path            <- here("data", "bulk_rnaseq", "compiled_RSEM_gene_counts.tsv")
gene_counts_rolandelli_path <- here("data", "bulk_rnaseq", "RSEM.genes.expected_count.all_samples.txt")
metadata_path               <- here("data", "metadata", "Tissue_Expression_Atlas_Metadata_Simplified.xlsx")
integrated_seurat_object_path <- here(
  "data", "processed",
  "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds"
)

# ------------------------------------------------------------------------------
# Load and format metadata (harmonize Ovaries spelling; append hemocytes)
# ------------------------------------------------------------------------------
metadata <- read_xlsx(metadata_path) %>%
  mutate(Tissue = str_extract(Description, "[A-Za-z]+")) %>%
  mutate(
    Description = case_when(
      Description == "Overies_R1" ~ "Ovaries_R1",
      Description == "Overies_R2" ~ "Ovaries_R2",
      Description == "Overies_R3" ~ "Ovaries_R3",
      Description == "Overies_R4" ~ "Ovaries_R4",
      TRUE ~ Description),
    Tissue = case_when(Tissue == "Overies" ~ "Ovaries", TRUE ~ Tissue)
  ) %>%
  mutate(
    Description = gsub("MalpighianTubules", "Malpighian_Tubules", Description),
    Description = gsub("SalivaryGlands",   "Salivary_Glands",   Description),
    Tissue      = gsub("MalpighianTubules", "Malpighian_Tubules", Tissue),
    Tissue      = gsub("SalivaryGlands",   "Salivary_Glands",   Tissue)
  ) %>%
  bind_rows(tibble(
    Sample_ID   = c("Unfed_Bulk1", "Unfed_Bulk2", "Unfed_Bulk3"),
    Description = c("Hemocytes_R1", "Hemocytes_R2", "Hemocytes_R3"),
    Tissue      = "Hemocytes"
  ))

# ------------------------------------------------------------------------------
# Assemble TPM matrix across both datasets (in-house atlas + Rolandelli et al.)
# ------------------------------------------------------------------------------
tpms_set1 <- read_tsv(tpms_path, show_col_types = FALSE) %>%
  pivot_longer(2:ncol(.), names_to = "Sample_ID", values_to = "TPM") %>%
  left_join(metadata, by = "Sample_ID") %>%
  dplyr::select(gene_id, Description, TPM) %>%
  pivot_wider(names_from = "Description", values_from = "TPM")

tpms_set2 <- read_tsv(tpms_rolandelli_path, show_col_types = FALSE) %>%
  dplyr::select(-GeneName, -Fed_Bulk1, -Fed_Bulk2, -Fed_Bulk3) %>%
  pivot_longer(2:ncol(.), names_to = "Sample_ID", values_to = "TPM") %>%
  left_join(metadata, by = "Sample_ID") %>%
  dplyr::select(gene_id, Description, TPM) %>%
  pivot_wider(names_from = "Description", values_from = "TPM")

tpms <- left_join(tpms_set1, tpms_set2, by = "gene_id")

# ------------------------------------------------------------------------------
# Differential expression helper (per-tissue vs all-others)
# ------------------------------------------------------------------------------
diff_exp_func <- function(metadata_for_de, counts_for_de, contrast) {
  counts_for_de <- counts_for_de[rowSums(counts_for_de) > 0, ]
  dds <- DESeqDataSetFromMatrix(countData = counts_for_de,
                                colData = metadata_for_de,
                                design = ~ Tissue)
  dds <- DESeq(dds)
  as.data.frame(results(dds, contrast = contrast, alpha = 0.1)) %>%
    rownames_to_column(var = "gene_id")
}

# ------------------------------------------------------------------------------
# Identify tissue-specific genes (Tau == 1 + DESeq2 enrichment)
# ------------------------------------------------------------------------------
identify_tissue_spec_genes <- function(samples) {
  gene_counts_set1 <- read_tsv(gene_counts_path, show_col_types = FALSE) %>%
    dplyr::rename("Ovaries_R1" = "Overies_R1", "Ovaries_R2" = "Overies_R2",
                  "Ovaries_R3" = "Overies_R3", "Ovaries_R4" = "Overies_R4")
  gene_counts_set2 <- read_tsv(gene_counts_rolandelli_path, show_col_types = FALSE) %>%
    dplyr::rename("Hemocytes_R1" = Unfed_Bulk1, "Hemocytes_R2" = Unfed_Bulk2,
                  "Hemocytes_R3" = Unfed_Bulk3) %>%
    dplyr::select(-GeneName)

  gene_counts <- left_join(gene_counts_set1, gene_counts_set2,
                           by = c("Gene_ID" = "gene_id")) %>%
    dplyr::select(Gene_ID, all_of(samples)) %>%
    distinct(Gene_ID, .keep_all = TRUE) %>%
    mutate_if(is.numeric, as.integer) %>%
    column_to_rownames(var = "Gene_ID")

  tissues <- unique(metadata$Tissue)
  diff_exp_results <- tibble()
  for (i in seq_along(tissues)) {
    temp_metadata <- metadata %>%
      filter(Description %in% samples) %>%
      mutate(Tissue = case_when(Tissue == tissues[i] ~ "Current_Tissue",
                                TRUE ~ "Other_Tissues")) %>%
      column_to_rownames(var = "Description")
    diff_exp_results <- bind_rows(
      diff_exp_results,
      diff_exp_func(temp_metadata, gene_counts,
                    c("Tissue", "Current_Tissue", "Other_Tissues")) %>%
        mutate(Tissue = tissues[i]) %>%
        dplyr::select(gene_id, Tissue, padj, log2FoldChange) %>%
        dplyr::rename("Diff_Exp_FDR" = padj, "Diff_Exp_log2FC" = log2FoldChange)
    )
  }

  mean_tpms <- tpms %>%
    pivot_longer(2:ncol(.), names_to = "Description", values_to = "TPM") %>%
    left_join(metadata, by = "Description") %>%
    filter(Description %in% samples) %>%
    group_by(gene_id, Tissue) %>%
    summarise(Mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Tissue, values_from = Mean_TPM)

  tpms %>%
    pivot_longer(2:ncol(.), names_to = "Description", values_to = "TPM") %>%
    left_join(metadata, by = "Description") %>%
    filter(Description %in% samples) %>%
    group_by(gene_id, Tissue) %>%
    summarise(Mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Tissue, values_from = Mean_TPM) %>%
    column_to_rownames(var = "gene_id") %>%
    tispec::log2Tran() %>% tispec::quantNorm() %>% calcTau() %>%
    rownames_to_column(var = "gene_id") %>%
    dplyr::select(-tau) %>%
    pivot_longer(2:ncol(.), names_to = "Tissue", values_to = "Tau") %>%
    filter(Tau == 1) %>%
    left_join(mean_tpms, by = "gene_id") %>%
    rowwise() %>%
    mutate(Max_Expression_Across_Tissues = max(
      Hemocytes, Malpighian_Tubules, Midgut, Other,
      Ovaries, Salivary_Glands, Synganglion, Trachia)) %>%
    dplyr::select(gene_id, Tissue, Max_Expression_Across_Tissues, Hemocytes:Trachia) %>%
    dplyr::rename("Gene_ID" = gene_id, "Tissue_Gene_is_Specific_To" = Tissue) %>%
    arrange(desc(Max_Expression_Across_Tissues)) %>%
    left_join(diff_exp_results,
              by = c("Gene_ID" = "gene_id",
                     "Tissue_Gene_is_Specific_To" = "Tissue"))
}

tissue_specific_genes_using_all_samples <- identify_tissue_spec_genes(
  samples = metadata %>% pull(Description)
)

# ------------------------------------------------------------------------------
## Panel A: Heatmap of tissue-specific markers across bulk tissues
# ------------------------------------------------------------------------------
samples <- metadata %>% pull(Description) %>% sort()

genes_to_plot <- tissue_specific_genes_using_all_samples %>%
  filter(Diff_Exp_log2FC > 1 & Diff_Exp_FDR < 0.01) %>%
  arrange(Tissue_Gene_is_Specific_To) %>%
  pull(Gene_ID)

scaled_numeric_matrix <- tpms %>%
  dplyr::select(gene_id, all_of(samples)) %>%
  filter(gene_id %in% genes_to_plot) %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix() %>% t() %>% scale() %>% t()
scaled_numeric_matrix <- scaled_numeric_matrix[genes_to_plot, ]

annotation <- tibble(Description = colnames(scaled_numeric_matrix)) %>%
  left_join(metadata, by = "Description") %>%
  dplyr::select(Description, Tissue) %>%
  column_to_rownames(var = "Description")

col_fun <- colorRamp2(c(-1, 0, 1), c("mediumblue", "white", "red3"))

panel_A <- Heatmap(
  scaled_numeric_matrix,
  col = col_fun(seq(floor(min(scaled_numeric_matrix)),
                    ceiling(max(scaled_numeric_matrix)), 1)),
  use_raster = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
  show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
  row_title = "Tissue-specific marker genes",
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_title = "Tissues", column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Z-score", title_gp = gpar(fontsize = 12, fontface = "bold"),
    legend_direction = "horizontal", title_position = "topcenter",
    at = c(-3, -2, -1, 0, 1, 2, 3), labels_gp = gpar(fontsize = 10)),
  top_annotation = HeatmapAnnotation(
    df = annotation,
    col = list(Tissue = c(
      Midgut = okabe(8)[1], Salivary_Glands = okabe(8)[2],
      Synganglion = okabe(8)[3], Malpighian_Tubules = okabe(8)[4],
      Ovaries = okabe(8)[5], Trachia = okabe(8)[6],
      Other = okabe(8)[7], Hemocytes = okabe(8)[8])),
    border = TRUE, simple_anno_size = unit(0.75, "cm"),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 12)))
)

draw(panel_A, heatmap_legend_side = "right",
     annotation_legend_side = "right", merge_legends = TRUE)

# ------------------------------------------------------------------------------
## Panel B: Heatmap of tissue-specific markers across 15 ISE6 clusters
# ------------------------------------------------------------------------------
integrated_seurat_object <- readRDS(integrated_seurat_object_path)
Idents(integrated_seurat_object) <- "SCT_snn_res.0.2"

integrated_seurat_object_scaled <- ScaleData(
  object = integrated_seurat_object, features = genes_to_plot
)

random_cells_150 <- Random_Cells_Downsample(
  seurat_object = integrated_seurat_object_scaled,
  num_cells = 150, allow_lower = TRUE
)

panel_B <- DoHeatmap(
  object = integrated_seurat_object_scaled, features = genes_to_plot,
  cells = random_cells_150, group.colors = cols25(16),
  angle = 0, hjust = 0.5
) +
  theme(axis.text.y = element_text(size = 0)) +
  guides(color = "none") +
  labs(y = "Tissue-specific marker genes", title = "Cluster Identity") +
  theme(axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5))

print(panel_B)
