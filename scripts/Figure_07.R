# ==============================================================================
# Figure_07.R
# GO enrichment + DEG heatmap comparing EARLY vs LATE infection (EME_d4 vs EME_d2)
# Panels: 7A (GO dot plot), 7B (DEGs within selected enriched GO terms heatmap)
# ==============================================================================

library(tidyverse)
library(readxl)
library(Seurat)
library(ComplexHeatmap)
library(edgeR)
library(circlize)
library(scales)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
diff_exp_go_results_path <- here("results", "go", "DiffExp_res0.20_pseudobulk_GO_Results.xlsx")
integrated_seurat_object_path <- here(
  "data", "processed",
  "iSE6_..._SCTransformed_RPCAintegrated_seurat_object.rds"
)
go_degs_of_interest_path <- here("data", "metadata", "GO_DEGs_of_interest.xlsx")

# ------------------------------------------------------------------------------
## Panel A: GO dot plot (EME_d4 vs EME_d2, All Cells)
# ------------------------------------------------------------------------------
plot_go_results <- function(go_results, legend_position, plot_title) {
  ggplot(go_results,
         aes(x = reorder(Term, Significant / Annotated),
             y = 100 * Significant / Annotated,
             size = Annotated, color = -log10(P_Value))) +
    geom_point() + coord_flip() +
    facet_grid(Ontology_Type ~ ., scales = "free", space = "free_y") +
    scale_color_gradient(low = "blue", high = "red") + theme_bw() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.position = legend_position,
          panel.border = element_rect(size = 3),
          strip.background = element_blank(),
          strip.text = element_text(size = 24),
          plot.margin = unit(c(1, 1, 1, 1.25), "cm")) +
    labs(x = NULL, y = "% of Annotated Genes\nFound Significant",
         size = "Annotated Genes\nper GO Term",
         color = "-log10(p-value)", title = plot_title)
}

EME_d4_vs_EME_d2_go_results <- read_xlsx(
  diff_exp_go_results_path, sheet = "EMEd4 vs EMEd2 GO Results"
) %>%
  filter(cluster == "All Cells")

print(
  plot_go_results(EME_d4_vs_EME_d2_go_results,
                  plot_title = "EME_d4 vs EME_d2",
                  legend_position = "right")
)

# ------------------------------------------------------------------------------
## Panel B: Heatmap of DEGs within selected enriched GO terms
# ------------------------------------------------------------------------------
integrated_seurat_object <- readRDS(integrated_seurat_object_path)

pseudobulk_integrated_seurat_object <- AggregateExpression(
  integrated_seurat_object, assays = "RNA", return.seurat = TRUE,
  group.by = c("Sample", "Sample_Group")
)

counts_for_de <- pseudobulk_integrated_seurat_object[["RNA"]]$counts
count_matrix  <- counts_for_de[rowSums(counts_for_de) > 0, ]
tmms <- cpm(calcNormFactors(DGEList(counts = count_matrix), method = "TMM")) %>%
  as.data.frame()

go_degs_of_interest_EME_d4_vs_EME_d2 <- read_xlsx(
  go_degs_of_interest_path, sheet = "GO_DEGs_EME_d4_vs_EME_d2"
) %>%
  filter(!str_detect(Product, "^LOC"))

metadata <- pseudobulk_integrated_seurat_object@meta.data %>%
  filter(Sample_Group %in% c("EME-d4", "EME-d2")) %>%
  mutate(Sample_Group = gsub("-", "_", Sample_Group))

unique_gos <- unique(go_degs_of_interest_EME_d4_vs_EME_d2$GO_ID)

scaled_numeric_matrix <- tmms %>%
  rownames_to_column(var = "gene_id") %>%
  dplyr::select(gene_id, all_of(rownames(metadata))) %>%
  filter(gene_id %in% filter(go_degs_of_interest_EME_d4_vs_EME_d2,
                             GO_ID == unique_gos[1])$Gene) %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix() %>% t() %>% scale() %>% t()

for (i in 2:length(unique_gos)) {
  tmp <- tmms %>%
    rownames_to_column(var = "gene_id") %>%
    dplyr::select(gene_id, all_of(rownames(metadata))) %>%
    filter(gene_id %in% filter(go_degs_of_interest_EME_d4_vs_EME_d2,
                               GO_ID == unique_gos[i])$Gene) %>%
    column_to_rownames(var = "gene_id") %>%
    as.matrix() %>% t() %>% scale() %>% t()
  scaled_numeric_matrix <- rbind(scaled_numeric_matrix, tmp)
}

annotation <- metadata %>% dplyr::select(Sample_Group)
col_fun <- colorRamp2(c(-2, 0, 2), c("mediumblue", "white", "red3"))

heatmap <- Heatmap(
  scaled_numeric_matrix,
  col = col_fun(seq(floor(min(scaled_numeric_matrix)),
                    ceiling(max(scaled_numeric_matrix)), 1)),
  use_raster = TRUE, show_row_dend = FALSE, show_row_names = TRUE,
  show_column_names = FALSE, cluster_columns = TRUE, cluster_rows = FALSE,
  split = go_degs_of_interest_EME_d4_vs_EME_d2$GO_ID,
  row_names_side = "right", gap = unit(0.2, "cm"), row_title_rot = 0,
  border = TRUE, row_names_gp = gpar(fontsize = 8),
  row_title_gp = gpar(fontsize = 10),
  row_labels = go_degs_of_interest_EME_d4_vs_EME_d2$Product,
  heatmap_width = unit(2.5, "npc"), width = unit(0.1, "npc"),
  heatmap_legend_param = list(
    title = "Z-score", title_gp = gpar(fontsize = 16, fontface = "bold"),
    direction = "horizontal", legend_width = unit(3, "cm"),
    at = c(-2, -1, 0, 1, 2), labels_gp = gpar(fontsize = 12),
    title_position = "topcenter"),
  top_annotation = HeatmapAnnotation(
    df = annotation,
    col = list(Sample_Group = c(`EME_d4` = hue_pal()(4)[2],
                                `EME_d2` = hue_pal()(4)[1])),
    border = TRUE, simple_anno_size = unit(1.5, "cm"),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 16)))
)

draw(heatmap, heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom", merge_legend = TRUE)
