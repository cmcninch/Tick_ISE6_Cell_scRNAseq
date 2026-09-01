# ==============================================================================
# scripts/pipeline/07_go_enrichment.R
# topGO enrichment (weight algorithm, Fisher test, p < 0.01) on pseudobulk DEGs
# for EMEd2_vs_UId2, EMEd4_vs_UId4, EMEd4_vs_EMEd2. Requires a topGO-compatible
# GO mapping file (build once via the commented block, from the VectorBase GAF).
# Output: results/go/DiffExp_res0.20_pseudobulk_GO_Results.xlsx
# ==============================================================================

library(tidyverse)
library(readxl)
library(openxlsx)
library(topGO)
library(GO.db)
library(here)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
go_gaf_file_path         <- here("data", "annotation_refs", "VectorBase_Reference_Files",
                                 "VectorBase-68_IscapularisPalLabHiFi_GO.gaf")
go_topGO_compatible_path <- here("data", "annotation_refs", "VectorBase_Reference_Files",
                                 "VectorBase-68_IscapularisPalLabHiFi_GO_topGOcompatible.txt")
locID_to_protein_path    <- here("data", "metadata", "LOC_ID_to_protein_info.txt")

degs_dir <- here("results", "degs")
EMEd2_vs_UId2_path    <- file.path(degs_dir, "EMEd2_vs_UId2_pseudobulk_res0.20.tsv")
EMEd4_vs_UId4_path    <- file.path(degs_dir, "EMEd4_vs_UId4_pseudobulk_res0.20.tsv")
EME_d4_vs_EME_d2_path <- file.path(degs_dir, "EME_d4_vs_EME_d2_pseudobulk_res0.20.tsv")

out_dir <- here("results", "go")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
diff_exp_go_results_path <- file.path(out_dir, "DiffExp_res0.20_pseudobulk_GO_Results.xlsx")

# ------------------------------------------------------------------------------
# GO enrichment function (returns annotated genes per term)
# ------------------------------------------------------------------------------
GO_enrichment_func <- function(myInterestingGenes, go_file_path,
                               target_ontologies, p_value_cutoff) {
  geneID2GO <- readMappings(file = go_file_path, sep = "\t", IDsep = ",")
  geneNames <- names(geneID2GO)
  geneList  <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames

  go_results <- tibble()
  for (j in seq_along(target_ontologies)) {
    go_data <- new(Class = "topGOdata", ontology = target_ontologies[j],
                   allGenes = geneList, nodeSize = 10,
                   annot = annFUN.gene2GO, gene2GO = geneID2GO)
    Fisher_Test <- runTest(go_data, algorithm = "weight", statistic = "fisher")
    temp_go_results <- GenTable(
      go_data, P_Value = Fisher_Test,
      topNodes = length(Fisher_Test@score[Fisher_Test@score <= p_value_cutoff]))
    colnames(temp_go_results) <- c("GO", "Term", "Annotated", "Significant",
                                   "Expected", "P_Value")

    AnnotatedGenes <- lapply(temp_go_results$GO, function(x)
      as.character(unlist(genesInTerm(object = go_data, whichGO = x))))
    names(AnnotatedGenes) <- temp_go_results$GO
    AnnotatedGenes_df <- data.frame(
      GO = names(AnnotatedGenes),
      Genes = sapply(AnnotatedGenes, function(x) paste(x, collapse = ", ")),
      stringsAsFactors = FALSE)

    temp_go_results <- temp_go_results %>%
      mutate(P_Value = as.numeric(sub('<', '', P_Value)),
             Ontology_Type = as.character(target_ontologies[j])) %>%
      left_join(AnnotatedGenes_df, by = "GO")

    go_results <- bind_rows(go_results, temp_go_results)
  }
  go_results
}

# ------------------------------------------------------------------------------
# (Build topGO-compatible GO file once from the VectorBase GAF -- run if missing)
# ------------------------------------------------------------------------------
# locID_to_protein <- read_tsv(locID_to_protein_path, show_col_types = FALSE) %>%
#   distinct(LOC_ID, `Gene ID`, .keep_all = FALSE)
# read_tsv(go_gaf_file_path, show_col_types = FALSE, skip = 1,
#          col_names = c("Extra_1","Gene ID","Extra_2","Extra_3","GO_ID","Extra_4",
#                        "Extra_5","Extra_6","Extra_7","Extra_8","Extra_9",
#                        "Extra_10","Extra_11")) %>%
#   dplyr::select(`Gene ID`, GO_ID) %>%
#   group_by(`Gene ID`) %>%
#   summarise(GO_ID = paste(GO_ID, collapse = ",")) %>%
#   left_join(locID_to_protein, by = "Gene ID") %>%
#   dplyr::select(LOC_ID, GO_ID) %>%
#   write_tsv(go_topGO_compatible_path, col_names = FALSE)

# ------------------------------------------------------------------------------
# Run GO enrichment per comparison, over all clusters present
# ------------------------------------------------------------------------------
run_go_over_clusters <- function(deg_path) {
  degs <- read_tsv(deg_path, show_col_types = FALSE)
  clusters <- unique(degs$Cluster)
  out <- tibble()
  for (i in seq_along(clusters)) {
    message(clusters[i])
    out <- bind_rows(out, GO_enrichment_func(
      myInterestingGenes = degs %>% filter(Cluster == clusters[i]) %>% pull(Gene),
      go_file_path = go_topGO_compatible_path,
      target_ontologies = c("BP", "CC", "MF"),
      p_value_cutoff = 0.01
    ) %>% mutate(cluster = clusters[i]))
  }
  out
}

EME_d2_vs_UI_d2_go_results  <- run_go_over_clusters(EMEd2_vs_UId2_path)
EME_d4_vs_UI_d4_go_results  <- run_go_over_clusters(EMEd4_vs_UId4_path)
EME_d4_vs_EME_d2_go_results <- run_go_over_clusters(EME_d4_vs_EME_d2_path)

write.xlsx(
  list(
    "EMEd2 vs UId2 GO Results"  = EME_d2_vs_UI_d2_go_results,
    "EMEd4 vs UId4 GO Results"  = EME_d4_vs_UI_d4_go_results,
    "EMEd4 vs EMEd2 GO Results" = EME_d4_vs_EME_d2_go_results
  ),
  firstRow = TRUE, colWidths = "auto", diff_exp_go_results_path
)
