# ==============================================================================
# scripts/pipeline/04_annotation_table.R
# Build consolidated gene annotation table (BLAST/FlyBase/FlyCellAtlas/
# OrthoFinder/past-publication cross-refs). Requires precomputed BLAST +
# OrthoFinder outputs (see confidence report items 10, 12).
# Output: data/metadata/consolidated_annotations.tsv
# ==============================================================================

library(Seurat)
library(tidyverse)
library(readxl)
library(openxlsx)
library(rtracklayer)
library(here)

# ------------------------------------------------------------------------------
# File paths (Drosophila references + OrthoFinder + VectorBase cross-refs)
# ------------------------------------------------------------------------------
dref <- here("data", "annotation_refs", "Drosophila_scRNAseq_references")
d_melanogaster_gtf_file_path            <- file.path(dref, "dmel-all-r6.59_just_polypeptides.gff")
d_melanogaster_FBgn_to_FBpp_path        <- file.path(dref, "fbgn_fbtr_fbpp_fb_2024_04.txt")
d_melanogaster_FBgn_to_gene_symbol_path <- file.path(dref, "fbgn_annotation_ID.txt")
d_melanogaster_tissue_markers_path      <- file.path(dref, "science.abk2432_table s2.xlsx")
d_melanogaster_MidBrain_atlas_path      <- file.path(dref, "elife-34550-fig1-data2-v2.csv")

d_melanogaster_blast_hits_path <- here(
  "results", "BLAST",
  "GCF_016920785.2_ASM1692078v2_translated_cds_BLASTP_against_Drosophila_melanogaster_translated_cds"
)

vbase <- here("data", "annotation_refs", "VectorBase_Reference_Files")
PalLabHiFi_gene_to_ISCP_protein_ID_path <- file.path(vbase, "VectorBase-68_IscapularisPalLabHiFi_GeneAliases.txt")
Saito_v_all_ixodida_bySaito_geneID_path <- here("results", "OrthoFinder", "Saito_v_all_annots_bySaito_geneID_Dec132024.tsv")

ioxedes_scapularis_gene_annotations_from_joe_path <- here(
  "data", "metadata", "Genes from tick publications 10.22.24.xlsx"
)
gtf_file_path <- here("data", "genome_files",
                      "GCF_016920785.2_ASM1692078v2_genomic_protein_coding_only.gtf")

pubs <- here("data", "annotation_refs", "Gene_Searches_For_Past_Publications_In_Vectorbase")
alberdi_2016_supp_table_2_path    <- file.path(pubs, "Alberdi_et al_Tissue-specific signatures--2016_SuppTable2.XLSX")
egekwu_2014_supp_table_1_path     <- file.path(pubs, "Egekwu_et al_Transcriptome of the--2014_VectorBase_Search_Results.xlsx")
mateos_2021_table_1_path          <- file.path(pubs, "Mateos-Hernandez_et al_Enlisting the ISE6--2021_VectorBase_Search_Results.xlsx")
oliver_2015_table_path            <- file.path(pubs, "Oliver_et al_An IS cell line--2015_Results_of_XP_Protein_ID_Search_In_VectorBase.xlsx")
rolandelli_2024_supp_table_6_path <- file.path(pubs, "Rolandelli_et al_Tick hemocytes--2024_Supp_Table6.xlsx")

out_path <- here("data", "metadata", "consolidated_annotations.tsv")

# ------------------------------------------------------------------------------
# Drosophila gene-info lookup (FBpp -> FBgn -> symbol)
# ------------------------------------------------------------------------------
d_melanogaster_FBgn_to_FBpp <- read_tsv(d_melanogaster_FBgn_to_FBpp_path, show_col_types = FALSE) %>%
  dplyr::select(FlyBase_FBgn, FlyBase_FBpp) %>%
  dplyr::rename("FlyBase_GeneID" = FlyBase_FBgn, "FlyBase_PolypeptideID" = FlyBase_FBpp)

d_melanogaster_FBgn_to_gene_symbol <- read_tsv(d_melanogaster_FBgn_to_gene_symbol_path, show_col_types = FALSE) %>%
  dplyr::select(primary_FBgn, gene_symbol) %>%
  dplyr::rename("FlyBase_GeneID" = primary_FBgn, "FlyBase_GeneName" = gene_symbol)

d_melanogaster_gene_info <- rtracklayer::import(d_melanogaster_gtf_file_path) %>%
  as.data.frame() %>%
  distinct(ID, Name) %>%
  dplyr::rename("FlyBase_PolypeptideID" = ID, "FlyBase_PolypetideName" = Name) %>%
  left_join(d_melanogaster_FBgn_to_FBpp, by = "FlyBase_PolypeptideID") %>%
  left_join(d_melanogaster_FBgn_to_gene_symbol, by = "FlyBase_GeneID") %>%
  dplyr::select(FlyBase_GeneID, FlyBase_GeneName, FlyBase_PolypeptideID, FlyBase_PolypetideName)

# ------------------------------------------------------------------------------
# Top BLAST hit per Ixodes gene (e-value < 1e-10, best e-value then %identity)
# ------------------------------------------------------------------------------
d_melanogaster_blast_hits <- read_tsv(
  d_melanogaster_blast_hits_path, show_col_types = FALSE,
  col_names = c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  filter(evalue < 1e-10) %>%
  group_by(qaccver) %>%
  filter(evalue == min(evalue)) %>%
  slice_max(order_by = pident, n = 1, with_ties = FALSE) %>%
  left_join(d_melanogaster_gene_info, by = c("saccver" = "FlyBase_PolypeptideID")) %>%
  dplyr::rename(
    "Ioxedes_scapularis_GeneID" = qaccver,
    "Top_Drosophila_BLAST_Hit_FlyBase_PolypeptideID" = saccver,
    "Top_Drosophila_BLAST_Hit_evalue" = evalue,
    "Top_Drosophila_BLAST_Hit_PercentIdentity" = pident) %>%
  ungroup() %>%
  dplyr::select(Ioxedes_scapularis_GeneID, Top_Drosophila_BLAST_Hit_FlyBase_PolypeptideID,
                FlyBase_GeneID, FlyBase_GeneName, FlyBase_PolypetideName,
                Top_Drosophila_BLAST_Hit_evalue, Top_Drosophila_BLAST_Hit_PercentIdentity)

# ------------------------------------------------------------------------------
# FlyCellAtlas + MidBrain atlas markers + Joe's manual annotations
# ------------------------------------------------------------------------------
FlyCellAtlas_tissue_markers <- read_xlsx(d_melanogaster_tissue_markers_path) %>%
  filter(!is.na(Tissue) & !is.na(Annotation)) %>%
  dplyr::select(Tissue, `Broad annotation`, Annotation, `Curators Markers`) %>%
  separate_longer_delim(cols = `Curators Markers`, delim = ", ") %>%
  dplyr::rename("FlyCellAtlas_TissueName" = Tissue,
                "FlyCellAtlas_BroadCellAnnotation" = `Broad annotation`,
                "FlyCellAtlas_SpecificCellAnnotation" = Annotation,
                "FlyCellAtlas_GeneName" = `Curators Markers`) %>%
  group_by(FlyCellAtlas_GeneName) %>%
  summarise(
    FlyCellAtlas_TissueName = paste(FlyCellAtlas_TissueName, collapse = "|||"),
    FlyCellAtlas_BroadCellAnnotation = paste(FlyCellAtlas_BroadCellAnnotation, collapse = "|||"),
    FlyCellAtlas_SpecificCellAnnotation = paste(FlyCellAtlas_SpecificCellAnnotation, collapse = "|||")
  ) %>%
  drop_na()

MidBrain_atlas_markers <- read_csv(d_melanogaster_MidBrain_atlas_path, show_col_types = FALSE) %>%
  dplyr::select(gene, cluster) %>%
  dplyr::rename("MidBrainCellAtlas_Marker_GeneName" = gene,
                "MidBrainCellAtlas_Marker_Cluster" = cluster) %>%
  group_by(MidBrainCellAtlas_Marker_GeneName) %>%
  summarise(MidBrainCellAtlas_Marker_Cluster = paste(MidBrainCellAtlas_Marker_Cluster, collapse = "|||"))

ioxedes_scapularis_gene_annotations_from_joe <- read_xlsx(
  ioxedes_scapularis_gene_annotations_from_joe_path, range = cell_cols("C:D")) %>%
  dplyr::rename("Joe_Annotation" = Anotation, "ID" = `HIFI Gene Name`)

consolidated_d_melanogaster_gene_info <- d_melanogaster_gene_info %>%
  full_join(d_melanogaster_blast_hits,
            by = c("FlyBase_PolypeptideID" = "Top_Drosophila_BLAST_Hit_FlyBase_PolypeptideID",
                   "FlyBase_GeneID", "FlyBase_GeneName", "FlyBase_PolypetideName")) %>%
  full_join(FlyCellAtlas_tissue_markers, by = c("FlyBase_GeneName" = "FlyCellAtlas_GeneName")) %>%
  full_join(MidBrain_atlas_markers, by = c("FlyBase_GeneName" = "MidBrainCellAtlas_Marker_GeneName")) %>%
  full_join(ioxedes_scapularis_gene_annotations_from_joe, by = c("Ioxedes_scapularis_GeneID" = "ID")) %>%
  dplyr::select(Ioxedes_scapularis_GeneID, Joe_Annotation, everything()) %>%
  filter(!is.na(Ioxedes_scapularis_GeneID))

# ------------------------------------------------------------------------------
# OrthoFinder cross-genome ortholog table (iSE6 / Wikel / PalLabHiFi + ixodida)
# ------------------------------------------------------------------------------
PalLabHiFi_gene_to_ISCP_protein_ID <- read_tsv(
  PalLabHiFi_gene_to_ISCP_protein_ID_path, show_col_types = FALSE,
  col_names = c("ISCP_ID", "LOC_ID"))

Saito_v_all_ixodida_bySaito_geneID <- read_tsv(
  Saito_v_all_ixodida_bySaito_geneID_path, show_col_types = FALSE) %>%
  mutate(ISCP_ID = gsub("^Saito-", "", Saito_geneID)) %>%
  left_join(PalLabHiFi_gene_to_ISCP_protein_ID, by = c("ISCP_ID")) %>%
  filter(!is.na(PalLabHiFi_geneID)) %>%
  dplyr::select(LOC_ID, ISCP_ID, PalLabHiFi_geneID:Rsan_gene_product) %>%
  dplyr::rename(
    "PalLabHiFi_ProteinID" = PalLabHiFi_geneID, "PalLabHiFi_ProteinProduct" = PalLabHiFi_gene_product,
    "ISE6_OrthologProteinID" = ISE6_geneID, "ISE6_OrthologProteinProduct" = ISE6_gene_product,
    "Wikel_OrthologProteinID" = Wikel_geneID, "Wikel_OrthologProteinProduct" = Wikel_gene_product,
    "Dsil_ProteinID" = Dsil_geneID, "Dsil_ProteinProduct" = Dsil_gene_product,
    "HaeL_ProteinID" = HaeL_geneID, "HaeL_ProteinProduct" = HaeL_gene_product,
    "Hyas_ProteinID" = Hyas_geneID, "Hyas_ProteinProduct" = Hyas_gene_product,
    "Rmic_ProteinID" = Rmic_geneID, "Rmic_ProteinProduct" = Rmic_gene_product,
    "Rsan_ProteinID" = Rsan_geneID, "Rsan_ProteinProduct" = Rsan_gene_product)

# ------------------------------------------------------------------------------
# Past-publication ID cross-references (LOC_ID <-> ortholog protein IDs)
# ------------------------------------------------------------------------------
locid_to_ortholog_protein_id <- Saito_v_all_ixodida_bySaito_geneID %>%
  dplyr::select(LOC_ID, PalLabHiFi_ProteinID, ISE6_OrthologProteinID, Wikel_OrthologProteinID) %>%
  pivot_longer(2:ncol(.), names_to = "Genome", values_to = "ID") %>%
  filter(!is.na(ID)) %>%
  mutate(Genome = gsub("_ProteinID", "", Genome),
         Genome = gsub("_OrthologProteinID", "", Genome))

alberdi_2016_supp_table_2 <- read_xlsx(alberdi_2016_supp_table_2_path) %>%
  dplyr::select(gene, test_id) %>% filter(gene != "-") %>%
  left_join(locid_to_ortholog_protein_id, by = c("gene" = "ID")) %>%
  filter(!is.na(LOC_ID)) %>%
  dplyr::rename("ID_In_Alberdi_2016_SuppTable2" = test_id) %>%
  dplyr::select(LOC_ID, ID_In_Alberdi_2016_SuppTable2)

egekwu_2014_supp_table_1 <- read_xlsx(egekwu_2014_supp_table_1_path) %>%
  dplyr::select(`Gene ID`, `Input ID`) %>%
  left_join(locid_to_ortholog_protein_id, by = c(`Gene ID` = "ID")) %>%
  filter(!is.na(LOC_ID)) %>%
  dplyr::rename("ID_In_Egekwu_2014_SuppTable1" = `Input ID`) %>%
  dplyr::select(LOC_ID, ID_In_Egekwu_2014_SuppTable1)

mateos_2021_table_1 <- read_xlsx(mateos_2021_table_1_path) %>%
  dplyr::select(`Gene ID`, `Input ID`) %>%
  left_join(locid_to_ortholog_protein_id, by = c(`Gene ID` = "ID")) %>%
  filter(!is.na(LOC_ID)) %>%
  dplyr::rename("ID_In_Mateos_2021_Table1" = `Input ID`) %>%
  dplyr::select(LOC_ID, ID_In_Mateos_2021_Table1)

oliver_2015_table <- read_xlsx(oliver_2015_table_path) %>%
  dplyr::select(`Gene ID`, `Input ID`) %>%
  left_join(locid_to_ortholog_protein_id, by = c(`Gene ID` = "ID")) %>%
  filter(!is.na(LOC_ID)) %>%
  dplyr::rename("ID_In_Oliver_2015_Table" = `Input ID`) %>%
  dplyr::select(LOC_ID, ID_In_Oliver_2015_Table)

rolandelli_2024_supp_table_6 <- bind_rows(
  lapply(c("Immune1 cluster", "Immune2 cluster", "Proliferative1 cluster",
           "Proliferative2 cluster", "Transitional cluster", "Gut cluster"),
         function(sheet_name) {
           read_xlsx(rolandelli_2024_supp_table_6_path, sheet = sheet_name, skip = 5) %>%
             dplyr::select(`LOC#`) %>% mutate(Cluster = sheet_name)
         })
) %>%
  group_by(`LOC#`) %>%
  summarise(Cluster = paste(Cluster, collapse = "|")) %>%
  dplyr::rename("LOC_ID" = `LOC#`, "Cluster_In_Rolandelli_2024_SuppTable6" = Cluster)

# ------------------------------------------------------------------------------
# All genes -> join everything -> write consolidated annotation table
# ------------------------------------------------------------------------------
gene_id_to_gene_name <- rtracklayer::import(gtf_file_path) %>%
  as.data.frame() %>%
  distinct(gene_id, product) %>%
  dplyr::rename("gene" = gene_id) %>%
  distinct(gene, .keep_all = TRUE)

gene_id_to_gene_name %>%
  left_join(Saito_v_all_ixodida_bySaito_geneID, by = c("gene" = "LOC_ID")) %>%
  left_join(consolidated_d_melanogaster_gene_info, by = c("gene" = "Ioxedes_scapularis_GeneID")) %>%
  left_join(alberdi_2016_supp_table_2, by = c("gene" = "LOC_ID")) %>%
  left_join(egekwu_2014_supp_table_1, by = c("gene" = "LOC_ID")) %>%
  left_join(mateos_2021_table_1, by = c("gene" = "LOC_ID")) %>%
  left_join(oliver_2015_table, by = c("gene" = "LOC_ID")) %>%
  left_join(rolandelli_2024_supp_table_6, by = c("gene" = "LOC_ID")) %>%
  write_tsv(out_path)
