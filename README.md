# Insights into tick-pathogen interactions – a single cell RNA sequencing approach of transcriptional changes during ehrlichial infection

Single-cell RNA-Seq–based analyses of transcriptional heterogeneity in the *Ixodes scapularis* embryonic ISE6 cell line and its response to *Ehrlichia muris eauclairensis* (EME) infection. This repository contains the data files and scripts used in analyses undertaken in Adegoke et al. (bioRxiv, 2026; https://doi.org/10.64898/2026.03.19.712879).

# data
This directory contains all data files used in analyses.
- `data/raw/` — raw inputs (FASTQ flowcell directories, where applicable)
- `data/genome_files/` — reference genome FASTA/GTF and CellRanger reference
- `data/processed/` — merged/integrated Seurat objects and other intermediate objects
- `data/bulk_rnaseq/` — bulk tissue-atlas TPM/count matrices (in-house + Rolandelli et al. 2024)
- `data/metadata/` — sample metadata and curated gene lists
- `data/annotation_refs/` — Drosophila, VectorBase, and past-publication reference files

Single-cell RNA-sequencing data and detailed metadata can be accessed through the NCBI Sequence Read Archive under BioProject accession **PRJNA1416822** (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1416822).

# scripts
- `scripts/pipeline/` — numbered upstream pipeline scripts (`00`–`07`) that build every processed intermediate. Run these in order before the figure scripts.
- `scripts/Figure_03.R` … `scripts/Figure_10.R` — one self-contained R script per publication figure, organized by figure panel. Figures 1–2 (microscopy and flow cytometry) are not code-generated.

# results
Intermediate and final outputs, including CellRanger matrices (`results/CellRanger/`), decontX estimates (`results/decontX/`), marker gene tables (`results/markers/`), differential-expression tables (`results/degs/`), GO enrichment tables (`results/go/`), and rendered figures (`results/figures/`).

# Pipeline execution order

The upstream pipeline must be run in sequence; the figure scripts consume its outputs.

```
00_cellranger.sh            (HPC / SLURM)  -> results/CellRanger/*_bc_matrices
01_decontX.R                               -> results/decontX/*_estimates.tsv
02_merge_qc_seurat.R                        -> data/processed/*_unfiltered_merged_seurat_object.rds
03_integration.R                            -> data/processed/*_RPCAintegrated_seurat_object.rds
04_annotation_table.R                       -> data/metadata/consolidated_annotations.tsv
05_marker_genes.R                           -> results/markers/marker_genes_uninfected_only_resolution_0.2.*
06_differential_expression.R                -> results/degs/*_pseudobulk_res0.20.tsv
07_go_enrichment.R                          -> results/go/DiffExp_res0.20_pseudobulk_GO_Results.xlsx
```

Dependency chain: `00 → 01 → 02 → 03`, then `03` feeds both `04 → 05` and `06 → 07`.

Figure-to-input mapping:
- **Fig 3, 9** — need `03` only
- **Fig 4** — needs `03` + `05`
- **Fig 5** — needs `03` + `06`
- **Fig 6, 7, 8** — need `03` + `06` + `07` (plus curated `GO_DEGs_of_interest.xlsx`)
- **Fig 10** — needs `03` + the bulk-atlas files in `data/bulk_rnaseq/`

# Setup and reproducibility

This project uses [`renv`](https://rstudio.github.io/renv/) to manage package dependencies. The `renv.lock`, `.Rprofile`, and `renv/activate.R` files pin the exact environment.

```bash
# 1. Clone the repository
git clone https://github.com/<owner>/<REPO_NAME>.git
cd <REPO_NAME>

# 2. Launch R in the project directory (auto-sources renv/activate.R)
R
```

```r
# 3. Restore the pinned package environment
renv::restore()
```

Then run the pipeline, followed by any figure script, e.g.:

```r
source("scripts/pipeline/01_decontX.R")
# ... 02–07 ...
source("scripts/Figure_03.R")
```

All R paths are relative (via the `here` package), assuming the repository root as the working directory. `00_cellranger.sh` runs on an HPC with SLURM and CellRanger 7.2.0.

# Citation
Adegoke A, Aspinwall J, McNinch C, Ho MCW, Miranda AX, Hoyt F, Nair V, Lack J, Saito TB. *Insights into tick-pathogen interactions – a single cell RNA sequencing approach of transcriptional changes during ehrlichial infection.* bioRxiv (2026). https://doi.org/10.64898/2026.03.19.712879
