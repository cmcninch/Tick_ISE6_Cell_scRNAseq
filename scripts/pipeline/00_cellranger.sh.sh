#!/bin/bash
# ==============================================================================
# scripts/pipeline/00_cellranger.sh
# Build CellRanger reference (protein-coding only) and run cellranger count
# Runs on HPC (SLURM). Produces filtered/raw feature-barcode matrices.
# ==============================================================================
set -euo pipefail

################################################################################
### Define paths (edit REPO_ROOT to your clone location) ###
################################################################################

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
genome_dir="${REPO_ROOT}/data/genome_files"
flowcell1_fastqs_path="${REPO_ROOT}/data/raw/Flowcell_AACCCK3M5"
flowcell2_fastqs_path="${REPO_ROOT}/data/raw/Flowcell_HT7V7DRX2"
out_dir="${REPO_ROOT}/results/CellRanger"
cell_ranger_reference="${genome_dir}/Ixodes_scapularis_NCBI_v2"

mkdir -p "${genome_dir}" "${out_dir}/logs"

################################################################################
### Download genome FASTA + GTF (Ixodes scapularis GCF_016920785.2) ###
################################################################################

cd "${genome_dir}"
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/920/785/GCF_016920785.2_ASM1692078v2/GCF_016920785.2_ASM1692078v2_genomic.fna.gz .
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/920/785/GCF_016920785.2_ASM1692078v2/GCF_016920785.2_ASM1692078v2_genomic.gtf.gz .
gunzip -f *.gz

module load cellranger/7.2.0-dntehee

################################################################################
### Filter GTF for protein-coding exons ###
################################################################################

grep "protein_coding" \
  "${genome_dir}/GCF_016920785.2_ASM1692078v2_genomic.gtf" | \
  grep -E -o "LOC[0-9]+" > "${genome_dir}/protein_coding_genes.txt"

grep -f "${genome_dir}/protein_coding_genes.txt" \
  "${genome_dir}/GCF_016920785.2_ASM1692078v2_genomic.gtf" | \
  grep -w "exon" > \
  "${genome_dir}/GCF_016920785.2_ASM1692078v2_genomic_protein_coding_only.gtf"

################################################################################
### Build CellRanger reference ###
################################################################################

cd "${genome_dir}"
cellranger mkref \
  --genome=Ixodes_scapularis_NCBI_v2 \
  --fasta="${genome_dir}/GCF_016920785.2_ASM1692078v2_genomic.fna" \
  --genes="${genome_dir}/GCF_016920785.2_ASM1692078v2_genomic_protein_coding_only.gtf"

################################################################################
### Symbolic links across flowcells ###
################################################################################

mkdir -p "${flowcell1_fastqs_path}/symbolic_links" "${flowcell2_fastqs_path}/symbolic_links"
ln -s "${flowcell1_fastqs_path}"/Sample*/*fastq.gz "${flowcell1_fastqs_path}/symbolic_links"
ln -s "${flowcell2_fastqs_path}"/Sample*/*fastq.gz "${flowcell2_fastqs_path}/symbolic_links"

################################################################################
### Run cellranger count (combine flowcells; spl4_Sept22 is single-flowcell) ###
################################################################################

cd "${out_dir}"
module purge
module load cellranger/7.2.0-dntehee

ls -1 "${flowcell1_fastqs_path}"/symbolic_links/*fastq.gz \
      "${flowcell2_fastqs_path}"/symbolic_links/*fastq.gz | \
  sed -E 's/_S[0-9]+_L[0-9]+_[IR][12]_001.fastq.gz//g' | \
  xargs -n1 basename | sort | uniq > "${out_dir}/samples_to_analyze.txt"

cat "${out_dir}/samples_to_analyze.txt" | grep -v "spl4_Sept22" | while read sample_name
do
  sbatch \
    --job-name=CellRanger_count_${sample_name} \
    --output=${out_dir}/logs/CellRanger_count_${sample_name}.out \
    --error=${out_dir}/logs/CellRanger_count_${sample_name}.error \
    --mem=8g --time=24:00:00 \
    --wrap="cellranger count \
      --id ${sample_name} --sample ${sample_name} \
      --transcriptome ${cell_ranger_reference} \
      --fastqs ${flowcell1_fastqs_path}/symbolic_links,${flowcell2_fastqs_path}/symbolic_links \
      --jobmode slurm --localcores 10 --localmem 100"
done

sbatch \
  --job-name=CellRanger_spl4_Sept22 \
  --output=${out_dir}/logs/CellRanger_count_spl4_Sept22.out \
  --error=${out_dir}/logs/CellRanger_count_spl4_Sept22.error \
  --mem=8g --time=24:00:00 \
  --wrap="cellranger count \
    --id spl4_Sept22 --sample spl4_Sept22 \
    --transcriptome ${cell_ranger_reference} \
    --fastqs ${flowcell1_fastqs_path}/symbolic_links \
    --jobmode slurm --localcores 10 --localmem 100"

################################################################################
### Collect filtered + raw feature-barcode matrices ###
################################################################################

mkdir -p "${out_dir}/filtered_feature_bc_matrices" "${out_dir}/raw_feature_bc_matrices"

for h5_file in "${out_dir}"/*/outs/filtered_feature_bc_matrix.h5; do
  sample_name=$(echo "${h5_file}" | xargs -n1 dirname | xargs -n1 dirname | xargs -n1 basename)
  cp "${h5_file}" "${out_dir}/filtered_feature_bc_matrices/${sample_name}_filtered_feature_bc_matrix.h5"
done

for h5_file in "${out_dir}"/*/outs/raw_feature_bc_matrix.h5; do
  sample_name=$(echo "${h5_file}" | xargs -n1 dirname | xargs -n1 dirname | xargs -n1 basename)
  cp "${h5_file}" "${out_dir}/raw_feature_bc_matrices/${sample_name}_raw_feature_bc_matrix.h5"
done
