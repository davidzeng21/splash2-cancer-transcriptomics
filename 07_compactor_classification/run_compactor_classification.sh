#!/bin/bash -l

#SBATCH -A naiss2025-22-738
#SBATCH -p main
#SBATCH -c 64
#SBATCH -t 4:00:00
#SBATCH -J classify_compactors
#SBATCH -o classify_compactors_%j.out
#SBATCH --mem=120G
#SBATCH --begin=now
#SBATCH --mail-type=ALL

ml PDC R rnastar/2.7.11b samtools/1.20 bedtools/2.31.0 bowtie2/2.5.4

REF_DIR=/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/ref_genomes
WORK_DIR=/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched
COMPACTORS_DIR=${WORK_DIR}/2026-02-20_min_50_samples/2026-02-20_min_50_samples_compactors
OUTPUT_DIR=${WORK_DIR}/compactor_classification/

mkdir -p ${OUTPUT_DIR}

Rscript ${WORK_DIR}/SPLASH2_compactor_classification.R \
    ${OUTPUT_DIR} \
    ${COMPACTORS_DIR}/after_correction.scores.top_effect_size_bin.tsv \
    ${REF_DIR}/STAR_index_v49 \
    ${REF_DIR}/Bowtie2_index/GRCh38 \
    ${REF_DIR}/GRCh38_gencode_v49_known_splice_sites.txt \
    ${REF_DIR}/GRCh38_gencode_v49_exon_coordinates.bed \
    ${REF_DIR}/GRCh38_gencode_v49_genes.bed \
    25
