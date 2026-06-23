#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J compactor_cptac
#SBATCH -o compactor_cptac_%j.out
#SBATCH -t 6:00:00
#SBATCH -p main

#SBATCH -n 1
#SBATCH -c 256
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

# Load required modules
ml PDC singularity

# Run Compactor with filtered anchors
srun -c 256 singularity exec -B /cfs/klemming ../splash compactors \
    fastq_list_test10.txt \
    extendor_threshold_sweep/anchors_superset.txt \
    compactors_output.tsv \
    --epsilon 0.01 \
    --num_kmers 2 \
    --kmer_len 24 \
    --num_threads 256 \
    --log compactors_run.log
