#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J splash_tcga
#SBATCH -o splash_tcga_%j.out
#SBATCH -t 12:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

# Load required modules
ml PDC singularity

# Run SPLASH with input file
# Using 128 threads (c-1) as per best practices
# Adjusted anchor_len and target_len for 50bp reads (default is 31+31=62bp for 101bp reads)
srun -c 127 singularity exec -B /cfs/klemming ../splash splash input_CPTAC.txt \
    --anchor_len 27 --target_len 27 --gap_len 0 --poly_ACGT_len 6 \
    --compactors_config compactors.json --kmc_use_RAM_only_mode --n_most_freq_targets 4 \
    --dump_sample_anchor_target_count_binary --outname_prefix 2025-12-04_test
