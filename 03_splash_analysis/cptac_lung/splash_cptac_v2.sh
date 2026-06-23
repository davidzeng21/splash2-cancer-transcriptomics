#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J splash_cptac
#SBATCH -o splash_cptac_%j.out
#SBATCH -t 24:00:00
#SBATCH -p main

#SBATCH -n 1
#SBATCH -c 256
#SBATCH --mem=440G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

# Load required modules
ml PDC singularity

# Run SPLASH with input file
srun -c 256 singularity exec -B /cfs/klemming ../splash splash input_CPTAC_107.txt \
    --anchor_len 25 --target_len 25 --gap_len 0 \
    --compactors_config compactors.json \
    --n_threads_stage_1 128 --n_threads_stage_1_internal 2 --n_threads_stage_2 256 \
    --kmc_use_RAM_only_mode --kmc_max_mem_GB 20 \
    --poly_ACGT_len 6 --n_most_freq_targets 4 --anchor_samples_threshold 100 \
    --dump_sample_anchor_target_count_binary --outname_prefix 2026-02-05_min_100

