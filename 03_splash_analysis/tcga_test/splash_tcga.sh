#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J splash_tcga
#SBATCH -o splash_tcga_%j.out
#SBATCH -t 3:00:00
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem=256G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

# Load required modules
ml PDC singularity

# Run SPLASH with input file
# Using 63 threads (c-1) as per best practices
# Adjusted anchor_len and target_len for 50bp reads (default is 31+31=62bp for 101bp reads)
srun -c 63 singularity exec -B /cfs/klemming ../splash splash input.txt --anchor_len 17 --target_len 17
