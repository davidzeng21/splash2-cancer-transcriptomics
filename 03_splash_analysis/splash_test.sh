#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J splash
#SBATCH -o splash_%j.out
#SBATCH -t 1:00:00
#SBATCH -p shared
#SBATCH -c 25
#SBATCH --mail-type=ALL
#SBATCH --begin=now

ml PDC singularity
srun -c 24 singularity exec -B /cfs/klemming splash splash test_run/input.txt
