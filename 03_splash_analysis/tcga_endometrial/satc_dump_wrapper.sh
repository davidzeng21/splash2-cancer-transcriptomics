#!/bin/bash
# Wrapper script to run satc_dump via singularity
# This avoids GLIBC dependency issues

# Load singularity if not already loaded
module load PDC singularity 2>/dev/null || true

# Path to singularity image
SPLASH_IMAGE="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/splash"

# Run satc_dump inside singularity container with all passed arguments
singularity exec -B /cfs/klemming "${SPLASH_IMAGE}" /usr/local/bin/satc_dump "$@"

