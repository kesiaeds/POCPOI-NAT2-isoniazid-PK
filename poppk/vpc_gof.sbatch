#!/bin/bash
#SBATCH --job-name=vpc_gof
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal
#SBATCH --output=output/vpc_%j.log
#SBATCH --error=output/vpc_%j.err

set -euo pipefail

ml R/4.3 gcc/12 cmake openblas/0.3.20

POCPOI_WORKDIR="${POCPOI_WORKDIR:-/scratch/users/$USER/POCPOI}"
cd "$POCPOI_WORKDIR"
export R_LIBS_USER="$POCPOI_WORKDIR/Rlib/4.3"

Rscript scripts/vpc_gof.R
