#!/bin/bash
#SBATCH --job-name=poppk_inh
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=normal
#SBATCH --output=output/poppk_%j.log
#SBATCH --error=output/poppk_%j.err

set -euo pipefail

ml R/4.3 gcc/12 cmake openblas/0.3.20

POCPOI_WORKDIR="${POCPOI_WORKDIR:-/scratch/users/$USER/POCPOI}"
cd "$POCPOI_WORKDIR"
export R_LIBS_USER="$POCPOI_WORKDIR/Rlib/4.3"

Rscript scripts/pk_analysis_combined.R
