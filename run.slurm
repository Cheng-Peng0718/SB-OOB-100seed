#!/bin/bash --login
#SBATCH -J oob_sb
#SBATCH --time=0-23:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=120G
#SBATCH --array=1-100
#SBATCH --constraint="agx|agg"
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

ulimit -c 0

echo "[INFO] Array task ${SLURM_ARRAY_TASK_ID} started on $(hostname) at $(date)"
echo "[INFO] CPUs: ${SLURM_CPUS_PER_TASK}"

module purge
module load R/4.4.1
cd ${SLURM_SUBMIT_DIR}

Rscript --vanilla Unified_parallel.R "${SLURM_ARRAY_TASK_ID}"

echo "[INFO] Task ${SLURM_ARRAY_TASK_ID} finished at $(date)"