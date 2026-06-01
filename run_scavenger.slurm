#!/bin/bash --login
#SBATCH -J oob_sb_scv
#SBATCH --partition=scavenger
#SBATCH --qos=scavenger
#SBATCH --requeue
#SBATCH --time=0-23:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=120G
#SBATCH --array=1-100
#SBATCH --constraint="agx|agg"
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

ulimit -c 0

SEED=${SLURM_ARRAY_TASK_ID}
SEED_DIR="${SLURM_SUBMIT_DIR}/results/seed_$(printf "%03d" ${SEED})"
DONE_MARK="${SEED_DIR}/exp5_raw.RData"   # 最后一个生成的文件，存在即代表整seed跑完

echo "[INFO] Array task ${SEED} started on $(hostname) at $(date)"
echo "[INFO] CPUs: ${SLURM_CPUS_PER_TASK}  Partition: ${SLURM_JOB_PARTITION}  QOS: ${SLURM_JOB_QOS}"

## ----------- 幂等性检查：seed 已完成则跳过 -----------
if [[ -f "${DONE_MARK}" ]]; then
  echo "[SKIP] seed ${SEED} already complete (${DONE_MARK} exists). Exiting cleanly."
  exit 0
fi

## ----------- 半完成清理：避免和并行 array job 互相覆盖损坏文件 -----------
## 如果 SEED_DIR 存在但 exp5_raw.RData 不存在，说明上次没跑完（可能被 preempt）
## 整目录删掉重跑，保证不会读到坏 RData
if [[ -d "${SEED_DIR}" ]]; then
  echo "[INFO] seed ${SEED} dir exists but incomplete, wiping for clean restart"
  rm -rf "${SEED_DIR}"
fi

module purge
module load R/4.4.1
cd ${SLURM_SUBMIT_DIR}

Rscript --vanilla Unified_parallel.R "${SEED}"
RC=$?

if [[ ${RC} -ne 0 ]]; then
  echo "[ERROR] Rscript exited with code ${RC}, likely preempted or crashed"
  ## 不要保留半成品 SEED_DIR，否则会迷惑后续重跑逻辑
  ## 但 SLURM --requeue 会自动重新调度，无需手动 resubmit
fi

echo "[INFO] Task ${SEED} finished at $(date) (rc=${RC})"
exit ${RC}
