#!/bin/bash

#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-16
#SBATCH -t 12:00:00
#SBATCH -C centos7
#SBATCH -p sched_any
#SBATCH -o outputs/analyze_array_%A_%a.out
#SBATCH -e outputs/analyze_array_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=nhirt@mit.edu

set -euo pipefail

mkdir -p outputs

RESULTS_ROOT="${1:-SlabDesignFactors/results/remote_results_full_sweep}"
RUN_NAME="${2:-full_sweep}"
COMPLETION_DIR="${3:-${RESULTS_ROOT}/completion}"
PARAMS_FILE="${4:-/home/nhirt/2024_Slab-Design-Factors/SlabDesignFactors/executable/params.txt}"
EXECUTABLE_PATH="${5:-SlabDesignFactors/executable/executable_analyze_sharded.jl}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" || -z "${SLURM_ARRAY_TASK_COUNT:-}" ]]; then
    echo "This script must be submitted as a Slurm array job."
    exit 1
fi

echo "Running shard ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_COUNT} with ${SLURM_CPUS_PER_TASK} threads."
echo "Results root: ${RESULTS_ROOT}"
echo "Run name: ${RUN_NAME}"
echo "Completion dir: ${COMPLETION_DIR}"
echo "Params file: ${PARAMS_FILE}"

julia --threads "${SLURM_CPUS_PER_TASK}" \
    "${EXECUTABLE_PATH}" \
    "${RESULTS_ROOT}" \
    "${RUN_NAME}" \
    "${COMPLETION_DIR}" \
    "${PARAMS_FILE}"
