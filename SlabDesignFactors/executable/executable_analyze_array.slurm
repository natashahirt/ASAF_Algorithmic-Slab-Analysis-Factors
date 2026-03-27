#!/bin/bash

#SBATCH -N 1
#SBATCH --cpus-per-task=8
# Explicit RAM (partition default OOM’d Julia+Gurobi shards). Tune after measuring MaxRSS:
#   sacct -j JOBID -o JobID,State,ReqMem,MaxRSS --units=G
# https://orcd-docs.mit.edu/running-jobs/requesting-resources/#memory
#SBATCH --mem=64G
#SBATCH --array=1-16
#SBATCH -t 12:00:00
#SBATCH -p mit_normal
# Direct sbatch: flat names under outputs/ (no extra subdir). Wrappers override with outputs/<JOB_LOG_ID>/...
#SBATCH -o outputs/slurm_analyze_%A_%a.out
#SBATCH -e outputs/slurm_analyze_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=nhirt@mit.edu

set -euo pipefail

module load community-modules
module load julia
module load gurobi

# Batch nodes: skip CairoMakie (avoids Glib_jll / system glib symbol mismatches like g_dir_unref).
export SLABDESIGN_SKIP_MAKIE="${SLABDESIGN_SKIP_MAKIE:-1}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
  REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
[[ -f "${REPO_ROOT}/env_gurobi.sh" ]] && source "${REPO_ROOT}/env_gurobi.sh"

mkdir -p outputs

RESULTS_ROOT="${1:-SlabDesignFactors/results/remote_results_full_sweep}"
RUN_NAME="${2:-full_sweep}"
COMPLETION_DIR="${3:-${RESULTS_ROOT}/completion}"
PARAMS_FILE="${4:-${REPO_ROOT}/SlabDesignFactors/executable/params.txt}"
EXECUTABLE_PATH="${5:-SlabDesignFactors/executable/executable_analyze_sharded.jl}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "This script must be submitted as a Slurm array job."
    exit 1
fi
# Some Slurm builds omit SLURM_ARRAY_TASK_COUNT; derive span from MIN/MAX when needed.
if [[ -z "${SLURM_ARRAY_TASK_COUNT:-}" && -n "${SLURM_ARRAY_TASK_MIN:-}" && -n "${SLURM_ARRAY_TASK_MAX:-}" ]]; then
    export SLURM_ARRAY_TASK_COUNT=$((SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1))
fi
if [[ -z "${SLURM_ARRAY_TASK_COUNT:-}" ]]; then
    echo "Could not determine array task count (set SLURM_ARRAY_TASK_COUNT or use Slurm that exports it)."
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
