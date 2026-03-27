#!/bin/bash

#SBATCH -N 1
#SBATCH --cpus-per-task=8
# https://orcd-docs.mit.edu/running-jobs/requesting-resources/#memory
#SBATCH --mem=64G
#SBATCH --array=1-6
#SBATCH -t 12:00:00
#SBATCH -p mit_normal
# Direct sbatch: flat names. Wrappers override with outputs/<JOB_LOG_ID>/...
#SBATCH -o outputs/slurm_experiments_%A_%a.out
#SBATCH -e outputs/slurm_experiments_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=nhirt@mit.edu

set -euo pipefail

module load community-modules
module load julia
module load gurobi

export SLABDESIGN_SKIP_MAKIE="${SLABDESIGN_SKIP_MAKIE:-1}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
  REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi
[[ -f "${REPO_ROOT}/env_gurobi.sh" ]] && source "${REPO_ROOT}/env_gurobi.sh"

mkdir -p outputs

RESULTS_PATH="${1:-SlabDesignFactors/results/remote_results_experiments}"
COMPLETION_DIR="${2:-${RESULTS_PATH}/completion}"
STUDY_LIST_CSV="${3:-max_depths,strip_resolution,constrained_inventory,nlp_solver_comparison,material_scenario_mc,validation_mip}"
EXECUTABLE_PATH="${4:-SlabDesignFactors/executable/executable_experiments.jl}"

IFS=',' read -r -a STUDIES <<< "$STUDY_LIST_CSV"
N_STUDIES="${#STUDIES[@]}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "This script must be submitted as a Slurm array job."
    exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID}"
if (( TASK_ID < 1 || TASK_ID > N_STUDIES )); then
    echo "Task index ${TASK_ID} is outside study list length ${N_STUDIES}; exiting."
    exit 0
fi

STUDY="${STUDIES[$((TASK_ID-1))]}"
STUDY="$(echo "$STUDY" | xargs)"
if [[ -z "$STUDY" ]]; then
    echo "Resolved empty study for task ${TASK_ID}; exiting."
    exit 1
fi

mkdir -p "$COMPLETION_DIR"
COMPLETION_FILE="${COMPLETION_DIR}/${STUDY}.complete"

echo "Running experiments study '${STUDY}' (task ${TASK_ID}/${N_STUDIES})"
echo "Results path: ${RESULTS_PATH}"
echo "Completion file: ${COMPLETION_FILE}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"

julia --threads "${SLURM_CPUS_PER_TASK}" \
    "${EXECUTABLE_PATH}" \
    "${RESULTS_PATH}" \
    "${COMPLETION_FILE}" \
    "${STUDY}"
