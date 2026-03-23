#!/bin/bash
#
# Single-node Slurm job: run the full experiment suite (or a subset) in one
# Julia process with threading. Walltime is capped at 12h; use executable_wrapper.bash
# to resubmit until experiments_complete.txt is written (resumable CSVs).
# For parallel-by-study runs, use executable_experiments_array.slurm instead.
#
# Usage (positional args to this batch script after sbatch):
#   sbatch executable.slurm [results_path] [completion_file] [studies_csv] [executable_path]
#
#   studies_csv  — optional; omit to run all default studies (see executable_experiments.jl).
#   executable_path — optional path to .jl file (default: executable_experiments.jl).
#
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH -p mit_normal
#SBATCH -o outputs/executable_%j.out
#SBATCH -e outputs/executable_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=nhirt@mit.edu

set -euo pipefail

module load julia
module load gurobi

mkdir -p outputs

RESULTS_PATH="${1:-SlabDesignFactors/results/remote_results_experiments}"
COMPLETION_FILE="${2:-${RESULTS_PATH}/experiments_complete.txt}"
STUDIES_CSV="${3:-}"
EXECUTABLE_PATH="${4:-SlabDesignFactors/executable/executable_experiments.jl}"

mkdir -p "$(dirname "$COMPLETION_FILE")"

echo "Results path: ${RESULTS_PATH}"
echo "Completion file: ${COMPLETION_FILE}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"
if [[ -n "$STUDIES_CSV" ]]; then
    echo "Studies: ${STUDIES_CSV}"
    julia --threads "${SLURM_CPUS_PER_TASK}" \
        "${EXECUTABLE_PATH}" \
        "${RESULTS_PATH}" \
        "${COMPLETION_FILE}" \
        "${STUDIES_CSV}"
else
    echo "Studies: (all defaults)"
    julia --threads "${SLURM_CPUS_PER_TASK}" \
        "${EXECUTABLE_PATH}" \
        "${RESULTS_PATH}" \
        "${COMPLETION_FILE}"
fi
