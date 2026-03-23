#!/bin/bash
#
# Submit single-node experiment jobs (executable.slurm) until the run finishes.
#
# The cluster enforces a 12:00:00 walltime on executable.slurm; jobs that hit the
# limit exit without writing the completion file. This wrapper resubmits after each
# job ends until COMPLETION_FILE appears (experiments are resumable via CSV done-sets).
#
# Override any variable below or pass extra args to sbatch after the script name.
#
# Example:
#   bash SlabDesignFactors/executable/executable_wrapper.bash
#
#   RESULTS_PATH=... COMPLETION_FILE=... STUDIES_CSV=max_depths bash .../executable_wrapper.bash
#

set -euo pipefail

SLURM_SCRIPT="${SLURM_SCRIPT:-SlabDesignFactors/executable/executable.slurm}"
RESULTS_PATH="${RESULTS_PATH:-SlabDesignFactors/results/remote_results_experiments}"
COMPLETION_FILE="${COMPLETION_FILE:-${RESULTS_PATH}/experiments_complete.txt}"
STUDIES_CSV="${STUDIES_CSV:-}"
EXECUTABLE_PATH="${EXECUTABLE_PATH:-SlabDesignFactors/executable/executable_experiments.jl}"
MAX_RESUBMISSIONS="${MAX_RESUBMISSIONS:-200}"
POLL_SECONDS="${POLL_SECONDS:-60}"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_FILE:-logs/executable_wrapper_${TIMESTAMP}.log}"

mkdir -p logs
exec > >(tee -i "$LOG_FILE")
exec 2>&1

wait_for_job_completion() {
    local job_id="$1"
    echo "$(date) - Waiting for job ${job_id} to leave the queue..."
    while [[ -n "$(squeue -j "$job_id" -h 2>/dev/null)" ]]; do
        sleep "$POLL_SECONDS"
    done
    echo "$(date) - Job ${job_id} no longer in squeue."
}

echo "SLURM_SCRIPT=${SLURM_SCRIPT}"
echo "RESULTS_PATH=${RESULTS_PATH}"
echo "COMPLETION_FILE=${COMPLETION_FILE}"
echo "STUDIES_CSV=${STUDIES_CSV:-<all default studies>}"
echo "EXECUTABLE_PATH=${EXECUTABLE_PATH}"
echo "MAX_RESUBMISSIONS=${MAX_RESUBMISSIONS}"
echo "POLL_SECONDS=${POLL_SECONDS}"
echo "$(date) - Log (wrapper): ${LOG_FILE}"

mkdir -p "$(dirname "$COMPLETION_FILE")"

attempt=0
while [[ "$attempt" -lt "$MAX_RESUBMISSIONS" ]]; do
    if [[ -f "$COMPLETION_FILE" ]]; then
        echo "$(date) - Completion file found: ${COMPLETION_FILE}"
        echo "$(date) - Done."
        exit 0
    fi

    attempt=$((attempt + 1))
    echo "$(date) - Submission attempt ${attempt}/${MAX_RESUBMISSIONS}"

    job_id=$(sbatch --parsable "$SLURM_SCRIPT" \
        "$RESULTS_PATH" \
        "$COMPLETION_FILE" \
        "${STUDIES_CSV}" \
        "$EXECUTABLE_PATH")

    echo "$(date) - Submitted job ${job_id}"
    wait_for_job_completion "$job_id"
done

echo "$(date) - Reached MAX_RESUBMISSIONS (${MAX_RESUBMISSIONS}) without completion file."
echo "$(date) - Expected: ${COMPLETION_FILE}"
echo "$(date) - Check Slurm logs under outputs/ and Julia output; increase MAX_RESUBMISSIONS if the suite is still progressing."
exit 1
