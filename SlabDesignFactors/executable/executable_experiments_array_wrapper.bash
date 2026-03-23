#!/bin/bash

set -euo pipefail

# Resubmission wrapper for experiments array runs.
# Parallelization is by study (one array task per study).

SLURM_SCRIPT="${SLURM_SCRIPT:-SlabDesignFactors/executable/executable_experiments_array.slurm}"
RESULTS_PATH="${RESULTS_PATH:-SlabDesignFactors/results/remote_results_experiments}"
COMPLETION_DIR="${COMPLETION_DIR:-${RESULTS_PATH}/completion}"
STUDY_LIST_CSV="${STUDY_LIST_CSV:-max_depths,strip_resolution,constrained_inventory,nlp_solver_comparison,material_scenario_mc}"
EXECUTABLE_PATH="${EXECUTABLE_PATH:-SlabDesignFactors/executable/executable_experiments.jl}"
MAX_RESUBMISSIONS="${MAX_RESUBMISSIONS:-20}"
POLL_SECONDS="${POLL_SECONDS:-60}"
LOG_FILE="${LOG_FILE:-logs/experiments_array_wrapper.log}"

mkdir -p logs
exec > >(tee -i "$LOG_FILE")
exec 2>&1

IFS=',' read -r -a STUDIES <<< "$STUDY_LIST_CSV"
N_STUDIES="${#STUDIES[@]}"

build_missing_index_list() {
    local missing=()
    local idx
    for ((idx=1; idx<=N_STUDIES; idx++)); do
        study="$(echo "${STUDIES[$((idx-1))]}" | xargs)"
        marker="${COMPLETION_DIR}/${study}.complete"
        if [[ ! -f "$marker" ]]; then
            missing+=("$idx")
        fi
    done
    if [[ "${#missing[@]}" -eq 0 ]]; then
        echo ""
    else
        (IFS=,; echo "${missing[*]}")
    fi
}

all_studies_complete() {
    local idx
    for ((idx=1; idx<=N_STUDIES; idx++)); do
        study="$(echo "${STUDIES[$((idx-1))]}" | xargs)"
        marker="${COMPLETION_DIR}/${study}.complete"
        if [[ ! -f "$marker" ]]; then
            return 1
        fi
    done
    return 0
}

submit_missing_studies() {
    local missing_indices="$1"
    echo "$(date) - Submitting study indices: ${missing_indices}"
    job_id=$(sbatch --parsable --array="${missing_indices}" \
        "$SLURM_SCRIPT" \
        "$RESULTS_PATH" \
        "$COMPLETION_DIR" \
        "$STUDY_LIST_CSV" \
        "$EXECUTABLE_PATH")
    echo "$(date) - Submitted job ${job_id}"
}

wait_for_job_completion() {
    local job_id="$1"
    echo "$(date) - Waiting for job ${job_id}..."
    while [[ -n "$(squeue -j "$job_id" -h)" ]]; do
        sleep "$POLL_SECONDS"
    done
    echo "$(date) - Job ${job_id} left queue."
}

mkdir -p "$COMPLETION_DIR"
echo "$(date) - Starting experiments array wrapper"
echo "SLURM_SCRIPT=${SLURM_SCRIPT}"
echo "RESULTS_PATH=${RESULTS_PATH}"
echo "COMPLETION_DIR=${COMPLETION_DIR}"
echo "STUDY_LIST_CSV=${STUDY_LIST_CSV}"
echo "EXECUTABLE_PATH=${EXECUTABLE_PATH}"
echo "MAX_RESUBMISSIONS=${MAX_RESUBMISSIONS}"

attempt=0
while [[ "$attempt" -lt "$MAX_RESUBMISSIONS" ]]; do
    if all_studies_complete; then
        echo "$(date) - All study completion files found."
        exit 0
    fi

    missing_indices=$(build_missing_index_list)
    if [[ -z "$missing_indices" ]]; then
        echo "$(date) - No missing indices found."
        exit 0
    fi

    attempt=$((attempt + 1))
    echo "$(date) - Attempt ${attempt}/${MAX_RESUBMISSIONS}"
    submit_missing_studies "$missing_indices"
    wait_for_job_completion "$job_id"
done

echo "$(date) - Reached max resubmissions (${MAX_RESUBMISSIONS}) with incomplete studies."
echo "$(date) - Check completion markers in ${COMPLETION_DIR}."
exit 1
