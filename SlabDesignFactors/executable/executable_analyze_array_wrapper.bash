#!/bin/bash

set -euo pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_REPO_ROOT="$(cd "${_SCRIPT_DIR}/../.." && pwd)"

# Resubmission wrapper for the sharded full-sweep runner.
# It submits only incomplete shards, waits for completion, and repeats until
# all shard completion markers are present.

SLURM_SCRIPT="${SLURM_SCRIPT:-SlabDesignFactors/executable/executable_analyze_array.slurm}"
RESULTS_ROOT="${RESULTS_ROOT:-SlabDesignFactors/results/remote_results_full_sweep}"
RUN_NAME="${RUN_NAME:-full_sweep}"
COMPLETION_DIR="${COMPLETION_DIR:-${RESULTS_ROOT}/completion}"
PARAMS_FILE="${PARAMS_FILE:-${_REPO_ROOT}/SlabDesignFactors/executable/params.txt}"
EXECUTABLE_PATH="${EXECUTABLE_PATH:-SlabDesignFactors/executable/executable_analyze_sharded.jl}"
TOTAL_SHARDS="${TOTAL_SHARDS:-16}"
MAX_RESUBMISSIONS="${MAX_RESUBMISSIONS:-20}"
POLL_SECONDS="${POLL_SECONDS:-60}"
MERGE_SCRIPT="${MERGE_SCRIPT:-SlabDesignFactors/executable/executable_merge_shards.jl}"
MERGED_DIR="${MERGED_DIR:-${RESULTS_ROOT}/merged}"
# One folder per wrapper invocation (override with JOB_LOG_ID).
JOB_LOG_ID="${JOB_LOG_ID:-analyze_array_$(date +%Y%m%d_%H%M%S)_$$}"
OUT_JOB_DIR="${_REPO_ROOT}/outputs/${JOB_LOG_ID}"
LOG_JOB_DIR="${_REPO_ROOT}/logs/${JOB_LOG_ID}"
mkdir -p "$OUT_JOB_DIR" "$LOG_JOB_DIR"
LOG_FILE="${LOG_FILE:-${LOG_JOB_DIR}/wrapper.log}"

exec > >(tee -i "$LOG_FILE")
exec 2>&1

completion_run_dir="${COMPLETION_DIR}/${RUN_NAME}"
run_completion_file="${completion_run_dir}/run_complete.txt"
merge_done_file="${MERGED_DIR}/merged_${RUN_NAME}.done"

run_merge_if_needed() {
    mkdir -p "$MERGED_DIR"
    if [[ -f "$merge_done_file" ]]; then
        echo "$(date) - Merge skipped: found ${merge_done_file} (delete it to force re-merge)."
        return 0
    fi
    echo "$(date) - Running merge → ${MERGED_DIR}"
    if [[ "${MERGE_SCRIPT}" = /* ]]; then
        _merge_jl="${MERGE_SCRIPT}"
    else
        _merge_jl="${_REPO_ROOT}/${MERGE_SCRIPT}"
    fi
    julia --project="${_REPO_ROOT}" "${_merge_jl}" "$RESULTS_ROOT" "$RUN_NAME" "$MERGED_DIR"
    touch "$merge_done_file"
    echo "$(date) - Merge complete; wrote ${merge_done_file}"
}

build_missing_shard_list() {
    local missing=()
    local i
    for i in $(seq 1 "$TOTAL_SHARDS"); do
        shard_suffix=$(printf "%03d" "$i")
        marker="${completion_run_dir}/shard_${shard_suffix}.complete"
        if [[ ! -f "$marker" ]]; then
            missing+=("$i")
        fi
    done
    if [[ "${#missing[@]}" -eq 0 ]]; then
        echo ""
    else
        (IFS=,; echo "${missing[*]}")
    fi
}

submit_missing_shards() {
    local missing_list="$1"
    echo "$(date) - Submitting missing shard list: ${missing_list}"
    job_id=$(sbatch --parsable --array="${missing_list}" \
        --output="${OUT_JOB_DIR}/slurm_%A_%a.out" \
        --error="${OUT_JOB_DIR}/slurm_%A_%a.err" \
        "$SLURM_SCRIPT" \
        "$RESULTS_ROOT" \
        "$RUN_NAME" \
        "$COMPLETION_DIR" \
        "$PARAMS_FILE" \
        "$EXECUTABLE_PATH")
    echo "$(date) - Submitted array job ${job_id} for shards ${missing_list}"
}

wait_for_job_completion() {
    local job_id="$1"
    echo "$(date) - Waiting for job ${job_id} to finish..."
    while [[ -n "$(squeue -j "$job_id" -h)" ]]; do
        sleep "$POLL_SECONDS"
    done
    echo "$(date) - Job ${job_id} no longer in queue."
}

echo "$(date) - Starting array wrapper."
echo "JOB_LOG_ID=${JOB_LOG_ID}"
echo "OUT_JOB_DIR=${OUT_JOB_DIR}"
echo "LOG_JOB_DIR=${LOG_JOB_DIR}"
echo "LOG_FILE=${LOG_FILE}"
echo "SLURM_SCRIPT=${SLURM_SCRIPT}"
echo "RESULTS_ROOT=${RESULTS_ROOT}"
echo "RUN_NAME=${RUN_NAME}"
echo "COMPLETION_DIR=${COMPLETION_DIR}"
echo "PARAMS_FILE=${PARAMS_FILE}"
echo "EXECUTABLE_PATH=${EXECUTABLE_PATH}"
echo "TOTAL_SHARDS=${TOTAL_SHARDS}"
echo "MAX_RESUBMISSIONS=${MAX_RESUBMISSIONS}"
echo "MERGED_DIR=${MERGED_DIR}"
echo "merge_done_file=${merge_done_file}"

attempt=0
while [[ "$attempt" -lt "$MAX_RESUBMISSIONS" ]]; do
    if [[ -f "$run_completion_file" ]]; then
        echo "$(date) - Run completion file found: ${run_completion_file}"
        echo "$(date) - All shards completed."
        run_merge_if_needed
        exit 0
    fi

    mkdir -p "$completion_run_dir"
    missing_list=$(build_missing_shard_list)
    if [[ -z "$missing_list" ]]; then
        echo "$(date) - No missing shard markers found."
        echo "$(date) - Creating run completion file."
        echo "all shards completed" > "$run_completion_file"
        run_merge_if_needed
        exit 0
    fi

    attempt=$((attempt + 1))
    echo "$(date) - Submission attempt ${attempt}/${MAX_RESUBMISSIONS}"

    submit_missing_shards "$missing_list"
    wait_for_job_completion "$job_id"

    # If all done after this cycle, merge and exit.
    if [[ -f "$run_completion_file" ]]; then
        echo "$(date) - Run completion file found after attempt ${attempt}."
        run_merge_if_needed
        exit 0
    fi
done

echo "$(date) - Maximum resubmissions reached (${MAX_RESUBMISSIONS})."
echo "$(date) - Please inspect shard markers in ${completion_run_dir} and rerun if needed."
exit 1
