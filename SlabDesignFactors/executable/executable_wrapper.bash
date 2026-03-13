#!/bin/bash

# Configuration
SLURM_SCRIPT="SlabDesignFactors/executable/executable.slurm"
MAX_RESUBMISSIONS=10
COUNTER=0
EMAIL="nhirt@mit.edu"
RESULTS_PATH="SlabDesignFactors/results/remote_results_nodeflection_noslabmin/"
COMPLETION_FILE="${RESULTS_PATH}experiments_complete.txt"
LOG_FILE="logs/slurm_monitor.log"
EXECUTABLE_PATH="SlabDesignFactors/executable/executable_experiments.jl"

# Log setup
exec > >(tee -i $LOG_FILE)
exec 2>&1

# Function to submit a job
submit_job() {
    echo "$(date) - Submitting job..."
    JOB_ID=$(sbatch --parsable "$SLURM_SCRIPT" "$EXECUTABLE_PATH" "$RESULTS_PATH" "$COMPLETION_FILE" 2>&1)
    if [ $? -ne 0 ]; then
        echo "$(date) - Failed to submit job: $JOB_ID"
        echo -e "Job submission failed:\n$JOB_ID" | mail -s "SLURM Job Submission Error" $EMAIL
        exit 1
    fi
    echo "$(date) - Job $JOB_ID submitted successfully."
    echo -e "Job $JOB_ID submitted.\nScript: $SLURM_SCRIPT" | mail -s "SLURM Job Submitted" $EMAIL
}

# Initial job submission
submit_job

# Monitoring loop
while [ $COUNTER -lt $MAX_RESUBMISSIONS ] && [ ! -f "$COMPLETION_FILE" ]; do
    echo "$(date) - Monitoring job $JOB_ID..."

    # Wait for job completion
    while squeue --job "$JOB_ID" > /dev/null 2>&1; do
        sleep 60  # Check job status every minute
    done

    # Check job state
    JOB_STATE=$(sacct -j "$JOB_ID" --format=State --noheader | head -n 1 | tr -d ' ')
    if [[ "$JOB_STATE" == "FAILED" || "$JOB_STATE" == "TIMEOUT" ]]; then
        echo "$(date) - Job $JOB_ID failed with state $JOB_STATE."
        FULL_ERROR=$(sacct -j "$JOB_ID" --format=JobID,JobName,State,ExitCode --noheader)
        echo -e "Job $JOB_ID failed with state $JOB_STATE:\n$FULL_ERROR" | mail -s "SLURM Job Failure" $EMAIL
    elif [ "$JOB_STATE" == "COMPLETED" ]; then
        echo "$(date) - Job $JOB_ID completed successfully."
        echo -e "Job $JOB_ID completed." | mail -s "SLURM Job Completed" $EMAIL
    else
        echo "$(date) - Unknown state: $JOB_STATE"
    fi

    # Check for completion file
    if [ -f "$COMPLETION_FILE" ]; then
        echo "$(date) - Completion file found. Ending monitoring."
        echo -e "Analysis completed. No further resubmissions." | mail -s "Analysis Completed" $EMAIL
        exit 0
    fi

    # Resubmit the job if the limit is not reached
    COUNTER=$((COUNTER + 1))
    if [ $COUNTER -lt $MAX_RESUBMISSIONS ]; then
        echo "$(date) - Resubmitting job ($COUNTER/$MAX_RESUBMISSIONS)..."
        submit_job
    else
        echo "$(date) - Maximum resubmissions reached."
        echo -e "Maximum resubmissions ($MAX_RESUBMISSIONS) reached. No further submissions." | mail -s "SLURM Resubmission Limit Reached" $EMAIL
        exit 0
    fi
done