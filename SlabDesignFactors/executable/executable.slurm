#!/bin/bash

#SBATCH -n 5
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -C centos7
#SBATCH -p sched_any
#SBATCH -o outputs/output_%j.txt
#SBATCH -e outputs/error_%j.txt
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=nhirt@mit.edu

# Set paths
EXECUTABLE_PATH="$1"
RESULTS_PATH="$2"  # Use the first argument as the results path
COMPLETION_FILE="$3"

# Ensure output directories exist
mkdir -p outputs
mkdir -p $RESULTS_PATH

echo "Running Julia script..."

# Run the Julia script
julia -t 3 $EXECUTABLE_PATH $RESULTS_PATH $COMPLETION_FILE

while [ ! -f "$COMPLETION_FILE" ]; do
    sleep 60  # Check every 60 seconds
done

echo "Experiments complete. Terminating SLURM job."

# Cancel the SLURM job
scancel $SLURM_JOB_ID
echo "SLURM job $SLURM_JOB_ID has been cancelled."

# Optionally, you can add a command to clean up or notify
# For example, send a notification email or clean up temporary files

echo "Script completed."