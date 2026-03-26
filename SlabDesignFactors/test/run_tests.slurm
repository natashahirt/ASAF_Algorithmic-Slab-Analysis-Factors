#!/bin/bash

#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH -p mit_normal
#SBATCH -o outputs/slurm_tests_%j.out
#SBATCH -e outputs/slurm_tests_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=nhirt@mit.edu

set -euo pipefail

module load community-modules
module load julia/1.10.4
module load gurobi

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
  REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi

cd "${REPO_ROOT}"
mkdir -p outputs
[[ -f "${REPO_ROOT}/env_gurobi.sh" ]] && source "${REPO_ROOT}/env_gurobi.sh"

# Cluster-safe Julia test settings.
export SLABDESIGN_SKIP_MAKIE="${SLABDESIGN_SKIP_MAKIE:-1}"
export SLABDESIGN_SKIP_REVISE="${SLABDESIGN_SKIP_REVISE:-1}"
export JULIA_NUM_PRECOMPILE_TASKS="${JULIA_NUM_PRECOMPILE_TASKS:-1}"

echo "Running SlabDesignFactors tests via Slurm batch node"
echo "Repo root: ${REPO_ROOT}"
echo "Threads: ${SLURM_CPUS_PER_TASK:-1}"

julia --project=. -e 'using Pkg; Pkg.instantiate()'

julia --threads "${SLURM_CPUS_PER_TASK:-1}" \
  --compiled-modules=no \
  --project=. \
  SlabDesignFactors/test/run.jl
