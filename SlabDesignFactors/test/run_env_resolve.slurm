#!/bin/bash
#
# Instantiate / resolve / precompile the Julia project on a compute node.
# Submit from the repository root:
#   cd /path/to/2024_Slab-Design-Factors
#   sbatch SlabDesignFactors/test/run_env_resolve.slurm
#
# Runs `Pkg.instantiate()` (+ precompile) like run_preflight.slurm — uses committed
# Manifest.toml, including path deps AsapOptim / AsapToolkit (see Project.toml [sources]).
#
# To regenerate Manifest.toml from Project.toml (after changing [compat], etc.):
#   export REBUILD_MANIFEST=1
# before sbatch, or add `export REBUILD_MANIFEST=1` below.
#
# Outputs: outputs/slurm_env_resolve_<jobid>.{out,err}

#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -p mit_normal
#SBATCH -o outputs/slurm_env_resolve_%j.out
#SBATCH -e outputs/slurm_env_resolve_%j.err
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

export SLABDESIGN_SKIP_MAKIE="${SLABDESIGN_SKIP_MAKIE:-1}"
export SLABDESIGN_SKIP_REVISE="${SLABDESIGN_SKIP_REVISE:-1}"
export JULIA_NUM_PRECOMPILE_TASKS="${JULIA_NUM_PRECOMPILE_TASKS:-8}"

# Do not inherit REBUILD_MANIFEST=1 from a previous shell session unless explicitly set for this job.
unset REBUILD_MANIFEST
# To force a new Manifest from Project.toml (develop path deps + resolve), uncomment:
# export REBUILD_MANIFEST=1

echo "Julia environment (instantiate + precompile; see run_env_resolve.jl)"
echo "Repo root: ${REPO_ROOT}"
echo "Julia: $(command -v julia)"
julia --version

exec julia --threads "${SLURM_CPUS_PER_TASK:-1}" \
  --project=. \
  SlabDesignFactors/test/run_env_resolve.jl
