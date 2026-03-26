#!/usr/bin/env bash
# Run tests on shared login / low-ulimit nodes where LLVM/lld fails during
# package precompile ("Resource temporarily unavailable" while linking Asap, etc.).
#
# --compiled-modules=no avoids building/loading native pkgimages (no lld); startup is slower.
# SLABDESIGN_SKIP_REVISE=1 avoids Revise precompile (same class of failure).
# JULIA_NUM_PRECOMPILE_TASKS=1 limits parallel precompile if you omit --compiled-modules=no.
#
# Use ONE Julia install per shell: do not mix `module load julia` with juliaup on PATH
# (mixed LD_LIBRARY_PATH can break linking).
#
# Usage (from repo root):
#   bash SlabDesignFactors/test/run_cluster_login.sh
#   # or:
#   chmod +x SlabDesignFactors/test/run_cluster_login.sh
#   ./SlabDesignFactors/test/run_cluster_login.sh

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

export SLABDESIGN_SKIP_REVISE="${SLABDESIGN_SKIP_REVISE:-1}"
export JULIA_NUM_PRECOMPILE_TASKS="${JULIA_NUM_PRECOMPILE_TASKS:-1}"

exec julia --compiled-modules=no --project=. SlabDesignFactors/test/run.jl
