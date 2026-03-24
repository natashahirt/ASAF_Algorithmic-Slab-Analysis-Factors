#!/usr/bin/env bash
# Source this file before Julia when using a repo-local license:
#   source env_gurobi.sh
# Slurm scripts in SlabDesignFactors/executable/ source it automatically.
_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
if [[ -f "${_root}/secrets/gurobi.lic" ]]; then
  export GRB_LICENSE_FILE="${_root}/secrets/gurobi.lic"
fi
