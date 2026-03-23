# Slab executable runners (cluster / Slurm)

This folder now includes two execution patterns:

1. `executable_experiments.jl` + `executable.slurm`: single-job threaded, resumable experiment suite (`max_depths`, `strip_resolution`, `constrained_inventory`, `nlp_solver_comparison`, `material_scenario_mc`).
2. `executable_analyze_sharded.jl` + `executable_analyze_array.slurm`: full-sweep, resumable, hybrid distributed+multicore runner for 12-hour walltime clusters.

## Experiments runner (`executable_experiments.jl`)

### Included studies

- `max_depths` (topology JSONs)
- `strip_resolution` sensitivity
- `constrained_inventory` sensitivity (resized + max-depth substudies)
- `nlp_solver_comparison` (solver benchmark table including `MIP` baseline plus NLP methods with common MIP warm-start)
- `material_scenario_mc` (concrete material scenarios × slab sizers with Monte Carlo on ECCs)

### NLP solver postprocess utility

After `nlp_solver_comparison.csv` is generated, compute per-geometry and per-solver deltas vs MIP:

`julia SlabDesignFactors/executable/executable_nlp_solver_summary.jl <results_path> [output_path]`

Example:

`julia SlabDesignFactors/executable/executable_nlp_solver_summary.jl SlabDesignFactors/results/remote_results_experiments SlabDesignFactors/results/remote_results_experiments`

Outputs:
- `nlp_solver_vs_mip_per_geometry.csv`
- `nlp_solver_vs_mip_summary.csv`

### Material scenario Monte Carlo study

`material_scenario_mc` runs structural re-analysis for each concrete material
scenario (normal-weight 4 ksi, normal-weight 5 ksi, lightweight 4 ksi,
lightweight 3 ksi) × slab sizer, then performs 10,000-sample Monte Carlo
simulations varying ECCs independently for steel, concrete, and rebar (and
jointly). Output columns include:

- Per-scenario structural metrics (`norm_mass_beams`, `ec_total`, `max_deflection_in`, etc.)
- Independent MC summaries per material (`mc_steel_*`, `mc_conc_*`, `mc_rebar_*`)
- Joint MC summary (`mc_joint_*`: mean, std, p5/p25/p50/p75/p95)
- Variance decomposition (`var_frac_steel`, `var_frac_conc`, `var_frac_rebar`)

Concrete scenarios are defined in `SlabDesignFactors/core/materials.jl` as
`ConcreteMaterial` structs with pre-converted SI and imperial properties. The
material's density propagates into both sizing (load generation) and
post-processing (mass/ECC), so lightweight scenarios produce genuinely different
structural designs.

### Usage

`julia SlabDesignFactors/executable/executable_experiments.jl <results_path> <completion_file> [studies_csv]`

Examples:

- Run all included studies (default):
  - `julia SlabDesignFactors/executable/executable_experiments.jl SlabDesignFactors/results/remote_results_experiments/ SlabDesignFactors/results/remote_results_experiments/experiments_complete.txt`
- Run only selected studies:
  - `julia SlabDesignFactors/executable/executable_experiments.jl SlabDesignFactors/results/remote_results_experiments/ SlabDesignFactors/results/remote_results_experiments/experiments_complete.txt max_depths,strip_resolution`
  - `julia SlabDesignFactors/executable/executable_experiments.jl SlabDesignFactors/results/remote_results_experiments/ SlabDesignFactors/results/remote_results_experiments/experiments_complete.txt nlp_solver_comparison`

### Slurm array (parallel by study)

- Slurm script: `SlabDesignFactors/executable/executable_experiments_array.slurm`
- Wrapper: `SlabDesignFactors/executable/executable_experiments_array_wrapper.bash`

Submit directly:

`sbatch SlabDesignFactors/executable/executable_experiments_array.slurm <results_path> <completion_dir> <study_list_csv> <executable_path>`

Example:

`sbatch SlabDesignFactors/executable/executable_experiments_array.slurm SlabDesignFactors/results/remote_results_experiments SlabDesignFactors/results/remote_results_experiments/completion max_depths,strip_resolution,constrained_inventory,nlp_solver_comparison,material_scenario_mc SlabDesignFactors/executable/executable_experiments.jl`

Wrapper (recommended for repeated 12-hour retries):

`bash SlabDesignFactors/executable/executable_experiments_array_wrapper.bash`

Optional overrides:

`RESULTS_PATH=SlabDesignFactors/results/remote_results_experiments COMPLETION_DIR=SlabDesignFactors/results/remote_results_experiments/completion STUDY_LIST_CSV=max_depths,strip_resolution,constrained_inventory,nlp_solver_comparison,material_scenario_mc MAX_RESUBMISSIONS=20 POLL_SECONDS=60 bash SlabDesignFactors/executable/executable_experiments_array_wrapper.bash`

## Full sweep runner (recommended on ORCD)

### What it does

- Builds the full config set from `params.txt` (geometry groups + slab combinations).
- Deterministically shards work by config index across array tasks.
- Uses Julia threading inside each array task.
- Writes per-shard CSVs (no cross-task file collisions).
- Resumes from existing shard CSVs by rebuilding a done-set at startup.
- Emits per-shard progress JSON and completion marker files.

### Output layout

Given:
- `RESULTS_ROOT=SlabDesignFactors/results/remote_results_full_sweep`
- `RUN_NAME=full_sweep`

The runner writes:
- `RESULTS_ROOT/RUN_NAME/shards/<results_name>_shard_<id>.csv`
- `RESULTS_ROOT/RUN_NAME/progress/progress_shard_<id>.json`
- `COMPLETION_DIR/RUN_NAME/shard_<id>.complete`
- `COMPLETION_DIR/RUN_NAME/run_complete.txt` (when all shards complete)

### Submit command

From repo root on ORCD:

`sbatch SlabDesignFactors/executable/executable_analyze_array.slurm <results_root> <run_name> <completion_dir> <params_file> <executable_path>`

All arguments are optional; defaults are defined in the `.slurm` file.

Example:

`sbatch SlabDesignFactors/executable/executable_analyze_array.slurm SlabDesignFactors/results/remote_results_full_sweep full_sweep SlabDesignFactors/results/remote_results_full_sweep/completion /home/nhirt/2024_Slab-Design-Factors/SlabDesignFactors/executable/params.txt SlabDesignFactors/executable/executable_analyze_sharded.jl`

### Wrapper for repeated 12-hour resubmissions

Use the wrapper to submit only incomplete shards until the run is complete:

`bash SlabDesignFactors/executable/executable_analyze_array_wrapper.bash`

Optional environment overrides:

`TOTAL_SHARDS=16 MAX_RESUBMISSIONS=20 POLL_SECONDS=60 RESULTS_ROOT=SlabDesignFactors/results/remote_results_full_sweep RUN_NAME=full_sweep COMPLETION_DIR=SlabDesignFactors/results/remote_results_full_sweep/completion PARAMS_FILE=/home/nhirt/2024_Slab-Design-Factors/SlabDesignFactors/executable/params.txt EXECUTABLE_PATH=SlabDesignFactors/executable/executable_analyze_sharded.jl bash SlabDesignFactors/executable/executable_analyze_array_wrapper.bash`

The wrapper:
- checks which shard completion markers are missing,
- submits only missing shard indices via `--array`,
- waits for that submission to finish,
- repeats until `run_complete.txt` exists or max retries are reached.

### Restart behavior after 12-hour timeout

Re-submit the same array command. Each task:
- reloads its shard CSV done-set,
- skips completed configs,
- continues from unfinished work.

The run is idempotent; reruns are expected and safe.

### Merge shards after completion

When `run_complete.txt` exists, merge shard outputs:

`julia SlabDesignFactors/executable/executable_merge_shards.jl <results_root> <run_name> <merged_dir>`

Example:

`julia SlabDesignFactors/executable/executable_merge_shards.jl SlabDesignFactors/results/remote_results_full_sweep full_sweep SlabDesignFactors/results/remote_results_full_sweep/merged`

## Notes on core usage

- `executable_analyze_array.slurm` currently uses `--cpus-per-task=8` and `--array=1-16`.
- Start with these defaults, then tune:
  - increase array size for more distributed parallelism,
  - increase `cpus-per-task` for more intra-node threading.
- The Julia runner sets BLAS threads to 1 to avoid nested oversubscription.
