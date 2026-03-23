# Slab executable runners (cluster / Slurm)

This folder now includes three execution patterns:

1. **`executable.slurm`** (single node): one Julia process runs **all** default studies (or an optional `studies_csv`) with threading. Wrapper: `executable_wrapper.bash`.
2. **`executable_experiments_array.slurm`**: Slurm **array** — one task per study (recommended for 12 h walltime / large suites). Wrapper: `executable_experiments_array_wrapper.bash`.
3. **`executable_analyze_array.slurm`**: sharded **full sweep** from `params.txt`. Wrapper: `executable_analyze_array_wrapper.bash` (optional auto-merge + `merged_<run>.done` guard).

All Slurm scripts in this folder use **`-p mit_normal`** and **`module load julia`** + **`module load gurobi`** (adjust module names on your cluster if needed).

## Experiments runner (`executable_experiments.jl`)

### Included studies

- `max_depths` (topology JSONs)
- `strip_resolution` sensitivity
- `constrained_inventory` sensitivity (resized + max-depth substudies, validation JSONs)
- `nlp_solver_comparison` (solver benchmark table including `MIP` baseline plus NLP methods with common MIP warm-start)
- `material_scenario_mc` (concrete material scenarios × slab sizers with Monte Carlo on ECCs)
- `validation_mip` (validation JSONs with drawn geometry, perimeter beams, holes, element IDs, building-type loads, unconstrained MIP sizing)

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

### Validation MIP study

`validation_mip` runs the three validation building geometries (office, school,
warehouse) with `drawn=true` so that `type_information` is available. It passes
`i_perimeter`, `i_holes`, and `element_ids` to the analysis and sizing params,
uses building-specific live loads and `beam_sizer=:discrete` with
`n_max_sections=0` (unconstrained MIP). Results are appended via
`append_results_to_csv`.

### Single-node Slurm (`executable.slurm`)

Submit one job that runs the experiment runner (all studies by default):

`sbatch SlabDesignFactors/executable/executable.slurm [results_path] [completion_file] [studies_csv] [executable_path]`

Or use the wrapper (timestamped log under `logs/`). The wrapper **resubmits** after each job finishes until `COMPLETION_FILE` exists (so 12h walltime limits do not stop a long resumable suite):

`bash SlabDesignFactors/executable/executable_wrapper.bash`

Environment overrides: `RESULTS_PATH`, `COMPLETION_FILE`, `STUDIES_CSV`, `EXECUTABLE_PATH`, `SLURM_SCRIPT`, plus `MAX_RESUBMISSIONS` (default 200), `POLL_SECONDS` (default 60).

### Usage (CLI)

`julia SlabDesignFactors/executable/executable_experiments.jl <results_path> <completion_file> [studies_csv]`

Examples:

- Run all included studies (default):
  - `julia SlabDesignFactors/executable/executable_experiments.jl SlabDesignFactors/results/remote_results_experiments/ SlabDesignFactors/results/remote_results_experiments/experiments_complete.txt`
- Run only selected studies:
  - `julia SlabDesignFactors/executable/executable_experiments.jl SlabDesignFactors/results/remote_results_experiments/ SlabDesignFactors/results/remote_results_experiments/experiments_complete.txt max_depths,strip_resolution`
  - `julia SlabDesignFactors/executable/executable_experiments.jl SlabDesignFactors/results/remote_results_experiments/ SlabDesignFactors/results/remote_results_experiments/experiments_complete.txt validation_mip`

### Slurm array (parallel by study)

- Slurm script: `SlabDesignFactors/executable/executable_experiments_array.slurm`
- Wrapper: `SlabDesignFactors/executable/executable_experiments_array_wrapper.bash`

Submit directly:

`sbatch SlabDesignFactors/executable/executable_experiments_array.slurm <results_path> <completion_dir> <study_list_csv> <executable_path>`

Example:

`sbatch SlabDesignFactors/executable/executable_experiments_array.slurm SlabDesignFactors/results/remote_results_experiments SlabDesignFactors/results/remote_results_experiments/completion max_depths,strip_resolution,constrained_inventory,nlp_solver_comparison,material_scenario_mc,validation_mip SlabDesignFactors/executable/executable_experiments.jl`

Wrapper (recommended for repeated 12-hour retries):

`bash SlabDesignFactors/executable/executable_experiments_array_wrapper.bash`

Optional overrides:

`RESULTS_PATH=SlabDesignFactors/results/remote_results_experiments COMPLETION_DIR=SlabDesignFactors/results/remote_results_experiments/completion STUDY_LIST_CSV=max_depths,strip_resolution,constrained_inventory,nlp_solver_comparison,material_scenario_mc,validation_mip MAX_RESUBMISSIONS=20 POLL_SECONDS=60 bash SlabDesignFactors/executable/executable_experiments_array_wrapper.bash`

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
- repeats until `run_complete.txt` exists or max retries are reached,
- automatically runs the merge step (`executable_merge_shards.jl`) once all shards complete,
- writes `MERGED_DIR/merged_<run_name>.done` after a successful merge and **skips** re-merging if that file exists (delete it to force a fresh merge).

Both wrapper scripts write timestamped log files (e.g. `logs/analyze_array_wrapper_20260323_140000.log`).

### Full-sweep load defaults (shared with `max_depths`)

Topology sweep load factors, SDL/live in ksi, strip spacing, and serviceability limit
are defined once in `SlabDesignFactors/core/full_sweep_defaults.jl` as `FULL_SWEEP_*`
and used by `executable_analyze_sharded.jl` and the topology experiments in
`executable_experiments.jl` (via `DEFAULT_*` aliases). Validation studies
(`constrained_inventory`, `validation_mip`) keep their own building-specific loads.

### Restart behavior after 12-hour timeout

Re-submit the same array command. Each task:
- reloads its shard CSV done-set,
- skips completed configs,
- continues from unfinished work.

The run is idempotent; reruns are expected and safe.

### Merge shards after completion

**Manual merge (required if you do not use `executable_analyze_array_wrapper.bash`):**
If you submit **`executable_analyze_array.slurm` directly** (or any custom orchestration), the wrapper’s automatic merge and `merged_<run_name>.done` guard do **not** run. After **all** shard markers exist and `run_complete.txt` is present, merge yourself:

`julia SlabDesignFactors/executable/executable_merge_shards.jl <results_root> <run_name> <merged_dir>`

Example:

`julia SlabDesignFactors/executable/executable_merge_shards.jl SlabDesignFactors/results/remote_results_full_sweep full_sweep SlabDesignFactors/results/remote_results_full_sweep/merged`

The analyze **wrapper** runs this step for you once the run is complete (unless `merged_<run_name>.done` already exists in `MERGED_DIR`).

## Notes on core usage

- The sharded full sweep now matches the `max_depths` experiment sizing policy
  (including `deflection_limit=true` and `slab_dead_load=0` / material-driven slab
  weight). Older sharded CSVs from runs with `deflection_limit=false` are not
  comparable without re-running.
- `executable.slurm`, `executable_experiments_array.slurm`, and `executable_analyze_array.slurm` use `-p mit_normal` and load `julia` + `gurobi`.
- `executable_analyze_array.slurm` currently uses `--cpus-per-task=8` and `--array=1-16`.
- `executable_experiments_array.slurm` uses `--array=1-6` (one task per study).
- Start with these defaults, then tune:
  - increase array size for more distributed parallelism,
  - increase `cpus-per-task` for more intra-node threading.
- The Julia runner sets BLAS threads to 1 to avoid nested oversubscription.
- `SlabDesignFactors.jl` imports `Gurobi`; ensure `module load gurobi` and a valid license on the cluster.
