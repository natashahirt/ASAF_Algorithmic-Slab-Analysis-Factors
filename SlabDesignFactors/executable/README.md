# Cluster runs (Slurm)

Submit **from the repository root** so `SLURM_SUBMIT_DIR`, `Project.toml`, and `secrets/gurobi.lic` resolve correctly.

All `.slurm` scripts here use `-p mit_normal`, then `module load community-modules`, `julia`, and `gurobi` (ORCD-style) — edit those lines if your site differs.

## Two independent pipelines

| Pipeline | What it runs | Default results directory | Slurm script | Wrapper (resubmit until done) |
|----------|----------------|---------------------------|--------------|--------------------------------|
| **Experiments** | Named studies (`max_depths`, `validation_mip`, …) | `SlabDesignFactors/results/remote_results_experiments/` | `executable_experiments_array.slurm` | `executable_experiments_array_wrapper.bash` |
| **Full sweep** | Cartesian product from `params.txt` | Wrapper default: `SlabDesignFactors/results/remote_results_full_sweep_<JOB_LOG_ID>/` (direct `sbatch`: `.../remote_results_full_sweep/`) | `executable_analyze_array.slurm` | `executable_analyze_array_wrapper.bash` |

They do **not** depend on each other.

### Experiments (array: one task per study)

```bash
cd /path/to/2024_Slab-Design-Factors
sbatch SlabDesignFactors/executable/executable_experiments_array.slurm
```

Optional args: `results_path`, `completion_dir`, `study1,study2,...`, `executable_experiments.jl`.

For a one-off local run (no Slurm), use the CLI under **Studies** below.

### Full sweep / analyze (array: one task per shard)

```bash
cd /path/to/2024_Slab-Design-Factors
sbatch SlabDesignFactors/executable/executable_analyze_array.slurm
```

Positional args (all optional; defaults match the table above): `results_root`, `run_name`, `completion_dir`, `params.txt`, `executable_analyze_sharded.jl`.  
Default `params.txt` is next to the sharded driver; geometry paths inside it are **relative to the repo root** (see comment in `params.txt`).

**Merge shards** (only if you did not use the analyze wrapper):  
`julia SlabDesignFactors/executable/executable_merge_shards.jl <results_root> <run_name> <merged_dir>`

---

## Studies (`executable_experiments.jl`)

- `max_depths` — topology JSONs  
- `strip_resolution`, `constrained_inventory`, `nlp_solver_comparison`, `material_scenario_mc`, `validation_mip` — see source docstring for detail  

CLI:

`julia SlabDesignFactors/executable/executable_experiments.jl <results_path> <completion_file> [studies_csv]`

**NLP vs MIP summary** (after `nlp_solver_comparison.csv` exists):

`julia SlabDesignFactors/executable/executable_nlp_solver_summary.jl <results_path> [output_path]`

---

## Full sweep layout

With `RESULTS_ROOT=.../remote_results_full_sweep` (or `.../remote_results_full_sweep_<JOB_LOG_ID>` from the analyze wrapper) and `RUN_NAME=full_sweep`:

- Shards: `RESULTS_ROOT/RUN_NAME/shards/shard_<id>/<name>.csv` (one folder per shard)  
- Progress: `RESULTS_ROOT/RUN_NAME/progress/progress_shard_<id>.json`  
- Done markers: `COMPLETION_DIR/RUN_NAME/shard_<id>.complete`, then `run_complete.txt`  

Re-submitting the same array job is safe (resume from CSVs). Shared load defaults live in `SlabDesignFactors/core/full_sweep_defaults.jl`.

---

## Wrappers

Both wrappers resubmit until the run finishes (12 h walltime chunks). Each **wrapper invocation** gets its own id (default `analyze_array_<timestamp>_<pid>` or `experiments_array_<timestamp>_<pid>`), with:

- `logs/<JOB_LOG_ID>/wrapper.log` — tee of everything printed on the terminal  
- `outputs/<JOB_LOG_ID>/slurm_%A_%a.out` / `.err` — Slurm stdout/stderr per array task (each `sbatch` in a resubmit loop uses the same folder so filenames stay unique by job id)

For the **analyze** wrapper, the default **`RESULTS_ROOT`** is `SlabDesignFactors/results/remote_results_full_sweep_<JOB_LOG_ID>/` so results share the same trace id as `logs/` and `outputs/`. Override with **`RESULTS_ROOT=...`** for a fixed path (e.g. resume into an existing tree).

Set **`JOB_LOG_ID=my_label`** to choose the id (and thus the default results folder suffix); set **`LOG_FILE`** to override the wrapper log path.

**Direct** `sbatch` (no wrapper) writes flat files: `outputs/slurm_analyze_%A_%a.out` or `outputs/slurm_experiments_%A_%a.out` (still under repo `outputs/`).

- Experiments: `bash SlabDesignFactors/executable/executable_experiments_array_wrapper.bash`  
- Analyze: `bash SlabDesignFactors/executable/executable_analyze_array_wrapper.bash` (also runs merge when complete; see script for `MERGED_DIR` / `merged_<run>.done`)

Other env overrides are printed at the start of each wrapper.

---

## Tunables

- `executable_analyze_array.slurm`: `--array=1-16`, `--cpus-per-task=8`  
- `executable_experiments_array.slurm`: `--array=1-6` (must cover the number of studies in `STUDY_LIST_CSV`)  
- Julia sets BLAS threads to 1 inside the runners to avoid oversubscription.
